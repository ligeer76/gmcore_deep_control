! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This is the main module of GMCORE, which provides initialization, run,
!   finalization subroutines.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
!   - Jianghao Li
! ==============================================================================

module gmcore_mod

  use mpi
  use flogger
  use string
  use const_mod
  use namelist_mod
  use latlon_parallel_mod
  use process_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use restart_mod
  use prepare_mod
  use block_mod
  use vert_coord_mod
  use time_schemes_mod
  use operators_mod
  use tracer_mod
  use interp_mod
  use gas_mod
  use adv_mod
  use formula_mod
  use pgf_mod
  use nh_mod
  use damp_mod
  use physics_mod
  use filter_mod
  use test_forcing_mod
  use perf_mod
  use regrid_mod
  use latlon_operators_mod

  implicit none

  private

  public gmcore_init
  public gmcore_init_stage0
  public gmcore_init_stage1
  public gmcore_init_stage2
  public gmcore_init_stage3
  public gmcore_run
  public gmcore_final

  public adv_accum_wind
  public block_type
  public dstate_type
  public dtend_type
  public proc

  procedure(space_operators_interface), pointer :: operators

contains

  subroutine gmcore_init(namelist_path, comm)

    character(*), intent(in) :: namelist_path
    integer, intent(in), optional :: comm

    call gmcore_init_stage0()
    call gmcore_init_stage1(namelist_path)
    call gmcore_init_stage2(namelist_path)
    call gmcore_init_stage3()

  end subroutine gmcore_init

  ! ============================================================================
  ! Description:
  !
  !   This subroutine runs some basic initialization procedures, including setup
  !   MPI environment and constants. Users can modify the constants after it.
  ! ============================================================================

  subroutine gmcore_init_stage0(comm)

    integer, intent(in), optional :: comm

    logical is_exist

    call log_init()
    call gas_mixture_init(planet)
    call const_init(planet)
    call time_scheme_init()
    call time_init(dt_dyn)
    call process_init(comm)

    if (use_async_io .and. nproc_io > 0 .and. nproc_io * proc_io_stride /= proc%np_io) then
      if (proc%is_root()) call log_error('Parameter nproc_io * proc_io_stride /= proc%np_io!')
    end if

#ifdef FC_IS_INTEL
    inquire(directory=gmcore_root, exist=is_exist)
#else
    inquire(file=gmcore_root, exist=is_exist)
#endif
    if (.not. is_exist) then
      if (proc%is_root()) call log_error('Cannot find the root directory of GMCORE!')
    end if

    if (proc%is_root()) then
      select case (planet)
      case ('mars')
        print *, ''
        print *, '  ██████            ███      ███                            '
        print *, ' ██    ██           ████    ████                            '
        print *, '██         ███████  ██ ██  ██ ██  ██████   ██ ████  ███████ '
        print *, '██  █████ ██     ██ ██  ████  ██       ██  ███     ██       '
        print *, '██     ██ ██     ██ ██   ██   ██  ███████  ██       ███████ '
        print *, ' ██   ███ ██     ██ ██        ██ ██    ██  ██            ░██'
        print *, '  █████    ███████  ██        ██  █████ ██ ██       ███████ '
      case default
        print *, ''
        print *, '  ██████  ███      ███   ██████    ██████   █████████  ██████████'
        print *, ' ██    ██ ████    ████  ██    ██  ██    ██  ██      ██ ██        '
        print *, '██        ██ ██  ██ ██ ██        ██      ██ ██      ██ ██        '
        print *, '██  █████ ██  ████  ██ ██        ██      ██ █████████  █████████ '
        print *, '██     ██ ██   ██   ██ ██        ██      ██ ██    ██   ██        '
        print *, ' ██   ███ ██        ██  ██    ██  ██    ██  ██     ██  ██        '
        print *, '  █████   ██        ██   ██████    ██████   ██      ██ ██████████'
      end select
    end if

  end subroutine gmcore_init_stage0

  ! ============================================================================
  ! Description:
  !
  !   This subroutine creates blocks which include dynamics mesh, state, tend
  !   objects. Users can prepare some surface static data, such as topography,
  !   after it, and add tracers.
  ! ============================================================================

  subroutine gmcore_init_stage1(namelist_path)

    character(*), intent(in) :: namelist_path

    call global_mesh%init_global(nlon, nlat, nlev, lon_hw=lon_hw, lat_hw=lat_hw, lev_hw=3)
    call vert_coord_init(namelist_path)
    call process_create_blocks()
    if (proc%is_model()) then
      associate (mesh => blocks(1)%mesh)
      min_lon = mesh%full_lon_deg(mesh%full_ims)
      max_lon = mesh%full_lon_deg(mesh%full_ime)
      min_lat = mesh%full_lat_deg(max(mesh%full_jms, 1))
      max_lat = mesh%full_lat_deg(min(mesh%full_jme, global_mesh%full_nlat))
      end associate
    end if
    call damp_init()
    call tracer_init_stage1()
    call physics_init_stage1(namelist_path)
    call history_init_stage1()

  end subroutine gmcore_init_stage1

  ! ============================================================================
  ! Description:
  !
  !   This subroutine initialize vertical coordinate which may depend on the
  !   topography. The tracers are also allocated, so users should already add
  !   the necessary tracer species.
  ! ============================================================================

  subroutine gmcore_init_stage2(namelist_path)

    character(*), intent(in) :: namelist_path

    integer iblk

    call physics_init_stage2(namelist_path)
    call pgf_init()
    call interp_init()
    call operators_init()
    call tracer_init_stage2()
    call adv_init()
    call restart_init_stage2()
    call history_init_stage2()

    operators => space_operators

    if (proc%is_root()) call print_namelist()

    do iblk = 1, size(blocks)
      blocks(iblk)%mesh%full_lev  = global_mesh%full_lev
      blocks(iblk)%mesh%half_lev  = global_mesh%half_lev
      blocks(iblk)%mesh%full_dlev = global_mesh%full_dlev
      blocks(iblk)%mesh%half_dlev = global_mesh%half_dlev
      call blocks(iblk)%static%init_stage2(blocks(iblk)%mesh)
    end do

  end subroutine gmcore_init_stage2

  ! ============================================================================
  ! Description:
  !
  !   This subroutine runs some initializations that need to be after initial
  !   conditions.
  ! ============================================================================

  subroutine gmcore_init_stage3()

    character(10) time_value, time_units
    real(r8) seconds

    time_value = split_string(print_interval, ' ', 1)
    time_units = split_string(print_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days', 'sol')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid print interval ' // trim(print_interval) // '!')
    end select
    call time_add_alert('print', seconds=seconds)

    call regrid_init()
    if (proc%is_model()) then
      if (.not. advection) call operators_prepare(blocks, old, dt_dyn)
      call physics_init_stage3()
    end if
    call history_init_stage3()

  end subroutine gmcore_init_stage3

  subroutine gmcore_run()

    integer iblk, itime

    if (proc%is_model()) then
      do iblk = 1, size(blocks)
        associate (block  => blocks(iblk)     , &
                   mesh   => blocks(iblk)%mesh, &
                   dstate => blocks(iblk)%dstate(old))
        if (baroclinic) then
          ! Ensure bottom gz_lev is the same as gzs.
          do itime = lbound(block%dstate, 1), ubound(block%dstate, 1)
            call block%dstate(itime)%gz_lev%copy(block%static%gzs, k=mesh%half_kde)
            call fill_halo(block%dstate(itime)%gz_lev)
          end do
        end if
        if (.not. advection) call wind_c2a_operator(block%dstate(old)%u_lon, block%dstate(old)%v_lat, block%aux%u, block%aux%v)
        call calc_div(block, dstate, .false.)
        if (baroclinic .and. .not. restart) call block%static%ref_ps%copy(dstate%mgs, with_halo=.true.)
        end associate
      end do
    end if

    if (proc%is_model()) then
      call prepare_run()
      call pgf_init_after_ic()
      call operators_prepare(blocks, old, dt_dyn)
      call adv_prepare(old)
      call diagnose(blocks, old)
      call regrid_run(old)
      if (proc%is_root()) call log_print_diag(curr_time%isoformat())
    end if
    call output(old)

    model_main_loop: do while (.not. time_is_finished())
      if (proc%is_model()) then
        ! ----------------------------------------------------------------------
        !                              Dynamical Core
        do iblk = 1, size(blocks)
          call time_integrator(operators, blocks(iblk), old, new, dt_dyn)
          call damp_run(blocks(iblk), blocks(iblk)%dstate(new), dt_dyn)
          call physics_update_after_dynamics(blocks(iblk), new, dt_dyn)
        end do
        ! ----------------------------------------------------------------------
        ! Advance to n+1 time level.
        ! NOTE: Time indices are swapped, e.g. new <=> old.
        call time_advance(dt_dyn)
        ! ----------------------------------------------------------------------
        !                            Tracer Advection
        call adv_run_tracers(new, old)
        ! ----------------------------------------------------------------------
        !                                Physics
        call test_forcing_run(dt_dyn, old)
        if (baroclinic) then
          do iblk = 1, size(blocks)
            call prepare_physics(blocks(iblk), old)
            call physics_run(blocks(iblk), old, dt_phys)
            call physics_update_after_physics(blocks(iblk), old, dt_phys)
          end do
        end if
        ! ----------------------------------------------------------------------
        if (time_is_alerted('print')) then
          call diagnose(blocks, old)
          if (proc%is_root()) call log_print_diag(curr_time%isoformat())
        end if
        call blocks(1)%accum(old)
      else
        call time_advance(dt_dyn)
      end if
      call output(old)
    end do model_main_loop

  end subroutine gmcore_run

  subroutine gmcore_final()

    call log_final()
    call time_final()
    call interp_final()
    call gas_mixture_final()
    call vert_coord_final()
    call tracer_final()
    call pgf_final()
    call adv_final()
    call damp_final()
    call history_final()
    call perf_final()
    call regrid_final()
    call process_final()

  end subroutine gmcore_final

  subroutine output(itime)

    integer, intent(in) :: itime

    logical, save :: first_call = .true.
    real(8), save :: time1, time2

    if (first_call .or. time_is_alerted('history_write')) then
      if (first_call) time1 = MPI_WTIME()
      time2 = MPI_WTIME()
      if (.not. first_call) then
        if (proc%is_root()) call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      if (output_h0) call history_write_h0(itime)
      if (output_h1) call history_write_h1(itime)
      if (output_h2) then
        call regrid_run(itime)
        call history_write_h2()
      end if
    end if
    if (time_is_alerted('restart_write') .or. (time_step == 0 .and. restart_interval /= 'N/A')) then
      call restart_write(itime)
    end if
    first_call = .false.

  end subroutine output

  subroutine diagnose(blocks, itime)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime

    integer i, j, k, m, iblk
    real(r8) tm, te, tpe, tpt, tqv, tqm, max_w
    real(r8) te_ke, te_ie, te_pe
    real(r8) area_cell, area_lon, area_lat, area_sfc, rsfc
    real(r8) g_cell, g_lon, g_lat, g_sfc
    real(r8), save :: tm_ref = 0, te_ref = 0, tpe_ref = 0
    real(r8), save :: tqv_ref = 0, tqm_ref = 0, tpt_ref = 0
    logical, save :: reset_ref = .true.

    tm    = 0
    te    = 0
    tpe   = 0
    tpt   = 0
    tqv   = 0
    tqm   = 0
    te_ke = 0
    te_ie = 0
    te_pe = 0
    do iblk = 1, size(blocks)
      associate (block   => blocks(iblk)                    , &
                 mesh    => blocks(iblk)%mesh               , &
                 gzs     => blocks(iblk)%static%gzs         , &
                 mgs     => blocks(iblk)%dstate(itime)%mgs  , &
                 dmg     => blocks(iblk)%dstate(itime)%dmg  , &
                 dmg_lon => blocks(iblk)%aux%dmg_lon        , &
                 dmg_lat => blocks(iblk)%aux%dmg_lat        , &
                 u_lon   => blocks(iblk)%dstate(itime)%u_lon, &
                 v_lat   => blocks(iblk)%dstate(itime)%v_lat, &
                 tv      => blocks(iblk)%aux%tv             , &
                 pt      => blocks(iblk)%dstate(itime)%pt   , &
                 pv_lon  => blocks(iblk)%aux%pv_lon         , &
                 pv_lat  => blocks(iblk)%aux%pv_lat         , &
#ifdef USE_DEEP_ATM
                 rdp     => blocks(iblk)%aux%rdp            , &
                 rdp_lon => blocks(iblk)%aux%rdp_lon        , &
                 rdp_lat => blocks(iblk)%aux%rdp_lat        , &
#endif
                 q       => tracers(iblk)%q                 )
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            area_cell = mesh%area_cell(j)
            g_cell = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) area_cell = area_cell * (rdp%d(i,j,k) / radius)**2
            if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_cell = gravity_from_radius(rdp%d(i,j,k))
#endif
            tm = tm + dmg%d(i,j,k) * area_cell / g_cell
            if (idx_qv > 0) tqv = tqv + dmg%d(i,j,k) * q%d(i,j,k,idx_qv) * area_cell / g_cell
            if (ntracers_water > 0) then
              do m = 1, ntracers
                if (is_water_tracer(m)) tqm = tqm + dmg%d(i,j,k) * q%d(i,j,k,m) * area_cell / g_cell
              end do
            end if
          end do
        end do
      end do

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            area_lon = mesh%area_lon(j)
            g_lon = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) area_lon = area_lon * (rdp_lon%d(i,j,k) / radius)**2
            if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_lon = gravity_from_radius(rdp_lon%d(i,j,k))
#endif
            te_ke = te_ke + dmg_lon%d(i,j,k) * 0.5_r8 * u_lon%d(i,j,k)**2 * area_lon * 2 / g_lon
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            area_lat = mesh%area_lat(j)
            g_lat = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) area_lat = area_lat * (rdp_lat%d(i,j,k) / radius)**2
            if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_lat = gravity_from_radius(rdp_lat%d(i,j,k))
#endif
            te_ke = te_ke + dmg_lat%d(i,j,k) * 0.5_r8 * v_lat%d(i,j,k)**2 * area_lat * 2 / g_lat
          end do
        end do
      end do
      if (baroclinic) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              area_cell = mesh%area_cell(j)
              g_cell = g
#ifdef USE_DEEP_ATM
              if (deepwater .and. use_mesh_change) area_cell = area_cell * (rdp%d(i,j,k) / radius)**2
              if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_cell = gravity_from_radius(rdp%d(i,j,k))
#endif
              te_ie = te_ie + dmg%d(i,j,k) * cpd * tv%d(i,j,k) * area_cell / g_cell
            end do
          end do
        end do
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            area_sfc = mesh%area_cell(j)
            g_sfc = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) then
              if (use_variable_gravity) then
                rsfc = radius_from_geopotential(gzs%d(i,j))
                g_sfc = gravity_from_radius(rsfc)
              else
                rsfc = radius + gzs%d(i,j) / g
              end if
              area_sfc = area_sfc * (rsfc / radius)**2
            end if
#endif
            te_pe = te_pe + gzs%d(i,j) * mgs%d(i,j) * area_sfc / g_sfc
          end do
        end do
      else
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            area_sfc = mesh%area_cell(j)
            g_sfc = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) then
              if (use_variable_gravity) then
                rsfc = radius_from_geopotential(gzs%d(i,j))
                g_sfc = gravity_from_radius(rsfc)
              else
                rsfc = radius + gzs%d(i,j) / g
              end if
              area_sfc = area_sfc * (rsfc / radius)**2
            end if
#endif
            te_pe = te_pe + (dmg%d(i,j,1)**2 * g * 0.5_r8 + dmg%d(i,j,1) * gzs%d(i,j)) * area_sfc / g_sfc
          end do
        end do
      end if

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            area_lon = mesh%area_lon(j)
            g_lon = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) area_lon = area_lon * (rdp_lon%d(i,j,k) / radius)**2
            if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_lon = gravity_from_radius(rdp_lon%d(i,j,k))
#endif
            tpe = tpe + dmg_lon%d(i,j,k) * pv_lon%d(i,j,k)**2 * 0.5_r8 * area_lon * 2 / g_lon
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            area_lat = mesh%area_lat(j)
            g_lat = g
#ifdef USE_DEEP_ATM
            if (deepwater .and. use_mesh_change) area_lat = area_lat * (rdp_lat%d(i,j,k) / radius)**2
            if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_lat = gravity_from_radius(rdp_lat%d(i,j,k))
#endif
            tpe = tpe + dmg_lat%d(i,j,k) * pv_lat%d(i,j,k)**2 * 0.5_r8 * area_lat * 2 / g_lat
          end do
        end do
      end do

      if (baroclinic) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              area_cell = mesh%area_cell(j)
              g_cell = g
#ifdef USE_DEEP_ATM
              if (deepwater .and. use_mesh_change) area_cell = area_cell * (rdp%d(i,j,k) / radius)**2
              if (deepwater .and. use_mesh_change .and. use_variable_gravity) g_cell = gravity_from_radius(rdp%d(i,j,k))
#endif
              tpt = tpt + dmg%d(i,j,k) * pt%d(i,j,k) * area_cell / g_cell
            end do
          end do
        end do
      end if
      end associate
    end do
                    tm    = global_sum(proc%comm_model, tm   )
    if (idx_qv > 0) tqv   = global_sum(proc%comm_model, tqv  )
    if (ntracers_water > 0) tqm = global_sum(proc%comm_model, tqm)
                    te_ke = global_sum(proc%comm_model, te_ke)
                    te_ie = global_sum(proc%comm_model, te_ie)
                    te_pe = global_sum(proc%comm_model, te_pe)
                    tpe   = global_sum(proc%comm_model, tpe  )
    if (baroclinic) tpt   = global_sum(proc%comm_model, tpt  )
    te = te_ke + te_ie + te_pe

    if (reset_ref .or. time_step == 0) then
      tm_ref  = tm
      tqv_ref = tqv
      tqm_ref = tqm
      tpt_ref = tpt
      te_ref  = te
      tpe_ref = tpe
      reset_ref = .false.
    end if

    do iblk = 1, size(blocks)
      blocks(iblk)%dstate(itime)%tm  = tm
      blocks(iblk)%dstate(itime)%tqv = tqv
      blocks(iblk)%dstate(itime)%tqm = tqm
      blocks(iblk)%dstate(itime)%tpt = tpt
      blocks(iblk)%dstate(itime)%te  = te
      blocks(iblk)%dstate(itime)%tpe = tpe
      blocks(iblk)%dstate(itime)%te_ke = te_ke
      blocks(iblk)%dstate(itime)%te_ie = te_ie
      blocks(iblk)%dstate(itime)%te_pe = te_pe
      blocks(iblk)%dstate(itime)%tm_relerr  = relative_drift(tm , tm_ref )
      blocks(iblk)%dstate(itime)%tqv_relerr = relative_drift(tqv, tqv_ref)
      blocks(iblk)%dstate(itime)%tqm_relerr = relative_drift(tqm, tqm_ref)
      blocks(iblk)%dstate(itime)%tpt_relerr = relative_drift(tpt, tpt_ref)
      blocks(iblk)%dstate(itime)%te_relerr  = relative_drift(te , te_ref )
      blocks(iblk)%dstate(itime)%tpe_relerr = relative_drift(tpe, tpe_ref)
    end do

    if (planet == 'mars') call log_add_diag('ls', curr_time%solar_longitude()*deg)
                          call log_add_diag('tm' , tm )
                          call log_add_diag('tm_relerr', relative_drift(tm, tm_ref))
    if (idx_qv > 0      ) call log_add_diag('tqv', tqv)
    if (idx_qv > 0      ) call log_add_diag('tqv_relerr', relative_drift(tqv, tqv_ref))
    if (ntracers_water > 0) call log_add_diag('tqm', tqm)
    if (ntracers_water > 0) call log_add_diag('tqm_relerr', relative_drift(tqm, tqm_ref))
    if (baroclinic      ) call log_add_diag('tpt', tpt)
    if (baroclinic      ) call log_add_diag('tpt_relerr', relative_drift(tpt, tpt_ref))
                          call log_add_diag('te' , te )
                          call log_add_diag('te_relerr', relative_drift(te, te_ref))
                          call log_add_diag('tpe', tpe)
                          call log_add_diag('tpe_relerr', relative_drift(tpe, tpe_ref))

  end subroutine diagnose

  pure elemental real(r8) function relative_drift(value, reference) result(res)

    real(r8), intent(in) :: value
    real(r8), intent(in) :: reference

    if (abs(reference) <= tiny(reference)) then
      if (abs(value) <= tiny(value)) then
        res = 0
      else
        res = sign(huge(value), value)
      end if
    else
      res = (value - reference) / reference
    end if

  end function relative_drift

  subroutine space_operators(block, old_dstate, star_dstate, new_dstate, dtend, dt, pass, substep)

    type(block_type ), intent(inout) :: block
    type(dstate_type), intent(in   ) :: old_dstate
    type(dstate_type), intent(inout) :: star_dstate
    type(dstate_type), intent(inout) :: new_dstate
    type(dtend_type ), intent(inout) :: dtend
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass
    integer, intent(in) :: substep

    integer i, j, k

    call dtend%reset_flags()

    associate (mesh => block%mesh)
    select case (pass)
    case (all_pass)
      call operators_prepare(block, star_dstate, dt, pass, substep)
      if (baroclinic) then
        call calc_grad_mf          (block, star_dstate)
        call calc_dmgsdt           (block, star_dstate, dtend)
        call calc_mfz              (block, star_dstate, dtend)
        call calc_grad_ptf         (block, star_dstate, dtend, dt)
        call calc_grad_ke          (block, star_dstate, dtend, dt)
        call pgf_run               (block, star_dstate, dtend)
        call calc_coriolis         (block, star_dstate, dtend, dt, substep)
        call calc_wedudlev_wedvdlev(block, star_dstate, dtend, dt)

        dtend%update_u   = .true.
        dtend%update_v   = .true.
        dtend%update_mgs = .true.
        dtend%update_pt  = .true.
      else
        call calc_grad_mf          (block, star_dstate)
        call calc_grad_ke          (block, star_dstate, dtend, dt)
        call pgf_run               (block, star_dstate, dtend)
        call calc_coriolis         (block, star_dstate, dtend, dt, substep)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              if (deepwater .and. use_mesh_change .and. use_variable_gravity) then
                dtend%dgzdt%d(i,j,k) = -block%aux%dmf%d(i,j,k) * gravity_from_geopotential(star_dstate%gz%d(i,j,k))
              else
                dtend%dgzdt%d(i,j,k) = -block%aux%dmf%d(i,j,k) * g
              end if
            end do
          end do
        end do

        dtend%update_u  = .true.
        dtend%update_v  = .true.
        dtend%update_gz = .true.
      end if
    case (forward_pass)
      call operators_prepare(block, star_dstate, dt, pass, substep)
      if (baroclinic) then
        call calc_grad_mf          (block, star_dstate)
        call calc_dmgsdt           (block, star_dstate, dtend)
        call calc_mfz              (block, star_dstate, dtend)
        call calc_grad_ptf         (block, star_dstate, dtend, dt)

        dtend%update_mgs = .true.
        dtend%update_pt  = .true.
      else
        call calc_grad_mf          (block, star_dstate)

        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              if (deepwater .and. use_mesh_change .and. use_variable_gravity) then
                dtend%dgzdt%d(i,j,k) = -block%aux%dmf%d(i,j,k) * gravity_from_geopotential(star_dstate%gz%d(i,j,k))
              else
                dtend%dgzdt%d(i,j,k) = -block%aux%dmf%d(i,j,k) * g
              end if
            end do
          end do
        end do

        dtend%update_gz = .true.
      end if
    case (backward_pass)
      if (nonhydrostatic) call nh_solve(block, old_dstate, star_dstate, new_dstate, dt)
      call operators_prepare(block, new_dstate, dt, pass, substep)
      if (baroclinic) then
        call calc_grad_ke          (block, star_dstate, dtend, dt)
        call pgf_run               (block, new_dstate , dtend)
        call calc_coriolis         (block, star_dstate, dtend, dt, substep)
#ifdef USE_DEEP_ATM
        if (deepwater .and. use_hor_nct) call calc_nct_coriolis(block, star_dstate, dtend, dt)
#endif
        call calc_wedudlev_wedvdlev(block, star_dstate, dtend, dt)

        dtend%update_u   = .true.
        dtend%update_v   = .true.
      else
        call calc_grad_ke          (block, star_dstate, dtend, dt)
        call pgf_run               (block, new_dstate , dtend)
        call calc_coriolis         (block, star_dstate, dtend, dt, substep)

        dtend%update_u  = .true.
        dtend%update_v  = .true.
      end if
    end select
    end associate

  end subroutine space_operators

  subroutine prepare_physics(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    call calc_grad_mf(block, block%dstate(itime))
    call calc_div    (block, block%dstate(itime), .false.)
    call calc_omg    (block, block%dstate(itime))
    call wind_c2a_operator(block%dstate(itime)%u_lon, block%dstate(itime)%v_lat, block%aux%u, block%aux%v)

  end subroutine prepare_physics

end module gmcore_mod
