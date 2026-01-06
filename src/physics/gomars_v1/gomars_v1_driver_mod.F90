! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module gomars_v1_driver_mod

  use datetime
  use const_mod
  use formula_mod
  use namelist_mod, only: restart, ptop
  use physics_parallel_mod
  use process_mod, only: proc
  use gomars_v1_namelist_mod
  use gomars_v1_objects_mod
  use gomars_v1_output_mod
  use gomars_v1_tracers_mod
  use gomars_v1_orbit_mod
  use gomars_v1_rad_mod
  use gomars_v1_pbl_mod
  use gomars_v1_lsm_mod
  use gomars_v1_co2cyc_mod
  use gomars_v1_dustcyc_mod
  use gomars_v1_mp_mod
  use gomars_v1_cnvadj_mod
  use gomars_v1_damp_mod
  use gomars_v1_utils_mod

  implicit none

  private

  public gomars_v1_init_stage2
  public gomars_v1_init_stage3
  public gomars_v1_run
  public gomars_v1_final
  public gomars_v1_d2p
  public gomars_v1_p2d
  public gomars_v1_add_output
  public gomars_v1_output
  public objects

contains

  ! ============================================================================
  ! Description:
  !
  !   This subroutine initialize vertical coordinate which may depend on the
  !   topography. The tracers are also allocated, so users should already add
  !   the necessary tracer species.
  ! ============================================================================

  subroutine gomars_v1_init_stage2(namelist_path, mesh, dt_adv, dt_phys, input_ngroup, model_root)

    character(*), intent(in) :: namelist_path
    type(physics_mesh_type), intent(in), target :: mesh(:)
    real(r8), intent(in) :: dt_adv
    real(r8), intent(in) :: dt_phys
    integer , intent(in) :: input_ngroup
    character(*), intent(in), optional :: model_root

    dt         = dt_phys
    dt_mp      = dt_phys / nsplit
    nlev       = mesh(1)%nlev
    nlayrad    = nlev + 1
    nlevrad    = nlev + 2
    ptrop      = ptop
    pstrat     = ptrop * 0.5_r8
    lnpstrat   = log(pstrat)
    pstratk    = (pstrat / p0)**rd_o_cpd
    tsat_strat = dewpoint_temperature_mars(pstrat)

    if (present(model_root)) data_root = trim(model_root) // '/data/mars'

    call gomars_v1_parse_namelist(namelist_path)
    call gomars_v1_tracers_init(dt_adv)
    call gomars_v1_objects_init(mesh)
    call gomars_v1_read_static_data()
    call gomars_v1_orbit_init      ()
    call gomars_v1_rad_init        ()
    call gomars_v1_lsm_init        ()
    call gomars_v1_pbl_init        ()
    call gomars_v1_mp_init         ()
    call gomars_v1_damp_init       ()

  end subroutine gomars_v1_init_stage2

  ! ============================================================================
  ! Description:
  !
  !   This subroutine runs some initializations that need to be after initial
  !   conditions.
  ! ============================================================================

  subroutine gomars_v1_init_stage3()

    integer iblk, icol, k, l, is, ig, n

    ! Check model top pressure. It must be lower than ptop in nasa_rad_mod.

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state, &
                 tend  => objects(iblk)%tend )
      ! Set some initial variables (from init1).
      if (.not. restart) then
        ! FIXME: Here is for cold run.
        do icol = 1, mesh%ncol
          state%tstrat    (icol) = state%t_bot(icol)
          state%co2ice_sfc(icol) = 0
          state%tg        (icol) = state%t_bot(icol)
        end do
      end if
      if (.not. restart) then
        do icol = 1, mesh%ncol
          state%irflx_sfc_dn(icol) = 1
          state%vsdif_sfc_dn(icol) = 0
          do is = 1, nspectv
            do ig = 1, ngauss
              state%detau(icol,is,ig) = 0.1_r8
            end do
          end do
        end do
      end if
      call gomars_v1_diags(state, tend)
      ! Initialize local time.
      do icol = 1, mesh%ncol
        state%local_time(icol) = mesh%lon(icol) / pi2 * 24
      end do
      end associate
    end do

  end subroutine gomars_v1_init_stage3

  subroutine gomars_v1_read_static_data()

    use fiona
    use latlon_interp_mod

    real(r8), allocatable :: lon(:), ilon(:)
    real(r8), allocatable :: lat(:), ilat(:)
    real(r8), allocatable :: dat(:,:)
    real(r8), allocatable :: col(:)
    real(r8) dlon, dlat
    integer i, j, k

    if (allocated(objects)) then
      associate (mesh => objects(1)%mesh, state => objects(1)%state)
      ! Surface albedo
      call fiona_open_dataset('alb', file_path=trim(data_root)//'/mgs_albedo.nc', mpi_comm=proc%comm_model)
      ! call fiona_open_dataset('alb', file_path=trim(data_root)//'/osu_alb_6x5_2011.nc', mpi_comm=proc%comm_model)
      call fiona_set_dim('alb', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('alb', 'lat', span=[-90, 90])
      call fiona_start_input('alb')
      call fiona_input_range('alb', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('alb', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('alb', 'albedo', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      ! call fiona_input_range('alb', 'alb', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('alb')
      call latlon_interp_bilinear_column(lon, lat, dat, mesh%lon, mesh%lat, state%alsp)
      call physics_pole_sum(mesh%lat, state%alsp)
      deallocate(lon, lat, dat)

      ! Surface thermal inertia
      call fiona_open_dataset('zin', file_path=trim(data_root)//'/osu_zin_6x5_2011.nc', mpi_comm=proc%comm_model)
      ! call fiona_open_dataset('zin', file_path=trim(data_root)//'/nasa_zin_2011.nc', mpi_comm=proc%comm_model)
      call fiona_set_dim('zin', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('zin', 'lat', span=[-90, 90])
      call fiona_start_input('zin')
      call fiona_input_range('zin', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('zin', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('zin', 'zin', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('zin')
      call latlon_interp_bilinear_column(lon, lat, dat, mesh%lon, mesh%lat, state%zin(:,1))
      do k = 2, nsoil
        state%zin(:,k) = state%zin(:,1)
      end do
      call physics_pole_sum(mesh%lat, state%zin)
      deallocate(lon, lat, dat)

      ! Flag of northern polar cap of water ice
      call fiona_open_dataset('npcflag', file_path=trim(data_root)//'/npcflag_osu_550.nc', mpi_comm=proc%comm_model)
      call fiona_set_dim('npcflag', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('npcflag', 'lat', span=[-90, 90])
      call fiona_start_input('npcflag')
      call fiona_input_range('npcflag', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('npcflag', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('npcflag', 'npcflag', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('npcflag')
      allocate(col(mesh%ncol))
      call latlon_interp_bilinear_column(lon, lat, dat, mesh%lon, mesh%lat, col)
      do i = 1, mesh%ncol
        state%npcflag(i) = col(i) > 0.5_r8
      end do
      deallocate(lon, lat, dat, col)

      ! Flag of GRS ground ice
      call fiona_open_dataset('gnd_ice', file_path=trim(data_root)//'/gnd_ice_map.nc', mpi_comm=proc%comm_model)
      call fiona_set_dim('gnd_ice', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('gnd_ice', 'lat', span=[-90, 90])
      call fiona_start_input('gnd_ice')
      call fiona_input_range('gnd_ice', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('gnd_ice', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('gnd_ice', 'gnd_ice', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('gnd_ice')
      allocate(col(mesh%ncol))
      call latlon_interp_bilinear_column(lon, lat, dat, mesh%lon, mesh%lat, col)
      do i = 1, mesh%ncol
        if (mesh%lat(i) > 0) then
          state%grsn(i) = col(i) > 0.5_r8
        else
          state%grss(i) = col(i) > 0.5_r8
        end if
      end do
      deallocate(lon, lat, dat, col)

      ! Zonal averaged ground temperature for cold run
      call fiona_open_dataset('zavgtg', file_path=trim(data_root)//'/nasa_zavgtg.nc', mpi_comm=proc%comm_model)
      call fiona_set_dim('zavgtg', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('zavgtg', 'lat', span=[-90, 90])
      call fiona_start_input('zavgtg')
      call fiona_input_range('zavgtg', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('zavgtg', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('zavgtg', 'zavgtg', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('zavgtg')
      call latlon_interp_bilinear_column(lon, lat, dat, mesh%lon, mesh%lat, state%zavgtg)
      call physics_pole_sum(mesh%lat, state%zavgtg)
      deallocate(lon, lat, dat)
      end associate
    end if

  end subroutine gomars_v1_read_static_data

  subroutine gomars_v1_run(time)

    type(datetime_type), intent(in) :: time

    integer iblk, icyc
    real(r8) ls

    ls = time%solar_longitude()
    time_of_day = time%time_of_day()
    call update_solar_decl_angle(ls)
    call update_toa_solar_flux(ls)

    do iblk = 1, size(objects)
      associate (mesh  => objects(iblk)%mesh , &
                 state => objects(iblk)%state, &
                 tend  => objects(iblk)%tend )
      call tend%reset()
      state%local_time = mod(state%local_time + (dt / mars_time_scale / 3600.0_r8), 24.0_r8)
      call gomars_v1_orbit_cosz     (state)
      call update_sfc_solar_flux    (state)
      call gomars_v1_lsm_run        (state, tend)
      call gomars_v1_orbit_cosz_avg (state)
      call gomars_v1_co2cyc_run     (state, tend)
      call interp_temperature       (state)
      call gomars_v1_rad_run        (state)
      call gomars_v1_pbl_run        (state)
      call interp_temperature       (state)
      do icyc = 1, nsplit
        call gomars_v1_dustcyc_run  (state, dt / nsplit)
        call gomars_v1_mp_run       (state, dt / nsplit)
      end do
      call gomars_v1_cnvadj_run     (state)
      call gomars_v1_damp_run       (state)
      call gomars_v1_diags          (state, tend)
      end associate
    end do

  end subroutine gomars_v1_run

  subroutine gomars_v1_final()

    call gomars_v1_objects_final()
    call gomars_v1_rad_final    ()
    call gomars_v1_lsm_final    ()
    call gomars_v1_pbl_final    ()
    call gomars_v1_mp_final     ()
    call gomars_v1_damp_final   ()

  end subroutine gomars_v1_final

  subroutine gomars_v1_d2p()

    integer iblk, i, k

    ! Calculate ts.

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      ! Change to use surface pressure as reference.
      do k = 1, mesh%nlev
        do i = 1, mesh%ncol
          state%pk(i,k) = (state%p(i,k) / state%ps(i))**rd_o_cpd
        end do
      end do
      do k = 1, mesh%nlev + 1
        do i = 1, mesh%ncol
          state%pk_lev(i,k) = (state%p_lev(i,k) / state%ps(i))**rd_o_cpd
        end do
      end do
      end associate
    end do

  end subroutine gomars_v1_d2p

  subroutine gomars_v1_p2d()

    integer iblk, i, k, m

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state, tend => objects(iblk)%tend)
      do k = 1, mesh%nlev
        do i = 1, mesh%ncol
          tend%dudt(i,k) = (state%u(i,k) - state%u_old(i,k)) / dt
          tend%dvdt(i,k) = (state%v(i,k) - state%v_old(i,k)) / dt
          tend%dtdt(i,k) = (state%t(i,k) - state%t_old(i,k)) / dt
        end do
      end do
      do m = 1, ntracers
        do k = 1, mesh%nlev
          do i = 1, mesh%ncol
            tend%dqdt(i,k,m) = (state%q(i,k,m) - state%q_old(i,k,m)) / dt
          end do
        end do
      end do
      tend%updated_ps = .true.
      tend%updated_u  = .true.
      tend%updated_v  = .true.
      tend%updated_t  = .true.
      tend%updated_q  = .true.
      end associate
    end do

  end subroutine gomars_v1_p2d

  subroutine gomars_v1_diags(state, tend)

    type(gomars_v1_state_type), intent(inout) :: state
    type(gomars_v1_tend_type), intent(in) :: tend

    integer i, k

    associate (mesh       => state%mesh      , &
               ps         => state%ps        , & ! in
               dp         => state%dp        , & ! in
               dpsdt      => tend %dpsdt     , & ! in
               co2ice_sfc => state%co2ice_sfc, & ! in
               q          => state%q         , & ! in
               qsfc       => state%qsfc      , & ! in
               tm_co2     => state%tm_co2    , & ! out
               tm_dst     => state%tm_dst    , & ! out
               tm_h2o     => state%tm_h2o    )   ! out
    tm_co2 = 0
    do i = 1, mesh%ncol
      tm_co2 = tm_co2 + mesh%area(i) * (ps(i) / g + dpsdt(i) / g * dt + co2ice_sfc(i))
    end do
    tm_co2 = physics_sum(tm_co2)
    tm_dst = 0
    do k = 1, mesh%nlev
      do i = 1, mesh%ncol
        tm_dst = tm_dst + mesh%area(i) * q(i,k,iMa_dst) * dp(i,k) / g
      end do
    end do
    do i = 1, mesh%ncol
      tm_dst = tm_dst + mesh%area(i) * qsfc(i,iMa_dst)
    end do
    tm_dst = physics_sum(tm_dst)
    tm_h2o = 0
    do k = 1, mesh%nlev
      do i = 1, mesh%ncol
        tm_h2o = tm_h2o + mesh%area(i) * q(i,k,iMa_vap) * dp(i,k) / g
      end do
    end do
    do i = 1, mesh%ncol
      tm_h2o = tm_h2o + mesh%area(i) * qsfc(i,iMa_vap)
    end do
    tm_h2o = physics_sum(tm_h2o)
    ! --------------------------------------------------------------------------
    ! T15
    call diag_t15(state)
    end associate

  end subroutine gomars_v1_diags

  subroutine diag_t15(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k, w
    real(r8) pw, tw, tavg, sum_wgt, sum_area
    real(r8), parameter :: wgt(62) = [                                        &
      0.03500_r8, 0.03500_r8, 0.03600_r8, 0.03600_r8, 0.03700_r8, 0.03700_r8, &
      0.03800_r8, 0.03900_r8, 0.04100_r8, 0.04300_r8, 0.04600_r8, 0.04900_r8, &
      0.05200_r8, 0.05600_r8, 0.06000_r8, 0.06500_r8, 0.07000_r8, 0.07600_r8, &
      0.08200_r8, 0.09000_r8, 0.09700_r8, 0.10500_r8, 0.11500_r8, 0.12500_r8, &
      0.13600_r8, 0.14700_r8, 0.15900_r8, 0.17200_r8, 0.18600_r8, 0.20100_r8, &
      0.21700_r8, 0.23400_r8, 0.25200_r8, 0.27000_r8, 0.28800_r8, 0.30800_r8, &
      0.32600_r8, 0.34300_r8, 0.36000_r8, 0.37200_r8, 0.38200_r8, 0.38800_r8, &
      0.39100_r8, 0.38800_r8, 0.37900_r8, 0.36400_r8, 0.34200_r8, 0.31000_r8, &
      0.27800_r8, 0.24600_r8, 0.21300_r8, 0.17800_r8, 0.14300_r8, 0.10900_r8, &
      0.08100_r8, 0.06100_r8, 0.04500_r8, 0.03200_r8, 0.02200_r8, 0.01500_r8, &
      0.01000_r8, 0.00700_r8]

    associate (mesh => state%mesh, &
               t    => state%t   , & ! in
               p    => state%p   , & ! in
               q    => state%q   , & ! in
               t15  => state%t15 )   ! out
    t15 = 0
    sum_wgt = sum(wgt)
    sum_area = 0
    do i = 1, mesh%ncol
      if (abs(mesh%lat(i)) <= 40 * rad) then
        tavg = 0
        do w = 1, size(wgt)
          pw = exp(w / 10.0_r8 - 4.6_r8) * 100
          if (p(i,1) >= pw) then
            tw = t(i,1)
          else if (pw >= p(i,mesh%nlev)) then
            tw = t(i,mesh%nlev)
          else
            do k = 1, mesh%nlev - 1
              if (p(i,k) <= pw .and. pw <= p(i,k+1)) exit
            end do
            tw = ((p(i,k+1) - pw) * t(i,k) + (pw - p(i,k)) * t(i,k+1)) / (p(i,k+1) - p(i,k))
          end if
          tavg = tavg + wgt(w) * (-0.0181075_r8 - 39.4312_r8 / (tw - 2353.76_r8 - 62644.0_r8 / (tw + 64.445_r8 + 84263.9_r8 / (tw - 185.333_r8))))
        end do
        tavg = tavg / sum_wgt
        t15 = t15 + mesh%area(i) * (881.042_r8 - 2.40183_r8 / (tavg + 0.00364298_r8 - 0.61044e-7_r8 / (tavg + 0.162965e-3_r8 - 0.113959e-8_r8 / (tavg + 0.228812e-4_r8))))
        sum_area = sum_area + mesh%area(i)
      end if
    end do
    t15 = physics_sum(t15) / physics_sum(sum_area)
    end associate

  end subroutine diag_t15

end module gomars_v1_driver_mod
