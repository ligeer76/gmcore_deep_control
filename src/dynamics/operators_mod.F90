module operators_mod

  use const_mod
  use perf_mod
  use vert_coord_mod
  use block_mod
  use latlon_parallel_mod
  use latlon_field_types_mod
  use latlon_operators_mod
  use process_mod, only: process_stop
  use formula_mod
  use namelist_mod
  use tracer_mod
  use pgf_mod
  use adv_mod
  use interp_mod
  use filter_mod

  implicit none

  private

  public operators_init
  public operators_prepare
  public calc_mg
  public calc_ph
  public calc_omg
  public calc_dmg
  public calc_t
  public calc_rhod
#ifdef USE_DEEP_ATM
  public calc_rdp
  public calc_nct_coriolis
  public deep_hamiton_modify_3d
#endif
  public calc_gz
  public calc_mfz
  public calc_div
  public calc_vor
  public calc_coriolis
  public calc_grad_ke
  public calc_grad_mf
  public calc_grad_ptf
  public calc_dmgsdt
  public calc_wedudlev_wedvdlev

  interface operators_prepare
    module procedure operators_prepare_1
    module procedure operators_prepare_2
  end interface operators_prepare

  interface
    subroutine interp_pv_interface(block, dstate, dt, substep)
      import block_type, dstate_type, r8
      type(block_type), intent(inout) :: block
      type(dstate_type), intent(inout) :: dstate
      real(r8), intent(in) :: dt
      integer, intent(in) :: substep
    end subroutine interp_pv_interface
  end interface

  procedure(interp_pv_interface), pointer :: interp_pv => null()

contains

  subroutine operators_init()

    select case (pv_adv_scheme)
    case ('midpoint')
      interp_pv => interp_pv_midpoint
    case ('upwind')
      interp_pv => interp_pv_upwind
    case ('weno')
      interp_pv => interp_pv_weno
    case default
      if (proc%is_root()) call log_error('Invalid pv_scheme ' // trim(pv_adv_scheme) // '!')
    end select

  end subroutine operators_init

  ! First time call before main model loop.

  subroutine operators_prepare_1(blocks, itime, dt)

    type(block_type), intent(inout) :: blocks(:)
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    integer iblk

    do iblk = 1, size(blocks)
      if (baroclinic    ) call calc_mg    (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_dmg                       (blocks(iblk), blocks(iblk)%dstate(itime))
      call tracer_calc_qm                 (blocks(iblk))
      if (baroclinic    ) call calc_ph    (blocks(iblk), blocks(iblk)%dstate(itime))
      if (nonhydrostatic .and. .not. restart) then
        if (sum(blocks(iblk)%dstate(itime)%p%d) == 0) then
          ! Set pressure to hydrostatic pressure in the initial condition.
          blocks(iblk)%dstate(itime)%p    %d = blocks(iblk)%dstate(itime)%ph    %d
          blocks(iblk)%dstate(itime)%p_lev%d = blocks(iblk)%dstate(itime)%ph_lev%d
        end if
      end if
      if (baroclinic    ) call calc_t     (blocks(iblk), blocks(iblk)%dstate(itime))
      !!!
      ! SA version calc_gz need tv and ph,in DA the r calc need to use gz, so we need calc gz in advance
      !!!
      if (baroclinic .and. (hydrostatic .or. init_hydrostatic_gz)) then
        call calc_gz                      (blocks(iblk), blocks(iblk)%dstate(itime))
      end if
      if (baroclinic    ) call calc_rhod  (blocks(iblk), blocks(iblk)%dstate(itime))
#ifdef USE_DEEP_ATM
      ! if (deepwater     ) call calc_rdp (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_rdp (blocks(iblk), blocks(iblk)%dstate(itime))
#endif
      call calc_mf                        (blocks(iblk), blocks(iblk)%dstate(itime), dt)
      call calc_ke                        (blocks(iblk), blocks(iblk)%dstate(itime),     total_substeps)
      call calc_pv                        (blocks(iblk), blocks(iblk)%dstate(itime),     total_substeps)
      call interp_pv                      (blocks(iblk), blocks(iblk)%dstate(itime), dt, total_substeps)
      ! if (baroclinic .and. (hydrostatic .or. init_hydrostatic_gz)) then
      !   call calc_gz                      (blocks(iblk), blocks(iblk)%dstate(itime))
      ! end if
      ! if (baroclinic    ) call calc_rhod  (blocks(iblk), blocks(iblk)%dstate(itime))
    end do

  end subroutine operators_prepare_1

  subroutine operators_prepare_2(block, dstate, dt, pass, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass
    integer, intent(in) :: substep

    select case (pass)
    case (all_pass)
      call calc_mf                        (block, dstate, dt)
      call calc_ke                        (block, dstate,     substep)
      call calc_pv                        (block, dstate,     substep)
      call interp_pv                      (block, dstate, dt, substep)
      if (baroclinic    ) call calc_t     (block, dstate)
      if (hydrostatic   ) call calc_gz    (block, dstate)
      if (baroclinic    ) call calc_rhod  (block, dstate)
#ifdef USE_DEEP_ATM
      ! if (deepwater     ) call calc_rdp (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_rdp (block, dstate)
#endif
    case (forward_pass)
      call calc_mf                        (block, dstate, dt)
      call calc_ke                        (block, dstate,     substep)
      call calc_pv                        (block, dstate,     substep)
      call interp_pv                      (block, dstate, dt, substep)
    case (backward_pass)
      if (baroclinic    ) call calc_t     (block, dstate)
      if (hydrostatic   ) call calc_gz    (block, dstate)
      if (baroclinic    ) call calc_rhod  (block, dstate)
#ifdef USE_DEEP_ATM
      ! if (deepwater     ) call calc_rdp (blocks(iblk), blocks(iblk)%dstate(itime))
      call calc_rdp (block, dstate)
#endif

    end select

  end subroutine operators_prepare_2

  subroutine calc_mg(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k, is, ie, js, je

    call perf_start('calc_mg')

    associate (mesh    => block%mesh    , &
               mgs     => dstate%mgs    , & ! in
               mg_lev  => dstate%mg_lev , & ! out
               mg      => dstate%mg     )   ! out
    is = mesh%full_ids - 1
    ie = mesh%full_ide + 1
    js = mesh%full_jds - merge(0, 1, mesh%has_south_pole())
    je = mesh%full_jde + merge(0, 1, mesh%has_north_pole())
    do k = mesh%half_kds, mesh%half_kde
      do j = js, je
        do i = is, ie
          mg_lev%d(i,j,k) = vert_coord_calc_mg_lev(k, mgs%d(i,j), block%static%ref_ps_perb%d(i,j))
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = js, je
        do i = is, ie
          mg%d(i,j,k) = 0.5_r8 * (mg_lev%d(i,j,k) + mg_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_mg')

  end subroutine calc_mg

  subroutine calc_ph(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k, is, ie, js, je

    call perf_start('calc_ph')

    associate (mesh    => block%mesh          , &
               mg      => dstate%mg           , & ! in
               mg_lev  => dstate%mg_lev       , & ! in
               dmg     => dstate%dmg          , & ! in
               qm      => tracers(block%id)%qm, & ! in
               ph_lev  => dstate%ph_lev       , & ! out
               ph      => dstate%ph           )   ! out
    is = mesh%full_ids - 1
    ie = mesh%full_ide + 1
    js = mesh%full_jds - merge(0, 1, mesh%has_south_pole())
    je = mesh%full_jde + merge(0, 1, mesh%has_north_pole())
    k = mesh%half_kds
    ph_lev%d(:,:,k) = mg_lev%d(:,:,k)
    do k = mesh%half_kds + 1, mesh%half_kde
      do j = js, je
        do i = is, ie
          ph_lev%d(i,j,k) = ph_lev%d(i,j,k-1) + dmg%d(i,j,k-1) * (1 + qm%d(i,j,k-1))
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = js, je
        do i = is, ie
          ph%d(i,j,k) = 0.5_r8 * (ph_lev%d(i,j,k) + ph_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_ph')

  end subroutine calc_ph

  subroutine calc_omg(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k
    real(r8) sum_dmf(block%mesh%full_ids:block%mesh%full_ide, &
                     block%mesh%full_jds:block%mesh%full_jde)
    real(r8) work   (block%mesh%full_ids:block%mesh%full_ide, &
                     block%mesh%full_kds:block%mesh%full_kde)
    real(r8) pole   (block%mesh%full_kds:block%mesh%full_kde)

    call perf_start('calc_omg')

    associate (mesh  => block%mesh   , &
               ph    => dstate%ph    , & ! in
               u_lon => dstate%u_lon , & ! in
               v_lat => dstate%v_lat , & ! in
               dmf   => block%aux%dmf, & ! in
               div   => block%aux%div, & ! in
               omg   => block%aux%omg)   ! out
    sum_dmf = 0
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          sum_dmf(i,j) = sum_dmf(i,j) + dmf%d(i,j,k)
          omg%d(i,j,k) = 0.5_r8 * ((                                              &
            u_lon%d(i  ,j,k) * (ph%d(i,j,k) + ph%d(i+1,j,k)) -                    &
            u_lon%d(i-1,j,k) * (ph%d(i,j,k) + ph%d(i-1,j,k))                      &
          ) * mesh%le_lon(j) + (                                                  &
            v_lat%d(i,j  ,k) * (ph%d(i,j,k) + ph%d(i,j+1,k)) * mesh%le_lat(j  ) - &
            v_lat%d(i,j-1,k) * (ph%d(i,j,k) + ph%d(i,j-1,k)) * mesh%le_lat(j-1)   &
          )) / mesh%area_cell(j) - ph%d(i,j,k) * div%d(i,j,k) - sum_dmf(i,j)
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = 0.5_r8 * v_lat%d(i,j,k) * (ph%d(i,j,k) + ph%d(i,j+1,k)) * mesh%le_lat(j)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / mesh%area_pole_cap
      do k = mesh%full_kds, mesh%full_kde
        i = mesh%full_ids
        sum_dmf(i,j) = sum_dmf(i,j) + dmf%d(i,j,k)
        omg%d(:,j,k) = pole(k) - ph%d(i,j,k) * div%d(i,j,k) - sum_dmf(i,j)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = 0.5_r8 * v_lat%d(i,j-1,k) * (ph%d(i,j,k) + ph%d(i,j-1,k)) * mesh%le_lat(j-1)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = -pole / mesh%area_pole_cap
      do k = mesh%full_kds, mesh%full_kde
        i = mesh%full_ids
        sum_dmf(i,j) = sum_dmf(i,j) + dmf%d(i,j,k)
        omg%d(:,j,k) = pole(k) - ph%d(i,j,k) * div%d(i,j,k) - sum_dmf(i,j)
      end do
    end if
    end associate

    call perf_stop('calc_omg')

  end subroutine calc_omg

  subroutine calc_t(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call wait_halo(dstate%pt)

    call perf_start('calc_t')

    associate (mesh => block%mesh         , &
               pt   => dstate%pt          , & ! in
               p    => dstate%p           , & ! in
               q    => tracers(block%id)%q, & ! in
               t    => block%aux%t        , & ! out
               tv   => block%aux%tv       )   ! out
    if (idx_qv > 0) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            t%d(i,j,k) = temperature(pt%d(i,j,k), p%d(i,j,k), q%d(i,j,k,idx_qv))
            tv%d(i,j,k) = virtual_temperature_from_modified_potential_temperature(pt%d(i,j,k), p%d(i,j,k)**rd_o_cpd, q%d(i,j,k,idx_qv))
          end do
        end do
      end do
    else
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            t%d(i,j,k) = temperature(pt%d(i,j,k), p%d(i,j,k), 0.0_r8)
            tv%d(i,j,k) = t%d(i,j,k)
          end do
        end do
      end do
    end if
    end associate

    call perf_stop('calc_t')

  end subroutine calc_t

  subroutine calc_rhod(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call wait_halo(dstate%gz_lev)

    call perf_start('calc_rhod')

    associate (mesh   => block%mesh    , &
               gz_lev => dstate%gz_lev , & ! in
               dmg    => dstate%dmg    , & ! in
               rhod   => block%aux%rhod)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          rhod%d(i,j,k) = dmg%d(i,j,k) / (gz_lev%d(i,j,k) - gz_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_rhod')

  end subroutine calc_rhod

  subroutine calc_mfz(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(in) :: dtend

    integer i, j, k

    call perf_start('calc_mfz')

    associate (mesh        => block%mesh       , &
               dmf         => block%aux%dmf    , & ! in
               dmgsdt      => dtend%dmgsdt     , & ! in
               mfz_lev     => block%aux%mfz_lev)   ! out
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          mfz_lev%d(i,j,k) = -vert_coord_calc_dmgdt_lev(k, dmgsdt%d(i,j)) - sum(dmf%d(i,j,1:k-1))
        end do
      end do
    end do
    call fill_halo(mfz_lev, west_halo=.false., south_halo=.false., async=.true.)
    end associate

    call perf_stop('calc_mfz')

  end subroutine calc_mfz

  subroutine calc_ke(block, dstate, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    integer, intent(in) :: substep

    integer i, j, k
    real(r8) ke_vtx(4)
    real(r8) work(block%mesh%full_ids:block%mesh%full_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)

    call perf_start('calc_ke')

    associate (mesh => block%mesh  , &
               u    => dstate%u_lon, & ! in
               v    => dstate%v_lat, & ! in
               ke   => block%aux%ke)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          ke%d(i,j,k) = 0.25_r8 * (             &
            u%d(i-1,j  ,k)**2 + u%d(i,j,k)**2 + &
            v%d(i  ,j-1,k)**2 + v%d(i,j,k)**2   &
          )
        end do
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v%d(i,j,k)**2
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          ke%d(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          work(i,k) = v%d(i,j-1,k)**2
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole / global_mesh%full_nlon
      do k = mesh%full_kds, mesh%full_kde
        do i = mesh%full_ids, mesh%full_ide
          ke%d(i,j,k) = pole(k)
        end do
      end do
    end if

    if (save_dyn_calc .and. substep < total_substeps) then
      call perf_stop('calc_ke')
      return
    end if

    if (ke_scheme == 2) then
      !
      !      ________u_________________u________
      !     |     i-1,j+1     |       i,j+1     |
      !     |                 |                 |
      !     |                 |                 |
      !     |        1        |        4        |
      !     v        o--------v--------o        v
      !  i-1,j    i-1,j      i,j      i,j    i+1,j
      !     |        |        |        |        |
      !     |        |        |        |        |
      !     |________u________|________u________|
      !     |     i-1,j      i,j      i,j       |
      !     |        |        |        |        |
      !     |        |        |        |        |
      !     |        |        |        |        |
      !     v        o--------v--------o        v
      !  i-1,j-1  i-1,j-1    i,j-1    i,j-1  i+1,j-1
      !     |        2        |        3        |
      !     |                 |                 |
      !     |________u________|________u________|
      !           i-1,j-1             i,j-1
      !
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            ke_vtx(1) = 0.25_r8 * (                   &
              v%d(i-1,j  ,k)**2 + v%d(i  ,j  ,k)**2 + &
              u%d(i-1,j  ,k)**2 + u%d(i-1,j+1,k)**2   &
            )
            ke_vtx(2) = 0.25_r8 * (                   &
              v%d(i-1,j-1,k)**2 + v%d(i  ,j-1,k)**2 + &
              u%d(i-1,j-1,k)**2 + u%d(i-1,j  ,k)**2   &
            )
            ke_vtx(3) = 0.25_r8 * (                   &
              v%d(i  ,j-1,k)**2 + v%d(i+1,j-1,k)**2 + &
              u%d(i  ,j-1,k)**2 + u%d(i  ,j  ,k)**2   &
            )
            ke_vtx(4) = 0.25_r8 * (                   &
              v%d(i  ,j  ,k)**2 + v%d(i+1,j  ,k)**2 + &
              u%d(i  ,j  ,k)**2 + u%d(i  ,j+1,k)**2   &
            )
            ke%d(i,j,k) = (1.0_r8 - ke_cell_wgt) * 0.25_r8 * ( &
              ke_vtx(1) + ke_vtx(4) + ke_vtx(2) + ke_vtx(3)    &
            ) + ke_cell_wgt * ke%d(i,j,k)
          end do
        end do
      end do
    end if
    end associate

    call perf_stop('calc_ke')

  end subroutine calc_ke

  subroutine calc_div(block, dstate, use_filter)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    logical, intent(in) :: use_filter

    integer i, j, k

    call wait_halo(dstate%u_lon)
    call wait_halo(dstate%v_lat)

    call perf_start('calc_div')

    associate (mesh  => block%mesh        , &
               u_lon => dstate%u_lon      , & ! in
               v_lat => dstate%v_lat      , & ! in
               div   => block%aux%div     , & ! out
               divx  => block%aux%g_3d_lon, & ! working array
               divy  => block%aux%g_3d_lat, & ! working array
               div2  => block%aux%div2    )   ! out
#ifdef USE_DEEP_ATM
    if (deepwater .and. use_mesh_change) then
      call wait_halo(block%aux%rdp_lon)
      call wait_halo(block%aux%rdp_lat)
  
      call div_operator(u_lon,block%aux%rdp_lon, v_lat,block%aux%rdp_lat, div, with_halo=.true.)
    else
      call div_operator(u_lon, v_lat, div, with_halo=.true.) 
    end if
    ! if (deepwater) call deep_hamiton_modify_3d(block,div,with_halo=.true.)
#else
    call div_operator(u_lon, v_lat, div, with_halo=.true.)
#endif
    if (div_damp_order == 4) then
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide + 1
            divx%d(i,j,k) = (div%d(i+1,j,k) - div%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = (div%d(i,j+1,k) - div%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      call div_operator(divx, divy, div2)
      if (use_filter) then
        call filter_run(block%big_filter, div2)
      end if
      call fill_halo(div2, west_halo=.false., south_halo=.false.)
    else
      if (use_filter) then
        call filter_run(block%small_filter, div)
        call fill_halo(div, west_halo=.false., south_halo=.false.)
      end if
    end if
    end associate

    call perf_stop('calc_div')

  end subroutine calc_div

  subroutine calc_gz(block, dstate)

    type(block_type), intent(in) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k

    call perf_start('calc_gz')

    associate (mesh   => block%mesh      , &
               tv     => block%aux%tv    , & ! in
               ph_lev => dstate%ph_lev   , & ! in
               gz_lev => dstate%gz_lev   , & ! out
               gz     => dstate%gz       )   ! out
    do k = mesh%half_kde - 1, mesh%half_kds, -1
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          gz_lev%d(i,j,k) = gz_lev%d(i,j,k+1) + rd * tv%d(i,j,k) * log(ph_lev%d(i,j,k+1) / ph_lev%d(i,j,k))
        end do
      end do
    end do
    ! For output
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
        do i = mesh%full_ids, mesh%full_ide + 1
          gz%d(i,j,k) = 0.5_r8 * (gz_lev%d(i,j,k) + gz_lev%d(i,j,k+1))
        end do
      end do
    end do
    end associate

    call perf_stop('calc_gz')

  end subroutine calc_gz

  subroutine calc_dmg(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    integer i, j, k, l, is, ie, js, je

    call perf_start('calc_dmg')

    associate (mesh    => block%mesh       , &
               mg      => dstate%mg        , & ! in
               mg_lev  => dstate%mg_lev    , & ! in
               gz      => dstate%gz        , & ! in
               gzs     => block%static%gzs , & ! in
               dmg     => dstate%dmg       , & ! out
               dmg_lon => block%aux%dmg_lon, & ! out
               dmg_lat => block%aux%dmg_lat, & ! out
               dmg_lev => dstate%dmg_lev   , & ! out
               dmg_vtx => block%aux%dmg_vtx)   ! out
    is = mesh%full_ids - 1
    ie = mesh%full_ide + 1
    js = mesh%full_jds - merge(0, 1, mesh%has_south_pole())
    je = mesh%full_jde + merge(0, 1, mesh%has_north_pole())
    if (baroclinic .or. advection) then
      do k = mesh%full_kds, mesh%full_kde
        do j = js, je
          do i = is, ie
            dmg%d(i,j,k) = mg_lev%d(i,j,k+1) - mg_lev%d(i,j,k)
            if (dmg%d(i,j,k) <= 0) then
              do l = mesh%half_kds, mesh%half_kde
                print *, l, mg_lev%d(i,j,l)
              end do
              print *, 'mgs(i,j) =', dstate%mgs%d(i,j)
              print *, mesh%full_lon_deg(i), '(', to_str(i), ')', mesh%full_lat_deg(j), '(', to_str(j), ')', k
              call log_warning('The dry-air weight levels are not monotonic!', __FILE__, __LINE__)
              call process_stop(1)
            end if
          end do
        end do
      end do
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = js, je
          do i = is, ie
            dmg_lev%d(i,j,k) = mg%d(i,j,k) - mg%d(i,j,k-1)
          end do
        end do
      end do
      ! Top boundary
      k = mesh%half_kds
      do j = js, je
        do i = is, ie
          dmg_lev%d(i,j,k) = mg%d(i,j,k) - mg_lev%d(i,j,k)
        end do
      end do
      ! Bottom boundary
      k = mesh%half_kde
      do j = js, je
        do i = is, ie
          dmg_lev%d(i,j,k) = mg_lev%d(i,j,k) - mg%d(i,j,k-1)
        end do
      end do
      call fill_halo(dmg_lev)
    else
      do j = js, je
        do i = is, ie
          if (deepwater .and. use_mesh_change .and. use_variable_gravity) then
            dmg%d(i,j,1) = height_from_geopotential(gz%d(i,j,1)) - height_from_geopotential(gzs%d(i,j))
          else
            dmg%d(i,j,1) = (gz%d(i,j,1) - gzs%d(i,j)) / g
          end if
        end do
      end do
    end if

    call average_run(dmg, dmg_lon); call fill_halo(dmg_lon, async=.true.)
    call average_run(dmg, dmg_lat); call fill_halo(dmg_lat, async=.true.)
    call interp_run (dmg, dmg_vtx)
    end associate

    call perf_stop('calc_dmg')

  end subroutine calc_dmg

  subroutine calc_mf(block, dstate, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt

    integer i, j, k
    call wait_halo(dstate%u_lon)
    call wait_halo(dstate%v_lat)
    call wait_halo(block%aux%dmg_lon)
    call wait_halo(block%aux%dmg_lat)

    call perf_start('calc_mf')

    associate (mesh    => block%mesh       , &
               dmg     => dstate%dmg       , & ! in
               dmg_lon => block%aux%dmg_lon, & ! in
               dmg_lat => block%aux%dmg_lat, & ! in
               u_lon   => dstate%u_lon     , & ! in
               v_lat   => dstate%v_lat     , & ! in
               u_lat   => block%aux%u_lat  , & ! out
               v_lon   => block%aux%v_lon  , & ! out
               mfx_lon => block%aux%mfx_lon, & ! out
               mfy_lat => block%aux%mfy_lat)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
        do i = mesh%half_ids - 1, mesh%half_ide
          mfx_lon%d(i,j,k) = dmg_lon%d(i,j,k) * u_lon%d(i,j,k)
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide + 1
          mfy_lat%d(i,j,k) = dmg_lat%d(i,j,k) * v_lat%d(i,j,k)
        end do
      end do
    end do

    call interp_run(u_lon, u_lat)
    call interp_run(v_lat, v_lon)
    end associate

    call perf_stop('calc_mf')

  end subroutine calc_mf

  subroutine calc_vor(block, dstate, with_halo)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    logical, intent(in), optional :: with_halo

    integer i, j, k
    real(r8) work(block%mesh%half_ids:block%mesh%half_ide,block%mesh%full_nlev)
    real(r8) pole(block%mesh%full_nlev)

    call wait_halo(dstate%u_lon)
    call wait_halo(dstate%v_lat)

    call perf_start('calc_vor')

    associate (mesh  => block%mesh     , &
               u_lon => dstate%u_lon   , & ! in
               v_lat => dstate%v_lat   , & ! in
               u_lat => block%aux%u_lat, & ! in
               vor   => block%aux%vor  )   ! out
#ifdef USE_DEEP_ATM
    ! if (deepwater) call deep_hamiton_modify_3d(block,vor,with_halo)
    if (deepwater .and. use_mesh_change) then
      call wait_halo(block%aux%rdp_lon)
      call wait_halo(block%aux%rdp_lat)
      call curl_operator_deep(u_lon,block%aux%rdp_lon, v_lat,block%aux%rdp_lat, vor, with_halo)
    else
      call curl_operator(u_lon, v_lat, vor, with_halo)
    end if 
#else
    call curl_operator(u_lon, v_lat, vor, with_halo)
#endif
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_jds
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            work(i,k) = u_lat%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = -pole * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids - 1, mesh%half_ide
            vor%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_jde
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids, mesh%half_ide
            work(i,k) = u_lat%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = mesh%full_kds, mesh%full_kde
          do i = mesh%half_ids - 1, mesh%half_ide
            vor%d(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
    end associate

    call perf_stop('calc_vor')

  end subroutine calc_vor

  subroutine calc_pv(block, dstate, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    integer, intent(in) :: substep

    integer i, j, k

    call perf_start('calc_pv')

    associate (mesh    => block%mesh       , &
               dmg_vtx => block%aux%dmg_vtx, & ! in
               vor     => block%aux%vor    , & ! in
               pv      => block%aux%pv     )   ! out
    call calc_vor(block, dstate)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide
          pv%d(i,j,k) = (vor%d(i,j,k) + mesh%f_lat(j)) / dmg_vtx%d(i,j,k)
        end do
      end do
    end do
    call fill_halo(pv)
    end associate

    call perf_stop('calc_pv')

  end subroutine calc_pv

  subroutine interp_pv_midpoint(block, dstate, dt, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    integer i, j, k

    call perf_start('interp_pv_midpoint')

    associate (mesh   => block%mesh      , &
               pv     => block%aux%pv    , & ! in
               pv_lon => block%aux%pv_lon, & ! out
               pv_lat => block%aux%pv_lat)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          pv_lat%d(i,j,k) = 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          pv_lon%d(i,j,k) = 0.5_r8 * (pv%d(i,j,k) + pv%d(i,j-1,k))
        end do
      end do
    end do
    call fill_halo(pv_lon, east_halo=.false., south_halo=.false., async=.true.)
    call fill_halo(pv_lat, west_halo=.false., north_halo=.false., async=.true.)
    end associate

    call perf_stop('interp_pv_midpoint')

  end subroutine interp_pv_midpoint

  subroutine interp_pv_upwind(block, dstate, dt, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    real(r8) b
    integer i, j, k

    call perf_start('interp_pv_upwind')

    associate (mesh   => block%mesh      , &
               un     => dstate%u_lon    , & ! in
               vn     => dstate%v_lat    , & ! in
               ut     => block%aux%u_lat , & ! in
               vt     => block%aux%v_lon , & ! in
               pv     => block%aux%pv    , & ! in
               pv_lon => block%aux%pv_lon, & ! out
               pv_lat => block%aux%pv_lat)   ! out
    select case (upwind_order_pv)
    case (1)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind1(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-1:j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind1(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-1:i,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    case (3)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind3(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-2:j+1,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b  = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind3(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-2:i+1,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    case (5)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind5(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-3:j+2,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind5(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-3:i+2,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    end select
    if (substep == total_substeps .or. .not. save_dyn_calc) then
      call fill_halo(pv_lon, east_halo=.false., south_halo=.false., async=.true.)
      call fill_halo(pv_lat, west_halo=.false., north_halo=.false., async=.true.)
    end if
    end associate

    call perf_stop('interp_pv_upwind')

  end subroutine interp_pv_upwind

  subroutine interp_pv_weno(block, dstate, dt, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    real(r8) b
    integer i, j, k

    call perf_start('interp_pv_weno')

    associate (mesh   => block%mesh      , &
               un     => dstate%u_lon    , & ! in
               vn     => dstate%v_lat    , & ! in
               ut     => block%aux%u_lat , & ! in
               vt     => block%aux%v_lon , & ! in
               pv     => block%aux%pv    , & ! in
               pv_lon => block%aux%pv_lon, & ! out
               pv_lat => block%aux%pv_lat)   ! out
    select case (weno_order_pv)
    case (1)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * upwind1(sign(1.0_r8, vt%d(i,j,k)), upwind_wgt_pv, pv%d(i,j-1:j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * upwind1(sign(1.0_r8, ut%d(i,j,k)), upwind_wgt_pv, pv%d(i-1:i,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    case (3)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * weno3(sign(1.0_r8, vt%d(i,j,k)), pv%d(i,j-2:j+1,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * weno3(sign(1.0_r8, ut%d(i,j,k)), pv%d(i-2:i+1,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    case (5)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids, mesh%half_ide
            b = abs(vt%d(i,j,k)) / (sqrt(un%d(i,j,k)**2 + vt%d(i,j,k)**2) + eps)
            pv_lon%d(i,j,k) = b * weno5(sign(1.0_r8, vt%d(i,j,k)), pv%d(i,j-3:j+2,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i,j-1,k) + pv%d(i,j,k))
          end do
        end do
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            b = abs(ut%d(i,j,k)) / (sqrt(ut%d(i,j,k)**2 + vn%d(i,j,k)**2) + eps)
            pv_lat%d(i,j,k) = b * weno5(sign(1.0_r8, ut%d(i,j,k)), pv%d(i-3:i+2,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv%d(i-1,j,k) + pv%d(i,j,k))
          end do
        end do
      end do
    end select
    if (substep == total_substeps .or. .not. save_dyn_calc) then
      call fill_halo(pv_lon, east_halo=.false., south_halo=.false., async=.true.)
      call fill_halo(pv_lat, west_halo=.false., north_halo=.false., async=.true.)
    end if
    end associate

    call perf_stop('interp_pv_weno')

  end subroutine interp_pv_weno

  subroutine calc_coriolis(block, dstate, dtend, dt, substep)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt
    integer, intent(in) :: substep

    real(r8) tmp
    integer i, j, k

    call wait_halo(block%aux%pv_lon)
    call wait_halo(block%aux%pv_lat)

    call perf_start('calc_coriolis')

    associate (mesh    => block%mesh       , &
               mfx_lon => block%aux%mfx_lon, & ! in
               mfy_lat => block%aux%mfy_lat, & ! in
               pv_lon  => block%aux%pv_lon , & ! in
               pv_lat  => block%aux%pv_lat , & ! in
               dudt    => dtend%dudt       , & ! out
               dvdt    => dtend%dvdt       )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp = - (                                                            &
            mesh%tg_wgt_lat(1,j) * (                                           &
              mfx_lon%d(i-1,j  ,k) * (pv_lat%d(i,j,k) + pv_lon%d(i-1,j  ,k)) + &
              mfx_lon%d(i  ,j  ,k) * (pv_lat%d(i,j,k) + pv_lon%d(i  ,j  ,k))   &
            ) +                                                                &
            mesh%tg_wgt_lat(2,j) * (                                           &
              mfx_lon%d(i-1,j+1,k) * (pv_lat%d(i,j,k) + pv_lon%d(i-1,j+1,k)) + &
              mfx_lon%d(i  ,j+1,k) * (pv_lat%d(i,j,k) + pv_lon%d(i  ,j+1,k))   &
            )                                                                  &
          ) * 0.5_r8
          dvdt%d(i,j,k) = dvdt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_coriolis%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          tmp = (                                                              &
            mesh%tg_wgt_lon(1,j) * (                                           &
              mfy_lat%d(i  ,j-1,k) * (pv_lon%d(i,j,k) + pv_lat%d(i  ,j-1,k)) + &
              mfy_lat%d(i+1,j-1,k) * (pv_lon%d(i,j,k) + pv_lat%d(i+1,j-1,k))   &
            ) +                                                                &
            mesh%tg_wgt_lon(2,j) * (                                           &
              mfy_lat%d(i  ,j  ,k) * (pv_lon%d(i,j,k) + pv_lat%d(i  ,j  ,k)) + &
              mfy_lat%d(i+1,j  ,k) * (pv_lon%d(i,j,k) + pv_lat%d(i+1,j  ,k))   &
            )                                                                  &
          ) * 0.5_r8
          dudt%d(i,j,k) = dudt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_coriolis%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

    call perf_stop('calc_coriolis')

  end subroutine calc_coriolis

#ifdef USE_DEEP_ATM
  subroutine calc_nct_coriolis(block, dstate, dtend, dt)  !!add by cui ,nct means nontraditional coriolis term
    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend

    type(latlon_field3d_type) :: w_lon, gz_lon, w_lat, gz_lat
    real(r8), intent(in) :: dt

    real(r8) tmp
    integer i, j, k
    call wait_halo(block%aux%rdp_lon)
    call wait_halo(block%aux%rdp_lat)

    call perf_start('calc_nct_coriolis')
    ! print*,"nct_hor"

    associate (mesh    => block%mesh       , &
               u       => dstate%u_lon     , &
               v       => dstate%v_lat     , &
               w       => dstate%w         , &
               w_lev   => dstate%w_lev     , &
              !  gz      => dstate1%gz        , &
               w_lon   => block%aux%w_lon  , &
               rdp_lon => block%aux%rdp_lon, &
               w_lat   => block%aux%w_lat  , &
               rdp_lat => block%aux%rdp_lat, &
               dudt    => dtend%dudt       , & ! out
               dvdt    => dtend%dvdt      )   ! out

    call interp_run(w_lev,w);call fill_halo(w, async=.true.)

    call interp_run(w,w_lon);call fill_halo(w_lon, async=.true.)
    call interp_run(w,w_lat);call fill_halo(w_lat, async=.true.)

    ! call interp_run(gz,gz_lon)
    ! call interp_run(gz,gz_lat)
    
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          ! tmp = (mesh%fd_lon(j) + u%d(i,j,k)/(gz_lon%d(i,j,k)/g + radius)) * w_lon%d(i,j,k)
          ! tmp = (mesh%fd_lon(j) ) * w_lon%d(i,j,k)
          ! tmp = u%d(i,j,k)/rdp_lon%d(i,j,k)* w_lon%d(i,j,k)
          tmp = (mesh%fd_lon(j) + u%d(i,j,k)/rdp_lon%d(i,j,k)) * w_lon%d(i,j,k)
          ! tmp = 0 !!!! for test the contribution
          ! tmp = (mesh%fd_lon(j) + u%d(i,j,k)/radius) * w_lon%d(i,j,k)
          ! print*,mesh%fd_lon(j),"1"
          ! print*,u%d(i,j,k),"2"

          ! print*,gz_lon%d(i,j,k),"3"
          ! print*,w_lon%d(i,j,k),"4"
          ! print*,radius
          ! if (mesh%full_lat_deg(j) .lt. 22.0 .and.mesh%full_lat_deg(j) .gt. 18.0 .and. &
          !     mesh%full_lon_deg(i) .lt. 22.0.and.mesh%full_lon_deg(i) .gt. 18.0) then
          !   ! print*,"2omegacos w ",mesh%fd_lon(j) * w_lon%d(i,j,k) 
          !   ! print*,"dudt", dudt%d(i,j,k)
          !   print*,"ratio", abs(tmp)/abs(dudt%d(i,j,k))
          !   print*,"mesh%full_lat_deg(j)",mesh%full_lat_deg(j),"k",k
          ! end if

          dudt%d(i,j,k) = dudt%d(i,j,k) - tmp
          ! print*,"tmp=",tmp,"mesh%fd_lon(j)",mesh%fd_lon(j)
          ! print*,"fw",mesh%fd_lon(j)* w_lon%d(i,j,k)

          ! dv%d(i,j,k) = dv%d(i,j,k) + (v%d(i,j,k)*g/gz_lon%d(i,j,k)) * w_lon%d(i,j,k)
! #ifdef USE_DEEP_ATM
          ! dtend%dudt_nct_coriolis%d(i,j,k) = -tmp
! #endif
        end do
      end do
    end do

    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          ! tmp = (mesh%fd_lat(j) + v%d(i,j,k)/(gz_lat%d(i,j,k)/g + radius)) * w_lat%d(i,j,k)
          tmp =  w_lat%d(i,j,k) * v%d(i,j,k)/rdp_lat%d(i,j,k)
          ! tmp = 0 !!!! for test the contribution
          ! tmp = w_lat%d(i,j,k) * v%d(i,j,k)/radius 

          dvdt%d(i,j,k) = dvdt%d(i,j,k) - tmp
        end do
      end do
    end do

    end associate
    call perf_stop('calc_nct_coriolis')
  end subroutine calc_nct_coriolis

  ! pure function integral_1rho_dmg(eta) result(res)
    

  subroutine calc_rdp(block, dstate)
    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    ! type(dtend_type), intent(inout) :: dtend
    ! real(r8), intent(in) :: dt

    !real integral_1rho_dmg
    integer i, j, k
    ! integer neval, ierr
    real(8) tmp_lev, tmp, r_h
    ! print*,"deepwater"
    call perf_start("calc_rdp")

    ! real(8) mgrho, abserr, integral_1rho_dmg

    associate(mesh    => block%mesh        , &
              rdp     => block%aux%rdp     , & ! out
              rdp_lev => block%aux%rdp_lev , & ! out
              rdp_lev_lon => block%aux%rdp_lev_lon , & ! out
              rdp_lev_lat => block%aux%rdp_lev_lat , & ! out
              rdp_lat => block%aux%rdp_lat , &
              rdp_lon => block%aux%rdp_lon , &
              rdp_vtx => block%aux%rdp_vtx , &
              gz      => dstate%gz         , &
              gz_lev  => dstate%gz_lev     , &
              ! rdp_lev => block%aux%rdp_lev , &              
              dmg     => dstate%dmg        , &
              dmg_lev     => dstate%dmg_lev        , &
              gzs     => block%static%gzs , & ! in
              ! mgs  => block%dstate%mgs , &
              rhod    => block%aux%rhod    )
    if (deepwater .and. use_mesh_change .and. use_variable_gravity) then
      do k = mesh%half_kds, mesh%half_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            rdp_lev%d(i,j,k) = radius_from_geopotential(gz_lev%d(i,j,k))
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide + 1
            rdp%d(i,j,k) = radius_from_geopotential(gz%d(i,j,k))
          end do
        end do
      end do
    else
      ! do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp_lev = 0.0_r8
          rdp_lev%d(i,j,mesh%half_kde) = radius + gzs%d(i,j) / g
          r_h = rdp_lev%d(i,j,mesh%half_kde)

          tmp = dstate%mgs%d(i,j) - dstate%mg%d(i,j,mesh%full_kde)
          rdp%d(i,j,mesh%full_kde) = (r_h**3 + 3*radius**2 * tmp / g)**(1.0_r8/3.0_r8)

          do k = mesh%half_kde - 1, mesh%half_kds, -1
            tmp_lev = tmp_lev + dmg%d(i,j,k) / rhod%d(i,j,k)
            rdp_lev%d(i,j,k) = (r_h**3 + 3*radius**2 * tmp_lev / g)**(1.0_r8/3.0_r8)
          end do
          do k = mesh%full_kde - 1, mesh%full_kds, -1
            tmp = tmp + dmg_lev%d(i,j,k+1) * (1 / rhod%d(i,j,k+1) + 1 / rhod%d(i,j,k)) / 2
            rdp%d(i,j,k) = (r_h**3 + 3*radius**2 * tmp / g)**(1.0_r8/3.0_r8)
          end do
        end do
      end do
    end if
    call fill_halo(rdp)
    call fill_halo(rdp_lev)
    call average_run(rdp, rdp_lon); call fill_halo(rdp_lon, async=.true.)
    call average_run(rdp, rdp_lat); call fill_halo(rdp_lat, async=.true.)
    call interp_run(rdp, rdp_vtx); call fill_halo(rdp_vtx, async=.true.)
    ! call interp_run(rdp, rdp_lev); call fill_halo(rdp_lev, async=.true.)
    call interp_run(rdp_lev, rdp_lev_lon); call fill_halo(rdp_lev_lon, async=.true.)
    call interp_run(rdp_lev, rdp_lev_lat); call fill_halo(rdp_lev_lat, async=.true.)
    ! if (proc%id_model == 1) then
    !   do k = mesh%full_kds, mesh%full_kde
    !     do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
    !       ! do i = mesh%half_ids - 1, mesh%half_ide
    !         print*,"before",rdp_lon%d(mesh%half_ids - 1,j,k),j,k
    !       ! end do
    !     end do
    !   end do
    !   end if

    ! do k = mesh%full_kds, mesh%full_kde
    !   do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
    !     do i = mesh%full_ids, mesh%full_ide
    !       integral_1rho_dmg = - 1/rhod(i,j,k)
    !       call qags(integral_1rho_dmg, mgs(i,j,k), mg(i,j,k),1.0d-12, 1.0d-3, mgrho,  abserr, neval, ierr)
    !       rdp = 3.0 * radius**2/g * mgrho 
    !       rdp = rdp
    
    end associate
    call perf_stop("calc_rdp")
  end subroutine calc_rdp

  subroutine deep_hamiton_modify_3d(block,fd,with_halo)
    type(block_type), intent(inout) :: block
    type(latlon_field3d_type), intent(inout)  :: fd
    logical, intent(in), optional :: with_halo
    ! type(latlon_field3d_type), intent(out) :: fd
    !
    integer i, j, k, is, ie, js, je, ks, ke
    real(r8) factor_r
    logical with_halo_opt
    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo
    ! ks = merge(mesh%full_kds, mesh%half_kds, fd%loc(1:3) /= 'lev')
    ! ke = merge(mesh%full_kde, mesh%half_kde, fd%loc(1:3) /= 'lev')
    select case (fd%loc)
    case("vtx")
      associate (mesh => fd%mesh ,&
                  rdp_vtx => block%aux%rdp_vtx)
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%half_ids, mesh%half_ide 
          do k = mesh%full_kds, mesh%full_kde
            factor_r = radius/rdp_vtx%d(i,j,k)
            fd%d(i,j,k) = fd%d(i,j,k)/factor_r
          end do
        end do
      end do
      end associate
    case("cell")
      associate (mesh => fd%mesh ,&
                  rdp => block%aux%rdp)
      is = mesh%full_ids
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = mesh%full_kds, mesh%full_kde
        do j = js, je
          do i = is, ie
            factor_r = radius/rdp%d(i,j,k)
            fd%d(i,j,k) = fd%d(i,j,k)/factor_r
          end do
        end do
      end do
      end associate
    case("lev")
      associate (mesh => fd%mesh ,&
                  rdp_lev => block%aux%rdp_lev)
      is = mesh%full_ids
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = mesh%half_kds, mesh%half_kde
        do j = js, je
          do i = is, ie
            factor_r = radius/rdp_lev%d(i,j,k)
            fd%d(i,j,k) = fd%d(i,j,k)/factor_r
          end do
        end do
      end do
      end associate
    case("lon")
      associate (mesh => fd%mesh ,&
        rdp_lon => block%aux%rdp_lon)
        is = mesh%half_ids
        ie = mesh%half_ide
        js = mesh%full_jds_no_pole
        je = mesh%full_jde_no_pole
        do k = mesh%full_kds, mesh%full_kde
          do j = js, je
            do i = is, ie
              factor_r = radius/rdp_lon%d(i,j,k)
              fd%d(i,j,k) = fd%d(i,j,k)/factor_r
            end do
          end do
        end do
      end associate
    case("lat")
      associate (mesh => fd%mesh ,&
        rdp_lat => block%aux%rdp_lat)
        is = mesh%full_ids
        ie = mesh%full_ide
        js = mesh%half_jds
        je = mesh%half_jde
        do k = mesh%full_kds, mesh%full_kde
          do j = js, je
            do i = is, ie
              factor_r = radius/rdp_lat%d(i,j,k)
              fd%d(i,j,k) = fd%d(i,j,k)/factor_r
            end do
          end do
        end do
      end associate
    end select
  end subroutine deep_hamiton_modify_3d
  !!!
#endif

  subroutine calc_grad_ke(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    real(r8) tmp
    integer i, j, k
#ifdef USE_DEEP_ATM
    real(r8) factor_r
    call wait_halo(block%aux%rdp_lon)
    call wait_halo(block%aux%rdp_lat)
#endif

    call perf_start('calc_grad_ke')

    associate (mesh => block%mesh  , &
               ke   => block%aux%ke, & ! in
#ifdef USE_DEEP_ATM
               rdp_lat => block%aux%rdp_lat, & ! in
               rdp_lon => block%aux%rdp_lon, & ! in
#endif
               dudt => dtend%dudt  , & ! out
               dvdt => dtend%dvdt  )   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
#ifdef USE_DEEP_ATM
          factor_r = merge(radius / rdp_lon%d(i,j,k), 1.0_r8, deepwater .and. use_mesh_change) 
          tmp = -(ke%d(i+1,j,k) - ke%d(i,j,k)) / mesh%de_lon(j) * factor_r
#else
          tmp = -(ke%d(i+1,j,k) - ke%d(i,j,k)) / mesh%de_lon(j)
#endif
          dudt%d(i,j,k) = dudt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_dkedx%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
#ifdef USE_DEEP_ATM
          factor_r = merge(radius / rdp_lat%d(i,j,k), 1.0_r8, deepwater .and. use_mesh_change) 
          tmp = -(ke%d(i,j+1,k) - ke%d(i,j,k)) / mesh%de_lat(j) * factor_r
#else          
          tmp = -(ke%d(i,j+1,k) - ke%d(i,j,k)) / mesh%de_lat(j)
#endif
          dvdt%d(i,j,k) = dvdt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_dkedy%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

    call perf_stop('calc_grad_ke')

  end subroutine calc_grad_ke

  subroutine calc_grad_mf(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate

    call perf_start('calc_grad_mf')

    associate (mesh    => block%mesh       , &
               mfx_lon => block%aux%mfx_lon, & ! in
               mfy_lat => block%aux%mfy_lat, & ! in
               dmf     => block%aux%dmf    )   ! out
#ifdef USE_DEEP_ATM
    if (deepwater .and. use_mesh_change) then
      call wait_halo(block%aux%rdp_lon)
      call wait_halo(block%aux%rdp_lat)
  
      call div_operator(mfx_lon,block%aux%rdp_lon, mfy_lat,block%aux%rdp_lat, dmf)
    else 
      call div_operator(mfx_lon, mfy_lat, dmf)
    end if 
    ! if (deepwater) call deep_hamiton_modify_3d(block,dmf)
#else
    call div_operator(mfx_lon, mfy_lat, dmf)
#endif
    end associate

    call perf_stop('calc_grad_mf')

  end subroutine calc_grad_mf

  subroutine calc_grad_ptf(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate
    type(dtend_type ), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    integer i, j, k

    call wait_halo(dstate%pt)

    call perf_start('calc_grad_ptf')

    associate (mesh    => block%filter_mesh      , &
               u_lon   => dstate%u_lon           , & ! in
               v_lat   => dstate%v_lat           , & ! in
               mfx_lon => block%aux%mfx_lon      , & ! in
               mfy_lat => block%aux%mfy_lat      , & ! in
               mfz_lev => block%aux%mfz_lev      , & ! in
#ifdef USE_DEEP_ATM
               rdp_lev => block%aux%rdp_lev      , & ! in
#endif
               dmg     => dstate%dmg             , & ! in
               pt      => dstate%pt              , & ! in
               ptfx    => block%adv_batch_pt%qmfx, & ! out
               ptfy    => block%adv_batch_pt%qmfy, & ! out
               ptfz    => block%adv_batch_pt%qmfz, & ! out
               dptdt   => dtend%dptdt            )   ! out
    call block%adv_batch_pt%set_wind(u=u_lon, v=v_lat, mfx=mfx_lon, mfy=mfy_lat, mfz=mfz_lev, m=dmg, dt=dt &
#ifdef USE_DEEP_ATM
      , rdp=block%aux%rdp, rdp_x=block%aux%rdp_lon, rdp_y=block%aux%rdp_lat, rdp_z=block%aux%rdp_lev &
#endif
      )
    call swift_prepare(block%adv_batch_pt, dt)
    call adv_calc_tracer_hflx(block%adv_batch_pt, pt, ptfx, ptfy, dt)
#ifdef USE_DEEP_ATM
    if (deepwater .and. use_mesh_change) then
      call wait_halo(block%aux%rdp_lon)
      call wait_halo(block%aux%rdp_lat)
  
      call div_operator(ptfx,block%aux%rdp_lon, ptfy,block%aux%rdp_lat, dptdt)
    else
      call div_operator(ptfx, ptfy, dptdt)
    end if
#else
    call div_operator(ptfx, ptfy, dptdt)
#endif
    call adv_fill_vhalo(pt, no_negvals=.true.)
    call adv_calc_tracer_vflx(block%adv_batch_pt, pt, ptfz, dt)
#ifdef USE_DEEP_ATM
    if (deepwater .and. use_mesh_change) call wait_halo(rdp_lev)
#endif
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
#ifdef USE_DEEP_ATM
          if (deepwater .and. use_mesh_change) then
            ! Match the deep vertical shell-area form already used by mass and tracer updates.
            dptdt%d(i,j,k) = -dptdt%d(i,j,k) - ( &
              ptfz%d(i,j,k+1) * (rdp_lev%d(i,j,k+1) / radius)**2 -        &
              ptfz%d(i,j,k  ) * (rdp_lev%d(i,j,k  ) / radius)**2 ) /       &
              (rdp_lev%d(i,j,k) / radius)**2
          else
            dptdt%d(i,j,k) = -dptdt%d(i,j,k) - (ptfz%d(i,j,k+1) - ptfz%d(i,j,k))
          end if
#else
          dptdt%d(i,j,k) = -dptdt%d(i,j,k) - (ptfz%d(i,j,k+1) - ptfz%d(i,j,k))
#endif
        end do
      end do
    end do
    end associate

    call perf_stop('calc_grad_ptf')

  end subroutine calc_grad_ptf

  subroutine calc_dmgsdt(block, dstate, dtend)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend

    integer i, j, k

    call perf_start('calc_dmgsdt')

    associate (mesh   => block%mesh   , &
               dmf    => block%aux%dmf, & ! in
               dmgsdt => dtend%dmgsdt )   ! out
    dmgsdt%d = 0
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dmgsdt%d(i,j) = dmgsdt%d(i,j) - dmf%d(i,j,k)
        end do
      end do
    end do
    end associate

    call perf_stop('calc_dmgsdt')

  end subroutine calc_dmgsdt

  subroutine calc_wedudlev_wedvdlev(block, dstate, dtend, dt)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(in) :: dstate
    type(dtend_type), intent(inout) :: dtend
    real(r8), intent(in) :: dt

    real(r8) tmp
    integer i, j, k

    ! Follow SB81 vertical advection discretization.

    call wait_halo(block%aux%mfz_lev)

    call perf_start('calc_wedudlev_wedvdlev')

    associate (mesh        => block%mesh           , &
               u           => dstate%u_lon         , & ! in
               v           => dstate%v_lat         , & ! in
               dmg_lon     => block%aux%dmg_lon    , & ! in
               dmg_lat     => block%aux%dmg_lat    , & ! in
               mfz_lev     => block%aux%mfz_lev    , & ! in
               mfz_lev_lon => block%aux%mfz_lev_lon, & ! out
               mfz_lev_lat => block%aux%mfz_lev_lat, & ! out
               dudt        => dtend%dudt           , & ! out
               dvdt        => dtend%dvdt           )   ! out
    call interp_run(mfz_lev, mfz_lev_lon)
    call interp_run(mfz_lev, mfz_lev_lat)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          tmp = -(                                                 &
            mfz_lev_lon%d(i,j,k+1) * (u%d(i,j,k+1) - u%d(i,j,k)) - &
            mfz_lev_lon%d(i,j,k  ) * (u%d(i,j,k-1) - u%d(i,j,k))   &
          ) / dmg_lon%d(i,j,k) / 2.0_r8
          dudt%d(i,j,k) = dudt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dudt_wedudeta%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          tmp = -(                                                 &
            mfz_lev_lat%d(i,j,k+1) * (v%d(i,j,k+1) - v%d(i,j,k)) - &
            mfz_lev_lat%d(i,j,k  ) * (v%d(i,j,k-1) - v%d(i,j,k))   &
          ) / dmg_lat%d(i,j,k) / 2.0_r8
          dvdt%d(i,j,k) = dvdt%d(i,j,k) + tmp
#ifdef OUTPUT_H1_DTEND
          dtend%dvdt_wedvdeta%d(i,j,k) = tmp
#endif
        end do
      end do
    end do
    end associate

    call perf_stop('calc_wedudlev_wedvdlev')

  end subroutine calc_wedudlev_wedvdlev

end module operators_mod
