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

    dt       = dt_phys
    dt_mp    = dt_phys / nsplit
    nlev     = mesh(1)%nlev
    nlayrad  = nlev + 1
    nlevrad  = nlev + 2
    ptrop    = ptop
    pstrat   = ptrop * 0.5_r8
    lnpstrat = log(pstrat)
    pstratk  = (pstrat / p0)**rd_o_cpd

    if (present(model_root)) data_root = trim(model_root) // '/data/mars'

    call gomars_v1_parse_namelist(namelist_path)
    call gomars_v1_tracers_init(dt_adv)
    call gomars_v1_objects_init(mesh)
    call gomars_v1_read_static_data()
    call gomars_v1_orbit_init()
    call gomars_v1_rad_init()
    call gomars_v1_lsm_init()
    call gomars_v1_pbl_init()
    call gomars_v1_mp_init()
    call gomars_v1_damp_init()

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
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      ! Set some initial variables (from init1).
      if (.not. restart) then
        ! FIXME: Here is for cold run.
        do icol = 1, mesh%ncol
          state%tstrat    (icol) = state%t_bot(icol)
          state%co2ice_sfc(icol) = 0
          state%tg        (icol) = state%t_bot(icol)
        end do
      end if
      call ini_optdst(qextv, qscatv, gv, qexti, qscati, gi, &
                      state%qxvdst, state%qsvdst, state%gvdst, &
                      state%qxidst, state%qsidst, state%gidst, &
                      state%qextrefdst)
      call ini_optcld(state%qxvcld, state%qsvcld, state%gvcld, &
                      state%qxicld, state%qsicld, state%gicld, &
                      state%qextrefcld, state%taurefcld)
      ! firstcomp3:
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
      call fiona_set_dim('alb', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('alb', 'lat', span=[-90, 90])
      call fiona_start_input('alb')
      call fiona_input_range('alb', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('alb', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('alb', 'albedo', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('alb')
      call latlon_interp_bilinear_column(lon, lat, dat, mesh%lon, mesh%lat, state%alsp)
      deallocate(lon, lat, dat)

      ! Surface thermal inertia
      call fiona_open_dataset('zin', file_path=trim(data_root)//'/nasa_zin_2011.nc', mpi_comm=proc%comm_model)
      call fiona_set_dim('zin', 'lon', span=[0, 360], cyclic=.true.)
      call fiona_set_dim('zin', 'lat', span=[-90, 90])
      call fiona_start_input('zin')
      call fiona_input_range('zin', 'lon', lon, coord_range=[mesh%min_lon, mesh%max_lon]); lon = lon * rad
      call fiona_input_range('zin', 'lat', lat, coord_range=[mesh%min_lat, mesh%max_lat]); lat = lat * rad
      call fiona_input_range('zin', 'zin', dat, coord_range_1=[mesh%min_lon, mesh%max_lon], coord_range_2=[mesh%min_lat, mesh%max_lat])
      call fiona_end_input('zin')
      allocate(ilon(size(lon)+1))
      allocate(ilat(size(lat)+1))
      dlon = lon(2) - lon(1) ! Assuming equidistant.
      do i = 1, size(lon)
        ilon(i  ) = lon(i) - 0.5_r8 * dlon
        ilon(i+1) = lon(i) + 0.5_r8 * dlon
      end do
      dlat = lat(2) - lat(1) ! Assuming equidistant.
      do j = 1, size(lat)
        ilat(j  ) = lat(j) - 0.5_r8 * dlat
        ilat(j+1) = lat(j) + 0.5_r8 * dlat
      end do
      call latlon_interp_fill_column(ilon, ilat, dat, mesh%lon, mesh%lat, state%zin(:,1))
      do k = 2, nsoil
        state%zin(:,k) = state%zin(:,1)
      end do
      call physics_pole_sum(mesh%lat, state%zin)
      deallocate(lon, lat, ilon, ilat, dat)

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
      deallocate(lon, lat, dat)
      end associate
    end if

  end subroutine gomars_v1_read_static_data

  subroutine gomars_v1_run(time)

    type(datetime_type), intent(in) :: time

    integer iblk, icol, k, l, m, substep
    real(r8) ls, nonlte, tmp
    real(r8) nfluxtopv, nfluxtopi, diffvt, albi

    ls = time%solar_longitude()
    time_of_day = time%time_of_day()
    call update_solar_decl_angle(ls)
    call update_solar(ls)

    ! NOTE: Old time step values of u, v, pt, q are already saved in state%u_old, etc.
    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      call gomars_v1_orbit_cosz     (state)
      call direct_solar_flux        (state)
      call gomars_v1_orbit_cosz_avg (state)
      call gomars_v1_lsm_run        (state)
      call interp_temperature       (state)
      call gomars_v1_pbl_run        (state)
      call interp_temperature       (state)
      call gomars_v1_cnvadj_run     (state)
      end associate
    end do

  end subroutine gomars_v1_run

  subroutine gomars_v1_final()

    call gomars_v1_objects_final()
    call gomars_v1_rad_final()
    call gomars_v1_lsm_final()
    call gomars_v1_pbl_final()
    call gomars_v1_mp_final()
    call gomars_v1_damp_final()

  end subroutine gomars_v1_final

  subroutine gomars_v1_d2p()

    integer iblk, i, k

    ! Calculate ts.

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
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
      do i = 1, mesh%ncol
        state%ps_old(i) = state%ps(i)
      end do
      end associate
    end do

  end subroutine gomars_v1_d2p

  subroutine gomars_v1_p2d()

    integer iblk, icol, k, m

    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state, tend => objects(iblk)%tend)
      do icol = 1, mesh%ncol
        tend%dpsdt(icol) = (state%ps(icol) - state%ps_old(icol)) / dt
      end do
      do k = 1, mesh%nlev
        do icol = 1, mesh%ncol
          tend%dudt(icol,k) = (state%u(icol,k) - state%u_old(icol,k)) / dt
          tend%dvdt(icol,k) = (state%v(icol,k) - state%v_old(icol,k)) / dt
          tend%dtdt(icol,k) = (state%t(icol,k) - state%t_old(icol,k)) / dt
        end do
      end do
      do m = 1, ntracers
        do k = 1, mesh%nlev
          do icol = 1, mesh%ncol
            tend%dqdt(icol,k,m) = (state%q(icol,k,m) - state%q_old(icol,k,m)) / dt
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

end module gomars_v1_driver_mod