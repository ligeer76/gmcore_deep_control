MODULE subinertial_wave_test_mod

    !=======================================================================
    !
    !  Subinertial internal-wave test for deep/shallow atmosphere comparison
    !  (Non-traditional Coriolis term effect, weak stratification waveguide)
    !
    !  Structure follows baroclinic_wave_test_mod (Ullrich et al. style):
    !    - subinertial_wave_test_init(namelist_path)
    !    - subinertial_wave_test_set_ic(block)
    !    - subinertial_wave_test_sample(...): point-wise analytic fields
    !
    !  Background:
    !    u=v=0, horizontally uniform, hydrostatic.
    !    T(p) constructed with weak stratification layer (p_weak_bot..p_weak_top)
    !    and stronger stratification above (smooth transition).
    !
    !  Perturbation:
    !    local warm/cold anomaly in T (or theta) confined in weak layer.
    !
    !=======================================================================
    
      IMPLICIT NONE
      private
    
      public :: subinertial_wave_test_init
      public :: subinertial_wave_test_set_ic
    
    !----------------------------
    ! Physical constants
    !----------------------------
      REAL(8), PARAMETER :: pi    = 3.141592653589793d0
      REAL(8), PARAMETER :: Rd    = 287.0d0
      REAL(8), PARAMETER :: cp    = 1004.5d0
      REAL(8), PARAMETER :: g     = 9.80616d0
      REAL(8), PARAMETER :: p0    = 100000.d0
      REAL(8), PARAMETER :: kappa = Rd/cp
      REAL(8), PARAMETER :: omega = 7.29212d-5
      REAL(8), PARAMETER :: a     = 6371220.d0
      REAL(8), PARAMETER :: deg2rad = pi/180.d0
    
    !----------------------------
    ! Namelist controls
    !----------------------------
      integer :: deep        = 1      ! 1 deep-atmosphere terms enabled (framework switch)
      integer :: add_pert    = 1      ! 1 add initial perturbation
      real(8) :: ptop        = 500.d0 ! Pa, model top (can be overwritten)
      ! Weak-layer bounds (Pa)
      real(8) :: p_weak_bot  = 90000.d0
      real(8) :: p_weak_top  = 70000.d0
      ! Reference temperatures
      real(8) :: T_surf      = 280.d0
      real(8) :: T_strat     = 220.d0
      ! Smooth transition controls in log-pressure coordinate
      real(8) :: p_trans     = 30000.d0   ! Pa, transition center (~300 hPa)
      real(8) :: trans_width = 0.35d0     ! nondim width in ln(p) space
    
      ! Perturbation parameters
      real(8) :: lonc_deg    = 30.d0
      real(8) :: latc_deg    = 30.d0
      real(8) :: Lh_km       = 500.d0
      real(8) :: dT_amp      = 1.0d0      ! K, amplitude of temperature anomaly
    
      namelist /subinertial_wave_control/ deep, add_pert, ptop,                  &
                                          p_weak_bot, p_weak_top,               &
                                          T_surf, T_strat, p_trans, trans_width,&
                                          lonc_deg, latc_deg, Lh_km, dT_amp
    
    CONTAINS
    
    !=======================================================================
      subroutine subinertial_wave_test_init(namelist_path)
    
        use namelist_mod
    
        character(*), intent(in) :: namelist_path
        integer :: ignore
    
        ! defaults are set above; optional override by namelist
        open(11, file=namelist_path, status='old')
        read(11, nml=subinertial_wave_control, iostat=ignore)
        close(11)
    
        if (ptop <= 0.d0) ptop = 500.d0
    
      end subroutine subinertial_wave_test_init
    
    !=======================================================================
      subroutine subinertial_wave_test_set_ic(block)
    
        use namelist_mod
        use block_mod
        use tracer_mod
        use operators_mod
        use formula_mod
        use latlon_parallel_mod
        use latlon_operators_mod
    
        type(block_type), intent(inout), target :: block
    
        integer :: i, j, k
        real(8) :: rho, thetav, z_dummy, q_dummy
        real(8) :: u_loc, v_loc, t_loc, pt_loc
    
        associate ( mesh   => block%mesh            , &
                    gzs    => block%static%gzs      , & ! out
                    mgs    => block%dstate(1)%mgs   , & ! out (surface dry pressure)
                    u      => block%aux%u           , & ! work
                    v      => block%aux%v           , & ! work
                    u_lon  => block%dstate(1)%u_lon , & ! out
                    v_lat  => block%dstate(1)%v_lat , & ! out
                    t      => block%aux%t           , & ! work
                    pt     => block%dstate(1)%pt    , & ! out
                    mg     => block%dstate(1)%mg    , & ! work (full-level pressure-like)
                    ph_lev => block%dstate(1)%ph_lev )   ! work
    
          ! flat surface
          gzs%d = 0.d0
          mgs%d = p0
    
          ! Let the core compute mg and geopotential levels first (as in baroclinic test)
          call calc_mg(block, block%dstate(1))
          call calc_ph(block, block%dstate(1))
    
          z_dummy = 0.d0
          q_dummy = 0.d0
    
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
    
                call subinertial_wave_test_sample(                         &
                  deep    = deep,                                          &
                  lon     = mesh%full_lon(i),                              &
                  lat     = mesh%full_lat(j),                              &
                  p       = mg%d(i,j,k),                                   &
                  z       = z_dummy,                                       &
                  zcoords = 0,                                             &
                  u       = u_loc,                                         &
                  v       = v_loc,                                         &
                  t       = t_loc,                                         &
                  pt      = pt_loc,                                        &
                  thetav  = thetav,                                        &
                  phis    = gzs%d(i,j),                                    &
                  ps      = mgs%d(i,j),                                    &
                  rho     = rho,                                           &
                  q       = q_dummy )
    
                u%d(i,j,k)  = u_loc
                v%d(i,j,k)  = v_loc
                t%d(i,j,k)  = t_loc
                pt%d(i,j,k) = pt_loc
    
              end do
            end do
          end do
    
          ! If your dynamical core expects u_lon/v_lat to be filled from aux u/v:
          u_lon%d = u%d
          v_lat%d = v%d

          call fill_halo(u_lon)
          call fill_halo(v_lat)
          call fill_halo(pt)
          call calc_gz(block,block%dstate(1))
        !   call fill_halo()
    
        end associate
    
      end subroutine subinertial_wave_test_set_ic
    
    !=======================================================================
    ! Point-wise analytic fields
    !=======================================================================
      subroutine subinertial_wave_test_sample( deep, lon, lat, p, z, zcoords, &
                                              u, v, t, pt, thetav, phis, ps, rho, q )
    
        integer, intent(in) :: deep
        real(8), intent(in) :: lon, lat
        real(8), intent(inout) :: p
        real(8), intent(inout) :: z
        integer, intent(in) :: zcoords
    
        real(8), intent(out) :: u, v, t, pt, thetav, phis, ps, rho, q
    
        real(8) :: Tbg, dT, lnpr, lnpr0, Strans
        real(8) :: lonc, latc, Lh, dist, cosd
        real(8) :: Sp, p_mid, p_halfwidth
    
        ! outputs defaults
        phis = 0.d0
        ps   = p0
        q    = 0.d0
    
        if (zcoords == 1) then
          ! If caller gives z, you would need a p(z) relation.
          ! For this test we assume pressure coordinate usage.
          ! Provide a very simple fallback:
          p = p0 * exp(-z/8000.d0)
        end if
    
        !----------------------------
        ! Background temperature Tbg(p)
        !----------------------------
        ! Smoothly transition from warm lower troposphere to cool upper levels
        ! in log-pressure coordinate.
        lnpr  = log(max(p, 1.d0))
        lnpr0 = log(max(p_trans, 1.d0))
    
        ! Strans ~ 0 at high p, ~1 at low p
        Strans = 0.5d0 * (1.d0 - tanh( (lnpr - lnpr0)/max(trans_width,1.d-6) ))
    
        Tbg = (1.d0 - Strans) * T_surf + Strans * T_strat
    
        !----------------------------
        ! Initial perturbation confined in weak layer + localized horizontally
        !----------------------------
        dT = 0.d0
        if (add_pert == 1) then
          lonc = lonc_deg * deg2rad
          latc = latc_deg * deg2rad
          Lh   = Lh_km * 1000.d0
    
          ! great-circle distance
          cosd = sin(latc)*sin(lat) + cos(latc)*cos(lat)*cos(lon-lonc)
          cosd = max(-1.d0, min(1.d0, cosd))
          dist = a * acos(cosd)
    
          ! vertical shape: sine in pressure within weak layer
          if ( (p <= p_weak_bot) .and. (p >= p_weak_top) ) then
            Sp = sin( pi * (p - p_weak_top) / max(p_weak_bot - p_weak_top,1.d0) )
          else
            Sp = 0.d0
          end if
    
          dT = dT_amp * exp( -(dist*dist)/(Lh*Lh) ) * Sp
        end if
    
        t  = Tbg + dT
        pt = t * (p0/max(p,1.d0))**kappa
        thetav = pt  ! dry
    
        !----------------------------
        ! Background wind = 0 (clean deep/shallow contrast)
        !----------------------------
        u = 0.d0
        v = 0.d0
    
        !----------------------------
        ! Density (ideal gas, dry)
        !----------------------------
        rho = max(p,1.d0) / (Rd * max(t, 50.d0))
    
      end subroutine subinertial_wave_test_sample
    
    END MODULE subinertial_wave_test_mod