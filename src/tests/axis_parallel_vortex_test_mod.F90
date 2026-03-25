module axis_parallel_vortex_test_mod
    !=======================================================================
    !
    !  Axis-parallel vortex (NCT-sensitive) initial condition
    !  - Dry air, hydrostatic, terrain-following pressure (mg ~ p) coordinate
    !  - Background: weakly stable theta(p)
    !  - Vortex: axisymmetric tangential wind at low latitude
    !  - Hydrostatic balance is enforced discretely via:
    !        dPhi = -Rd * T * d ln p   (dry)
    !
    !  Usage pattern follows baroclinic_wave_test_mod / supercell_test_mod:
    !    - set gzs, mgs
    !    - calc_mg -> full-level mg (pressure)
    !    - vert_coord_calc_mg_lev -> half-level mg_lev (pressure)
    !    - integrate hydrostatic to get z_lev (m), then convert to gz_lev (m2/s2)
    !    - compute u/v/t/pt and u_lon/v_lat halos
    !
    !=======================================================================
    
      use mpi
      use string
      use flogger
      use math_mod
      use const_mod, only: r8, pi, rad, rd, g, omega, a => radius
      use block_mod
      use formula_mod,   only: modified_potential_temperature
      use vert_coord_mod
      use namelist_mod
      use tracer_mod
      use latlon_parallel_mod
      use operators_mod

      implicit none
      private
    
      public axis_parallel_vortex_test_init
      public axis_parallel_vortex_test_set_ic
    
      !--------------------------
      ! Background / vortex params
      !--------------------------
      real(r8), parameter :: p0     = 100000.0_r8  ! dry surface pressure (Pa)
      real(r8), parameter :: theta0 = 300.0_r8     ! reference theta (K)
      real(r8), parameter :: eps_th = 0.01_r8      ! weak stability in theta(p)
      real(r8), parameter :: p_ref  = 100000.0_r8  ! reference pressure for theta
    
      ! Vortex center (low latitude is important for NCT sensitivity)
      real(r8), parameter :: cen_lat = 10.0_r8 * rad
      real(r8), parameter :: cen_lon =  0.0_r8 * rad
    
      ! Tangential wind profile parameters
      real(r8), parameter :: Rm    = 200000.0_r8   ! radius of max wind (m)
      real(r8), parameter :: Vmax  = 30.0_r8       ! max tangential wind (m/s)
    
    contains
    
      subroutine axis_parallel_vortex_test_init()
        ! Nothing to read for now; keep consistent with other test modules.
        ! Ensure no moisture tracer is required.
      end subroutine axis_parallel_vortex_test_init
    
    
      pure real(r8) function theta_bg(p) result(th)
        real(r8), intent(in) :: p
        ! Weakly stable theta(p): theta = theta0 * (1 + eps * ln(p_ref/p))
        th = theta0 * (1.0_r8 + eps_th * log(p_ref / max(p, 1.0_r8)))
      end function theta_bg
    
    
      pure real(r8) function temp_bg(p) result(t)
        real(r8), intent(in) :: p
        real(r8), parameter :: kappa = 2.0_r8/7.0_r8  ! Rd/cp for dry air
        t = theta_bg(p) * (max(p, 1.0_r8) / p_ref)**kappa
      end function temp_bg
    
    
      pure real(r8) function vtheta_profile(r) result(vt)
        real(r8), intent(in) :: r
        ! Smooth single-peak profile:
        ! vtheta = Vmax * (r/Rm) * exp(1 - r/Rm)
        if (r <= 0.0_r8) then
          vt = 0.0_r8
        else
          vt = Vmax * (r / Rm) * exp(1.0_r8 - r / Rm)
        end if
      end function vtheta_profile
    
    
      pure subroutine gc_distance_and_bearing(lon, lat, lon0, lat0, r, bearing)
        ! Great-circle distance r (m) and initial bearing (rad) from center -> point.
        ! bearing is measured clockwise from north.
        real(r8), intent(in)  :: lon, lat, lon0, lat0
        real(r8), intent(out) :: r, bearing
        real(r8) dlon, cosds, ds
    
        dlon = lon - lon0
    
        cosds = sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(dlon)
        cosds = max(-1.0_r8, min(1.0_r8, cosds))
        ds    = acos(cosds)
        r     = a * ds
    
        ! initial bearing formula
        bearing = atan2( sin(dlon)*cos(lat), &
                         cos(lat0)*sin(lat) - sin(lat0)*cos(lat)*cos(dlon) )
      end subroutine gc_distance_and_bearing
    
    
      subroutine axis_parallel_vortex_test_set_ic(block)
        type(block_type), intent(inout), target :: block
    
        integer :: i, j, k
        real(r8) :: r, brg, dir, vt, p_half_up, p_half_dn, p_mid
        real(r8) :: Tmid, dz, ztmp
    
        associate (mesh   => block%mesh            , &
                   lon    => block%mesh%full_lon   , &
                   lat    => block%mesh%full_lat   , &
                   gzs    => block%static%gzs      , &
                   mgs    => block%dstate(1)%mgs   , &
                   mg     => block%dstate(1)%mg    , &
                   mg_lev => block%dstate(1)%mg_lev, &
                   p      => block%dstate(1)%p     , &
                   p_lev  => block%dstate(1)%p_lev , &
                   z      => block%dstate(1)%gz    , &
                   z_lev  => block%dstate(1)%gz_lev, &
                   u      => block%aux      %u     , &
                   u_lon  => block%dstate(1)%u_lon , &
                   v      => block%aux      %v     , &
                   v_lat  => block%dstate(1)%v_lat , &
                   t      => block%aux      %t     , &
                   pt     => block%dstate(1)%pt    )
    
          !--------------------------
          ! 0) flat surface, constant ps
          !--------------------------
          gzs%d = 0.0_r8
          mgs%d = p0
          call fill_halo(mgs)
    
          !--------------------------
          ! 1) compute full-level pressure (mg) from vertical coordinate definition
          !--------------------------
          call calc_mg(block, block%dstate(1))
          call fill_halo(mg)
    
          ! Copy mg -> p for convenience / consistency (some parts of code expect p)
          p%d = mg%d
          call fill_halo(p)
    
          !--------------------------
          ! 2) compute half-level pressure mg_lev and p_lev
          !--------------------------
          do k = mesh%half_kds, mesh%half_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                mg_lev%d(i,j,k) = vert_coord_calc_mg_lev(k, mgs%d(i,j))
              end do
            end do
          end do
          call fill_halo(mg_lev)
    
          p_lev%d = mg_lev%d
          call fill_halo(p_lev)
    
          !--------------------------
          ! 3) Hydrostatic integration for height z_lev (meters)
          !
          ! Index convention in your other tests:
          !   half level k=half_kds is TOP, k=half_kde is BOTTOM (surface)
          ! and z_lev(bottom) is set to 0.
          !
          ! Discrete hydrostatic (dry):
          !   Phi(k) - Phi(k+1) = Rd * T_layer * ln(p(k+1)/p(k))
          ! with k top->bottom; here we integrate upward from surface:
          !   z_lev(bottom)=0, then for k=bottom-1..top:
          !     dz = (Rd*T_layer/g) * ln(p_dn/p_up)
          !--------------------------
          ! First set bottom half level height = 0
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              z_lev%d(i,j,mesh%half_kde) = 0.0_r8
            end do
          end do
    
          do k = mesh%half_kde-1, mesh%half_kds, -1
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                p_half_up = p_lev%d(i,j,k)     ! smaller (upper)
                p_half_dn = p_lev%d(i,j,k+1)   ! larger (lower)
    
                ! Use geometric mean pressure for layer midpoint
                p_mid = sqrt(max(p_half_up,1.0_r8) * max(p_half_dn,1.0_r8))
                Tmid  = temp_bg(p_mid)
    
                dz = (rd * Tmid / g) * log(p_half_dn / max(p_half_up,1.0_r8))
                z_lev%d(i,j,k) = z_lev%d(i,j,k+1) + dz
              end do
            end do
          end do
    
          ! full-level geopotential height z (meters) is average of adjacent half levels
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                ztmp = 0.5_r8 * (z_lev%d(i,j,k) + z_lev%d(i,j,k+1))
                ! We'll store gz (m^2/s^2) later by multiplying g, following your style.
                z%d(i,j,k) = ztmp
              end do
            end do
          end do
    
          call fill_halo(z_lev)
          call fill_halo(z)
    
          !--------------------------
          ! 4) Thermodynamics: T(p), pt from modified_potential_temperature (dry qv=0)
          !--------------------------
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                t%d(i,j,k)  = temp_bg(p%d(i,j,k))
                pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), p%d(i,j,k), 0.0_r8)
              end do
            end do
          end do
          call fill_halo(t)
          call fill_halo(pt)
    
          !--------------------------
          ! 5) Winds: axisymmetric tangential wind around (cen_lon, cen_lat)
          !    Convert tangential speed to (u east, v north) using bearing.
          !    Cyclonic (NH) means counterclockwise: direction = bearing - pi/2.
          !--------------------------
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%full_ids, mesh%full_ide
                call gc_distance_and_bearing(lon(i), lat(j), cen_lon, cen_lat, r, brg)
                vt  = vtheta_profile(r)
                dir = brg - 0.5_r8*pi      ! from north, clockwise
    
                ! (u east, v north)
                u%d(i,j,k) = vt * sin(dir)
                v%d(i,j,k) = vt * cos(dir)
              end do
            end do
          end do
          call fill_halo(u)
          call fill_halo(v)
    
          !--------------------------
          ! 6) Staggered winds u_lon / v_lat
          !--------------------------
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              do i = mesh%half_ids, mesh%half_ide
                u_lon%d(i,j,k) = 0.5_r8 * (u%d(i,j,k) + u%d(i+1,j,k))
              end do
            end do
          end do
          call fill_halo(u_lon)
    
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%half_jds, mesh%half_jde
              do i = mesh%full_ids, mesh%full_ide
                v_lat%d(i,j,k) = 0.5_r8 * (v%d(i,j,k) + v%d(i,j+1,k))
              end do
            end do
          end do
          call fill_halo(v_lat)
    
          !--------------------------
          ! 7) Convert stored heights (m) to geopotential (m^2/s^2) like other tests
          !--------------------------
          z_lev%d = z_lev%d * g
          z%d     = z%d     * g
          call fill_halo(z_lev)
          call fill_halo(z)
    
        end associate
    
        if (proc%is_root()) call log_notice('Axis-parallel vortex IC set (dry hydrostatic, weakly stable background).')
    
      end subroutine axis_parallel_vortex_test_set_ic
    
    end module axis_parallel_vortex_test_mod