! module critical_wave_test_mod
!   use const_mod, only: r8
!   implicit none

!   private

!   public critical_wave_test_init
!   public critical_wave_test_set_ic

!   REAL(8), PARAMETER ::               &
!   a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
!   Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
!   g     = 9.80616d0,             & ! Gravity (m s^2)
!   cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
!   Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
!   Rvap  = 461.5d0,               & ! Ideal gas constnat for water vapor
!   Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
!   pi    = 3.14159265358979d0,    & ! pi
!   p0    = 100000.0d0,            & ! surface pressure (Pa)
!   kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
!   omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
!   deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

!   REAL(8), PARAMETER ::          &
! !   N_freq = 1.0d-4,               & 
!     N_freq = 6.0d-5,               & 
!   theta0 = 288.0_r8

!   REAL(8), PARAMETER ::          &
! !   pertu0     = 0.5d0      ,      & ! SF Perturbation wind velocity (m/s)
! !   pertr      = 1.d0/6.d0  ,      & ! SF Perturbation radius (Earth radii)
! !   pertup     = 1.0d0      ,      & ! Exp. perturbation wind velocity (m/s)
! !   pertexpr   = 0.1d0      ,      & ! Exp. perturbation radius (Earth radii)
!   lon_c    = pi/9.d0    ,      & ! Perturbation longitude
!   lat_c    = pi/7.d0,     & ! Perturbation latitude
!   z_c      = 5000.d0   ,      & ! Perturbation height cap
!   sigma_h  = 3.0d5      ,      & ! 水平尺度 300km
!   sigma_v  = 3000.0_r8  ,      &  ! 垂直尺度 3km
!   k_x      = 2.0_r8 * pi / 4.0d6 ,&! 水平波数
!   m_z      = 2.0_r8 * pi / 5.0d3 ,&! 垂直波数
!   pertsigma  = 500000.d0  , & ! perturbation sigma
!   dxepsilon  = 1.d-5               ! Small value for numerical derivatives

! !   integer :: deep = 0
!   integer :: pert = 0

! !   namelist /critical_wave_control/ deep
!   namelist /critical_wave_control/ pert

! CONTAINS
!   subroutine critical_wave_test_init(namelist_path)
!     use namelist_mod
    
!     character(*), intent(in) :: namelist_path

!     integer ignore

!     open(11, file=namelist_path, status='old')
!     read(11, nml=critical_wave_control, iostat=ignore)
!     close(11)

!     ! ptop = 226.0d0
!   end subroutine critical_wave_test_init

!   subroutine critical_wave_test_set_ic(block)

!     use namelist_mod
!     use block_mod
!     use operators_mod
!     use formula_mod
!     use latlon_parallel_mod
!     use latlon_operators_mod

!     type(block_type), intent(inout), target :: block

!     integer i, j, k
!     real(8) theta, t, rho, z
!     real(8) beta, term, theta_z, theta_pert
!     real(r8) :: dist_h, k_x, m_z

!     beta  = (N_freq**2) / g

!     associate (mesh   => block%mesh            , &
!         mgs    => block%dstate(1)%mgs   , & ! out
!         gzs    => block%static%gzs      , & ! out
!         lon    => block%mesh%full_lon   , &
!         lat    => block%mesh%full_lat   , &
!         u      => block%aux%u           , & ! work array
!         v      => block%aux%v           , & ! work array
!         u_lon  => block%dstate(1)%u_lon , & ! out
!         v_lat  => block%dstate(1)%v_lat , & ! out
!         t      => block%aux%t           , & ! work array
!         pt     => block%dstate(1)%pt    , & ! out
!         mg     => block%dstate(1)%mg    , & ! work array
!         gz     => block%dstate(1)%gz    , & ! work array
!         gz_lev => block%dstate(1)%gz_lev, & ! work array
!         ph     => block%dstate(1)%ph    , & ! work array
!         ph_lev => block%dstate(1)%ph_lev)   ! work array

!       gzs%d = 0.0_r8
!       mgs%d = p0

!     call calc_mg(block, block%dstate(1))
!     call calc_dmg(block,block%dstate(1))
!     call calc_ph(block, block%dstate(1))
!     ! print*,ph%d(30,30,11)
!     ! stop 999

!     do k = mesh%full_kds, mesh%full_kde
!         do j = mesh%full_jds, mesh%full_jde
!             do i = mesh%full_ids, mesh%full_ide
!                 term = (kappa * g**2) / (rd * theta0 * N_freq**2)
!                 gz%d(i,j,k) = - (1.0_r8 / beta) * log(1.0_r8 - (1.0_r8 - (ph%d(i,j,k)/p0)**kappa) / term)
!                 pt%d(i,j,k) = theta0 * exp(beta * gz%d(i,j,k))
!             end do
!         end do
!     end do

!     do k = mesh%half_kds, mesh%half_kde
!         do j = mesh%full_jds, mesh%full_jde
!             do i = mesh%full_ids, mesh%full_ide
!                 term = (kappa * g**2) / (rd * theta0 * N_freq**2)
!                 gz_lev%d(i,j,k) = - (1.0_r8 / beta) * log(1.0_r8 - (1.0_r8 - (ph_lev%d(i,j,k)/p0)**kappa) / term)
!                 ! pt%d(i,j,k) = theta0 * exp(beta * gz%d(i,j,k))
!             end do
!         end do
!     end do

!     u_lon%d = 0.0_r8
!     v_lat%d = 0.0_r8
!     if (pert .eq. 0) then
!         theta_pert = 0.0
!     else 
!         do k = mesh%full_kds, mesh%full_kde
!             do j = mesh%full_jds, mesh%full_jde
!                 do i = mesh%full_ids, mesh%full_ide
!                     dist_h = (a/scale_X) * acos(max(-1.0_r8, min(1.0_r8, &
!                      sin(lat_c)*sin(lat(j)) + &
!                      cos(lat_c)*cos(lat(j))*cos(lon(i)-lon_c))))

!                     ! 位温扰动波包
!                     theta_pert = 0.5_r8 * exp(- (dist_h/sigma_h)**2 ) * &
!                          exp(- ((gz%d(i,j,k)-z_c)/sigma_v)**2 ) * &
!                          cos(k_x * dist_h + m_z * (gz%d(i,j,k) - z_c))
!                     pt%d(i,j,k) = pt%d(i,j,k) + theta_pert 
!                 end do
!             end do
!         end do
!     end if

!     call fill_halo(u_lon)
!     call fill_halo(v_lat)
!     call fill_halo(pt)
!     gz%d = gz%d * g
!     gz_lev%d = gz_lev%d * g

!     call fill_halo(gz)
!     call fill_halo(gz_lev)

!     end associate
!   end subroutine critical_wave_test_set_ic

! end module critical_wave_test_mod

module critical_wave_test_mod
    !
    ! Critical-latitude internal-wave packet for testing TA vs NCT
    ! (Gerkema 2008 Fig.11-like frequency-band argument).
    !
    ! Background: dry, hydrostatic-pressure coordinate, u=v=0, constant N.
    ! Default case here is aimed at the user-requested N = 2 Omega regime:
    !   N = 2 Omega, omega_wave = 1.5 Omega, lat0 = 20 deg.
    ! For this case the traditional inertial latitude is about 48.6 deg,
    ! while the full-Coriolis critical latitude from Gerkema's Eq. (22)
    ! is about 64.1 deg, giving a clear poleward-extension target.
    !
    ! Two perturbation types:
    !   pert_type=1: theta-only (pt-only) Gaussian packet
    !   pert_type=2: polarized packet (u,v,w,pt) using local TA polarization
    !                at the source latitude so shallow/deep runs start from
    !                the same packet and diverge dynamically afterward.
    !
    use const_mod, only : r8
    implicit none
    private
    public :: critical_wave_test_init
    public :: critical_wave_test_set_ic
    public :: critical_wave_test_apply_forcing
    ! public :: critical_wave_test_set_params
  
    ! ---------------- constants ----------------
    real(8), parameter :: pi = 3.1415926535897932384626433832795d0
    real(8), parameter :: deg2rad = pi/180.d0
    real(8), parameter :: a_earth = 6371220.0d0
    real(8), parameter :: Rd   = 287.0d0
    real(8), parameter :: cp   = 1004.5d0
    real(8), parameter :: g    = 9.80616d0
    real(8), parameter :: p0   = 100000.0d0
    real(8), parameter :: kappa = Rd/cp
    real(8), parameter :: Omega_earth = 7.292115d-5
  
    ! ---------------- user controls ----------------
    integer :: pert_type = 0          ! 0 none, 1 theta-only, 2 polarized u/v/w/pt
    integer :: use_small_planet = 0
    real(8) :: scale_X = 1.0d0
    integer :: forcing_mode = 1       ! 0 none, 1 smooth thermal wavemaker
    integer :: forcing_zonal_uniform = 1
  
    ! Background N and reference surface theta
    real(8) :: N_freq  = 2.0d0*Omega_earth
    real(8) :: theta_s = 288.0d0
    real(8) :: N_ratio = -1.0d0

    ! Wave design
    real(8) :: omega_wave = 1.5d0*Omega_earth
    real(8) :: omega_wave_ratio = -1.0d0
    real(8) :: lon0 = 180.0d0*deg2rad
    real(8) :: lat0 = 20.0d0*deg2rad
    real(8) :: lon0_deg = -999.0d0
    real(8) :: lat0_deg = -999.0d0

    ! Packet geometry
    real(8) :: H_top = 200000.d0   ! reference depth used to set m=pi/H_top
    real(8) :: z0    = 12000.d0    ! packet center height (m)
    real(8) :: Lx    = 2000.d3     ! zonal envelope width (m)
    real(8) :: Ly    = 1500.d3     ! meridional envelope width (m)
    real(8) :: Lz    = 6000.d0     ! vertical envelope width (m)
  
    integer :: ky_from_disp = 1    ! 1: compute ky from dispersion; 0: use lambda_y
    real(8) :: lambda_y = 2500.d3  ! used if ky_from_disp=0
    real(8) :: kx = 0.0d0          ! keep 0 for mostly meridional propagation
  
    ! Amplitudes
    real(8) :: Atheta = 0.1d0      ! K (pert_type=1)
    real(8) :: Aw     = 0.05d0     ! m/s reference w amplitude used in polarization (pert_type=2)
    real(8) :: forcing_pt_amp = 2.0d-4 ! K/s peak potential-temperature forcing
    real(8) :: forcing_ncycles = 3.0d0
    real(8) :: forcing_phase = 0.0d0
    real(8) :: forcing_elapsed = 0.0d0
  
    ! Numeric safety
    real(8), parameter :: eps = 1.d-12
  
    namelist /critical_wave_control/ &
      pert_type, use_small_planet, scale_X, forcing_mode, forcing_zonal_uniform, &
      N_freq, theta_s, N_ratio, omega_wave, omega_wave_ratio, lon0, lon0_deg, lat0, lat0_deg, &
      H_top, z0, Lx, Ly, Lz, ky_from_disp, lambda_y, kx, &
      Atheta, Aw, forcing_pt_amp, forcing_ncycles, forcing_phase
  
  contains
  
    ! subroutine critical_wave_test_set_params()
    !   ! Optional: set global constants in your framework (radius, omega) if you use such.
    !   use const_mod, only : radius, omega
    !   if (use_small_planet == 1) then
    !     radius = radius / max(scale_X, 1.d0)
    !   end if
    !   omega = Omega_earth
    ! end subroutine critical_wave_test_set_params
  
    subroutine critical_wave_test_init(namelist_path)
      character(*), intent(in) :: namelist_path
      integer :: ios
      open(11, file=namelist_path, status='old', action='read')
      read(11, nml=critical_wave_control, iostat=ios)
      close(11)
      call apply_derived_inputs()
      forcing_elapsed = 0.0d0
      call check_case_inputs()
      call print_case_summary()
    end subroutine critical_wave_test_init

    subroutine critical_wave_test_set_ic(block)
      use block_mod
      use operators_mod
      use formula_mod
      use interp_mod
      use latlon_parallel_mod
      implicit none
      type(block_type), intent(inout), target :: block
  
      integer :: i, j, k
      real(8) :: beta, term
      real(8) :: rr
      real(8) :: f0, m, ky
      real(8) :: dlon, x, y, zz, amp, phase
      real(8) :: theta0_loc
  
      beta = (N_freq*N_freq) / g
      rr = a_earth
      if (use_small_planet == 1) rr = a_earth / max(scale_X, 1.d0)
  
      associate (mesh   => block%mesh, &
        mgs    => block%dstate(1)%mgs, &
        gzs    => block%static%gzs, &
        lon    => block%mesh%full_lon, &
        lat    => block%mesh%full_lat, &
        u_lon  => block%dstate(1)%u_lon, &
        v_lat  => block%dstate(1)%v_lat, &
        pt     => block%dstate(1)%pt, &
        w      => block%dstate(1)%w, &
        w_lev  => block%dstate(1)%w_lev, &
        gz     => block%dstate(1)%gz, &
        gz_lev => block%dstate(1)%gz_lev, &
        ph     => block%dstate(1)%ph, &
        ph_lev => block%dstate(1)%ph_lev )
  
        ! ---------------- surface ----------------
        gzs%d = 0.0_r8
        mgs%d = p0
  
        ! ---------------- pressure fields ----------------
        call calc_mg(block, block%dstate(1))
        call calc_dmg(block, block%dstate(1))
        call calc_ph(block, block%dstate(1))
  
        ! ---------------- background: constant N, u=v=0 ----------------
        ! z(p) inversion consistent with earlier constant-N derivation:
        ! z = -(1/beta) ln [ 1 - (1 - (p/p0)^kappa)/term ]
        term = (kappa * g*g) / (Rd * theta_s * (N_freq*N_freq))
  
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              gz%d(i,j,k) = - (1.0_r8 / beta) * log( max(eps, 1.0_r8 - (1.0_r8 - (ph%d(i,j,k)/p0)**kappa) / term ) )
              pt%d(i,j,k) = theta_s * exp( beta * gz%d(i,j,k) )
            end do
          end do
        end do
  
        do k = mesh%half_kds, mesh%half_kde
          do j = mesh%full_jds, mesh%full_jde
            do i = mesh%full_ids, mesh%full_ide
              gz_lev%d(i,j,k) = - (1.0_r8 / beta) * log( max(eps, 1.0_r8 - (1.0_r8 - (ph_lev%d(i,j,k)/p0)**kappa) / term ) )
            end do
          end do
        end do

        u_lon%d = 0.0_r8
        v_lat%d = 0.0_r8
        w_lev%d = 0.0_r8
        w%d     = 0.0_r8
  
        ! ---------------- wave numbers ----------------
        f0 = 2.0d0 * Omega_earth * sin(lat0)
        m  = pi / max(H_top, 1.d0)
  
        if (ky_from_disp == 1) then
          ! TA dispersion for kx≈0:
          ! omega^2 = (N^2 ky^2 + f0^2 m^2)/(ky^2 + m^2)
          ! => ky^2 = m^2 (omega^2 - f0^2)/(N^2 - omega^2)
          ky = m * sqrt( max(0.d0, (omega_wave*omega_wave - f0*f0) / max(eps, (N_freq*N_freq - omega_wave*omega_wave)) ) )
          ! safety: prevent absurdly huge wavelength
          if (ky < 2.d0*pi/1.0d7) ky = 2.d0*pi/1.0d7
        else
          ky = 2.0d0*pi / max(lambda_y, 1.d0)
        end if
  
        ! ---------------- perturbation ----------------
        select case (pert_type)
  
        case (1)
          ! theta-only packet:
          ! pt += Atheta * exp(-(x/Lx)^2-(y/Ly)^2-((z-z0)/Lz)^2) * cos(kx x + ky y + m(z-z0))
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              y = rr * (lat(j) - lat0)
              do i = mesh%full_ids, mesh%full_ide
                dlon = lon(i) - lon0
                if (dlon >  pi) dlon = dlon - 2.d0*pi
                if (dlon < -pi) dlon = dlon + 2.d0*pi
                x  = rr * cos(lat0) * dlon
  
                zz = gz%d(i,j,k)
                amp   = exp( - (x/Lx)**2 - (y/Ly)**2 - ((zz - z0)/Lz)**2 )
                phase = kx*x + ky*y + m*(zz - z0)
                pt%d(i,j,k) = pt%d(i,j,k) + Atheta * amp * cos(phase)
              end do
            end do
          end do
  
        case (2)
          ! Polarized packet (TA, local f0 at lat0, kx≈0):
          ! choose w = Aw * amp * cos(phase)
          ! v = -(m/ky) * Aw * amp * cos(phase)
          ! u =  (f0/omega) * (m/ky) * Aw * amp * sin(phase)
          ! theta' = - theta_bg * N^2/(g omega) * Aw * amp * sin(phase)
          do k = mesh%full_kds, mesh%full_kde
            do j = mesh%full_jds, mesh%full_jde
              y = rr * (lat(j) - lat0)
              do i = mesh%full_ids, mesh%full_ide
                dlon = lon(i) - lon0
                if (dlon >  pi) dlon = dlon - 2.d0*pi
                if (dlon < -pi) dlon = dlon + 2.d0*pi
                x  = rr * cos(lat0) * dlon
  
                zz    = gz%d(i,j,k)
                amp   = exp( - (x/Lx)**2 - (y/Ly)**2 - ((zz - z0)/Lz)**2 )
                phase = kx*x + ky*y + m*(zz - z0)
  
                theta0_loc = pt%d(i,j,k)
  
                v_lat%d(i,j,k) = v_lat%d(i,j,k) - (m/max(ky,eps)) * Aw * amp * cos(phase)
                u_lon%d(i,j,k) = u_lon%d(i,j,k) + (f0/max(omega_wave,eps)) * (m/max(ky,eps)) * Aw * amp * sin(phase)
  
                pt%d(i,j,k) = pt%d(i,j,k) - (theta0_loc*(N_freq*N_freq)/(g*max(omega_wave,eps))) * Aw * amp * sin(phase)
              end do
            end do
          end do

          do k = mesh%half_kds, mesh%half_kde
            do j = mesh%full_jds, mesh%full_jde
              y = rr * (lat(j) - lat0)
              do i = mesh%full_ids, mesh%full_ide
                dlon = lon(i) - lon0
                if (dlon >  pi) dlon = dlon - 2.d0*pi
                if (dlon < -pi) dlon = dlon + 2.d0*pi
                x  = rr * cos(lat0) * dlon

                zz    = gz_lev%d(i,j,k)
                amp   = exp( - (x/Lx)**2 - (y/Ly)**2 - ((zz - z0)/Lz)**2 )
                phase = kx*x + ky*y + m*(zz - z0)

                w_lev%d(i,j,k) = Aw * amp * cos(phase)
              end do
            end do
          end do

        case default
          ! do nothing
        end select
  
        call fill_halo(u_lon)
        call fill_halo(v_lat)
        call fill_halo(pt)
        call fill_halo(w_lev)
        call interp_run(w_lev, w)
        call fill_halo(w)

        ! convert z(m) -> gz(m^2/s^2)
        gz%d     = gz%d     * g
        gz_lev%d = gz_lev%d * g
  
        call fill_halo(gz)
        call fill_halo(gz_lev)

      end associate
    end subroutine critical_wave_test_set_ic

    subroutine print_case_summary()
      real(8) :: phi_i_deg, phi_c_deg
      real(8) :: phi_i_rad, phi_c_rad
      real(8) :: f0, m, ky, lambda_y_eff, lambda_z_eff

      phi_i_rad = traditional_inertial_latitude(omega_wave, Omega_earth)
      phi_c_rad = gerkema_critical_latitude(omega_wave, N_freq, Omega_earth)
      phi_i_deg = phi_i_rad / deg2rad
      phi_c_deg = phi_c_rad / deg2rad

      f0 = 2.0d0 * Omega_earth * sin(lat0)
      m  = pi / max(H_top, 1.d0)
      if (ky_from_disp == 1) then
        ky = m * sqrt( max(0.d0, (omega_wave*omega_wave - f0*f0) / &
          max(eps, (N_freq*N_freq - omega_wave*omega_wave)) ) )
      else
        ky = 2.0d0*pi / max(lambda_y, 1.d0)
      end if

      if (ky > eps) then
        lambda_y_eff = 2.0d0*pi / ky
      else
        lambda_y_eff = -1.0d0
      end if
      lambda_z_eff = 2.0d0*pi / max(m, eps)

      write(*,'(A)') 'critical_wave_test_mod: Gerkema Fig.11-inspired case summary'
      write(*,'(A,F12.6,A,F12.6,A,F10.3)') '  N/Omega=', N_freq/Omega_earth, &
        '  omega_wave/Omega=', omega_wave/Omega_earth, '  lat0(deg)=', lat0/deg2rad
      write(*,'(A,ES12.4,A,ES12.4)') '  N_freq(s^-1)=', N_freq, '  omega_wave(s^-1)=', omega_wave
      write(*,'(A,F10.3,A,F10.3,A,F10.3)') '  phi_i(deg)=', phi_i_deg, &
        '  phi_c_full(deg)=', phi_c_deg, '  delta(deg)=', phi_c_deg - phi_i_deg
      write(*,'(A,ES12.4,A,ES12.4)') '  lambda_y(m)=', lambda_y_eff, '  lambda_z(m)=', lambda_z_eff
      write(*,'(A,I3,A,I3,A,ES12.4,A,F8.3)') '  pert_type=', pert_type, &
        '  forcing_mode=', forcing_mode, '  forcing_pt_amp(K/s)=', forcing_pt_amp, &
        '  forcing_ncycles=', forcing_ncycles
    end subroutine print_case_summary

    subroutine apply_derived_inputs()
      if (N_ratio > 0.0d0) then
        N_freq = N_ratio * Omega_earth
      end if
      if (omega_wave_ratio > 0.0d0) then
        omega_wave = omega_wave_ratio * Omega_earth
      end if
      if (lon0_deg > -900.0d0) then
        lon0 = lon0_deg * deg2rad
      end if
      if (lat0_deg > -900.0d0) then
        lat0 = lat0_deg * deg2rad
      end if
    end subroutine apply_derived_inputs

    subroutine check_case_inputs()
      if (N_freq > 100.0d0*Omega_earth) then
        write(*,'(A)') 'critical_wave_test_mod error: N_freq is much larger than Omega.'
        write(*,'(A)') '  Fortran namelist does not evaluate expressions like 2.0d0*7.292115d-5.'
        write(*,'(A)') '  Use one of the following instead:'
        write(*,'(A)') '    N_ratio = 2.0d0'
        write(*,'(A)') '  or'
        write(*,'(A)') '    N_freq = 1.4584230d-4'
        stop 1
      end if
      if (omega_wave > 100.0d0*Omega_earth) then
        write(*,'(A)') 'critical_wave_test_mod error: omega_wave is much larger than Omega.'
        write(*,'(A)') '  Use omega_wave_ratio = 1.5d0 or omega_wave = 1.09381725d-4'
        stop 1
      end if
      if (omega_wave > N_freq) then
        write(*,'(A)') 'critical_wave_test_mod error: omega_wave > N_freq, the packet is not an IGW.'
        stop 1
      end if
      if (omega_wave <= 2.0d0*Omega_earth*sin(lat0)) then
        write(*,'(A)') 'critical_wave_test_mod error: omega_wave is not superinertial at the source latitude.'
        stop 1
      end if
      if (forcing_ncycles < 0.0d0) then
        write(*,'(A)') 'critical_wave_test_mod error: forcing_ncycles must be non-negative.'
        stop 1
      end if
    end subroutine check_case_inputs

    subroutine critical_wave_test_apply_forcing(block, dt, dstate)
      use block_mod
      use latlon_parallel_mod

      type(block_type), intent(in) :: block
      real(8), intent(in) :: dt
      type(dstate_type), intent(inout) :: dstate

      integer :: i, j, k
      real(8) :: rr, tmid, period, burst_time, window_t, temporal
      real(8) :: dlon, x, y, z_m, spatial

      if (forcing_mode == 0) return
      if (forcing_pt_amp == 0.0d0) return

      rr = a_earth
      if (use_small_planet == 1) rr = a_earth / max(scale_X, 1.d0)

      tmid = forcing_elapsed + 0.5d0 * dt
      period = 2.0d0 * pi / max(omega_wave, eps)
      burst_time = max(forcing_ncycles, eps) * period

      if (tmid >= burst_time) then
        forcing_elapsed = forcing_elapsed + dt
        return
      end if

      window_t = sin(pi * tmid / burst_time)**2
      temporal = sin(omega_wave * tmid + forcing_phase)

      associate (mesh => block%mesh, &
                 lon  => block%mesh%full_lon, &
                 lat  => block%mesh%full_lat, &
                 pt   => dstate%pt, &
                 gz   => dstate%gz)
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          y = rr * (lat(j) - lat0)
          do i = mesh%full_ids, mesh%full_ide
            z_m = gz%d(i,j,k) / g
            if (forcing_zonal_uniform == 1) then
              spatial = exp( - (y/Ly)**2 - ((z_m - z0)/Lz)**2 )
            else
              dlon = lon(i) - lon0
              if (dlon >  pi) dlon = dlon - 2.d0*pi
              if (dlon < -pi) dlon = dlon + 2.d0*pi
              x = rr * cos(lat0) * dlon
              spatial = exp( - (x/Lx)**2 - (y/Ly)**2 - ((z_m - z0)/Lz)**2 )
            end if
            pt%d(i,j,k) = pt%d(i,j,k) + dt * forcing_pt_amp * spatial * window_t * temporal
          end do
        end do
      end do
      call fill_halo(pt, async=.true.)
      end associate

      forcing_elapsed = forcing_elapsed + dt
    end subroutine critical_wave_test_apply_forcing

    pure real(8) function traditional_inertial_latitude(omega_wave_in, omega_planet)
      real(8), intent(in) :: omega_wave_in, omega_planet
      traditional_inertial_latitude = asin(clamp_unit(omega_wave_in / max(2.0d0*omega_planet, eps)))
    end function traditional_inertial_latitude

    pure real(8) function gerkema_critical_latitude(omega_wave_in, n_freq_in, omega_planet)
      real(8), intent(in) :: omega_wave_in, n_freq_in, omega_planet
      real(8) :: sin2_phi

      sin2_phi = (omega_wave_in*omega_wave_in * &
        (n_freq_in*n_freq_in + 4.0d0*omega_planet*omega_planet - omega_wave_in*omega_wave_in)) / &
        max(eps, 4.0d0*omega_planet*omega_planet*n_freq_in*n_freq_in)
      gerkema_critical_latitude = asin(sqrt(clamp_unit(sin2_phi)))
    end function gerkema_critical_latitude

    pure real(8) function clamp_unit(x)
      real(8), intent(in) :: x
      clamp_unit = min(1.0d0, max(0.0d0, x))
    end function clamp_unit

  end module critical_wave_test_mod
