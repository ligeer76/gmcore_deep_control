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
    ! Critical-latitude internal-wave packet for testing TA vs NCT (Gerkema 2008 Fig.11(b)-like).
    !
    ! Background: dry, hydrostatic-pressure coordinate, u=v=0, constant N = Omega.
    ! Forcing frequency: omega_wave = 0.5 * Omega (fixed).
    ! Center latitude: lat0 ~ 10deg.
    !
    ! Two perturbation types:
    !   pert_type=1: theta-only (pt-only) Gaussian packet
    !   pert_type=2: polarized packet (u,v,pt) using TA f-plane polarization at f0=2Ωsin(lat0)
    !
    use const_mod, only : r8
    implicit none
    private
    public :: critical_wave_test_init
    public :: critical_wave_test_set_ic
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
    integer :: pert_type = 2          ! 1 theta-only, 2 polarized u/v/pt
    integer :: use_small_planet = 0
    real(8) :: scale_X = 1.0d0
  
    ! Background N and reference surface theta
    real(8) :: N_freq  = Omega_earth      ! <-- your request: N = Omega
    real(8) :: theta_s = 288.0d0
  
    ! Wave design
    real(8) :: omega_wave = 0.5d0*Omega_earth  ! <-- your request: omega = 0.5 Omega
    real(8) :: lon0 = 0.0d0
    real(8) :: lat0 = 10.0d0*deg2rad           ! <-- your request: 10 deg
  
    ! Packet geometry
    real(8) :: H_top = 20000.d0    ! model top for m=pi/H_top
    real(8) :: z0    = 5000.d0     ! packet center height (m)
    real(8) :: Lx    = 1200.d3     ! zonal envelope width (m)
    real(8) :: Ly    = 1200.d3     ! meridional envelope width (m)
    real(8) :: Lz    = 4000.d0     ! vertical envelope width (m)
  
    integer :: ky_from_disp = 1    ! 1: compute ky from dispersion; 0: use lambda_y
    real(8) :: lambda_y = 2500.d3  ! used if ky_from_disp=0
    real(8) :: kx = 0.0d0          ! keep 0 for mostly meridional propagation
  
    ! Amplitudes
    real(8) :: Atheta = 0.2d0      ! K (pert_type=1)
    real(8) :: Aw     = 0.2d0      ! m/s reference w amplitude used in polarization (pert_type=2)
  
    ! Numeric safety
    real(8) :: eps = 1.d-12
  
    namelist /critical_wave_control/ &
      pert_type, use_small_planet, scale_X, &
      N_freq, theta_s, omega_wave, lon0, lat0, &
      H_top, z0, Lx, Ly, Lz, ky_from_disp, lambda_y, kx, &
      Atheta, Aw
  
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
    end subroutine critical_wave_test_init
  
    subroutine critical_wave_test_set_ic(block)
      use block_mod
      use operators_mod
      use formula_mod
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
        gz     => block%dstate(1)%gz, &
        gz_lev => block%dstate(1)%gz_lev, &
        ph     => block%dstate(1)%ph, &
        ph_lev => block%dstate(1)%ph_lev )
  
        ! ---------------- surface ----------------
        gzs%d = 0.0_r8
        mgs%d = p0
  
        ! ---------------- pressure fields ----------------
        call calc_mg(block, block%dstate(1))
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
  
        case default
          ! do nothing
        end select
  
        call fill_halo(u_lon)
        call fill_halo(v_lat)
        call fill_halo(pt)
  
        ! convert z(m) -> gz(m^2/s^2)
        gz%d     = gz%d     * g
        gz_lev%d = gz_lev%d * g
  
        call fill_halo(gz)
        call fill_halo(gz_lev)
  
      end associate
    end subroutine critical_wave_test_set_ic
  
  end module critical_wave_test_mod