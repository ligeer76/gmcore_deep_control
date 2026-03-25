module swg_wave_test_mod
  use const_mod, only: r8
  use namelist_mod

  implicit none

  public swg_wave_test_init
  public swg_wave_test_set_ic

  REAL(8), PARAMETER ::               &
  a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
  Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
  g     = 9.80616d0,             & ! Gravity (m s^2)
  cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
  Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
  Rvap  = 461.5d0,               & ! Ideal gas constnat for water vapor
  Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
  pi    = 3.14159265358979d0,    & ! pi
  p0    = 100000.0d0,            & ! surface pressure (Pa)
  kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
  omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
  deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

  REAL(8), PARAMETER ::               &
  pertu0     = 0.5d0      ,      & ! SF Perturbation wind velocity (m/s)
  pertr      = 1.d0/6.d0  ,      & ! SF Perturbation radius (Earth radii)
  pertup     = 1.0d0      ,      & ! Exp. perturbation wind velocity (m/s)
  pertexpr   = 0.1d0      ,      & ! Exp. perturbation radius (Earth radii)
  pertlon    = 2*pi/9.d0    ,      & ! Perturbation longitude (20E)
  pertlat    = pi/36.d0    ,      & ! Perturbation latitude (20N)
  pertz      = 15000.d0   ,      & ! Perturbation height cap (悬空在 15km)
  dxepsilon  = 1.d-5               ! Small value for numerical derivatives

  REAL(8), PARAMETER ::               &
  moistqlat  = 2.d0*pi/9.d0,     & ! Humidity latitudinal width
  moistqp    = 34000.d0,         & ! Humidity vertical pressure width
  moisttr    = 0.1d0,            & ! Vertical cut-off pressure for humidity
  moistqs    = 1.d-12,           & ! Humidity above cut-off
  moistq0    = 0.018d0,          & ! Maximum specific humidity
  moistqr    = 0.9d0,            & ! Maximum saturation ratio
  moisteps   = 0.622d0,          & ! Ratio of gas constants
  moistT0    = 273.16d0,         & ! Reference temperature (K)
  moistE0Ast = 610.78d0            ! Saturation vapor pressure at T0 (Pa) 


  REAL(8) ::                     &
  T0     = 300.0_r8,             &
  N_freq = 2.0_r8 * omega ,      & 
  theta0 = 300.0_r8

  ! === 惯性重力波 (pertt=2) 控制参数 ===
  REAL(8) :: pert_v0      = 0.01_r8    ! 波包振幅 (m/s)
  REAL(8) :: pert_Ly      = 500000.0_r8! 水平波长 (m)
  REAL(8) :: pert_sigma_x = 500000.0_r8! 纬向高斯包络半宽 (m)
  REAL(8) :: pert_sigma_y = 500000.0_r8! 经向高斯包络半宽 (m)
  REAL(8) :: pert_sigma_z = 2000.0_r8  ! 垂直高斯包络半宽 (m)
  
  ! 【极其关键】重置为地球自转频率，以测试 NTA 是否能穿透 30度 转折纬度
  REAL(8) :: target_omega = 2.5*omega     

  integer :: deep  = 0
  integer :: pertt = 0

  namelist /swg_wave_control/ N_freq, T0, deep, pertt
  namelist /swg_wave_control/ pert_v0, pert_Ly, pert_sigma_x, pert_sigma_y, pert_sigma_z, target_omega

contains
subroutine swg_wave_test_init(namelist_path)
  use namelist_mod
  
  character(*), intent(in) :: namelist_path
  integer ignore

  open(11, file=namelist_path, status='old')
  read(11, nml=swg_wave_control, iostat=ignore)
  close(11)

  ptop = 22600.0d0
end subroutine swg_wave_test_init

subroutine swg_wave_test_set_ic(block)

  use namelist_mod
  use block_mod
  use operators_mod
  use formula_mod
  use latlon_parallel_mod
  use latlon_operators_mod

  type(block_type), intent(inout), target :: block

  integer i, j, k
  real(8) theta, t, rho, z
  real(8) beta, term, theta_z, theta_pert

  ! ==========================================
  ! 新增 NTA 极化关系计算相关变量
  ! ==========================================
  real(r8) :: f_c, f_tilde_c, A_coeff, B_coeff, C_coeff, discriminant, tan_theta, ky, mz
  real(r8) :: lat_val, lon_val, z_val, dist_x, dist_y, dist_z, env, phase
  real(r8) :: u_pert, v_pert, w_pert, pt_pert
  real(r8), parameter :: W0 = 1.0e-3_r8  ! 强制指定的 w 振幅 (1 mm/s)，保障在线性范围内

  beta  = (N_freq**2) / g

  associate (mesh   => block%mesh            , &
      mgs    => block%dstate(1)%mgs   , & ! out
      gzs    => block%static%gzs      , & ! out
      lon    => block%mesh%full_lon   , &
      lat    => block%mesh%full_lat   , &
      u      => block%aux%u           , & ! work array
      v      => block%aux%v           , & ! work array
      u_lon  => block%dstate(1)%u_lon , & ! out
      v_lat  => block%dstate(1)%v_lat , & ! out
      t      => block%aux%t           , & ! work array
      pt     => block%dstate(1)%pt    , & ! out
      mg     => block%dstate(1)%mg    , & ! work array
      gz     => block%dstate(1)%gz    , & ! work array
      gz_lev => block%dstate(1)%gz_lev, & ! work array
      ph     => block%dstate(1)%ph    , & ! work array
      ph_lev => block%dstate(1)%ph_lev, & ! work array
      w_lev  => block%dstate(1)%w_lev )   ! 【新增】非静力半层垂直速度 out

    gzs%d = 0.0_r8
    mgs%d = p0

  call calc_mg(block, block%dstate(1))

  z = 0
  ! 1. 计算背景场 (完全静止的大气)
  do k = mesh%full_kds, mesh%full_kde
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        call swg_wave_test(  &
          deep   =deep            , &
          pertt  =pertt           , &
          X      =scale_X         , &
          lon    =mesh%full_lon(i), &
          lat    =mesh%full_lat(j), &
          p      =mg%d(i,j,k)     , &
          z      =z               , &
          zcoords=0               , &
          u      =u%d(i,j,k)      , &
          v      =v%d(i,j,k)      , &
          t      =t%d(i,j,k)      , &
          thetav = pt%d(i,j,k)    , &
          phis   =gzs%d(i,j)      , &
          ps     =mgs%d(i,j)      , &
          rho    =rho)
      end do
    end do
  end do

  ! 初始化所有风场为 0
  ! u_lon%d = 0.0_r8
  ! v_lat%d = 0.0_r8
  w_lev%d = 0.0_r8

  ! ======================================================================
  ! 2. 注入极其严密的 3D NTA 初值波包 (IVP)
  ! ======================================================================
  if (pertt .eq. 2) then
    ! 计算波源中心的科氏参数
    f_c       = 2.0_r8 * omega * sin(pertlat)
    f_tilde_c = 2.0_r8 * omega * cos(pertlat)

    ! Gerkema (2008) 频散关系求解 tan(theta) = m/k_y
    A_coeff = target_omega**2 - f_c**2
    B_coeff = - 2.0_r8 * f_c * f_tilde_c
    C_coeff = target_omega**2 - N_freq**2 - f_tilde_c**2
    discriminant = B_coeff**2 - 4.0_r8 * A_coeff * C_coeff

    if (discriminant >= 0.0_r8) then
      tan_theta = (-B_coeff + sqrt(discriminant)) / (2.0_r8 * A_coeff)
    else
      print *, 'ERROR: NTA Dispersion Error. Target frequency forbidden at perturbation latitude!'
      stop
    end if

    ky = 2.0_r8 * pi / pert_Ly
    mz = ky * tan_theta

    ! -------------------------------------------------------------
    ! 2a. 对全层预报变量 (u, v, pt) 施加 NTA 极化响应
    ! -------------------------------------------------------------
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          lat_val = mesh%full_lat(j)
          lon_val = mesh%full_lon(i)
          z_val   = gz%d(i,j,k) / g
          
          dist_x = a * cos(lat_val) * (lon_val - pertlon)
          dist_y = a * (lat_val - pertlat)
          dist_z = z_val - pertz
          
          ! 真正的 3D 悬空高斯包络
          env = exp( - (dist_x**2)/(2.0_r8 * pert_sigma_x**2) &
                     - (dist_y**2)/(2.0_r8 * pert_sigma_y**2) &
                     - (dist_z**2)/(2.0_r8 * pert_sigma_z**2) )
                     
          phase = ky * dist_y + mz * dist_z

          ! 绝对严密的 NTA 极化关系
          w_pert  = W0 * env * cos(phase)
          v_pert  = - W0 * tan_theta * env * cos(phase)
          u_pert  = W0 * (f_c * tan_theta + f_tilde_c) / target_omega * env * sin(phase)
          ! 热力学关系 (利用刚刚从背景场得到的 pt%d)
          pt_pert = W0 * (pt%d(i,j,k) * N_freq**2) / (g * target_omega) * env * sin(phase)

          u%d(i,j,k)  = u%d(i,j,k)  + u_pert
          v%d(i,j,k)  = v%d(i,j,k)  + v_pert
          pt%d(i,j,k) = pt%d(i,j,k) + pt_pert
        end do
      end do
    end do

    ! -------------------------------------------------------------
    ! 2b. 对半层非静力垂直风场 (w_lev) 施加对应的垂直扰动
    ! -------------------------------------------------------------
    do k = mesh%half_kds, mesh%half_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          lat_val = mesh%full_lat(j)
          lon_val = mesh%full_lon(i)
          z_val   = gz_lev%d(i,j,k) / g
          
          dist_x = a * cos(lat_val) * (lon_val - pertlon)
          dist_y = a * (lat_val - pertlat)
          dist_z = z_val - pertz
          
          env = exp( - (dist_x**2)/(2.0_r8 * pert_sigma_x**2) &
                     - (dist_y**2)/(2.0_r8 * pert_sigma_y**2) &
                     - (dist_z**2)/(2.0_r8 * pert_sigma_z**2) )
                     
          phase = ky * dist_y + mz * dist_z

          ! w 风场赋予
          w_pert = W0 * env * cos(phase)
          w_lev%d(i,j,k) = w_lev%d(i,j,k) + w_pert
        end do
      end do
    end do

  end if
  ! ======================================================================

  call fill_halo(u)
  call fill_halo(v)
  
  call wind_a2c_operator(u, v, u_lon, v_lat)

  call fill_halo(u_lon)
  call fill_halo(v_lat)
  call fill_halo(pt)
  call fill_halo(w_lev) ! 确保半层垂直速度完成 Halo 交换

  gz%d = gz%d * g
  call fill_halo(gz)

  end associate
end subroutine swg_wave_test_set_ic

!=======================================================================
!    Generate the swg instability initial conditions (Background Only)
!=======================================================================
SUBROUTINE swg_wave_test(deep,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho) &
  BIND(c, name = "swg_wave_test")

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: deep, pertt
  REAL(8), INTENT(IN)  :: lon, lat, X
  REAL(8), INTENT(INOUT) :: p, z
  INTEGER, INTENT(IN) :: zcoords
  REAL(8), INTENT(OUT) :: u, v, t, thetav, phis, ps, rho

  if (zcoords .eq. 1) then
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)
  else
    CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t)
  end if

  ! 强制保持大气绝对静止
  ps = p0
  u = 0.d0
  v = 0.1d0
  phis = 0.d0
  rho = p / (Rd * t)
  thetav = theta0 * exp(N_freq**2/g * z)

END SUBROUTINE swg_wave_test

!-----------------------------------------------------------------------
! Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
SUBROUTINE evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)

  INTEGER, INTENT(IN)  :: deep
  REAL(8), INTENT(IN)  :: X, lon, lat, z
  REAL(8), INTENT(OUT) :: p, t

  REAL(8) :: TT, NNG

  TT  = g**2 / (cp * N_freq**2)
  NNG = N_freq **2 / g

  t = TT + (T0 - TT ) * exp(z / g * N_freq**2)
  
  ! 严密计算离散静力气压
  p = p0 * exp( - cp / Rd * (NNG * z - log(abs(TT+(T0-TT)*exp(NNG*z))) + log(T0)) )

END SUBROUTINE evaluate_pressure_temperature

!-----------------------------------------------------------------------
! Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
SUBROUTINE evaluate_z_temperature(deep, X, lon, lat, p, z, t)
  
  INTEGER, INTENT(IN)  :: deep
  REAL(8), INTENT(IN)  :: X, lon, lat, p
  REAL(8), INTENT(OUT) :: z, t
  INTEGER :: ix
  REAL(8) :: z0, z1, z2, p0, p1, p2

  z0 = 0.d0
  z1 = 10000.d0
  CALL evaluate_pressure_temperature(deep, X, lon, lat, z0, p0, t)
  CALL evaluate_pressure_temperature(deep, X, lon, lat, z1, p1, t)

  DO ix = 1, 100
    z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z2, p2, t)
    IF (ABS((p2 - p)/p) .lt. 1.0d-8) EXIT
    z0 = z1
    p0 = p1
    z1 = z2
    p1 = p2
  END DO

  z = z2
  CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p0, t)

END SUBROUTINE evaluate_z_temperature

end module swg_wave_test_mod