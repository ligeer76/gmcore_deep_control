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
    pertlon    = pi/9.d0    ,      & ! Perturbation longitude
    pertlat    = pi/9.d0,     & ! Perturbation latitude
    pertz      = 10000.d0   ,      & ! Perturbation height cap
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
    N_freq = 2.0 * omega,          &
    ! N_freq = 1.0d-2,               & 
    ! N_freq = 6.0d-5,               & 
    theta0 = 300.0_r8

    ! === 新增：惯性重力波 (pertt=2) 控制参数 ===
    !!!! chatgpt
    REAL(r8) :: packet_v0       = 0.15_r8          ! v 振幅 (m/s)
    REAL(r8) :: packet_Ly       = 8.0e5_r8         ! 经向载波波长 (m)
    REAL(r8) :: packet_sigma_x  = 3.0e6_r8         ! zonal envelope 半宽 (m)
    REAL(r8) :: packet_sigma_y  = 3.0e5_r8         ! meridional envelope 半宽 (m)
    REAL(r8) :: packet_sigma_z  = 4.0e3_r8         ! vertical envelope 半宽 (m)
    REAL(r8) :: packet_lon0     = 20.0_r8 * deg2rad
    REAL(r8) :: packet_lat0     = 10.0_r8 * deg2rad
    REAL(r8) :: packet_z0       = 10000.0_r8       ! 波包垂直中心 (m)
    REAL(r8) :: packet_omega    = 1.0_r8 * omega   ! 目标频率 (rad/s)
    INTEGER  :: packet_branch   = +1               ! +1/-1: 两个 m/ky 根
    
    namelist /swg_wave_control/ packet_v0, packet_Ly, packet_sigma_x, packet_sigma_y, &
                                packet_sigma_z, packet_lon0, packet_lat0, packet_z0,  &
                                packet_omega, packet_branch
    !!!! chatgpt
    REAL(8) :: pert_v0      = 2.1_r8    ! 波包振幅 (m/s)
    REAL(8) :: pert_Ly      = 100000.0_r8! 水平波长 (m)
    REAL(8) :: pert_sigma_x = 200000.0_r8! 【新增】纬向高斯包络半宽 (m)，控制经度范围
    REAL(8) :: pert_sigma_y = 200000.0_r8! 经向高斯包络半宽 (m)
    REAL(8) :: pert_sigma_z = 50000.0_r8  ! 垂直高斯包络半宽 (m)
    ! REAL(8) :: target_omega = 0.8d-4     ! 【新增】目标波频率 (rad/s)，建议0.8d-4以测试60度转折
    REAL(8) :: target_omega = 1.0 * omega     ! 【新增】目标波频率 (rad/s)，建议0.8d-4以测试60度转折

    integer :: deep  = 0
    integer :: pertt = 0

    namelist /swg_wave_control/ N_freq
    namelist /swg_wave_control/ T0
    namelist /swg_wave_control/ deep,pertt
    ! namelist /swg_wave_control/ pert_v0, pert_Ly, pert_sigma_y, pert_sigma_z
    namelist /swg_wave_control/ pert_v0, pert_Ly, pert_sigma_x, pert_sigma_y, pert_sigma_z, target_omega

contains
  subroutine swg_wave_test_init(namelist_path)
    use namelist_mod
    
    character(*), intent(in) :: namelist_path

    integer ignore

    open(11, file=namelist_path, status='old')
    read(11, nml=swg_wave_control, iostat=ignore)
    close(11)

    ptop = 226.0d0
  end subroutine swg_wave_test_init

  subroutine swg_wave_test_set_ic(block)

    use namelist_mod
    use block_mod
    use operators_mod
    use formula_mod
    use latlon_parallel_mod
    use latlon_operators_mod
    use interp_mod

    type(block_type), intent(inout), target :: block

    ! integer i, j, k
    real(8) theta, t, rho, z
    real(8) beta, term, theta_z, theta_pert
    real(r8) :: dist_h, k_x, m_z

    integer :: i, j, k
    real(r8) :: ky, mz
    real(r8) :: z_here, t_here, z_half
    real(r8) :: up, vp, wp, ptp, tp
    real(r8), allocatable :: zfull(:)
    real(r8) :: dz1

    beta  = (N_freq**2) / g

    associate (mesh   => block%mesh            , &
        mgs    => block%dstate(1)%mgs   , & ! out
        gzs    => block%static%gzs      , & ! out
        lon    => block%mesh%full_lon   , &
        lat    => block%mesh%full_lat   , &
        u      => block%aux%u           , & ! work array
        v      => block%aux%v           , & ! work array
        w      => block%dstate(1)%w     , & ! chatgpt 
        w_lev  => block%dstate(1)%w_lev , & ! chatgpt
        u_lon  => block%dstate(1)%u_lon , & ! out
        v_lat  => block%dstate(1)%v_lat , & ! out
        t      => block%aux%t           , & ! work array
        pt     => block%dstate(1)%pt    , & ! out
        mg     => block%dstate(1)%mg    , & ! work array
        gz     => block%dstate(1)%gz    , & ! work array
        gz_lev => block%dstate(1)%gz_lev, & ! work array
        ph     => block%dstate(1)%ph    , & ! work array
        ph_lev => block%dstate(1)%ph_lev)   ! work array

      gzs%d = 0.0_r8
      mgs%d = p0

    call calc_mg(block, block%dstate(1))
    call fill_halo(mg)
    ! call calc_dmg(block,block%dstate(1))
    ! call calc_ph(block, block%dstate(1))
    ! print*,ph%d(30,30,11)
    ! stop 999
    z = 0
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
            ! thetav =thetav          , &
            thetav = pt%d(i,j,k)    , &
            phis   =gzs%d(i,j)      , &
            ps     =mgs%d(i,j)      , &
            rho    =rho)
        end do
      end do
    end do
    ! do k = mesh%full_kds, mesh%full_kde
    !   do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
    !     do i = mesh%full_ids, mesh%full_ide
    !     ! print*,u%d(i,j,k),"u",i,j,k
    !     ! print*,v%d(i,j,k),"v"
    !     print*,pt%d(i,j,k),"pt",i,j,k
    !     print*,t%d(i,j,k),"t"              
    !     end do
    !   end do
    ! end do
    ! stop 999

    ! print*,u%d(i,j,k),"u",i,j,k
    ! print*,v%d(i,j,k),"v"
    ! print*,pt%d(1,1,1),"pt"
    ! print*,t%d(1,1,1),"t"              
    ! stop 999                      

    ! do k = mesh%full_kds, mesh%full_kde
    !     do j = mesh%full_jds, mesh%full_jde
    !         do i = mesh%full_ids, mesh%full_ide
    !             term = (kappa * g**2) / (rd * theta0 * N_freq**2)
    !             gz%d(i,j,k) = - (1.0_r8 / beta) * log(1.0_r8 - (1.0_r8 - (ph%d(i,j,k)/p0)**kappa) / term)
    !             pt%d(i,j,k) = theta0 * exp(beta * gz%d(i,j,k))
    !         end do
    !     end do
    ! end do

    ! do k = mesh%half_kds, mesh%half_kde
    !     do j = mesh%full_jds, mesh%full_jde
    !         do i = mesh%full_ids, mesh%full_ide
    !             term = (kappa * g**2) / (rd * theta0 * N_freq**2)
    !             gz_lev%d(i,j,k) = - (1.0_r8 / beta) * log(1.0_r8 - (1.0_r8 - (ph_lev%d(i,j,k)/p0)**kappa) / term)
    !             ! pt%d(i,j,k) = theta0 * exp(beta * gz%d(i,j,k))
    !         end do
    !     end do
    ! end do

    ! u_lon%d = 0.0_r8
    ! v_lat%d = 0.0_r8
    ! if (pert .eq. 0) then
    !     theta_pert = 0.0
    ! else 
    !     do k = mesh%full_kds, mesh%full_kde
    !         do j = mesh%full_jds, mesh%full_jde
    !             do i = mesh%full_ids, mesh%full_ide
    !                 dist_h = (a/scale_X) * acos(max(-1.0_r8, min(1.0_r8, &
    !                  sin(lat_c)*sin(lat(j)) + &
    !                  cos(lat_c)*cos(lat(j))*cos(lon(i)-lon_c))))

    !                 ! 位温扰动波包
    !                 theta_pert = 0.5_r8 * exp(- (dist_h/sigma_h)**2 ) * &
    !                      exp(- ((gz%d(i,j,k)-z_c)/sigma_v)**2 ) * &
    !                      cos(k_x * dist_h + m_z * (gz%d(i,j,k) - z_c))
    !                 pt%d(i,j,k) = pt%d(i,j,k) + theta_pert 
    !             end do
    !         end do
    !     end do
    ! end if
    call fill_halo(u)
    call fill_halo(v)
    ! call fill_halo(u_lon)
    ! call fill_halo(v_lat)
    ! call fill_halo(pt)
    !     print*,u%d(i,j,k),"u",i,j,k
    ! print*,v%d(i,j,k),"v"
    ! print*,pt%d(1,1,1),"pt"
    ! print*,t%d(1,1,1),"t"              
    ! stop 999                      

    ! call fill_halo(pt)
    ! call fill_halo(w_lev)

    ! call interp_run(w_lev, w)
    ! call fill_halo(w)

    if (pertt .eq. 3) then
      call add_free_igw_packet_ic(block)
    !   call packet_get_wavenumbers(ky, mz)
  
    !   allocate(zfull(mesh%full_kds:mesh%full_kde))
  
    !   do j = mesh%full_jds, mesh%full_jde
    !     do i = mesh%full_ids, mesh%full_ide
  
    !       ! ---- full levels: u, v, pt, t ----
    !       do k = mesh%full_kds, mesh%full_kde
    ! ! print*,u%d(i,j,k),"u",i,j,k
    ! ! print*,v%d(i,j,k),"v"
    ! ! print*,pt%d(1,1,1),"pt"
    ! ! print*,t%d(1,1,1),"t"              
    ! ! stop 999                      
      
    !         call evaluate_z_temperature(deep, scale_X, mesh%full_lon(i), mesh%full_lat(j), &
    !                                     mg%d(i,j,k), z_here, t_here)
    !         zfull(k) = z_here
  
    !         call packet_polarization_full(mesh%full_lon(i), mesh%full_lat(j), z_here,        &
    !                                       pt%d(i,j,k), t%d(i,j,k), ky, mz, up, vp, ptp, tp)
    ! !   print*,u%d(i,j,k),"u",i,j,k
    ! ! print*,v%d(i,j,k),"v"
    ! ! print*,pt%d(1,1,1),"pt"
    ! ! print*,t%d(1,1,1),"t"              
    ! ! stop 999                      

    !         u%d(i,j,k)  = u%d(i,j,k)  + up
    !         v%d(i,j,k)  = v%d(i,j,k)  + vp
    !         pt%d(i,j,k) = pt%d(i,j,k) + ptp
    !         t%d(i,j,k)  = t%d(i,j,k)  + tp
    ! print*,u%d(i,j,k),up,"u",i,j,k
    ! print*,v%d(i,j,k),vp,"v"
    ! print*,pt%d(1,1,1),"pt"
    ! print*,t%d(1,1,1),"t"              
    ! stop 999                      
    !       end do
  
    !       ! ---- half levels: w_lev ----
    !       do k = mesh%half_kds, mesh%half_kde
    !         if (k .eq. mesh%half_kds) then
    !           dz1 = zfull(mesh%full_kds+1) - zfull(mesh%full_kds)
    !           z_half = max(0.0_r8, zfull(mesh%full_kds) - 0.5_r8*dz1)
    !         elseif (k .eq. mesh%half_kde) then
    !           dz1 = zfull(mesh%full_kde) - zfull(mesh%full_kde-1)
    !           z_half = zfull(mesh%full_kde) + 0.5_r8*dz1
    !         else
    !           z_half = 0.5_r8 * ( zfull(k) + zfull(k-1) )
    !         end if
  
    !         call packet_vertical_velocity(mesh%full_lon(i), mesh%full_lat(j), z_half, ky, mz, wp)
    !         w_lev%d(i,j,k) = wp
    !       end do
  
      !   end do
      ! end do
    end if
    
    call wind_a2c_operator(u, v, u_lon, v_lat)

    call fill_halo(u_lon)
    call fill_halo(v_lat)
    call fill_halo(pt)
    ! call fill_halo(pt)
    call fill_halo(w_lev)

    call interp_run(w_lev, w)
    call fill_halo(w)

    ! gz%d = gz%d * g
    ! gz_lev%d = gz_lev%d * g

    ! call fill_halo(gz)
    ! call fill_halo(gz_lev)
    ! init_hydrostatic_gz = .true.

    end associate
  end subroutine swg_wave_test_set_ic

!=======================================================================
!    Generate the swg instability initial conditions
!=======================================================================
  SUBROUTINE swg_wave_test(deep,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho) &
    BIND(c, name = "swg_wave_test")
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
                deep,       & ! Deep (1) or Shallow (0) test case
                ! moist,      & ! Moist (1) or Dry (0) test case
                pertt         ! Perturbation type

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Earth scaling parameter

    REAL(8), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa)
                z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                thetav,     & ! Virtual potential temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constH, constC, scaledZ, inttau2, rratio
    REAL(8) :: inttermU, bigU, rcoslat, omegarcoslat
    REAL(8) :: eta, qratio, qnum, qden

    !------------------------------------------------
    !   Pressure and temperature
    !------------------------------------------------
    if (zcoords .eq. 1) then
      CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)
    else
      CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    end if

    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = a / X
    omegaref = omega * X

    ! T0 = 0.5d0 * (T0E + T0P)

    ! constH = Rd * T0 / g

    ! constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)

    ! scaledZ = z / (B * constH)

    ! inttau2 = constC * z * exp(- scaledZ**2)

    ! ! radius ratio
    ! if (deep .eq. 0) then
    !   rratio = 1.d0
    ! else
    !   rratio = (z + aref) / aref;
    ! end if

    !-----------------------------------------------------
    !   Initialize surface pressure
    !-----------------------------------------------------
    ps = p0

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    ! inttermU = (rratio * cos(lat))**(K - 1.d0) - (rratio * cos(lat))**(K + 1.d0)
    ! bigU = g / aref * K * inttau2 * inttermU * t
    ! if (deep .eq. 0) then
    !   rcoslat = aref * cos(lat)
    ! else
    !   rcoslat = (z + aref) * cos(lat)
    ! end if

    ! omegarcoslat = omegaref * rcoslat
    
    ! u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    ! u = 0.01d0 * cos(lat)

    ! if (lat == -pi/2.0 .or. lat == pi/2.0) then
    !     u = 0.0
    ! end if
    ! end do
    u = 0.d0
    v = 0.d0

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    ! Exponential type
    if (pertt .eq. 0) then
      u = u + evaluate_exponential(lon, lat, z)

    ! Stream function type
    elseif (pertt .eq. 1) then
      u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
          ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
          - evaluate_streamfunction(lon, lat - dxepsilon, z))

      v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
          ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
          - evaluate_streamfunction(lon - dxepsilon, lat, z))
    ! === 新增：激发频率为2*Omega的非传统惯性重力波 ===
    elseif (pertt .eq. 2) then
      v = v + evaluate_igw_v(lon, lat, z)
      elseif (pertt .eq. 3) then
        ! free packet 在 swg_wave_test_set_ic 里统一加入
      end if
    ! end if

    !-----------------------------------------------------
    !   Initialize surface geopotential
    !-----------------------------------------------------
    phis = 0.d0

    !-----------------------------------------------------
    !   Initialize density
    !-----------------------------------------------------
    rho = p / (Rd * t)

    !-----------------------------------------------------
    !   Initialize specific humidity
    !-----------------------------------------------------
    ! if (moist .eq. 1) then
    !   eta = p/p0

    !   if (eta .gt. moisttr) then
    !     q = moistq0 * exp(- (lat/moistqlat)**4)          &
    !                 * exp(- ((eta-1.d0)*p0/moistqp)**2)
    !   else
    !     q = moistqs
    !   end if

    !   ! Convert virtual temperature to temperature
    !   t = t / (1.d0 + Mvap * q)

    ! else
    !   q = 0.d0
    ! end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature
    !-----------------------------------------------------
    ! thetav = t * (1.d0 + 0.61d0 * q) * (p0 / p)**(Rd / cp)
    thetav = theta0 * exp(N_freq**2/g *z)
    ! print*,"theta0",theta0,lat,lon
    ! print*,"N_freq",N_freq,lat,lon
    ! print*,"z",z,lat,lon
    ! print*,"g",g,lat,lon


  END SUBROUTINE swg_wave_test

!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    REAL(8), INTENT(OUT) :: &
                p,          & ! Pressure (Pa)
                t             ! Temperature (K)

    REAL(8) :: aref, omegaref
    REAL(8) :: constA, constB, constC, constH, scaledZ
    REAL(8) :: tau1, tau2, inttau1, inttau2
    REAL(8) :: rratio, inttermT
    REAL(8) :: q
    REAL(8) :: TT, NNG

    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / X
    omegaref = omega * X

    ! T0 = 0.5d0 * (T0E + T0P)
    ! constA = 1.d0 / lapse
    ! constB = (T0 - T0P) / (T0 * T0P)
    ! constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    ! constH = Rd * T0 / g

    ! scaledZ = z / (B * constH)

    ! !--------------------------------------------
    ! !    tau values
    ! !--------------------------------------------
    ! tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
    !      + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    ! tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    ! inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
    !         + constB * z * exp(- scaledZ**2)
    ! inttau2 = constC * z * exp(- scaledZ**2)

    ! !--------------------------------------------
    ! !    radius ratio
    ! !--------------------------------------------
    ! if (deep .eq. 0) then
    !   rratio = 1.d0
    ! else
    !   rratio = (z + aref) / aref;
    ! end if

    ! !--------------------------------------------
    ! !    interior term on temperature expression
    ! !--------------------------------------------
    ! inttermT = (rratio * cos(lat))**K &
    !          - K / (K + 2.d0) * (rratio * cos(lat))**(K + 2.d0)

    TT  = g**2 / (cp * N_freq**2)
    ! print*,"tt=",TT,"g=",g,"cp=",cp,"N_freq",N_freq
    NNG = N_freq **2 /g
    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    ! t = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))
    t = TT + (T0 - TT ) * exp(z /g * N_freq**2)
    ! print*,t,z,TT,g
    
    q = - cp / Rd * (NNG * z - log(abs(TT+(T0-TT)*exp(NNG*z)))+log(T0))
    ! print*,q
    
    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    ! p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))
    p = p0 * exp(q )

  END SUBROUTINE evaluate_pressure_temperature

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    
    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    REAL(8), INTENT(OUT) :: &
                z,          & ! Altitude (m)
                t             ! Temperature (K)

    INTEGER :: ix

    REAL(8) :: z0, z1, z2
    REAL(8) :: p0, p1, p2

    z0 = 0.d0
    z1 = 10000.d0
    ! print*,"condition: p=",p
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z0, p0, t)
    ! print*,"z0=",z0,"p0=",p0,"lat=",lat
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z1, p1, t)
    ! print*,"z1=",z1,"p1=",p1,"lat=",lat

    ! stop 999


    DO ix = 1, 100
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
      ! print*,"z2=",z2

      CALL evaluate_pressure_temperature(deep, X, lon, lat, z2, p2, t)
      ! print*,"p2=",p2
    
      IF (ABS((p2 - p)/p) .lt. 1.0d-6) THEN
        EXIT
      END IF

      z0 = z1
      p0 = p1

      z1 = z2
      p1 = p2
    END DO

    z = z2

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p0, t)


  END SUBROUTINE evaluate_z_temperature

!-----------------------------------------------------------------------
!    Exponential perturbation function
!-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_exponential(lon, lat, z)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(8) :: greatcircler, perttaper

    ! Great circle distance
    greatcircler = 1.d0 / pertexpr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Zonal velocity perturbation
    if (greatcircler < 1.d0) then
      evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
      evaluate_exponential = 0.d0
    end if

  END FUNCTION evaluate_exponential

!-----------------------------------------------------------------------
!    Stream function perturbation function
!-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_streamfunction(lon, lat, z)

    REAL(8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    REAL(8) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.d0 / pertr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
      cospert = cos(0.5d0 * pi * greatcircler)
    else
      cospert = 0.d0
    end if

    evaluate_streamfunction = &
        (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction

  !-----------------------------------------------------------------------
! Inertia-Gravity Wave perturbation function (target freq = 2*Omega)
!-----------------------------------------------------------------------
  ! REAL(8) FUNCTION evaluate_igw_v(lon, lat, z)

  !   REAL(8), INTENT(IN)  :: lon, lat, z

  !   REAL(8) :: ky, mz, dist_y, dist_z, env, phase
  !   REAL(8) :: target_ratio

  !   ! 根据色散关系反解，确定垂直波数与水平波数的比值
  !   ! 使用负号来指定相位的传播倾角方向
  !   target_ratio = (sin(pertlat) - N_freq / (2.d0 * omega)) / cos(pertlat)

  !   ky = 2.d0 * pi / pert_Ly
  !   mz = ky * target_ratio

  !   dist_y = a * (lat - pertlat)
  !   dist_z = z - pertz

  !   ! 高斯包络限制扰动范围
  !   env = exp( - (dist_y**2)/(2.d0 * pert_sigma_y**2) &
  !              - (dist_z**2)/(2.d0 * pert_sigma_z**2) )

  !   ! 余弦载波
  !   phase = ky * dist_y + mz * dist_z

  !   evaluate_igw_v = pert_v0 * cos(phase) * env

  ! END FUNCTION evaluate_igw_v
  !-----------------------------------------------------------------------
! Inertia-Gravity Wave perturbation function (Dynamic Frequency & Lon-Bounded)
!-----------------------------------------------------------------------
  ! REAL(8) FUNCTION evaluate_igw_v(lon, lat, z)

  !   REAL(8), INTENT(IN)  :: lon, lat, z

  !   REAL(8) :: ky, mz, dist_x, dist_y, dist_z, env, phase
  !   REAL(8) :: A_coeff, B_coeff, C_coeff, discriminant, tan_theta
  !   REAL(8) :: f_c, f_tilde_c

  !   ! 计算波源中心的传统与非传统科氏参数
  !   f_c       = 2.d0 * omega * sin(pertlat)
  !   f_tilde_c = 2.d0 * omega * cos(pertlat)

  !   ! 通过 Gerkema (2008) 频散关系构造一元二次方程求解 tan(theta) = m/k_y
  !   ! A*x^2 + B*x + C = 0
  !   A_coeff = target_omega**2 - f_c**2
  !   B_coeff = - 2.d0 * f_c * f_tilde_c
  !   C_coeff = target_omega**2 - N_freq**2 - f_tilde_c**2

  !   discriminant = B_coeff**2 - 4.d0 * A_coeff * C_coeff

  !   ! 确保在该纬度下该频率的波可以存在
  !   if (discriminant >= 0.d0) then
  !     ! 取其中一个根决定相位倾斜方向 (向下或向上)
  !     tan_theta = (-B_coeff + sqrt(discriminant)) / (2.d0 * A_coeff)
  !   else
  !     ! 理论上超出频散允许范围，设为0以防止崩溃
  !     tan_theta = 0.d0
  !   end if

  !   ky = 2.d0 * pi / pert_Ly
  !   mz = ky * tan_theta

  !   ! 计算三维空间距离
  !   dist_x = a * cos(lat) * (lon - pertlon)  ! 经度方向距离
  !   dist_y = a * (lat - pertlat)             ! 纬度方向距离
  !   dist_z = z - pertz                       ! 垂直方向距离

  !   ! 真正的 3D 局域化高斯包络
  !   env = exp( - (dist_x**2)/(2.d0 * pert_sigma_x**2) &
  !              - (dist_y**2)/(2.d0 * pert_sigma_y**2) &
  !              - (dist_z**2)/(2.d0 * pert_sigma_z**2) )

  !   ! 余弦载波 (纯经向-垂直二维相位传播)
  !   phase = ky * dist_y + mz * dist_z

  !   evaluate_igw_v = pert_v0 * cos(phase) * env

  ! END FUNCTION evaluate_igw_v
  !-----------------------------------------------------------------------
! 极简版：纯经向自由传播波包 (Zero-Brainer Meridional Wave Packet)
!-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_igw_v(lon, lat, z)

    REAL(8), INTENT(IN)  :: lon, lat, z
    REAL(8) :: ky, mz, dist_x, dist_y, dist_z, env, phase
    REAL(8) :: center_lon, center_lat, center_z

    ! --- 1. 设定固定的波源中心 (例如: 20°E, 20°N, 15km) ---
    center_lon = 20.d0 * pi / 180.d0
    center_lat = 20.d0 * pi / 180.d0
    center_z   = 15000.d0

    ! --- 2. 极简波数设定 (完全丢掉科氏力和频散公式) ---
    ky = 2.d0 * pi / 500000.d0  ! 经向波长设为 500 km
    mz = 2.d0 * pi / 10000.d0   ! 垂直波长设为 10 km

    ! --- 3. 计算纯几何距离 ---
    dist_x = a * cos(lat) * (lon - center_lon)
    dist_y = a * (lat - center_lat)
    dist_z = z - center_z

    ! --- 4. 极简的三维高斯包络 ---
    env = exp( - (dist_x**2)/(2.d0 * 200000.d0**2) &  ! 纬向衰减
               - (dist_y**2)/(2.d0 * 200000.d0**2) &  ! 经向衰减
               - (dist_z**2)/(2.d0 * 3000.d0**2) )    ! 垂直衰减

    ! --- 5. 纯经向-垂直的二维相位 ---
    phase = ky * dist_y + mz * dist_z

    ! --- 6. 暴力赋值 (振幅 1.0 m/s) ---
    evaluate_igw_v = 1.0_r8 * cos(phase) * env

  END FUNCTION evaluate_igw_v

  subroutine add_free_igw_packet_ic(block)
    use block_mod
    use operators_mod
    use latlon_parallel_mod
    use interp_mod
  
    type(block_type), intent(inout), target :: block
  
    integer :: i, j, k
    real(r8) :: ky, mz
    real(r8) :: z_here, t_here, z_half
    real(r8) :: up, vp, wp, ptp, tp
    real(r8), allocatable :: zfull(:)
    real(r8) :: dz1
  
    associate (mesh   => block%mesh            , &
               u      => block%aux%u           , &
               v      => block%aux%v           , &
               t      => block%aux%t           , &
               pt     => block%dstate(1)%pt    , &
               mg     => block%dstate(1)%mg    , &
               w      => block%dstate(1)%w     , &
               w_lev  => block%dstate(1)%w_lev )
  
      call packet_get_wavenumbers(ky, mz)
  
      allocate(zfull(mesh%full_kds:mesh%full_kde))
  
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
  
          ! ---- full levels: u, v, pt, t ----
          do k = mesh%full_kds, mesh%full_kde
            call evaluate_z_temperature(deep, scale_X, mesh%full_lon(i), mesh%full_lat(j), &
                                        mg%d(i,j,k), z_here, t_here)
            zfull(k) = z_here
  
            call packet_polarization_full(mesh%full_lon(i), mesh%full_lat(j), z_here,        &
                                          pt%d(i,j,k), t%d(i,j,k), ky, mz, up, vp, ptp, tp)
  
            u%d(i,j,k)  = u%d(i,j,k)  + up
            v%d(i,j,k)  = v%d(i,j,k)  + vp
            pt%d(i,j,k) = pt%d(i,j,k) + ptp
            t%d(i,j,k)  = t%d(i,j,k)  + tp
            ! print*,u%d(i,j,k),"u",i,j,k
            ! print*,v%d(i,j,k),"v"
            ! print*,pt%d(i,j,k),"pt"
            ! print*,t%d(i,j,k),"t"                                    
          end do
  
          ! ---- half levels: w_lev ----
          do k = mesh%half_kds, mesh%half_kde
            if (k .eq. mesh%half_kds) then
              dz1 = zfull(mesh%full_kds+1) - zfull(mesh%full_kds)
              z_half = max(0.0_r8, zfull(mesh%full_kds) - 0.5_r8*dz1)
            elseif (k .eq. mesh%half_kde) then
              dz1 = zfull(mesh%full_kde) - zfull(mesh%full_kde-1)
              z_half = zfull(mesh%full_kde) + 0.5_r8*dz1
            else
              z_half = 0.5_r8 * ( zfull(k) + zfull(k-1) )
            end if
  
            call packet_vertical_velocity(mesh%full_lon(i), mesh%full_lat(j), z_half, ky, mz, wp)
            w_lev%d(i,j,k) = wp
          end do
  
        end do
      end do
  
      call fill_halo(u)
      call fill_halo(v)
      call fill_halo(t)
      call fill_halo(pt)
      call fill_halo(w_lev)
  
      call interp_run(w_lev, w)
      call fill_halo(w)
  
      deallocate(zfull)
    end associate
  
  end subroutine add_free_igw_packet_ic
  
  subroutine packet_get_wavenumbers(ky, mz)
    real(r8), intent(out) :: ky, mz
  
    real(r8) :: f0, ft0, Acoef, Bcoef, Ccoef, disc, root_r
    real(r8), parameter :: eps = 1.0e-14_r8
  
    f0  = 2.0_r8 * omega * sin(packet_lat0)
    ft0 = 2.0_r8 * omega * cos(packet_lat0)
  
    Acoef = packet_omega**2 - f0**2
    Bcoef = -2.0_r8 * f0 * ft0
    Ccoef = packet_omega**2 - N_freq**2 - ft0**2
  
    disc = Bcoef**2 - 4.0_r8*Acoef*Ccoef
    if (disc < 0.0_r8) then
      write(*,*) 'ERROR(packet_get_wavenumbers): discriminant < 0, packet_omega not allowed.'
      stop 911
    end if
  
    ky = 2.0_r8*pi / packet_Ly
  
    if (abs(Acoef) < eps) then
      write(*,*) 'ERROR(packet_get_wavenumbers): Acoef ~ 0, choose another packet_omega or source latitude.'
      stop 912
    end if
  
    root_r = (-Bcoef + real(packet_branch, r8) * sqrt(disc)) / (2.0_r8*Acoef)
    mz = ky * root_r
  
    write(*,'(A,1PE12.4,A,1PE12.4,A,1PE12.4)') 'packet ky=', ky, '  mz=', mz, '  m/ky=', root_r
  end subroutine packet_get_wavenumbers
  
  subroutine packet_polarization_full(lon, lat, z, pt_bg, t_bg, ky, mz, up, vp, ptp, tp)
    real(r8), intent(in)  :: lon, lat, z, pt_bg, t_bg, ky, mz
    real(r8), intent(out) :: up, vp, ptp, tp
  
    real(r8) :: env, phase, dlon, x, y, zz
    real(r8) :: f0, ft0, alpha_u, alpha_th
    real(r8) :: carrier_c, carrier_s
  
    f0  = 2.0_r8 * omega * sin(packet_lat0)
    ft0 = 2.0_r8 * omega * cos(packet_lat0)
  
    dlon = wrapped_lon_diff(lon, packet_lon0)
    x  = a * cos(packet_lat0) * dlon
    y  = a * (lat - packet_lat0)
    zz = z - packet_z0
  
    env = exp( -0.5_r8 * ( (x/packet_sigma_x)**2 + (y/packet_sigma_y)**2 + (zz/packet_sigma_z)**2 ) )
    phase = ky*y + mz*zz
  
    carrier_c = cos(phase)
    carrier_s = sin(phase)
  
    ! v' : in phase with cos(phase)
    vp = packet_v0 * env * carrier_c
  
    ! u' : 与 v' 正交 90 度，来自 -iω u - f v + f_tilde w = 0
    alpha_u = ( f0 + ft0*(ky/mz) ) / packet_omega
    up = - alpha_u * packet_v0 * env * carrier_s
  
    ! theta' : 由 b' = g theta'/theta_bg, 以及 b' = -(N^2 * ky / (ω m)) * v0 * sin(phase)
    alpha_th = (N_freq**2 * ky) / (packet_omega * mz * g)
    ptp = - pt_bg * alpha_th * packet_v0 * env * carrier_s
  
    ! 近似在固定 pressure 下令 T'/T ≈ theta'/theta
    tp = t_bg * (ptp / max(pt_bg, 1.0e-12_r8))
  end subroutine packet_polarization_full
  
  subroutine packet_vertical_velocity(lon, lat, z, ky, mz, wp)
    real(r8), intent(in)  :: lon, lat, z, ky, mz
    real(r8), intent(out) :: wp
  
    real(r8) :: env, phase, dlon, x, y, zz
  
    dlon = wrapped_lon_diff(lon, packet_lon0)
    x  = a * cos(packet_lat0) * dlon
    y  = a * (lat - packet_lat0)
    zz = z - packet_z0
  
    env = exp( -0.5_r8 * ( (x/packet_sigma_x)**2 + (y/packet_sigma_y)**2 + (zz/packet_sigma_z)**2 ) )
    phase = ky*y + mz*zz
  
    ! continuity: ky*v + mz*w = 0  (kx=0)
    wp = - (ky/mz) * packet_v0 * env * cos(phase)
  end subroutine packet_vertical_velocity
  
  real(r8) function wrapped_lon_diff(lon, lon0)
    real(r8), intent(in) :: lon, lon0
    wrapped_lon_diff = lon - lon0
    if (wrapped_lon_diff >  pi) wrapped_lon_diff = wrapped_lon_diff - 2.0_r8*pi
    if (wrapped_lon_diff < -pi) wrapped_lon_diff = wrapped_lon_diff + 2.0_r8*pi
  end function wrapped_lon_diff

end module swg_wave_test_mod
