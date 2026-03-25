module weak_stratified_wave_test_mod

    !===============================================================================
    ! 弱层结大气亚惯性内波传播与赤道捕获测试算例
    !
    ! 本算例用于检验大气动力框架对非传统科氏力项（即 Coriolis 参数的余弦分量）的响应。
    ! 基本态为静止大气，温度廓线设计为低层（0-8 km）接近干绝热递减率（弱层结，N≈0），
    ! 高层（12 km 以上）为等温层（强层结）。初始时刻在特定位置添加局地热力扰动，
    ! 以激发亚惯性内波。在深大气（包含非传统科氏力）和浅大气（传统近似）两组实验中，
    ! 波动的经向传播行为将出现显著差异：深大气下波动可穿越惯性纬度进入极区并被捕获在弱层结层，
    ! 浅大气下波动在惯性纬度处反射，无法进入高纬。
    !
    ! 参考：Gerkema et al. (2008) Reviews of Geophysics, 46, RG2004.
    !
    ! 接口设计参考了 baroclinic_wave_test_mod 等现有模块。
    !
    ! 作者：根据用户要求生成
    !===============================================================================
    
      use const_mod, only: r8
      use flogger
      use string
      use tracer_mod
      use block_mod
      use formula_mod
      use vert_coord_mod
      use namelist_mod
      use process_mod
      use latlon_parallel_mod
      use operators_mod
      use latlon_operators_mod
      use mpi
    
      implicit none
    
      private
    
      public weak_stratified_wave_test_init
      public weak_stratified_wave_test_set_ic
    
      REAL(8), PARAMETER ::               &
      a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
      rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
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
    
  real(r8), parameter :: default_Ts    = 300.0_r8      ! 地面温度 (K)
  real(r8), parameter :: default_Tstrat= 220.0_r8      ! 平流层温度 (K)
  real(r8), parameter :: default_z1    = 8000.0_r8     ! 中性层顶高度 (m)
  real(r8), parameter :: default_z2    = 12000.0_r8    ! 等温层底高度 (m)
  real(r8), parameter :: default_gamma_d = g / cp      ! 干绝热递减率 (K/m)

  ! 热力扰动参数
  real(r8), parameter :: default_pert_dtheta = 3.0_r8  ! 扰动幅度 (K)
  real(r8), parameter :: default_pert_lonc   = 30.0_r8 ! 扰动中心经度 (度)
  real(r8), parameter :: default_pert_latc   = 30.0_r8 ! 扰动中心纬度 (度)
  real(r8), parameter :: default_pert_zc     = 4000.0_r8   ! 扰动中心高度 (m)
  real(r8), parameter :: default_pert_rh     = 500000.0_r8 ! 水平半宽 (m)
  real(r8), parameter :: default_pert_rz     = 1500.0_r8   ! 垂直半宽 (m)

  ! 模块内全局变量（由 namelist 设置）
  real(r8) :: Ts      = default_Ts
  real(r8) :: Tstrat  = default_Tstrat
  real(r8) :: z1      = default_z1
  real(r8) :: z2      = default_z2
  real(r8) :: gamma_d = default_gamma_d

  real(r8) :: pert_dtheta = default_pert_dtheta
  real(r8) :: pert_lonc   = default_pert_lonc
  real(r8) :: pert_latc   = default_pert_latc
  real(r8) :: pert_zc     = default_pert_zc
  real(r8) :: pert_rh     = default_pert_rh
  real(r8) :: pert_rz     = default_pert_rz

  ! 控制选项
  integer :: moist = 0   ! 本算例为干大气，始终为0
  integer :: deep  = 0   ! 该参数仅用于接口兼容，实际由模型配置决定
  integer :: pertt = 1   ! 是否添加热力扰动 (1:添加, 0:不添加)

  namelist /weak_stratified_wave_control/ Ts, Tstrat, z1, z2, gamma_d, &
       pert_dtheta, pert_lonc, pert_latc, pert_zc, pert_rh, pert_rz, pertt

  ! 用于迭代求解高度的容差和最大迭代次数
  real(r8), parameter :: eps_iter = 1.0e-12_r8
  integer,  parameter :: max_iter = 100

contains

  !=============================================================================
  ! 初始化子程序：读取 namelist 并注册 tracer（如有必要，此处无 tracer）
  !=============================================================================
  subroutine weak_stratified_wave_test_init(namelist_path)

    character(*), intent(in) :: namelist_path

    integer :: ignore

    open(11, file=namelist_path, status='old', action='read')
    read(11, nml=weak_stratified_wave_control, iostat=ignore)
    close(11)

    ! 本算例不使用 tracer，但若框架要求可在此注册（保持兼容）
    ! call tracer_add('dry', dt_adv, ...) 等

  end subroutine weak_stratified_wave_test_init

  !=============================================================================
  ! 设置整个块的初始条件（调用采样子程序填充每个格点）
  !=============================================================================
  subroutine weak_stratified_wave_test_set_ic(block)

    use string
    use block_mod

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) :: thetav, rho, qv
    real(r8) :: time1, time2

    time1 = MPI_Wtime()

    associate (mesh   => block%mesh            , &
               lon    => block%mesh%full_lon   , &
               lat    => block%mesh%full_lat   , &
               mgs    => block%dstate(1)%mgs   , & ! 表面干空气压力
               mg     => block%dstate(1)%mg    , & ! 各层干空气压力
               mg_lev => block%dstate(1)%mg_lev, & ! 半层干空气压力
               u      => block%aux%u           , &
               v      => block%aux%v           , &
               u_lon  => block%dstate(1)%u_lon , &
               v_lat  => block%dstate(1)%v_lat , &
               t      => block%aux%t           , &
               pt     => block%dstate(1)%pt    , &
               gz     => block%dstate(1)%gz    , &
               gz_lev => block%dstate(1)%gz_lev, &
               q      => tracers(block%id)%q   )
    ! 表面气压设为常数 p0
    mgs%d = p0
    call fill_halo(mgs)

    ! 计算各层干空气压力（根据垂直坐标定义）
    call calc_mg(block, block%dstate(1))
    print*,"1"
    ! print*,i,j
    ! print*,ptop, lon(1), lat(1), pertt, scale_X

    ! 计算各层位势高度（给定压力求高度）
    ! 注意：本算例基本态由温度廓线定义，气压与高度一一对应。
    ! 我们采用迭代方法求解高度。
    ! 先计算模型顶的高度（假定模型顶压力 ptop 已知，这里取 ptop = 100 Pa 作为示例）
    ! 实际中 ptop 应由 namelist 传入，此处简化使用固定值
    ! block%ptop = 100.0_r8
    do j = mesh%full_jds, mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        print*,i,j,lon(i),lat(j),pertt,scale_X
        gz_lev%d(i,j,1) = get_top_height(ptop, lon(i), lat(j), pertt, scale_X)
      end do
    end do
    ! 删除了原代码中错误的 fill_halo(gz_lev, 1) 调用，后续统一填充
    print*,"1.1"


    ! 对每一层，从上一层的压力/高度通过迭代得到当前层高度
    do k = mesh%half_kds+1, mesh%half_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          gz_lev%d(i,j,k) = get_height_from_p( &
               mg_lev%d(i,j,k), mg_lev%d(i,j,k-1), gz_lev%d(i,j,k-1), &
               lon(i), lat(j), pertt, scale_X)
        end do
      end do
    end do
    print*,"1.2"

    ! 底层高度应为0
    where (gz_lev%d(:,:,mesh%half_kde) /= 0) gz_lev%d(:,:,mesh%half_kde) = 0

    print*,"1.3"

    ! 全层高度取相邻半层的平均
    do k = mesh%full_kds, mesh%full_kde
      gz%d(:,:,k) = 0.5_r8 * (gz_lev%d(:,:,k) + gz_lev%d(:,:,k+1))
    end do
    print*,"2"
    ! 调用采样子程序填充各变量
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          call weak_stratified_wave_test( &
               deep  = deep,       &   ! 该参数仅用于接口兼容
               moist = moist,      &
               pertt = pertt,      &
               X     = scale_X,    &
               lon   = lon(i),     &
               lat   = lat(j),     &
               p     = mg%d(i,j,k),&
               z     = gz%d(i,j,k),&
               zcoords = 1,        &   ! 此处已知高度，求其他量
               u     = u%d(i,j,k), &
               v     = v%d(i,j,k), &
               t     = t%d(i,j,k), &
               thetav = thetav,    &
               phis  = block%static%gzs%d(i,j), &
               ps    = mgs%d(i,j), &
               rho   = rho,        &
               q     = qv)
          ! 本算例为干大气，q 保持0，但框架可能需要，故赋值
          q%d(i,j,k,1) = qv
          ! 计算位温
          pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k),q%d(i,j,k,1))
        end do
      end do
    end do
    print*,"3"

    ! 填充边界
    call fill_halo(u)
    call fill_halo(v)
    call fill_halo(t)
    call fill_halo(pt)
    gz_lev%d = gz_lev%d * g
    gz%d = gz%d * g

    call fill_halo(gz)
    call fill_halo(gz_lev)

    ! 将 A 网格风转换为 C 网格
    call wind_a2c_operator(u, v, u_lon, v_lat)
    call fill_halo(u_lon)
    call fill_halo(v_lat)

    end associate

    time2 = MPI_Wtime()
    if (proc%is_root()) call log_notice('weak_stratified_wave_test_set_ic time: ' // to_str(time2-time1,5) // ' s')

  end subroutine weak_stratified_wave_test_set_ic

  !=============================================================================
  ! 采样子程序（核心接口）
  ! 输入：经纬度、高度或压力、坐标类型标志；输出所有变量
  !=============================================================================
  subroutine weak_stratified_wave_test(deep, moist, pertt, X, lon, lat, p, z, zcoords, &
       u, v, t, thetav, phis, ps, rho, q)

    integer,  intent(in)    :: deep      ! 深/浅标志（仅用于接口兼容）
    integer,  intent(in)    :: moist     ! 湿/干标志（本算例忽略）
    integer,  intent(in)    :: pertt     ! 是否添加热力扰动 (1:添加, 0:不添加)
    real(r8), intent(in)    :: X         ! 地球缩放因子
    real(r8), intent(in)    :: lon       ! 经度 (rad)
    real(r8), intent(in)    :: lat       ! 纬度 (rad)
    real(r8), intent(inout) :: p         ! 压力 (Pa)
    real(r8), intent(inout) :: z         ! 高度 (m)
    integer,  intent(in)    :: zcoords   ! 1: 输入为高度，输出压力；0: 输入为压力，输出高度
    real(r8), intent(out)   :: u         ! 纬向风 (m/s)
    real(r8), intent(out)   :: v         ! 经向风 (m/s)
    real(r8), intent(out)   :: t         ! 温度 (K)
    real(r8), intent(out)   :: thetav    ! 虚位温 (K)
    real(r8), intent(out)   :: phis      ! 地表位势 (m^2/s^2)
    real(r8), intent(out)   :: ps        ! 地表压力 (Pa)
    real(r8), intent(out)   :: rho       ! 密度 (kg/m^3)
    real(r8), intent(out)   :: q         ! 水汽混合比 (kg/kg)

    real(r8) :: t0     ! 未扰动温度
    real(r8) :: z_local, p_local

    ! 基本态无风
    u = 0.0_r8
    v = 0.0_r8

    ! 地表位势为0
    phis = 0.0_r8
    ! 地表压力为常数 p0
    ps = p0

    ! 根据坐标类型计算温度和压力/高度
    if (zcoords == 1) then
      ! 已知高度 z，计算压力 p 和未扰动温度 t0
      call evaluate_pressure_temperature(z, p, t0, X)
      z_local = z
    else
      ! 已知压力 p，计算高度 z 和未扰动温度 t0
      call evaluate_z_temperature(p, z, t0, X)
      z_local = z
    end if

    ! 初始温度设为未扰动温度
    t = t0

    ! 如果需要添加热力扰动，则叠加
    if (pertt == 1) then
      t = t + thermal_perturbation(lon, lat, z_local, X)
    end if

    ! 计算密度
    rho = p / (rd * t)
    ! 水汽混合比恒为0
    q = 0.0_r8

    ! 计算虚位温（干空气，虚位温 = 位温）
    thetav = modified_potential_temperature(t, p,q)


  end subroutine weak_stratified_wave_test

  !=============================================================================
  ! 辅助函数：给定高度 z，计算压力 p 和未扰动温度 t0
  ! 修改：将 p 声明为 INTENT(OUT) 以匹配调用处的用法
  !=============================================================================
  recursive subroutine evaluate_pressure_temperature(z, p, t, X)

    real(r8), intent(in)  :: z
    real(r8), intent(out) :: p        ! 修改为 OUT
    real(r8), intent(out) :: t
    real(r8), intent(in)  :: X   ! 缩放因子（此处未使用，因为基本态不依赖半径）

    real(r8) :: T1, p1, b, a, exponent, zz

    ! 分段计算温度
    if (z <= z1) then
      ! 中性层（干绝热）
      t = Ts - gamma_d * z
      p = p0 * (t / Ts)**(cp / rd)
    else if (z >= z2) then
      ! 等温层
      t = Tstrat
      ! 先计算 z2 处的压力
      call evaluate_pressure_temperature(z2, p1, t, X)   ! 递归调用，但需注意避免无限递归
      ! 等温层
      p = p1 * exp( -g / (rd * Tstrat) * (z - z2) )
    else
      ! 过渡层
      t = Ts - gamma_d * z1 + (Tstrat - (Ts - gamma_d * z1)) * (z - z1) / (z2 - z1)
      ! 计算过渡层底压力 p1
      T1 = Ts - gamma_d * z1
      p1 = p0 * (T1 / Ts)**(cp / rd)
      ! 线性温度：T(z) = a + b*z，其中 a = T1 - b*z1, b = (Tstrat - T1)/(z2 - z1)
      b = (Tstrat - T1) / (z2 - z1)
      a = T1 - b * z1
      ! 积分 dp/p = -g/(rd T) dz 得 p = p1 * ( (a + b*z)/(a + b*z1) )^(-g/(rd*b))
      exponent = -g / (rd * b)
      p = p1 * ( (a + b*z) / (a + b*z1) )**exponent
    end if

  end subroutine evaluate_pressure_temperature

  !=============================================================================
  ! 辅助函数：给定压力 p，求解高度 z 和未扰动温度 t0（使用牛顿迭代）
  !=============================================================================
  subroutine evaluate_z_temperature(p, z, t, X)

    real(r8), intent(in)  :: p
    real(r8), intent(out) :: z, t
    real(r8), intent(in)  :: X

    real(r8) :: z0, z1, z2, p0, p1, p2, t0
    integer  :: iter
    real(r8) :: p_local


    ! 初始猜测区间 [0, 50000 m]
    z0 = 0.0_r8
    z1 = 50000.0_r8
    call evaluate_pressure_temperature(z0, p0, t0, X)
    call evaluate_pressure_temperature(z1, p1, t0, X)

    do iter = 1, max_iter
      ! 线性插值求新高度
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
      call evaluate_pressure_temperature(z2, p2, t0, X)
      if (abs((p2 - p)/p) < eps_iter) exit
      z0 = z1; p0 = p1
      z1 = z2; p1 = p2
    end do

    if (iter == max_iter+1) then
      call log_error('evaluate_z_temperature: failed to converge', __FILE__, __LINE__)
    end if

    z = z2
    ! call evaluate_pressure_temperature(z, p, t, X)   ! 获得温度，这里 p 是输出，但调用时 p 已为输入，需注意
    ! 实际上这里 t 已经获得，p 是输入，所以应该调用 evaluate_pressure_temperature 时传入 z 得到 t 和 p（但 p 可能被覆盖）
    ! 由于 p 是输入且不需要输出，我们只关心 t。但 evaluate_pressure_temperature 会输出 p，覆盖了输入值。
    ! 为避免覆盖，我们可以用局部变量存储，或者重新设计。这里简单重新调用，但 p 的值会被修改，不过函数结束后 p 不会再使用。
    ! 也可以直接使用前面计算的 t0，但 t0 是未扰动的？实际上 t0 在迭代中每次调用都返回对应高度的温度，但最后一次调用后 t0 就是所求温度。
    ! 所以这里调用 evaluate_pressure_temperature 时，p 被覆盖，但函数返回后 p 值可能改变，不过我们不关心，因为 p 是输入。
    ! 为安全起见，我们使用一个局部变量接收 p 的输出。
    ! real(r8) :: p_local
    call evaluate_pressure_temperature(z, p_local, t, X)

  end subroutine evaluate_z_temperature

  !=============================================================================
  ! 辅助函数：计算模型顶高度（给定模型顶压力 ptop）
  !=============================================================================
  function get_top_height(ptop, lon, lat, pertt, X) result(z_top)

    real(r8), intent(in) :: ptop, lon, lat, X
    integer,  intent(in) :: pertt
    real(r8) :: z_top

    real(r8) :: t_unused, p_unused

    ! 采用与 evaluate_z_temperature 相同的迭代，但初始猜测可设为较大值
    real(r8) :: z0, z1, z2, p0, p1, p2
    integer  :: iter

    z0 = 0.0_r8
    z1 = 50000.0_r8
    call evaluate_pressure_temperature(z0, p0, t_unused, X)
    call evaluate_pressure_temperature(z1, p1, t_unused, X)

    do iter = 1, max_iter
      z2 = z1 - (p1 - ptop) * (z1 - z0) / (p1 - p0)
      call evaluate_pressure_temperature(z2, p2, t_unused, X)
      if (abs((p2 - ptop)/ptop) < eps_iter) exit
      z0 = z1; p0 = p1
      z1 = z2; p1 = p2
    end do

    if (iter == max_iter+1) then
      call log_error('get_top_height: failed to converge', __FILE__, __LINE__)
    end if

    z_top = z2

  end function get_top_height

  !=============================================================================
  ! 辅助函数：给定压力 p，已知上一层的压力 p_abv 和高度 z_abv，求当前层高度
  ! 使用割线法从上一层的值开始迭代
  !=============================================================================
  function get_height_from_p(p, p_abv, z_abv, lon, lat, pertt, X) result(z)

    real(r8), intent(in) :: p, p_abv, z_abv, lon, lat, X
    integer,  intent(in) :: pertt
    real(r8) :: z

    real(r8) :: z0, z1, z2, p0, p1, p2, t_unused
    integer  :: iter

    z0 = z_abv
    z1 = z_abv + 1000.0_r8   ! 向上延伸1km作为初始猜测
    call evaluate_pressure_temperature(z0, p0, t_unused, X)
    call evaluate_pressure_temperature(z1, p1, t_unused, X)

    do iter = 1, max_iter
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
      call evaluate_pressure_temperature(z2, p2, t_unused, X)
      if (abs((p2 - p)/p) < eps_iter) exit
      z0 = z1; p0 = p1
      z1 = z2; p1 = p2
    end do

    if (iter == max_iter+1) then
      call log_error('get_height_from_p: failed to converge', __FILE__, __LINE__)
    end if

    z = z2

  end function get_height_from_p

  !=============================================================================
  ! 热力扰动函数
  ! 返回在给定位置 (lon, lat, z) 的温度增量 (K)
  ! 扰动形状为 cos^2 型，在距离中心小于1的区域内有效
  !=============================================================================
  function thermal_perturbation(lon, lat, z, X) result(dT)

    real(r8), intent(in) :: lon, lat, z, X
    real(r8) :: dT

    real(r8) :: aref, gr, rtheta, lon_c, lat_c, z_c

    aref = a / X
    lon_c = pert_lonc * pi / 180.0_r8   ! 转换为弧度
    lat_c = pert_latc * pi / 180.0_r8
    z_c   = pert_zc

    ! 大圆距离
    gr = aref * acos( sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon - lon_c) )

    ! 归一化距离
    rtheta = sqrt( (gr / pert_rh)**2 + ((z - z_c) / pert_rz)**2 )

    if (rtheta <= 1.0_r8) then
      dT = pert_dtheta * cos(0.5_r8 * pi * rtheta)**2
    else
      dT = 0.0_r8
    end if

  end function thermal_perturbation

end module weak_stratified_wave_test_mod