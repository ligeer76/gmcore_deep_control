module equatorial_wave_test_mod

! =============================================================
!  Equatorial Wave Test Case for Deep Atmosphere Verification
!
!  Function: Sets up an equatorial perturbation (warm bubble)
!            to trigger Kelvin and Rossby waves.
!
!  Features: 1. Deep atmosphere support (r = a + z)
!            2. Full Coriolis terms (Traditional + Non-traditional)
!            3. Variable gravity g(r) = g0 * (a/r)^2
! =============================================================

  use mpi
  use string
  use flogger
  use math_mod
  use const_mod, only: r8, pi, omega, rd, g, cpd, a => radius
  use block_mod
  use formula_mod
  use vert_coord_mod
  use namelist_mod
  use tracer_mod
  use latlon_parallel_mod

  implicit none

  private

  public equatorial_wave_test_init
  public equatorial_wave_test_set_ic

  ! ===========================================================
  ! Test Case Parameters
  ! ===========================================================
  real(r8), parameter ::   &
    T_bg     = 240.0d0   , & ! 背景等温大气温度 (K)
    p0       = 100000.0d0, & ! 地表参考压力 (Pa)
    delta_T  = 2.0d0     , & ! 扰动振幅 (K)
    lon_c    = 180.0d0   , & ! 扰动中心经度 (deg)
    lat_c    = 0.0d0     , & ! 扰动中心纬度 (deg)
    z_c      = 10000.0d0 , & ! 扰动中心高度 (m)
    alpha_h  = 1.5d6     , & ! 水平衰减尺度 (m)
    alpha_v  = 5000.0d0      ! 垂直衰减尺度 (m)

  integer :: deep = 1        ! 1 = 深大气, 0 = 浅水近似
  ! real(r8) :: scale_X = 1.0d0 ! 地球缩放因子 (用于放大深大气效应)

  namelist /equatorial_wave_control/ deep, scale_X

contains

  subroutine equatorial_wave_test_init(namelist_path)
    character(*), intent(in) :: namelist_path
    integer :: ignore

    ! 读取参数
    open(11, file=namelist_path, status='old')
    read(11, nml=equatorial_wave_control, iostat=ignore)
    close(11)

    if (proc%is_root()) then
      call log_notice("Equatorial Wave Test Initialized")
      if (deep == 1) then
        call log_notice("Mode: Deep Atmosphere (Full Coriolis + Variable Gravity)")
      else
        call log_notice("Mode: Shallow Atmosphere (Traditional Approximation)")
      end if
    end if
  end subroutine equatorial_wave_test_init

  subroutine equatorial_wave_test_set_ic(block)
    type(block_type), intent(inout), target :: block
    integer :: i, j, k
    real(r8) :: p, t, u, v, w, rho, phis

    associate (mesh => block%mesh,             &
               lon  => block%mesh%full_lon,    &
               lat  => block%mesh%full_lat,    &
               p_d  => block%dstate(1)%p,      &
               z_d  => block%dstate(1)%gz,     &
               t_d  => block%aux%t,            &
               u_d  => block%aux%u,            &
               v_d  => block%aux%v,            &
               pt_d => block%dstate(1)%pt)

      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            
            ! 假设 z 坐标已经由框架预先设定（或在此计算）
            ! DCMIP 惯例：z_d 在此处通常存储位势高度，除以 g 得到物理高度
            call equatorial_wave_test(lon(i), lat(j), p_d%d(i,j,k), &
                                     z_d%d(i,j,k)/g, u_d%d(i,j,k),  &
                                     v_d%d(i,j,k), t_d%d(i,j,k))
            
            ! 计算位温
            pt_d%d(i,j,k) = t_d%d(i,j,k) * (p0 / p_d%d(i,j,k))**(rd/cpd)
          end do
        end do
      end do
      
      ! 此处省略 Halo 填充和 C-网格风场转换逻辑（参考 baroclinic 算例）
    end associate
  end subroutine equatorial_wave_test_set_ic

  subroutine equatorial_wave_test(lon, lat, p, z, u, v, t)
    real(r8), intent(in)  :: lon, lat, z
    real(r8), intent(out) :: p, u, v, t

    real(r8) :: aref, r, gr, dist_h, dist_v, t_pert
    real(r8) :: lon_c_rad, lat_c_rad

    ! 1. 缩放参数
    aref = a / scale_X
    lon_c_rad = lon_c * (pi/180.0d0)
    lat_c_rad = lat_c * (pi/180.0d0)

    ! 2. 确定半径 r 和重力加速度 g
    if (deep == 1) then
      r = aref + z
    else
      r = aref
    end if

    ! 3. 计算静力平衡背景场
    ! 深大气等温大气的解析压高公式: p = p0 * exp( (g0*a/Rd*T) * (a/r - 1) )
    if (deep == 1) then
      p = p0 * exp( (g * aref) / (rd * T_bg) * (aref / r - 1.0d0) )
    else
      p = p0 * exp( - (g * z) / (rd * T_bg) )
    end if

    ! 4. 设置初始风场（静止背景）
    u = 0.0d0
    v = 0.0d0

    ! 5. 添加赤道热泡扰动
    ! 计算球面大圆距离
    gr = aref * acos(sin(lat_c_rad)*sin(lat) + cos(lat_c_rad)*cos(lat)*cos(lon - lon_c_rad))
    
    dist_h = (gr / alpha_h)**2
    dist_v = ((z - z_c) / alpha_v)**2
    
    t_pert = delta_T * exp(-dist_h) * exp(-dist_v)
    t = T_bg + t_pert

  end subroutine equatorial_wave_test

end module equatorial_wave_test_mod