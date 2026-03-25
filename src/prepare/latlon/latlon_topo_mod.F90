! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module latlon_topo_mod

  use const_mod
  use namelist_mod
  use topo_reader_mod
  use block_mod
  use latlon_parallel_mod
  use latlon_field_types_mod
  use latlon_operators_mod
  use latlon_interp_mod
  use math_mod
  use filter_mod

  implicit none

  private

  public latlon_topo_regrid
  public latlon_topo_smooth

contains

  subroutine fill_grid(lon1, lon2, lat1, lat2, gzs, std, lnd, cnt)

    real(r8), intent(in) :: lon1
    real(r8), intent(in) :: lon2
    real(r8), intent(in) :: lat1
    real(r8), intent(in) :: lat2
    real(r8), intent(inout) :: gzs
    real(r8), intent(inout) :: std
    real(r8), intent(inout) :: lnd
    integer , intent(inout) :: cnt

    integer is, ie, js, je, i, j

    do is = 1, size(topo_lon)
      if (lon1 <= topo_lon(is)) exit
    end do
    do ie = size(topo_lon), 1, -1
      if (lon2 >= topo_lon(ie)) exit
    end do
    do js = 1, size(topo_lat)
      if (lat1 <= topo_lat(js)) exit
    end do
    do je = size(topo_lat), 1, -1
      if (lat2 >= topo_lat(je)) exit
    end do

    select case (planet)
    case ('earth')
      do j = js, je
        do i = is, ie
          if (topo_lon(i) < lon1 .or. topo_lon(i) > lon2 .or. topo_lat(j) < lat1 .or. topo_lat(j) > lat2) then
            stop 999
          end if
          if (allocated(topo_mask)) then
            if (topo_mask(i,j) == 1) then
              gzs = gzs + topo_gzs(i,j)
              std = std + topo_gzs(i,j)**2
            end if
          else
            if (topo_gzs(i,j) > 0) then
              gzs = gzs + topo_gzs(i,j)
              std = std + topo_gzs(i,j)**2
            end if
          end if
        end do
      end do
      lnd = lnd + count(topo_gzs(is:ie,js:je) > 0)
    case ('mars')
      gzs = sum(topo_gzs(is:ie,js:je))
      std = sum(topo_gzs(is:ie,js:je)**2)
      lnd = (ie - is + 1) * (je - js + 1)
    end select
    cnt = cnt + (ie - is + 1) * (je - js + 1)

  end subroutine fill_grid

  subroutine latlon_topo_regrid(block)

    type(block_type), intent(inout) :: block

    real(r8) min_lon, max_lon, min_lat, max_lat
    real(r8) lon1, lon2, lat1, lat2, pole_gzs, pole_std, pole_lnd
    integer i, j, pole_n
    integer n(block%mesh%full_ids:block%mesh%full_ide)

    associate (mesh  => block%filter_mesh    , &
               lnd   => block%static%landmask, &
               gzs   => block%static%gzs     , &
               std   => block%static%zs_std  , &
               dzsdx => block%static%dzsdx   , &
               dzsdy => block%static%dzsdy   )
    select case (topo_regrid_type)
    case ('fill')
      lnd%d = 0; gzs%d = 0; std%d = 0
      do j = mesh%full_jds, mesh%full_jde
        lat1 = mesh%half_lat_deg(j-1); lat1 = merge(lat1, -90.0_r8, lat1 /= inf)
        lat2 = mesh%half_lat_deg(j  ); lat2 = merge(lat2,  90.0_r8, lat2 /= inf)
        n = 0
        do i = mesh%full_ids, mesh%full_ide
          lon1 = mesh%half_lon_deg(i-1)
          lon2 = mesh%half_lon_deg(i  )
          call fill_grid(lon1, lon2, lat1, lat2, gzs%d(i,j), std%d(i,j), lnd%d(i,j), n(i))
          if (.not. mesh%is_pole(j)) then
            gzs%d(i,j) = gzs%d(i,j) / n(i)
            std%d(i,j) = (std%d(i,j) - 2 * gzs%d(i,j)**2 * n(i) + gzs%d(i,j)**2) / n(i) / g
            lnd%d(i,j) = lnd%d(i,j) / n(i)
          end if
        end do
        if (mesh%is_pole(j)) then
          call zonal_sum(proc%zonal_circle, gzs%d(mesh%full_ids:mesh%full_ide,j), pole_gzs)
          call zonal_sum(proc%zonal_circle, std%d(mesh%full_ids:mesh%full_ide,j), pole_std)
          call zonal_sum(proc%zonal_circle, lnd%d(mesh%full_ids:mesh%full_ide,j), pole_lnd)
          call zonal_sum(proc%zonal_circle, n, pole_n)
          gzs%d(mesh%full_ids:mesh%full_ide,j) = pole_gzs / pole_n
          std%d(mesh%full_ids:mesh%full_ide,j) = (pole_std - 2 * pole_gzs**2 * pole_n + pole_gzs**2) / pole_n / g
          lnd%d(mesh%full_ids:mesh%full_ide,j) = pole_lnd / pole_n
        end if
      end do
    case ('bilin')
      call latlon_interp_bilinear_cell(topo_lon, topo_lat, topo_gzs, mesh, gzs%d)
    end select
    call fill_halo(gzs)
    call fill_halo(std)
    call fill_halo(lnd)
    call calc_zs_slope(gzs, dzsdx, dzsdy)
    end associate

  end subroutine latlon_topo_regrid

  subroutine latlon_topo_smooth(block)

    type(block_type), intent(inout) :: block

    real(r8) tmp

    if (.not. use_zs_grad_filter .and. .not. use_zs_zonal_filter) return
    if (proc%is_root()) call log_notice('Filter topography.')

    associate (lnd    => block%static%landmask, &
               gzs    => block%static%gzs     , &
               dzsdx  => block%static%dzsdx   , &
               dzsdy  => block%static%dzsdy   )
    if (use_zs_zonal_filter) then
      call zs_zonal_filter(block, gzs, dzsdx, dzsdy)
    end if
    if (use_zs_grad_filter) then
      tmp = global_max(proc%comm_model, maxval(gzs%d / g))
      if (proc%is_root()) call log_notice('Maximum zs is ' // to_str(tmp, 'F8.1') // ' m.')
      tmp = global_max(proc%comm_model, max(dzsdx%absmax(), dzsdy%absmax()))
      if (proc%is_root()) call log_notice('Maximum topography slope angle before grad filter is ' // to_str(atan(tmp) * deg, 'F7.2') // ' deg.')
      call zs_grad_filter(block, lnd, gzs, dzsdx, dzsdy)
      tmp = global_max(proc%comm_model, maxval(gzs%d / g))
      if (proc%is_root()) call log_notice('Maximum zs is ' // to_str(tmp, 'F8.1') // ' m.')
      tmp = global_max(proc%comm_model, max(dzsdx%absmax(), dzsdy%absmax()))
      if (proc%is_root()) call log_notice('Maximum topography slope angle after grad filter is ' // to_str(atan(tmp) * deg, 'F7.2') // ' deg.')
    end if
    end associate

  end subroutine latlon_topo_smooth

  subroutine calc_zs_slope(gzs, dzsdx, dzsdy)

    type(latlon_field2d_type), intent(inout) :: gzs
    type(latlon_field2d_type), intent(inout) :: dzsdx
    type(latlon_field2d_type), intent(inout) :: dzsdy

    call grad_operator(gzs, dzsdx, dzsdy)
    dzsdx%d = dzsdx%d / g
    dzsdy%d = dzsdy%d / g
    call fill_halo(dzsdx)
    call fill_halo(dzsdy)

  end subroutine calc_zs_slope

  subroutine zs_zonal_filter(block, gzs, dzsdx, dzsdy)

    type(block_type), intent(in) :: block
    type(latlon_field2d_type), intent(inout) :: gzs
    type(latlon_field2d_type), intent(inout) :: dzsdx
    type(latlon_field2d_type), intent(inout) :: dzsdy

    integer j, cyc

    associate (mesh => block%filter_mesh, halo => block%filter_halo)
    do cyc = 1, zs_zonal_filter_cycles
      call filter_run(block%big_filter, gzs)
    end do
    call fill_halo(gzs)
    call calc_zs_slope(gzs, dzsdx, dzsdy)
    end associate

  end subroutine zs_zonal_filter

  subroutine zs_grad_filter(block, lnd, gzs, dzsdx, dzsdy) !fixme deep here grad div 

    type(block_type), intent(inout) :: block
    type(latlon_field2d_type), intent(in) :: lnd
    type(latlon_field2d_type), intent(inout) :: gzs
    type(latlon_field2d_type), intent(inout) :: dzsdx
    type(latlon_field2d_type), intent(inout) :: dzsdy

    real(r8) c0, max_slope
    integer i, j, k, cyc

    ! Do Laplacian damping on target grid with terrain slope larger than a threshold.

    associate (mesh => gzs%mesh       , &
               f    => block%aux%g_2d , &
               fx   => block%aux%fx_2d, &
               fy   => block%aux%fy_2d, &
               df   => block%aux%df_2d)
    c0 = (-1)**(zs_grad_filter_order / 2) * zs_grad_filter_coef
    do cyc = 1, zs_grad_filter_cycles
      call f%copy(gzs, with_halo=.true.)
      do k = 1, (zs_grad_filter_order - 2) / 2
        call grad_operator(f, fx, fy, with_halo=.true.)
        call div_operator(fx, fy, f)
        call fill_halo(f)
      end do
      call grad_operator(f, fx, fy, with_halo=.true.)
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          if (abs(dzsdx%d(i,j)) < topo_max_slope) then
            fx%d(i,j) = abs(dzsdx%d(i,j)) / topo_max_slope * fx%d(i,j)
          end if
          fx%d(i,j) = c0 * max(0.0_r8, min(lnd%d(i,j), lnd%d(i+1,j))) * fx%d(i,j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          if (abs(dzsdy%d(i,j)) < topo_max_slope) then
            fy%d(i,j) = abs(dzsdy%d(i,j)) / topo_max_slope * fy%d(i,j)
          end if
          fy%d(i,j) = c0 * max(0.0_r8, min(lnd%d(i,j), lnd%d(i,j+1))) * fy%d(i,j)
        end do
      end do
      if (zs_grad_filter_order > 2) then
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j) = fx%d(i,j) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j) * (gzs%d(i+1,j) - gzs%d(i,j))))
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j) = fy%d(i,j) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j) * (gzs%d(i,j+1) - gzs%d(i,j))))
          end do
        end do
      end if
      call div_operator(fx, fy, df)
      call filter_run(block%big_filter, df)
      call gzs%sub(df)
      call fill_halo(gzs)
      call calc_zs_slope(gzs, dzsdx, dzsdy)
    end do
    end associate

  end subroutine zs_grad_filter

end module latlon_topo_mod
