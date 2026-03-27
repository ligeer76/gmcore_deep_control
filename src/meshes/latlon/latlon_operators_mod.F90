! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   This module implements basic differential operators on the lat-lon meshes.
!
! History:
!
!   20240304: Initial creation.
!   20250623: Add grad_operator.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module latlon_operators_mod

  use const_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_parallel_mod
  use perf_mod

  implicit none

  private

  public divx_operator
  public divy_operator
  public divx_operator_deep
  public divy_operator_deep
  public div_operator
  public curl_operator
  public curl_operator_deep
  public grad_operator
  public wind_c2a_operator
  public wind_a2c_operator

  interface div_operator
    module procedure div_operator_2d
    module procedure div_operator_3d
    module procedure div_operator_3d_deep
  end interface div_operator

  interface grad_operator
    module procedure grad_operator_2d
    module procedure grad_operator_3d
  end interface grad_operator

contains

  subroutine divx_operator(fx, divx, with_halo)

    type(latlon_field3d_type), intent(inout) :: fx
    type(latlon_field3d_type), intent(inout) :: divx
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)

    call perf_start('divx_operator')

    associate (mesh => divx%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, divx%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, divx%loc(1:3) /= 'lev')
    select case (divx%loc)
    case ('cell', 'lev')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divx%d(i,j,k) = (fx%d(i,j,k) - fx%d(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) divx%d(:,mesh%full_jds,:) = 0
      if (mesh%has_north_pole()) divx%d(:,mesh%full_jde,:) = 0
    case ('vtx')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divx%d(i,j,k) = (fx%d(i+1,j,k) - fx%d(i,j,k)) * mesh%de_lat(j) / mesh%area_vtx(j)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('divx_operator')

  end subroutine divx_operator

  subroutine divy_operator(fy, divy, with_halo)

    type(latlon_field3d_type), intent(inout) :: fy
    type(latlon_field3d_type), intent(inout) :: divy
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    real(r8) work(divy%mesh%full_ids:divy%mesh%full_ide,divy%nlev+1)
    real(r8) pole(divy%nlev+1)
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fy)

    call perf_start('divy_operator')

    associate (mesh => divy%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, divy%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, divy%loc(1:3) /= 'lev')
    select case (divy%loc)
    case ('cell', 'lev')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divy%d(i,j,k) = (                    &
              fy%d(i,j  ,k) * mesh%le_lat(j  ) - &
              fy%d(i,j-1,k) * mesh%le_lat(j-1)   &
            ) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = is, ie
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = -pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('vtx')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divy%d(i,j,k) = (                    &
              fy%d(i,j+1,k) * mesh%de_lon(j+1) - &
              fy%d(i,j  ,k) * mesh%de_lon(j  )   &
            ) / mesh%area_vtx(j)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('divy_operator')

  end subroutine divy_operator

  subroutine divx_operator_deep(fx, rx, divx, with_halo)

    type(latlon_field3d_type), intent(inout) :: fx
    type(latlon_field3d_type), intent(inout) :: rx
    type(latlon_field3d_type), intent(inout) :: divx
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)
    call wait_halo(rx)

    call perf_start('divx_operator_deep')

    associate (mesh => divx%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, divx%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, divx%loc(1:3) /= 'lev')
    select case (divx%loc)
    case ('cell', 'lev')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divx%d(i,j,k) = (fx%d(i,j,k) * radius / rx%d(i,j,k) - fx%d(i-1,j,k) * radius / rx%d(i-1,j,k)) * &
                            mesh%le_lon(j) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) divx%d(:,mesh%full_jds,:) = 0
      if (mesh%has_north_pole()) divx%d(:,mesh%full_jde,:) = 0
    case ('vtx')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divx%d(i,j,k) = (fx%d(i+1,j,k) * radius / rx%d(i+1,j,k) - fx%d(i,j,k) * radius / rx%d(i,j,k)) * &
                            mesh%de_lat(j) / mesh%area_vtx(j)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('divx_operator_deep')

  end subroutine divx_operator_deep

  subroutine divy_operator_deep(fy, ry, divy, with_halo)

    type(latlon_field3d_type), intent(inout) :: fy
    type(latlon_field3d_type), intent(inout) :: ry
    type(latlon_field3d_type), intent(inout) :: divy
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    real(r8) work(divy%mesh%full_ids:divy%mesh%full_ide,divy%nlev+1)
    real(r8) pole(divy%nlev+1)
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fy)
    call wait_halo(ry)

    call perf_start('divy_operator_deep')

    associate (mesh => divy%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, divy%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, divy%loc(1:3) /= 'lev')
    select case (divy%loc)
    case ('cell', 'lev')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divy%d(i,j,k) = (fy%d(i,j  ,k) * radius / ry%d(i,j  ,k) * mesh%le_lat(j  ) - &
                             fy%d(i,j-1,k) * radius / ry%d(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j,k) * radius / ry%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = is, ie
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j-1,k) * radius / ry%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = -pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            divy%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('vtx')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            divy%d(i,j,k) = (fy%d(i,j+1,k) * radius / ry%d(i,j+1,k) * mesh%de_lon(j+1) - &
                             fy%d(i,j  ,k) * radius / ry%d(i,j  ,k) * mesh%de_lon(j  )) / mesh%area_vtx(j)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('divy_operator_deep')

  end subroutine divy_operator_deep

  subroutine div_operator_2d(fx, fy, div, with_halo)

    type(latlon_field2d_type), intent(inout) :: fx
    type(latlon_field2d_type), intent(inout) :: fy
    type(latlon_field2d_type), intent(inout) :: div
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    real(r8) work(div%mesh%full_ids:div%mesh%full_ide)
    real(r8) pole
    integer i, j, is, ie, js, je

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)
    call wait_halo(fy)

    call perf_start('div_operator_2d')

    associate (mesh => div%mesh)
    is = mesh%full_ids - merge(1, 0, with_halo_opt)
    ie = mesh%full_ide + merge(1, 0, with_halo_opt)
    js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
    je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
    do j = js, je
      do i = is, ie
        div%d(i,j) = ((                    &
          fx%d(i,j) - fx%d(i-1,j)          &
        ) * mesh%le_lon(j) + (             &
          fy%d(i,j  ) * mesh%le_lat(j  ) - &
          fy%d(i,j-1) * mesh%le_lat(j-1)   &
        )) / mesh%area_cell(j)
      end do
    end do
    if (mesh%has_south_pole()) then
      j = mesh%full_jds
      do i = mesh%full_ids, mesh%full_ide
        work(i) = fy%d(i,j)
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%area_pole_cap
      do i = mesh%full_ids, mesh%full_ide
        div%d(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_jde
      do i = mesh%full_ids, mesh%full_ide
        work(i) = fy%d(i,j-1)
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = -pole * mesh%le_lat(j-1) / global_mesh%area_pole_cap
      do i = mesh%full_ids, mesh%full_ide
        div%d(i,j) = pole
      end do
    end if
    end associate

    call perf_stop('div_operator_2d')

  end subroutine div_operator_2d

  subroutine div_operator_3d(fx, fy, div, with_halo)

    type(latlon_field3d_type), intent(inout) :: fx
    type(latlon_field3d_type), intent(inout) :: fy
    type(latlon_field3d_type), intent(inout) :: div
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    real(r8) work(div%mesh%full_ids:div%mesh%full_ide,div%nlev+1)
    real(r8) pole(div%nlev+1)
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)
    call wait_halo(fy)

    call perf_start('div_operator_3d')

    associate (mesh => div%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, div%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, div%loc(1:3) /= 'lev')
    select case (div%loc)
    case ('cell', 'lev')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            div%d(i,j,k) = ((                    &
              fx%d(i,j,k) - fx%d(i-1,j,k)        &
            ) * mesh%le_lon(j) + (               &
              fy%d(i,j  ,k) * mesh%le_lat(j  ) - &
              fy%d(i,j-1,k) * mesh%le_lat(j-1)   &
            )) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            div%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = -pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            div%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('lon')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            div%d(i,j,k) = ((                    &
              fx%d(i+1,j,k) - fx%d(i,j,k)        &
            ) * mesh%le_lon(j) + (               &
              fy%d(i,j  ,k) * mesh%le_lat(j  ) - &
              fy%d(i,j-1,k) * mesh%le_lat(j-1)   &
            )) / (2 * mesh%area_lon(j))
          end do
        end do
      end do
      ! No need to handle poles.
    case ('lat')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            div%d(i,j,k) = ((                    &
              fx%d(i,j,k) - fx%d(i-1,j,k)        &
            ) * mesh%de_lat(j) + (               &
              fy%d(i,j+1,k) * mesh%de_lon(j+1) - &
              fy%d(i,j  ,k) * mesh%de_lon(j)     &
            )) / (2 * mesh%area_lat(j))
          end do
        end do
      end do
    case default
      call log_error('div_operator_3d: unsupported location ' // trim(div%loc) // '!', __FILE__, __LINE__, pid=proc%id_model)
    end select
    end associate

    call perf_stop('div_operator_3d')

  end subroutine div_operator_3d

  subroutine div_operator_3d_deep(fx, rx, fy, ry, div, with_halo)

    type(latlon_field3d_type), intent(inout) :: fx
    type(latlon_field3d_type), intent(inout) :: rx
    type(latlon_field3d_type), intent(inout) :: fy
    type(latlon_field3d_type), intent(inout) :: ry
    type(latlon_field3d_type), intent(inout) :: div
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    real(r8) work(div%mesh%full_ids:div%mesh%full_ide,div%nlev+1)
    real(r8) pole(div%nlev+1)
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)
    call wait_halo(fy)

    call perf_start('div_operator_3d_deep')

    associate (mesh => div%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, div%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, div%loc(1:3) /= 'lev')
    select case (div%loc)
    case ('cell', 'lev')
      is = mesh%full_ids
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            div%d(i,j,k) = ((                    &
              fx%d(i,j,k) * radius / rx%d(i,j,k) - fx%d(i-1,j,k) * radius / rx%d(i-1,j,k)        &
            ) * mesh%le_lon(j) + (               &
              fy%d(i,j  ,k) * radius / ry%d(i,j  ,k) * mesh%le_lat(j  ) - &
              fy%d(i,j-1,k) * radius / ry%d(i,j-1,k) * mesh%le_lat(j-1)   &
            )) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_jds
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j,k) * radius / ry%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = pole(ks:ke) * mesh%le_lat(j) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            div%d(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_jde
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            work(i,k) = fy%d(i,j-1,k) * radius / ry%d(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work(:,ks:ke), pole(ks:ke))
        pole(ks:ke) = -pole(ks:ke) * mesh%le_lat(j-1) / global_mesh%area_pole_cap
        do k = ks, ke
          do i = mesh%full_ids, mesh%full_ide
            div%d(i,j,k) = pole(k)
          end do
        end do
      end if
    case ('lon')
      is = mesh%half_ids
      ie = mesh%half_ide
      js = mesh%full_jds_no_pole
      je = mesh%full_jde_no_pole
      do k = ks, ke
        do j = js, je
          do i = is, ie
            div%d(i,j,k) = ((                    &
              fx%d(i+1,j,k) * radius / rx%d(i+1,j,k) - fx%d(i,j,k) * radius / rx%d(i,j,k)     &
            ) * mesh%le_lon(j) + (               &
              fy%d(i,j  ,k) * radius / ry%d(i,j  ,k) * mesh%le_lat(j  ) - &
              fy%d(i,j-1,k) * radius / ry%d(i,j-1,k) * mesh%le_lat(j-1)   &
            )) / (2 * mesh%area_lon(j))
          end do
        end do
      end do
      ! No need to handle poles.
    case ('lat')
      is = mesh%full_ids
      ie = mesh%full_ide
      js = mesh%half_jds
      je = mesh%half_jde
      do k = ks, ke
        do j = js, je
          do i = is, ie
            div%d(i,j,k) = ((                    &
              fx%d(i,j,k) * radius / rx%d(i,j,k)- fx%d(i-1,j,k) * radius / rx%d(i-1,j,k)       &
            ) * mesh%de_lat(j) + (               &
              fy%d(i,j+1,k) * radius / ry%d(i,j+1,k) * mesh%de_lon(j+1) - &
              fy%d(i,j  ,k) * radius / ry%d(i,j  ,k) * mesh%de_lon(j)     &
            )) / (2 * mesh%area_lat(j))
          end do
        end do
      end do
    case default
      call log_error('div_operator_3d: unsupported location ' // trim(div%loc) // '!', __FILE__, __LINE__, pid=proc%id_model)
    end select
    end associate

    call perf_stop('div_operator_3d_deep')

  end subroutine div_operator_3d_deep

  subroutine curl_operator(fx, fy, curl, with_halo)

    type(latlon_field3d_type), intent(inout) :: fx
    type(latlon_field3d_type), intent(inout) :: fy
    type(latlon_field3d_type), intent(inout) :: curl
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)
    call wait_halo(fy)

    call perf_start('curl_operator')

    associate (mesh => curl%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, curl%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, curl%loc(1:3) /= 'lev')
    is = mesh%half_ids - merge(1, 0, with_halo_opt)
    ie = mesh%half_ide
    js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
    je = mesh%half_jde
    do k = ks, ke
      do j = js, je
        do i = is, ie
          curl%d(i,j,k) = (                                                     &
            fx%d(i  ,j,k) * mesh%de_lon(j) - fx%d(i,j+1,k) * mesh%de_lon(j+1) + &
            fy%d(i+1,j,k) * mesh%de_lat(j) - fy%d(i,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do
    end associate

    call perf_stop('curl_operator')

  end subroutine curl_operator

  subroutine curl_operator_deep(fx,rx, fy,ry, curl, with_halo)

    type(latlon_field3d_type), intent(inout) :: fx
    type(latlon_field3d_type), intent(inout) :: rx
    type(latlon_field3d_type), intent(inout) :: fy
    type(latlon_field3d_type), intent(inout) :: ry
    type(latlon_field3d_type), intent(inout) :: curl
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(fx)
    call wait_halo(fy)

    call perf_start('curl_operator')

    associate (mesh => curl%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, curl%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, curl%loc(1:3) /= 'lev')
    is = mesh%half_ids - merge(1, 0, with_halo_opt)
    ie = mesh%half_ide
    js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
    je = mesh%half_jde
    do k = ks, ke
      do j = js, je
        do i = is, ie
          curl%d(i,j,k) = (                                                     &
            fx%d(i  ,j,k) * mesh%de_lon(j)* radius / rx%d(i,j,k) - fx%d(i,j+1,k) * mesh%de_lon(j+1)* radius / rx%d(i,j+1,k) + &
            fy%d(i+1,j,k) * mesh%de_lat(j)* radius / ry%d(i+1,j,k) - fy%d(i,j  ,k) * mesh%de_lat(j  ) * radius / ry%d(i,j,k)  &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do
    end associate

    call perf_stop('curl_operator')

  end subroutine curl_operator_deep

  subroutine grad_operator_2d(f, gradx, grady, with_halo)

    type(latlon_field2d_type), intent(inout) :: f
    type(latlon_field2d_type), intent(inout) :: gradx
    type(latlon_field2d_type), intent(inout) :: grady
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(f)

    call perf_start('grad_operator_2d')

    associate (mesh => f%mesh)
    select case (f%loc)
    case ('cell')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do j = js, je
        do i = is, ie
          gradx%d(i,j) = (f%d(i+1,j) - f%d(i,j)) / mesh%de_lon(j)
        end do
      end do
      is = mesh%full_ids
      ie = mesh%full_ide
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde
      do j = js, je
        do i = is, ie
          grady%d(i,j) = (f%d(i,j+1) - f%d(i,j)) / mesh%de_lat(j)
        end do
      end do
    end select
    end associate

    call perf_stop('grad_operator_2d')

  end subroutine grad_operator_2d

  subroutine grad_operator_3d(f, gradx, grady, with_halo)

    type(latlon_field3d_type), intent(inout) :: f
    type(latlon_field3d_type), intent(inout) :: gradx
    type(latlon_field3d_type), intent(inout) :: grady
    logical, intent(in), optional :: with_halo

    logical with_halo_opt
    integer i, j, k, is, ie, js, je, ks, ke

    with_halo_opt = .false.; if (present(with_halo)) with_halo_opt = with_halo

    call wait_halo(f)

    call perf_start('grad_operator_3d')

    associate (mesh => f%mesh)
    ks = merge(mesh%full_kds, mesh%half_kds, f%loc(1:3) /= 'lev')
    ke = merge(mesh%full_kde, mesh%half_kde, f%loc(1:3) /= 'lev')
    select case (f%loc)
    case ('cell', 'lev')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            gradx%d(i,j,k) = (f%d(i+1,j,k) - f%d(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      is = mesh%full_ids
      ie = mesh%full_ide
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde
      do k = ks, ke
        do j = js, je
          do i = is, ie
            grady%d(i,j,k) = (f%d(i,j+1,k) - f%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    case ('lon')
      is = mesh%full_ids - merge(1, 0, with_halo_opt)
      ie = mesh%full_ide + merge(1, 0, with_halo_opt)
      js = mesh%full_jds_no_pole - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            gradx%d(i,j,k) = (f%d(i,j,k) - f%d(i-1,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      is = mesh%half_ids
      ie = mesh%half_ide
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde
      do k = ks, ke
        do j = js, je
          do i = is, ie
            grady%d(i,j,k) = (f%d(i,j+1,k) - f%d(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
    case ('lat')
      is = mesh%half_ids - merge(1, 0, with_halo_opt)
      ie = mesh%half_ide + merge(1, 0, with_halo_opt)
      js = mesh%half_jds - merge(1, 0, with_halo_opt .and. .not. mesh%has_south_pole())
      je = mesh%half_jde + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            gradx%d(i,j,k) = (f%d(i+1,j,k) - f%d(i,j,k)) / mesh%le_lat(j)
          end do
        end do
      end do
      is = mesh%full_ids
      ie = mesh%full_ide
      js = mesh%full_jds_no_pole
      je = mesh%full_jde_no_pole + merge(1, 0, with_halo_opt .and. .not. mesh%has_north_pole())
      do k = ks, ke
        do j = js, je
          do i = is, ie
            grady%d(i,j,k) = (f%d(i,j,k) - f%d(i,j-1,k)) / mesh%le_lon(j)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('grad_operator_3d')

  end subroutine grad_operator_3d

  subroutine wind_c2a_operator(u_lon, v_lat, u, v)

    type(latlon_field3d_type), intent(in   ) :: u_lon
    type(latlon_field3d_type), intent(in   ) :: v_lat
    type(latlon_field3d_type), intent(inout) :: u
    type(latlon_field3d_type), intent(inout) :: v

    integer i, j, k

    associate (mesh => u_lon%mesh)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%full_ids, mesh%full_ide
          u%d(i,j,k) = 0.5_r8 * (u_lon%d(i,j,k) + u_lon%d(i-1,j,k))
          v%d(i,j,k) = 0.5_r8 * (v_lat%d(i,j,k) + v_lat%d(i,j-1,k))
        end do
      end do
    end do
    end associate

  end subroutine wind_c2a_operator

  subroutine wind_a2c_operator(u, v, u_lon, v_lat)

    type(latlon_field3d_type), intent(in   ) :: u
    type(latlon_field3d_type), intent(in   ) :: v
    type(latlon_field3d_type), intent(inout) :: u_lon
    type(latlon_field3d_type), intent(inout) :: v_lat

    integer i, j, k

    associate (mesh => u%mesh)
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids, mesh%half_ide
          u_lon%d(i,j,k) = 0.5_r8 * (u%d(i,j,k) + u%d(i+1,j,k))
        end do
      end do
      do j = mesh%half_jds, mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          v_lat%d(i,j,k) = 0.5_r8 * (v%d(i,j,k) + v%d(i,j+1,k))
        end do
      end do
    end do
    end associate

  end subroutine wind_a2c_operator

end module latlon_operators_mod
