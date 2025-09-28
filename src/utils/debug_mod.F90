module debug_mod

  use, intrinsic :: ieee_arithmetic
  use flogger
  use const_mod
  use namelist_mod
  use latlon_mesh_mod
  use block_mod
  use latlon_parallel_mod
  use latlon_field_types_mod

  private

  public print_min_max
  public is_inf
  public is_greater_than
  public is_less_than

  interface print_min_max
    module procedure print_min_max_3d
  end interface print_min_max

  interface is_inf
    module procedure is_inf_r4
    module procedure is_inf_r8
  end interface is_inf

contains

  subroutine print_min_max_3d(field)

    type(latlon_field3d_type), intent(in) :: field

    integer is, ie, js, je, ks, ke
    real(r8) min_val, max_val

    is = merge(field%mesh%full_ids, field%mesh%half_ids, field%full_lon)
    ie = merge(field%mesh%full_ide, field%mesh%half_ide, field%full_lon)
    js = merge(field%mesh%full_jds, field%mesh%half_jds, field%full_lat)
    je = merge(field%mesh%full_jde, field%mesh%half_jde, field%full_lat)
    ks = merge(field%mesh%full_kds, field%mesh%half_kds, field%full_lev)
    ke = merge(field%mesh%full_kde, field%mesh%half_kde, field%full_lev)

    min_val = global_min(proc%comm_model, minval(field%d(is:ie,js:je,ks:ke)))
    max_val = global_max(proc%comm_model, maxval(field%d(is:ie,js:je,ks:ke)))

    if (proc%is_root()) write(6, *) trim(field%name), min_val, max_val

  end subroutine print_min_max_3d

  logical function is_inf_r4(x) result(res)

    real(4), intent(in) :: x

    res = .not. ieee_is_finite(x)

  end function is_inf_r4

  logical function is_inf_r8(x) result(res)

    real(8), intent(in) :: x

    res = .not. ieee_is_finite(x)

  end function is_inf_r8

  logical function is_greater_than(x, y) result(res)

    real(8), intent(in) :: x, y

    res = (x - y) / (x + eps) > 1.0d-10

  end function is_greater_than

  logical function is_less_than(x, y) result(res)

    real(8), intent(in) :: x, y

    res = (y - x) / (x + eps) > 1.0d-10

  end function is_less_than

end module debug_mod
