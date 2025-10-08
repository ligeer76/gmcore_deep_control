! ==============================================================================
! This file is part of GoMars since 2023.
!
! GoMars is a Martian general circulation model developed in Institute of
! Atmospheric Physics (IAP), Chinese Academy of Sciences (CAS).
!
! GMCORE is a dynamical core for atmospheric model used in GoMars.
!
! GoMars and GMCORE are distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module gomars_v1_output_mod

  use fiona
  use gomars_v1_objects_mod

  implicit none

  private

  public gomars_v1_add_output
  public gomars_v1_output

contains

  subroutine gomars_v1_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    character(3) :: dims(2) = ['lon', 'lat']

    call fiona_add_var(tag, 'alsp'   , long_name='Surface albedo'           , units='', dim_names=dims, dtype=dtype)
    call fiona_add_var(tag, 'zin'    , long_name='Surface thermal inertia'  , units='', dim_names=dims, dtype=dtype)

  end subroutine gomars_v1_add_output

  subroutine gomars_v1_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk

    associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
    call fiona_output(tag, 'alsp'   , reshape(state%alsp   , mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'zin'    , reshape(state%zin    , mesh%cell_count_2d(1:2)), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    end associate

  end subroutine gomars_v1_output

end module gomars_v1_output_mod