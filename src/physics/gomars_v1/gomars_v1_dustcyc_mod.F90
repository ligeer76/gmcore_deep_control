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

module gomars_v1_dustcyc_mod

  use gomars_v1_const_mod
  use gomars_v1_namelist_mod
  use gomars_v1_types_mod

  implicit none

  private

  public gomars_v1_dustcyc_run

contains

  subroutine gomars_v1_dustcyc_run(state)

    type(gomars_v1_state_type), intent(inout) :: state

  end subroutine gomars_v1_dustcyc_run

end module gomars_v1_dustcyc_mod