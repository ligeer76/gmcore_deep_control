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

module gomars_v1_damp_mod

  use vert_coord_mod
  use gomars_v1_const_mod
  use gomars_v1_types_mod

  implicit none

  private

  public gomars_v1_damp_init
  public gomars_v1_damp_run
  public gomars_v1_damp_final

  ! Number of top levels for Rayleigh friction
  integer , parameter :: lray       = 3
  integer , parameter :: k1         = 4
  real(r8), parameter :: alfray     = 2
  ! Relaxation time in days
  real(r8), parameter :: trefr      = 0.5_r8

  real(r8), allocatable :: rayk(:)

contains

  subroutine gomars_v1_damp_init()

    integer k
    real(r8) pl(nlev)

    call gomars_v1_damp_final()

    ! Set a reference pressure profile.
    do k = 1, nlev
      pl(k) = vert_coord_calc_mg(k, psl)
    end do

    allocate(rayk(nlev))

    rayk = 0

    do k = 1, lray
      rayk(k) = (pl(k1) / pl(k))**alfray / (trefr * mars_sol_seconds)
    end do

  end subroutine gomars_v1_damp_init

  subroutine gomars_v1_damp_final()

    if (allocated(rayk)) deallocate(rayk)

  end subroutine gomars_v1_damp_final

  subroutine gomars_v1_damp_run(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k
    real(r8) rkei, rkef

    associate (mesh => state%mesh, &
               u    => state%u   , & ! inout
               v    => state%v   , & ! inout
               t    => state%t   )   ! inout
    do i = 1, mesh%ncol
      do k = 1, mesh%nlev
        if (rayk(k) > 0.0_r8) then
          rkei = 0.5_r8 * (u(i,k)**2 + v(i,k)**2)
          u(i,k) = u(i,k) - rayk(k) * u(i,k) * dt
          v(i,k) = v(i,k) - rayk(k) * v(i,k) * dt
          rkef = 0.5_r8 * (u(i,k)**2 + v(i,k)**2)
          t(i,k) = t(i,k) + (rkei - rkef) / cpd
        end if
      end do      
    end do
    end associate

  end subroutine gomars_v1_damp_run

end module gomars_v1_damp_mod