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

module gomars_v1_lsm_mod

  use gomars_v1_const_mod
  use gomars_v1_objects_mod

  implicit none

  private

  public gomars_v1_lsm_init
  public gomars_v1_lsm_final
  public sthick, sdepth

  real(r8), public, parameter :: gidn  = 0.0545_r8
  real(r8), public, parameter :: gids  = 0.0805_r8
  real(r8), public, parameter :: factl = 0.25_r8
  real(r8), public, parameter :: factm = 1.2_r8
  real(r8), public, parameter :: skind = 0.06_r8

  real(r8), allocatable :: sthick(:)
  real(r8), allocatable :: sdepth(:)

contains

  subroutine gomars_v1_lsm_init()

    integer i, k, iblk

    call gomars_v1_lsm_final()

    allocate(sthick(2*nsoil+1))
    allocate(sdepth(2*nsoil+1))

    sthick(2) = factl * skind
    do k = 4, 2 * nsoil, 2
      sthick(k) = sthick(k-2) * factm
    end do

    sdepth(1) = 0
    do k = 3, 2 * nsoil + 1, 2
      sdepth(k) = sdepth(k-2) + sthick(k-1)
    end do

    do k = 2, 2 * nsoil, 2
      sdepth(k) = sdepth(k-1) + 0.5_r8 * sthick(k)
    end do

    do k = 3, 2 * nsoil - 1, 2
      sthick(k) = 0.5_r8 * (sthick(k-1) + sthick(k+1))
    end do

    ! Initialize soil model.
    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      do i = 1, mesh%ncol
        if (state%spcflag(i)) then
          ! Southern hemisphere
          do k = 1, nsoil
            if (sdepth(2*k-1) > gids) then
              state%zin    (i,  k) = 2236.995_r8
              state%rhosoil(i,  k) = 1781.99_r8
              state%cpsoil (i,  k) = 1404.09_r8
              state%scond  (i,2*k) = state%zin(i,k)**2 / (state%rhosoil(i,k) * state%cpsoil(i,k))
            end if
          end do
          state%scond(i,1) = state%scond(i,2)
          do k = 3, 2 * nsoil - 1, 2
            state%scond(i,k) = 0.5_r8 * (state%scond(i,k-1) + state%scond(i,k+1))
          end do
        else if (state%npcflag(i)) then
          ! Northern hemisphere
          do k = 1, nsoil
            if (sdepth(2*k-1) > gidn) then
              state%zin    (i,  k) = 1100.00_r8
              state%rhosoil(i,  k) = 1781.99_r8
              state%cpsoil (i,  k) = 1404.09_r8
              state%scond  (i,2*k) = state%zin(i,k)**2 / (state%rhosoil(i,k) * state%cpsoil(i,k))
            end if
          end do
          state%scond(i,1) = state%scond(i,2)
          do k = 3, 2 * nsoil - 1, 2
            state%scond(i,k) = 0.5_r8 * (state%scond(i,k-1) + state%scond(i,k+1))
          end do
        end if

        ! Soil temperature for cold run varying from 170K to 200K.
        do k = 1, 2 * nsoil + 1
          state%stemp(i,k) = state%zavgtg(i)
        end do
      end do
      end associate
    end do

  end subroutine gomars_v1_lsm_init

  subroutine gomars_v1_lsm_final()

    if (allocated(sthick)) deallocate(sthick)
    if (allocated(sdepth)) deallocate(sdepth)

  end subroutine gomars_v1_lsm_final

end module gomars_v1_lsm_mod