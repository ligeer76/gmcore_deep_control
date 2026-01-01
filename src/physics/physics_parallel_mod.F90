! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module physics_parallel_mod

  use mpi
  use const_mod
  use process_mod

  implicit none

  private

  public physics_pole_sum
  public physics_sum

  interface physics_pole_sum
    module procedure physics_pole_sum_0d
    module procedure physics_pole_sum_1d
  end interface physics_pole_sum

  interface physics_sum
    module procedure physics_sum_0d_r8
    module procedure physics_sum_1d_r8
  end interface physics_sum

contains

  subroutine physics_pole_sum_0d(lat, x)

    real(r8), intent(in   ) :: lat(:)
    real(r8), intent(inout) :: x  (:)

    integer icol, nlon, ierr
    integer n_south, i_south(size(x))
    integer n_north, i_north(size(x))
    real(r8) sum_south, val_south
    real(r8) sum_north, val_north

    n_south = 0
    n_north = 0
    do icol = 1, size(x)
      if (lat(icol) < 0 .and. abs(lat(icol) + pi05) < 1.0e-12) then
        n_south = n_south + 1
        i_south(n_south) = icol
      else if (lat(icol) > 0 .and. abs(lat(icol) - pi05) < 1.0e-12) then
        n_north = n_north + 1
        i_north(n_north) = icol
      end if
    end do

    sum_south = 0
    if (n_south > 0) then
      do icol = 1, n_south
        sum_south = sum_south + x(i_south(icol))
      end do
    end if
    call MPI_ALLREDUCE(sum_south, val_south, 1, MPI_DOUBLE, MPI_SUM, proc%comm_model, ierr)
    call MPI_ALLREDUCE(n_south, nlon, 1, MPI_INT, MPI_SUM, proc%comm_model, ierr)
    val_south = val_south / nlon
    sum_north = 0
    if (n_north > 0) then
      do icol = 1, n_north
        sum_north = sum_north + x(i_north(icol))
      end do
    end if
    call MPI_ALLREDUCE(sum_north, val_north, 1, MPI_DOUBLE, MPI_SUM, proc%comm_model, ierr)
    val_north = val_north / nlon

    if (n_south > 0) then
      do icol = 1, n_south
        x(i_south(icol)) = val_south
      end do
    end if
    if (n_north > 0) then
      do icol = 1, n_north
        x(i_north(icol)) = val_north
      end do
    end if

  end subroutine physics_pole_sum_0d

  subroutine physics_pole_sum_1d(lat, x)

    real(r8), intent(in   ) :: lat(:  )
    real(r8), intent(inout) :: x  (:,:)
    real(r8) res(size(x, 2))

    integer icol, nlon, ierr
    integer n_south, i_south(size(x, 1))
    integer n_north, i_north(size(x, 1))
    real(r8) sum_south(size(x, 2)), val_south(size(x, 2))
    real(r8) sum_north(size(x, 2)), val_north(size(x, 2))

    n_south = 0
    n_north = 0
    do icol = 1, size(x, 1)
      if (lat(icol) < 0 .and. abs(lat(icol) + pi05) < 1.0e-12) then
        n_south = n_south + 1
        i_south(n_south) = icol
      else if (lat(icol) > 0 .and. abs(lat(icol) - pi05) < 1.0e-12) then
        n_north = n_north + 1
        i_north(n_north) = icol
      end if
    end do

    sum_south(:) = 0
    if (n_south > 0) then
      do icol = 1, n_south
        sum_south(:) = sum_south(:) + x(i_south(icol),:)
      end do
    end if
    call MPI_ALLREDUCE(sum_south, val_south, size(x, 2), MPI_DOUBLE, MPI_SUM, proc%comm_model, ierr)
    call MPI_ALLREDUCE(n_south, nlon, 1, MPI_INT, MPI_SUM, proc%comm_model, ierr)
    val_south(:) = val_south(:) / nlon
    sum_north(:) = 0
    if (n_north > 0) then
      do icol = 1, n_north
        sum_north(:) = sum_north(:) + x(i_north(icol),:)
      end do
    end if
    call MPI_ALLREDUCE(sum_north, val_north, size(x, 2), MPI_DOUBLE, MPI_SUM, proc%comm_model, ierr)
    val_north(:) = val_north(:) / nlon

    if (n_south > 0) then
      do icol = 1, n_south
        x(i_south(icol),:) = val_south(:)
      end do
    end if
    if (n_north > 0) then
      do icol = 1, n_north
        x(i_north(icol),:) = val_north(:)
      end do
    end if

  end subroutine physics_pole_sum_1d

  real(8) function physics_sum_0d_r8(x) result(res)

    real(8), intent(in) :: x

    integer ierr

    call MPI_ALLREDUCE(x, res, 1, MPI_DOUBLE, MPI_SUM, proc%comm_model, ierr)

  end function physics_sum_0d_r8

  real(8) function physics_sum_1d_r8(x)

    real(8), intent(in) :: x(:)
    real(8) res(size(x))

    integer k, ierr

    call MPI_ALLREDUCE(x, res, size(x), MPI_DOUBLE, MPI_SUM, proc%comm_model, ierr)

  end function physics_sum_1d_r8

end module physics_parallel_mod
