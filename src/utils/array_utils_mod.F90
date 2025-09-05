module array_utils_mod

  use const_mod

  implicit none

  private

  public reverse_array

  interface reverse_array
    module procedure reverse_array_1d
    module procedure reverse_array_3d
  end interface

contains

  subroutine reverse_array_1d(x)

    real(r8), intent(inout) :: x(:)

    integer i, n
    real(r8), allocatable :: tmp(:)

    n = size(x)
    allocate(tmp(n))
    tmp = x
    do i = 1, n
      x(i) = tmp(n-i+1)
    end do
    deallocate(tmp)

  end subroutine reverse_array_1d

  subroutine reverse_array_3d(x, axis)

    real(r8), intent(inout) :: x(:,:,:)
    integer, intent(in) :: axis

    integer i, j, k, n1, n2, n3
    real(r8), allocatable :: tmp(:,:,:)

    n1 = size(x, 1)
    n2 = size(x, 2)
    n3 = size(x, 3)

    allocate(tmp(n1,n2,n3))
    tmp = x
    select case (axis)
    case (1)
      do i = 1, n1
        do j = 1, n2
          do k = 1, n3
            x(i,j,k) = tmp(n1-i+1,j,k)
          end do
        end do
      end do
    case (2)
      do i = 1, n1
        do j = 1, n2
          do k = 1, n3
            x(i,j,k) = tmp(i,n2-j+1,k)
          end do
        end do
      end do
    case (3)
      do i = 1, n1
        do j = 1, n2
          do k = 1, n3
            x(i,j,k) = tmp(i,j,n3-k+1)
          end do
        end do
      end do
    end select
    deallocate(tmp)

  end subroutine reverse_array_3d

end module array_utils_mod