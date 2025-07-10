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
!   This module implements Laplacian damping at variable order.
!
! History:
!
!   - 20250624: Use grad_operator and div_operator to simplify codes.
!
! Authors:
!
!   - Li Dong <dongli@lasg.iap.ac.cn>
! ==============================================================================

module laplace_damp_mod

  ! ∂q       n+1    2n
  ! -- = (-1)    K ∇ q
  ! ∂t

  use const_mod
  use math_mod
  use namelist_mod
  use latlon_mesh_mod, only: global_mesh
  use latlon_field_types_mod
  use latlon_parallel_mod
  use latlon_operators_mod
  use block_mod
  use filter_mod
  use perf_mod

  implicit none

  private

  public laplace_damp_init
  public laplace_damp_final
  public laplace_damp_run
  public sponge_layer_run
  public laplace_damp

  interface laplace_damp
    module procedure laplace_damp_2d
    module procedure laplace_damp_3d
  end interface laplace_damp

  real(r8), allocatable, dimension(:), target :: sponge_layer

contains

  subroutine laplace_damp_init()

    integer k

    call laplace_damp_final()

    allocate(sponge_layer(global_mesh%full_nlev))

    do k = global_mesh%full_kds, global_mesh%full_kde
      sponge_layer(k) = exp_two_values(1.0_r8, 0.0_r8, 1.0_r8, real(sponge_layer_k0, r8), real(k, r8))
    end do

  end subroutine laplace_damp_init

  subroutine laplace_damp_final()

    if (allocated(sponge_layer)) deallocate(sponge_layer)

  end subroutine laplace_damp_final

  subroutine laplace_damp_run(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    call laplace_damp(block, dstate%u_lon, laplace_damp_order, laplace_damp_coef, block%aux%dudt_damp)
    call laplace_damp(block, dstate%v_lat, laplace_damp_order, laplace_damp_coef, block%aux%dvdt_damp)
    if (nonhydrostatic) then
      call laplace_damp(block, dstate%w_lev, laplace_damp_order, laplace_damp_coef, block%aux%dwdt_damp)
    end if

  end subroutine laplace_damp_run

  subroutine sponge_layer_run(block, dstate)

    type(block_type), intent(inout) :: block
    type(dstate_type), intent(inout) :: dstate

    call laplace_damp(block, dstate%u_lon, 2, sponge_layer_coef, block%aux%dudt_damp, use_sponge_layer=.true.)
    call laplace_damp(block, dstate%v_lat, 2, sponge_layer_coef, block%aux%dvdt_damp, use_sponge_layer=.true.)
    if (nonhydrostatic) then
      call laplace_damp(block, dstate%w_lev, 2, sponge_layer_coef, block%aux%dwdt_damp, use_sponge_layer=.true.)
    end if

  end subroutine sponge_layer_run

  subroutine laplace_damp_2d(block, f, order, coef)

    type(block_type), intent(inout) :: block
    type(latlon_field2d_type), intent(inout) :: f
    integer, intent(in) :: order
    real(r8), intent(in) :: coef

    real(r8) work(f%mesh%full_ids:f%mesh%full_ide), pole
    real(r8) c0
    integer i, j, l

    call wait_halo(f)

    call perf_start('laplace_damp_2d')

    c0 = (-1)**(order / 2 + 1) * coef

    select case (f%loc)
    case ('cell')
      associate (mesh => block%mesh     , &
                 g    => block%aux%g_2d , &
                 fx   => block%aux%fx_2d, &
                 fy   => block%aux%fy_2d, &
                 df   => block%aux%df_2d)
      call g%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        call grad_operator(g, fx, fy, with_halo=.true.)
        call div_operator(fx, fy, g)
        call fill_halo(g)
      end do
      call grad_operator(g, fx, fy, with_halo=.true.)
      if (order > 2) then
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j) = fx%d(i,j) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j) * (f%d(i+1,j) - f%d(i,j))))
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j) = fy%d(i,j) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j) * (f%d(i,j+1) - f%d(i,j))))
          end do
        end do
      end if
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          fx%d(i,j) = c0 * fx%d(i,j)
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          fy%d(i,j) = c0 * fy%d(i,j)
        end do
      end do
      call div_operator(fx, fy, df)
      call filter_run(block%small_filter, df)
      call f%add(df)
      end associate
    end select

    call fill_halo(f)

    call perf_stop('laplace_damp_2d')

  end subroutine laplace_damp_2d

  subroutine laplace_damp_3d(block, f, order, coef, df, use_sponge_layer)

    type(block_type), intent(inout) :: block
    type(latlon_field3d_type), intent(inout) :: f
    integer, intent(in) :: order
    real(r8), intent(in) :: coef
    type(latlon_field3d_type), intent(inout) :: df
    logical, intent(in), optional :: use_sponge_layer

    logical use_sponge_layer_opt
    real(r8) work(f%mesh%full_ids:f%mesh%full_ide,f%mesh%full_nlev), pole(f%mesh%full_nlev)
    real(r8) c0
    integer i, j, k, l

    use_sponge_layer_opt = .false.; if (present(use_sponge_layer)) use_sponge_layer_opt = use_sponge_layer

    call wait_halo(f)

    call perf_start('laplace_damp_3d')

    c0 = (-1)**(order / 2 + 1) * coef

    select case (f%loc)
    case ('cell')
      associate (mesh => block%mesh     , &
                 g    => block%aux%g_3d , &
                 fx   => block%aux%fx_3d, &
                 fy   => block%aux%fy_3d)
      call g%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        call grad_operator(g, fx, fy, with_halo=.true.)
        call div_operator(fx, fy, g)
        call fill_halo(g)
      end do
      call grad_operator(g, fx, fy, with_halo=.true.)
      if (order > 2) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i+1,j,k) - f%d(i,j,k))))
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j+1,k) - f%d(i,j,k))))
            end do
          end do
        end do
      end if
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = c0 * fx%d(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = c0 * fy%d(i,j,k)
          end do
        end do
      end do
      call div_operator(fx, fy, df)
      call filter_run(block%small_filter, df)
      call f%add(df)
      end associate
    case ('lon')
      associate (mesh => block%mesh         , &
                 g    => block%aux%g_3d_lon , &
                 fx   => block%aux%fx_3d_lon, &
                 fy   => block%aux%fy_3d_lon)
      call g%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        call grad_operator(g, fx, fy, with_halo=.true.)
        call div_operator(fx, fy, g)
        call fill_halo(g)
      end do
      call grad_operator(g, fx, fy, with_halo=.true.)
      if (order > 2) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%full_ids, mesh%full_ide + 1
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i,j,k) - f%d(i-1,j,k))))
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%half_ids, mesh%half_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j+1,k) - f%d(i,j,k))))
            end do
          end do
        end do
      end if
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%full_ids, mesh%full_ide + 1
            fx%d(i,j,k) = c0 * fx%d(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%half_ids, mesh%half_ide
            fy%d(i,j,k) = c0 * fy%d(i,j,k)
          end do
        end do
      end do
      call div_operator(fx, fy, df)
      call filter_run(block%small_filter, df)
      call f%add(df)
      end associate
    case ('lat')
      associate (mesh => block%mesh         , &
                 g    => block%aux%g_3d_lat , &
                 fx   => block%aux%fx_3d_lat, &
                 fy   => block%aux%fy_3d_lat)
      call g%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        call grad_operator(g, fx, fy, with_halo=.true.)
        call div_operator(fx, fy, g)
        call fill_halo(g)
      end do
      call grad_operator(g, fx, fy, with_halo=.true.)
      if (order > 2) then
        do k = mesh%full_kds, mesh%full_kde
          do j = mesh%half_jds, mesh%half_jde
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i+1,j,k) - f%d(i,j,k))))
            end do
          end do
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j,k) - f%d(i,j-1,k))))
            end do
          end do
        end do
      end if
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = c0 * fx%d(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = c0 * fy%d(i,j,k)
          end do
        end do
      end do
      call div_operator(fx, fy, df)
      call filter_run(block%small_filter, df)
      call f%add(df)
      end associate
    case ('lev')
      associate (mesh => block%mesh         , &
                 g    => block%aux%g_3d_lev , &
                 fx   => block%aux%fx_3d_lev, &
                 fy   => block%aux%fy_3d_lev)
      call g%copy(f, with_halo=.true.)
      do l = 1, (order - 2) / 2
        call grad_operator(g, fx, fy, with_halo=.true.)
        call div_operator(fx, fy, g)
        call fill_halo(g)
      end do
      call grad_operator(g, fx, fy, with_halo=.true.)
      if (order > 2) then
        do k = mesh%half_kds, mesh%half_kde
          do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
            do i = mesh%half_ids - 1, mesh%half_ide
              fx%d(i,j,k) = fx%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fx%d(i,j,k) * (f%d(i+1,j,k) - f%d(i,j,k))))
            end do
          end do
          do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
            do i = mesh%full_ids, mesh%full_ide
              fy%d(i,j,k) = fy%d(i,j,k) * max(0.0_r8, sign(1.0_r8, -fy%d(i,j,k) * (f%d(i,j+1,k) - f%d(i,j,k))))
            end do
          end do
        end do
      end if
      do k = mesh%half_kds, mesh%half_kde
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            fx%d(i,j,k) = c0 * fx%d(i,j,k)
          end do
        end do
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            fy%d(i,j,k) = c0 * fy%d(i,j,k)
          end do
        end do
      end do
      call div_operator(fx, fy, df)
      call filter_run(block%small_filter, df)
      call f%add(df)
      end associate
    end select

    call fill_halo(f, async=.true.)

    call perf_stop('laplace_damp_3d')

  end subroutine laplace_damp_3d

end module laplace_damp_mod
