! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================

module mars_vert_coord_mod

  use flogger
  use const_mod
  use latlon_mesh_mod
  use latlon_parallel_types_mod

  implicit none

  private

  public mars_vert_coord_emars28
  public mars_vert_coord_nasa24

contains

  subroutine mars_vert_coord_emars28(p0, ptop, hyai, hybi)

    real(r8), intent(out) :: p0
    real(r8), intent(out) :: ptop
    real(r8), intent(out) :: hyai(29)
    real(r8), intent(out) :: hybi(29)

    if (global_mesh%full_nlev /= 28 .and. proc%is_root()) then
      call log_error('nlev should be 28 in namelist!')
    end if

    hyai = [          &
      0.02_r8       , & !  1
      0.05738127_r8 , & !  2
      0.1958398_r8  , & !  3
      0.5922958_r8  , & !  4
      1.566023_r8   , & !  5
      2.445497_r8   , & !  6
      2.768375_r8   , & !  7
      2.885169_r8   , & !  8
      2.917223_r8   , & !  9
      2.908704_r8   , & ! 10
      2.859894_r8   , & ! 11
      2.768765_r8   , & ! 12
      2.632705_r8   , & ! 13
      2.450922_r8   , & ! 14
      2.226681_r8   , & ! 15
      1.968468_r8   , & ! 16
      1.689483_r8   , & ! 17
      1.405581_r8   , & ! 18
      1.132426_r8   , & ! 19
      0.8828918_r8  , & ! 20
      0.6654847_r8  , & ! 21
      0.4840102_r8  , & ! 22
      0.3382412_r8  , & ! 23
      0.225107_r8   , & ! 24
      0.1399572_r8  , & ! 25
      0.07761155_r8 , & ! 26
      0.0330855_r8  , & ! 27
      0.002_r8      , & ! 28
      0.0_r8          & ! 29
    ]
    hybi = [          &
      0.0_r8        , & !  1
      0.0_r8        , & !  2
      0.0_r8        , & !  3
      0.0_r8        , & !  4
      0.0_r8        , & !  5
      0.001936639_r8, & !  6
      0.007441913_r8, & !  7
      0.01622727_r8 , & !  8
      0.02707519_r8 , & !  9
      0.043641_r8   , & ! 10
      0.0681068_r8  , & ! 11
      0.1028024_r8  , & ! 12
      0.1497195_r8  , & ! 13
      0.2098713_r8  , & ! 14
      0.2827023_r8  , & ! 15
      0.3658161_r8  , & ! 16
      0.4552023_r8  , & ! 17
      0.545936_r8   , & ! 18
      0.6331097_r8  , & ! 19
      0.7126763_r8  , & ! 20
      0.7819615_r8  , & ! 21
      0.8397753_r8  , & ! 22
      0.8862035_r8  , & ! 23
      0.9222317_r8  , & ! 24
      0.9493454_r8  , & ! 25
      0.9691962_r8  , & ! 26
      0.9833726_r8  , & ! 27
      0.9932694_r8  , & ! 28
      1.0_r8          & ! 29
    ]

    p0 = 701
    hyai = hyai / p0
    ptop = p0 * hyai(1)

  end subroutine mars_vert_coord_emars28

  subroutine mars_vert_coord_nasa24(ptop, sigi, sig)

    real(r8), intent(out) :: ptop
    real(r8), intent(out) :: sigi(25)
    real(r8), intent(out) :: sig(24)

    sigi = [                                                    &
      0.0, 0.0001237, 0.0003423, 0.00073, 0.0014177, 0.0026376, &
      0.0047998, 0.008657, 0.0143857, 0.0238144, 0.0393587,     &
      0.0649849, 0.1072458, 0.1768876, 0.2917842, 0.4448475,    &
      0.6020744, 0.7372935, 0.8585347, 0.9304322, 0.974,        &
      0.988, 0.996, 0.999, 1.0                                  &
    ]
    sig = [                                                     &
      0.00006185, 0.000233, 0.00053615, 0.00107385, 0.00202765, &
      0.0037187, 0.0067284, 0.01152135, 0.01910005, 0.03158655, &
      0.0521718, 0.08611535, 0.1420667, 0.2343359, 0.36831585,  &
      0.52346095, 0.66968395, 0.7979141, 0.89448345, 0.9522161, &
      0.981, 0.992, 0.9975, 0.9995                              &
    ]
    ptop = 0.08 ! Pa

  end subroutine mars_vert_coord_nasa24

end module mars_vert_coord_mod
