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

module gomars_v1_utils_mod

  use vert_coord_mod
  use gomars_v1_const_mod
  use gomars_v1_types_mod

  implicit none

  private

  public interp_temperature

contains

  subroutine interp_temperature(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k
    real(r8) ptstrat, dlnp1, dlnp2, slope

    !       Pressure             Temperature
    ! ..... pstrat               tstrat
    ! ----- p_lev(1) or ptrop    t_lev(1)
    ! ===== p    (1)             t    (1)
    ! ----- p_lev(2)             t_lev(2)
    ! ===== p    (2)             t    (2)
    !   .
    !   .
    !   .
    ! ----- p_lev(nlev)          t_lev(nlev)
    ! ===== p    (nlev)          t    (nlev)
    ! ----- p_lev(nlev+1) or ps  t_lev(nlev+1)
    ! /////                      tg
    
    associate (mesh    => state%mesh   , &
               tstrat  => state%tstrat , & ! in
               ps      => state%ps     , & ! in
               p       => state%p      , & ! in
               p_lev   => state%p_lev  , & ! in
               pk      => state%pk     , & ! in
               pk_lev  => state%pk_lev , & ! in
               lnp     => state%lnp    , & ! in
               lnp_lev => state%lnp_lev, & ! in
               t       => state%t      , & ! in
               t_lev   => state%t_lev  , & ! out
               pt      => state%pt     , & ! out
               pt_lev  => state%pt_lev , & ! out
               z       => state%z      , & ! out
               dz      => state%dz     , & ! out
               z_lev   => state%z_lev  , & ! out
               dz_lev  => state%dz_lev )   ! out
    do i = 1, mesh%ncol
      ! Calculate potential temperature on stratosphere level.
      ptstrat = tstrat(i) * (ps(i) / pstrat)**rd_o_cpd

      ! Calculate potential temperature on full levels.
      do k = 1, mesh%nlev
        pt(i,k) = t(i,k) / pk(i,k)
      end do

      ! Interpolate temperature on half levels.
      do k = 2, mesh%nlev - 1
        slope = (t(i,k) - t(i,k-1)) / (lnp(i,k) - lnp(i,k-1))
        t_lev(i,k) = t(i,k-1) + slope * (lnp_lev(i,k) - lnp(i,k-1))
      end do

      ! Interpolate temperature on model top level.
      k = 1
      slope = (t(i,k) - tstrat(i)) / (lnp(i,k) - lnpstrat)
      t_lev(i,k) = tstrat(i) + slope * (lnp_lev(i,k) - lnpstrat)

      ! Extrapolate from above two model levels.
      k = mesh%nlev
      slope = (t(i,k-1) - t(i,k-2)) / (lnp(i,k-1) - lnp(i,k-2))
      t_lev(i,k) = t(i,k-1) + slope * (lnp_lev(i,k) - lnp(i,k-1))

      ! Extrapolate temperature on model bottom level.
      k = mesh%nlev + 1
      slope = (t(i,k-1) - t(i,k-2)) / (lnp(i,k-1) - lnp(i,k-2))
      t_lev(i,k) = t(i,k-1) + slope * (lnp_lev(i,k) - lnp(i,k-1))

      ! Calculate potential temperature on half levels.
      do k = 1, nlev + 1
        pt_lev(i,k) = t_lev(i,k) / pk_lev(i,k)
      end do

      ! Update geopotential height thickness.
      do k = nlev, 1, -1
        z_lev(i,k) = z_lev(i,k+1) + rd * t(i,k) * log(p_lev(i,k+1) / p_lev(i,k)) / g
      end do
      do k = 1, nlev
        z(i,k) = 0.5_r8 * (z_lev(i,k) + z_lev(i,k+1))
        dz(i,k) = z_lev(i,k) - z_lev(i,k+1)
      end do
      k = 1
      dz_lev(i,k) = z_lev(i,k) - z(i,k)
      k = nlev + 1
      dz_lev(i,k) = z(i,k-1) - z_lev(i,k)
      do k = 2, nlev
        dz_lev(i,k) = z(i,k-1) - z(i,k)
      end do
    end do
    end associate

  end subroutine interp_temperature

end module gomars_v1_utils_mod
