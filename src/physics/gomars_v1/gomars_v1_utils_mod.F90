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
  public update_pressure

contains

  subroutine interp_temperature(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k
    real(r8) ptstrat, dlnp1, dlnp2, slope

    associate (mesh    => state%mesh   , &
               tstrat  => state%tstrat , & ! in
               ps      => state%ps     , & ! in
               pk      => state%pk     , & ! in
               pk_lev  => state%pk_lev , & ! in
               lnp     => state%lnp    , & ! in
               lnp_lev => state%lnp_lev, & ! in
               t       => state%t      , & ! in
               t_lev   => state%t_lev  , & ! out
               pt      => state%pt     , & ! out
               pt_lev  => state%pt_lev )   ! out
    do i = 1, mesh%ncol
      ! Calculate potential temperature on stratosphere level.
      ptstrat = tstrat(i) * (ps(i) / pstrat)**rd_o_cpd

      ! Calculate potential temperature on full levels.
      do k = 1, mesh%nlev
        pt(i,k) = t(i,k) / pk(i,k)
      end do

      ! Interpolate temperature on half levels.
      do k = 2, mesh%nlev - 1
        dlnp1 =            lnp_lev(i,k) - lnp(i,k-1)
        dlnp2 = lnp(i,k) - lnp_lev(i,k)
        slope = (pt(i,k) - pt(i,k-1)) / (dlnp1 + dlnp2)
        pt_lev (i,k) = pt(i,k-1) + slope * dlnp1
      end do

      ! Extrapolate temperature on model top level.
      k = 1
      dlnp1 =            lnp_lev(i,k) - lnpstrat
      dlnp2 = lnp(i,k) - lnp_lev(i,k)
      slope = (pt(i,k) - ptstrat) / (dlnp1 + dlnp2)
      pt_lev(i,k) = ptstrat + slope * dlnp1

      ! Extrapolate from above two model levels.
      k = mesh%nlev
      slope = (pt(i,k-1) - pt(i,k-2)) / (lnp(i,k-1) - lnp(i,k-2))
      pt_lev(i,k) = pt(i,k-1) + slope * (lnp_lev(i,k) - lnp(i,k-1))

      ! Extrapolate temperature on model bottom level.
      k = mesh%nlev + 1
      slope = (pt(i,k-1) - pt(i,k-2)) / (lnp(i,k-1) - lnp(i,k-2))
      pt_lev(i,k) = pt(i,k-1) + slope * (lnp_lev(i,k) - lnp(i,k-1))

      ! Calculate potential temperature on half levels.
      do k = 1, nlev + 1
        t_lev(i,k) = pt_lev(i,k) * pk_lev(i,k)
      end do
    end do
    end associate

  end subroutine interp_temperature

  subroutine update_pressure(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k

    associate (mesh   => state%mesh  , &
               ps     => state%ps    , & ! in
               p      => state%p     , & ! out
               p_lev  => state%p_lev , & ! out
               dp_dry => state%dp_dry, & ! out
               pk     => state%pk    , & ! out
               pk_lev => state%pk_lev)   ! out
    do i = 1, mesh%ncol
      do k = 1, mesh%nlev
        p (i,k) = vert_coord_calc_mg(k, ps(i))
        pk(i,k) = (p(i,k) / ps(i))**rd_o_cpd
      end do
      do k = 1, mesh%nlev + 1
        p_lev (i,k) = vert_coord_calc_mg_lev(k, ps(i))
        pk_lev(i,k) = (p_lev(i,k) / ps(i))**rd_o_cpd
      end do
      do k = 1, mesh%nlev
        dp_dry(i,k) = p_lev(i,k+1) - p_lev(i,k)
      end do
    end do
    end associate

  end subroutine update_pressure

end module gomars_v1_utils_mod