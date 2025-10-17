subroutine potemp( &
  tstrat         , &
  p              , &
  p_lev          , &
  lnp            , &
  lnp_lev        , &
  t              , &
  t_lev          , &
  pt             , &
  pt_lev         )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use gomars_v1_const_mod

  implicit none

  real(r8), intent(in ) :: tstrat
  real(r8), intent(in ) :: p      (  nlev  )
  real(r8), intent(in ) :: p_lev  (  nlev+1)
  real(r8), intent(in ) :: lnp    (  nlev  )
  real(r8), intent(in ) :: lnp_lev(  nlev+1)
  real(r8), intent(in ) :: t      (  nlev  )
  real(r8), intent(out) :: t_lev  (  nlev+1)
  real(r8), intent(out) :: pt     (  nlev  )
  real(r8), intent(out) :: pt_lev (  nlev+1)

  integer k
  real(r8) ps, dlnp1, dlnp2, slope, tsl

  ps = p_lev(nlev+1)

  ! Interpolate temperature on half levels.
  do k = 2, nlev
    dlnp1 =          lnp_lev(k) - lnp(k-1)
    dlnp2 = lnp(k) - lnp_lev(k)
    slope = (t(k) - t(k-1)) / (dlnp1 + dlnp2)
    t_lev (k) = t(k-1) + slope * dlnp1
  end do

  ! Extrapolate temperature on model top level.
  k = 1
  dlnp1 =          lnp_lev(k) - lnpstrat
  dlnp2 = lnp(k) - lnp_lev(k)
  slope = (t(k) - tstrat) / (dlnp1 + dlnp2)
  t_lev(k) = tstrat + slope * dlnp1

  ! Extrapolate temperature on model bottom level.
  k = nlev + 1
  dlnp1 =            lnp_lev(k-1) - lnp(k-2)
  dlnp2 = lnp(k-1) - lnp_lev(k-1)
  slope = (t(k-1) - t(k-2)) / (dlnp1 + dlnp2)
  t_lev(k) = t(k-1) + slope * (lnp_lev(k) - lnp(k-1))

  ! Calculate potential temperature on full levels.
  do k = 1, nlev
    pt(k) = t(k) * (ps / p(k))**rd_o_cpd
  end do
  ! Calculate potential temperature on half levels.
  do k = 1, nlev + 1
    pt_lev(k) = t_lev(k) * (ps / p_lev(k))**rd_o_cpd
  end do

end subroutine potemp