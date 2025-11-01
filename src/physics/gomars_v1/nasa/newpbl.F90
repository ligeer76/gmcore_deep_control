subroutine newpbl( &
  z0             , &
  tg             , &
  qrad           , &
  ps             , &
  ts             , &
  npcflag        , &
  rho            , &
  u              , &
  v              , &
  pt             , &
  pt_lev         , &
  q              , &
  z              , &
  dz             , &
  z_lev          , &
  dz_lev         , &
  shr2           , &
  ri             , &
  km             , &
  kh             , &
  ustar          , &
  tstar          , &
  taux           , &
  tauy           , &
  ht_pbl         , &
  rhouch         , &
  tm_sfc         , &
  h2osub_sfc     )

  ! Legacy Mars GCM v24
  ! Mars Climate Modeling Center
  ! NASA Ames Research Center

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_pbl_mod

  implicit none

  real(r8), intent(in   ) :: z0                       ! Roughness length (m)
  real(r8), intent(in   ) :: tg                       ! Ground temperature (K)
  real(r8), intent(in   ) :: qrad   (0:nlev  )        ! Radiative heating rate on full levels including TOA (K s-1)
  real(r8), intent(in   ) :: ps                       ! Surface pressure (Pa)
  real(r8), intent(in   ) :: ts                       ! Surface air temperature (K)
  logical , intent(in   ) :: npcflag
  real(r8), intent(inout) :: rho    (  nlev  )
  real(r8), intent(inout) :: u      (  nlev  )
  real(r8), intent(inout) :: v      (  nlev  )
  real(r8), intent(inout) :: pt     (  nlev  )
  real(r8), intent(in   ) :: pt_lev (  nlev+1)
  real(r8), intent(inout) :: q      (  nlev,ntracers)
  real(r8), intent(in   ) :: z      (  nlev  )
  real(r8), intent(inout) :: dz     (  nlev  )
  real(r8), intent(in   ) :: z_lev  (  nlev+1)
  real(r8), intent(inout) :: dz_lev (  nlev+1)
  real(r8), intent(inout) :: shr2   (  nlev+1)
  real(r8), intent(  out) :: ri     (  nlev+1)
  real(r8), intent(  out) :: km     (  nlev+1)
  real(r8), intent(  out) :: kh     (  nlev+1)
  real(r8), intent(  out) :: ustar
  real(r8), intent(  out) :: tstar
  real(r8), intent(  out) :: taux                     ! Wind stress in x direction (N m-2)
  real(r8), intent(  out) :: tauy                     ! Wind stress in y direction (N m-2)
  real(r8), intent(  out) :: ht_pbl                   ! Heat rate at the surface (K s-1)
  real(r8), intent(  out) :: rhouch
  real(r8), intent(inout) :: tm_sfc(ntracers)         ! Tracer mass on the ground (kg)
  real(r8), intent(inout) :: h2osub_sfc               ! Water ice upward sublimation flux at the surface (kg m-2 s-1)

  integer i, k
  real(r8) h    (nlev  ) ! Height of full levels from ground
  real(r8) h_lev(nlev+1) ! Height of half levels from ground
  real(r8) rhos
  real(r8) cdm
  real(r8) cdh
  real(r8) alpha
  real(r8) qsat
  real(r8) coef
  real(r8) kdf(nlev+1,4)
  real(r8) var(nlev  ,4)
  real(r8) bnd(       4)
  real(r8) rhs(nlev  ,4)

  ! rho = 1

  h = z - z_lev(nlev+1)
  h_lev = z_lev - z_lev(nlev+1)

  call eddycoef(h_lev, dz_lev, u, v, pt, pt_lev, q, shr2, ri, km, kh)
  call bndcond(u, v, pt, tg, h, z0, cdm, cdh, ustar, tstar)

  ! Diagnose the surface stress and heat fluxes.
  rhos   =  dry_air_density(ts, ps)
  alpha  =  atan2(v(nlev), u(nlev))
  taux   =  rhos * ustar**2 * cos(alpha)
  tauy   =  rhos * ustar**2 * sin(alpha)
  ht_pbl = -rhos * cpd * ustar * tstar
  rhouch =  rhos * cpd * ustar * cdh

  ! Reduce the sublimation flux by a coefficient. It is a tunable parameter to
  ! avoid the formation of low-lying clouds in summer above the north permanent cap.
  qsat = water_vapor_saturation_mixing_ratio_mars(tg, ps)
  coef = merge(1.0_r8, 1.0_r8, qsat > q(nlev,iMa_vap))

  i = 0
  ! Zonal wind
  i = i + 1
  kdf(:   ,i) = km
  bnd(     i) = dt / dz(nlev) * sqrt(cdm) * ustar
  rhs(:   ,i) = 0
  var(:   ,i) = rho * u
  ! Meridional wind
  i = i + 1
  kdf(:   ,i) = km
  bnd(     i) = dt / dz(nlev) * sqrt(cdm) * ustar
  rhs(:   ,i) = 0
  var(:   ,i) = rho * v
  ! Potential temperature
  i = i + 1
  kdf(:   ,i) = kh
  bnd(     i) = dt / dz(nlev) * cdh * ustar
  rhs(:   ,i) = rho * qrad(1:nlev) * dt
  rhs(nlev,i) = rhs(nlev,i) + dt / dz(nlev) * rho(nlev) * cdh * ustar * tg
  var(:   ,i) = rho * pt
  ! Water vapor
  i = i + 1
  kdf(:   ,i) = kh
  bnd(     i) = dt / dz(nlev) * ustar * cdh * coef
  rhs(:   ,i) = 0
  rhs(nlev,i) = h2osub_sfc + dt / dz(nlev) * rho(nlev) * ustar * cdh * coef
  var(:   ,i) = rho * q(:,iMa_vap)

  do i = 1, 4
    call pbl_solve(kdf(:,i), bnd(i), rhs(:,i), var(:,i))
  end do

  u            = var(:,1) / rho
  v            = var(:,2) / rho
  pt           = var(:,3) / rho
  q(:,iMa_vap) = var(:,4) / rho

  ! Update water ice budget on the surface.
  tm_sfc(iMa_vap) = tm_sfc(iMa_vap) - h2osub_sfc
  if (.not. npcflag .and. tm_sfc(iMa_vap) < 0) then
    tm_sfc(iMa_vap) = 0
  end if

contains

  subroutine pbl_solve(kdf, bnd, rhs, var)

    use math_mod

    real(r8), intent(in   ) :: kdf(nlev+1)
    real(r8), intent(in   ) :: bnd         ! Lower boundary condition for coefficient c
    real(r8), intent(in   ) :: rhs(nlev  ) ! Right hand side
    real(r8), intent(inout) :: var(nlev  )

    integer k
    real(r8) a(nlev), b(nlev), c(nlev), d(nlev)

    ! Setup the tridiagonal matrix.
    a(1) = 0
    do k = 2, nlev
      a(k) = -dt * kdf(k  ) / dz(k) / dz_lev(k  )
    end do
    do k = 1, nlev - 1
      c(k) = -dt * kdf(k+1) / dz(k) / dz_lev(k+1)
    end do
    c(nlev) = 0
    do k = 1, nlev - 1
      b(k) = 1 - a(k) - c(k)
    end do
    b(nlev) = 1 - a(nlev) + bnd
    d = var + rhs

    call tridiag_thomas(a, b, c, d, var)

  end subroutine pbl_solve

end subroutine newpbl
