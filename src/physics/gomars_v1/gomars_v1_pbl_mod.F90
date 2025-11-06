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

module gomars_v1_pbl_mod

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_types_mod
  use gomars_v1_tracers_mod

  implicit none

  private

  public gomars_v1_pbl_init
  public gomars_v1_pbl_final
  public gomars_v1_pbl_run

  real(r8), parameter :: z00    = 0.01_r8
  ! Maximum mixing length (m)
  real(r8), parameter :: ml0    = 150.0_r8
  ! Critical Richardson number above which turbulence is suppressed
  real(r8), parameter :: ric    = 0.195_r8
  ! Mixing parameters
  real(r8), parameter :: Sm     = 0.393_r8
  real(r8), parameter :: Sh     = 0.493_r8
  real(r8), parameter :: sqrtGM = sqrt(0.153_r8)

contains

  subroutine gomars_v1_pbl_init()

  end subroutine gomars_v1_pbl_init

  subroutine gomars_v1_pbl_final()

  end subroutine gomars_v1_pbl_final

  subroutine gomars_v1_pbl_run(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k
    real(r8) h(nlev), h_lev(nlev+1)
    real(r8) ml2, dudz, dvdz, dptdz, shr
    real(r8) km0, kh0, kmin
    real(r8) lnzz, dpt, wsp, rib, fm, fh
    real(r8) rhos, alpha, qsat, coef
    real(r8) kdf(nlev), rhs(nlev), bnd, var(nlev)

    associate (mesh       => state%mesh      , &
               zs         => state%zs        , & ! in
               ps         => state%ps        , & ! in
               ts         => state%ts        , & ! in
               tg         => state%tg        , & ! in
               npcflag    => state%npcflag   , & ! in
               t          => state%t         , & ! inout
               pt         => state%pt        , & ! inout
               pt_lev     => state%pt_lev    , & ! in
               rho        => state%rhod      , & ! in
               p          => state%p         , & ! in
               pk         => state%pk        , & ! in
               z          => state%z         , & ! in
               dz         => state%dz        , & ! in
               z_lev      => state%z_lev     , & ! in
               dz_lev     => state%dz_lev    , & ! in
               z0         => state%z0        , & ! in
               u          => state%u         , & ! inout
               v          => state%v         , & ! inout
               q          => state%q         , & ! inout
               qrad       => state%qrad      , & ! in
               co2ice_sfc => state%co2ice_sfc, & ! in
               h2osub_sfc => state%h2osub_sfc, & ! in
               shr2       => state%shr2      , & ! inout
               ri         => state%ri        , & ! out
               km         => state%km        , & ! out
               kh         => state%kh        , & ! out
               cdm        => state%cdm       , & ! out
               cdh        => state%cdh       , & ! out
               ustar      => state%ustar     , & ! out
               tstar      => state%tstar     , & ! out
               taux       => state%taux      , & ! out
               tauy       => state%tauy      , & ! out
               ht_pbl     => state%ht_pbl    , & ! out
               rhouch     => state%rhouch    , & ! out
               qsfc       => state%qsfc      )   ! inout
    do i = 1, mesh%ncol
      h     = z    (i,:) - zs(i)
      h_lev = z_lev(i,:) - zs(i)
      ! ------------------------------------------------------------------------
      ! Set surface roughness length.
      z0(i) = z00
      if (co2ice_sfc(i) > 0) z0(i) = 1.0e-4_r8
      if (qsfc(i,iMa_vap) > 100 .or. npcflag(i)) z0(i) = 1.0e-4_r8
      ! ------------------------------------------------------------------------
      ! Calculate eddy mixing coefficients.
      do k = 2, mesh%nlev ! Loop on half levels excluding top and bottom.
        ! Calculate mixing length.
        ml2 = (ml0 * ka * h_lev(k) / (ml0 + ka * h_lev(k)))**2
        ! Calculate gradient Richardson number.
        dudz    = -(u (i,k) - u (i,k-1)) / dz_lev(i,k)
        dvdz    = -(v (i,k) - v (i,k-1)) / dz_lev(i,k)
        dptdz   = -(pt(i,k) - pt(i,k-1)) / dz_lev(i,k)
        ! Smooth the wind shear,
        shr2(i,k) = shr2(i,k) - (shr2(i,k) - dudz**2 - dvdz**2) * dt / 1.0e4_r8
        shr       = sqrt(shr2(i,k))
        ri  (i,k) = g / pt_lev(i,k) * dptdz / (shr2(i,k) + 1.0e-9_r8)
        ! Calculate neutral eddy coefficients.
        km0 = Sm * ml2 * shr / sqrtGM
        kh0 = Sh * ml2 * shr / sqrtGM
        ! Calculate eddy mixing coefficients.
        if (ri(i,k) <= 0) then
          km(i,k) = km0 * (1 - 15 * ri(i,k))**0.25_r8
          kh(i,k) = kh0 * (1 - 15 * ri(i,k))**0.50_r8
        else
          km(i,k) = km0 * (1 - ri(i,k) / ric)
          kh(i,k) = kh0 * (1 - ri(i,k) / ric)
        end if
        ! Limit the coefficients.
        kmin = merge(0.1_r8, 0.001_r8, h_lev(k) < 300)
        km(i,k) = max(km(i,k), kmin)
        kh(i,k) = max(kh(i,k), kmin)
      end do
      ! ------------------------------------------------------------------------
      lnzz = log(h(nlev) / z0(i))
      ! For now, use psi at the bottom full level. Later interpolate to the surface.
      dpt = pt(i,nlev) - tg(i)
      ! Calculate the bulk Richardson number.
      wsp = sqrt(u(i,nlev)**2 + v(i,nlev)**2)
      rib = g * dpt * h(nlev) / (pt(i,nlev) * wsp**2 + 1.0e-9_r8)
      ! Calculate the stability functions.
      if (rib >= 0) then
        fm = 1.0_r8 / (1 + 10 * rib / sqrt(1 + 5 * rib))
        fh = 1.0_r8 / (1 + 15 * rib / sqrt(1 + 5 * rib))
      else
        fm = sqrt(1 - 16 * rib)
        fh = sqrt(1 - 64 * rib)
      end if
      ! Calculate the drag coefficients.
      cdm(i) = fm * (ka / lnzz)**2
      cdh(i) = sqrt(fh) * ka / lnzz
      ! Calculate the friction velocity.
      ustar(i) = sqrt(cdm(i)) * wsp
      ! Calculate the temperature scale.
      tstar(i) = cdh(i) * dpt
      ! ------------------------------------------------------------------------
      ! Diagnose the surface stress and heat fluxes.
      rhos      =  dry_air_density(ts(i), ps(i))
      alpha     =  atan2(v(i,nlev), u(i,nlev))
      taux  (i) =  rhos * ustar(i)**2 * cos(alpha)
      tauy  (i) =  rhos * ustar(i)**2 * sin(alpha)
      ht_pbl(i) = -rhos * cpd * ustar(i) * tstar(i)
      rhouch(i) =  rhos * cpd * ustar(i) * cdh  (i)
      ! ------------------------------------------------------------------------
      ! Reduce the sublimation flux by a coefficient. It is a tunable parameter to
      ! avoid the formation of low-lying clouds in summer above the north permanent cap.
      qsat = water_vapor_saturation_mixing_ratio_mars(tg(i), ps(i))
      coef = merge(1.0_r8, 1.0_r8, qsat > q(i,nlev,iMa_vap))
      ! ------------------------------------------------------------------------
      ! Solve the implicit diffusion equations.
      ! - Zonal wind
      kdf    = km(i,:)
      bnd    = dt / dz(i,nlev) * sqrt(cdm(i)) * ustar(i)
      rhs    = 0
      var    = rho(i,:) * u(i,:)
      call pbl_solve(dz(i,:), dz_lev(i,:), kdf, bnd, rhs, var)
      u(i,:) = var / rho(i,:)
      ! - Meridional wind
      kdf    = km(i,:)
      bnd    = dt / dz(i,nlev) * sqrt(cdm(i)) * ustar(i)
      rhs    = 0
      var    = rho(i,:) * v(i,:)
      call pbl_solve(dz(i,:), dz_lev(i,:), kdf, bnd, rhs, var)
      v(i,:) = var / rho(i,:)
      ! - Potential temperature
      kdf    = kh(i,:)
      bnd    = dt / dz(i,nlev) * cdh(i) * ustar(i)
      rhs    = rho(i,:) * qrad(i,:) * dt
      var    = rho(i,:) * pt(i,:)
      rhs(nlev) = rhs(nlev) + dt / dz(i,nlev) * rho(i,nlev) * cdh(i) * ustar(i) * tg(i)
      call pbl_solve(dz(i,:), dz_lev(i,:), kdf, bnd, rhs, var)
      pt(i,:) = var / rho(i,:)
      ! - Water vapor
      kdf    = kh(i,:)
      bnd    = dt / dz(i,nlev) * ustar(i) * cdh(i) * coef
      rhs    = 0
      var    = rho(i,:) * q(i,:,iMa_vap)
      rhs(nlev) = h2osub_sfc(i) + dt / dz(i,nlev) * rho(i,nlev) * ustar(i) * cdh(i) * coef
      call pbl_solve(dz(i,:), dz_lev(i,:), kdf, bnd, rhs, var)
      q(i,:,iMa_vap) = var / rho(i,:)
      ! ------------------------------------------------------------------------
      ! Update water ice budget on the surface.
      qsfc(i,iMa_vap) = qsfc(i,iMa_vap) - h2osub_sfc(i)
      if (.not. npcflag(i) .and. qsfc(i,iMa_vap) < 0) then
        qsfc(i,iMa_vap) = 0
      end if
      ! ------------------------------------------------------------------------
      do k = 1, mesh%nlev
        t(i,k) = pt(i,k) * pk(i,k)
      end do
    end do
    end associate

  end subroutine gomars_v1_pbl_run

  subroutine pbl_solve(dz, dz_lev, kdf, bnd, rhs, var)

    use math_mod

    real(r8), intent(in   ) :: dz    (nlev  )
    real(r8), intent(in   ) :: dz_lev(nlev+1)
    real(r8), intent(in   ) :: kdf   (nlev+1)
    real(r8), intent(in   ) :: bnd            ! Lower boundary condition for coefficient c
    real(r8), intent(in   ) :: rhs   (nlev  ) ! Right hand side
    real(r8), intent(inout) :: var   (nlev  )

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

end module gomars_v1_pbl_mod