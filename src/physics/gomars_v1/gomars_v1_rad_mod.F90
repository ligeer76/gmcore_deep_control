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

module gomars_v1_rad_mod

  use math_mod
  use gomars_v1_const_mod
  use gomars_v1_types_mod
  use gomars_v1_namelist_mod
  use gomars_v1_orbit_mod
  use gomars_v1_tracers_mod
  use gomars_v1_mp_mod

  implicit none

  private

  public gomars_v1_rad_init
  public gomars_v1_rad_final
  public gomars_v1_rad_run
  public update_toa_solar_flux
  public update_sfc_solar_flux

  real(r8), parameter :: ubari = 0.5_r8
  ! If optical depth is less than this value, place the Gauss-point into the
  ! "zero" channel.
  real(r8), parameter :: tlimits = 1.0e-3_r8
  real(r8), parameter :: tlimiti = 5.0e-3_r8
  real(r8), parameter :: maxexp  = 35.0_r8
  real(r8), parameter :: taumax  = 35.0_r8
  ! Rayleigh scattering reference pressure in hPa.
  real(r8), parameter :: pray0   = 9.423e+4_r8

  real(r8), allocatable, dimension(:) :: wnoi
  real(r8), allocatable, dimension(:) :: dwni
  real(r8), allocatable, dimension(:) :: wavei
  real(r8), allocatable, dimension(:) :: wnov
  real(r8), allocatable, dimension(:) :: dwnv
  real(r8), allocatable, dimension(:) :: wavev
  real(r8), allocatable, dimension(:) :: solar
  real(r8), allocatable, dimension(:) :: tauray

  real(r8), allocatable, dimension(:,:,:,:,:) :: co2i
  real(r8), allocatable, dimension(:,:,:,:,:) :: co2v

  real(r8), allocatable, dimension(:) :: fzeroi
  real(r8), allocatable, dimension(:) :: fzerov

  real(r8) qextref, ptop
  real(r8), allocatable, dimension(:) :: wv
  real(r8), allocatable, dimension(:) :: wi

  real(r8), allocatable, dimension(:,:) :: planckir

  real(r8), allocatable, dimension(:) :: pfgasref ! log(p [hPa])

  real(r8), parameter :: solar_1au(nspectv) = [ &
    12.7d0, 24.2d0, 54.6d0, 145.9d0, 354.9d0, 657.5d0, 106.3d0]

  ! Visible spectral boundaries (cm-1)
  real(r8), parameter :: bwnv(nspectv+1) = [ &
     2222.22d0, & ! 4.500 microns
     3087.37d0, & ! 3.239 microns
     4030.63d0, & ! 2.481 microns
     5370.57d0, & ! 1.862 microns
     7651.11d0, & ! 1.307 microns
    12500.00d0, & ! 0.800 microns
    25000.00d0, & ! 0.400 microns
    41666.67d0]   ! 0.240 microns

  ! Infrared spectral boundaries (cm-1)
  real(r8), parameter :: bwni(nspecti+1) = [ &
      10.000d0, & ! 1000.0 microns
     166.667d0, & !   60.0 microns
     416.667d0, & !   24.0 microns
     833.333d0, & !   12.0 microns
    1250.000d0, & !    8.0 microns
    2222.222d0]   !    4.5 microns

  real(r8), parameter :: gweight(ngauss) = [ &
    4.8083554740d-02, 1.0563099137d-01,      &
    1.4901065679d-01, 1.7227479710d-01,      &
    1.7227479710d-01, 1.4901065679d-01,      &
    1.0563099137d-01, 4.8083554740d-02,      &
    2.5307134073d-03, 5.5595258613d-03,      &
    7.8426661469d-03, 9.0670945845d-03,      &
    9.0670945845d-03, 7.8426661469d-03,      &
    5.5595258613d-03, 2.5307134073d-03, 0.0d0]

  ! Reference volume mixing ratio of CO2
  real(r8), parameter :: wrefco2(nrefh2o) = [      &
    9.999999d-1, 9.99999d-1, 9.9999d-1, 9.999d-1,  &
    9.99d-1, 9.9d-1, 9.0d-1, 8.0d-1, 7.0d-1, 6.0d-1]

  ! Reference volume ratio of H2O
  real(r8), parameter :: wrefh2o(nrefh2o) = [       &
    1.0d-7, 1.0d-6, 1.0d-5, 1.0d-4, 1.0d-3, 1.0d-2, &
    1.0d-1, 2.0d-1, 3.0d-1, 4.0d-1]

  real(r8), parameter :: qextv(nspectv) = [ &
    1.834d0, &
    2.296d0, &
    2.672d0, &
    2.829d0, &
    2.698d0, &
    2.452d0, &
    2.261d0]
  real(r8) :: qscatv(nspectv) = [ &
    1.695d0, &
    2.031d0, &
    2.583d0, &
    2.744d0, &
    2.626d0, &
    2.225d0, &
    1.525d0]
  real(r8), parameter :: gv(nspectv) = [ &
    0.551d0, &
    0.640d0, &
    0.661d0, &
    0.678d0, &
    0.690d0, &
    0.743d0, &
    0.868d0]

  real(r8), parameter :: qexti(nspecti) = [ &
    0.008d0, &
    0.262d0, &
    0.491d0, &
    1.017d0, &
    0.444d0]
  real(r8) :: qscati(nspecti) = [ &
    0.001d0, &
    0.037d0, &
    0.122d0, &
    0.351d0, &
    0.336d0]
  real(r8), parameter :: gi(nspecti) = [ &
    0.004d0, &
    0.030d0, &
    0.095d0, &
    0.214d0, &
    0.316d0]

  real(r8), parameter :: pgasref(npref) = [ & ! hPa
    1.0d-6, &
    1.0d-5, &
    1.0d-4, &
    1.0d-3, &
    1.0d-2, &
    1.0d-1, &
    1.0d0 , &
    1.0d+1, &
    1.0d+2, &
    1.0d+3, &
    1.0d+4]
  real(r8), parameter :: tgasref(ntref) = [ &
     50.0d0, &
    100.0d0, &
    150.0d0, &
    200.0d0, &
    250.0d0, &
    300.0d0, &
    350.0d0]

contains

  subroutine gomars_v1_rad_init()

    integer is, it, i
    real(r8) a, b, ans, y, bpa, bma

    ! c1 and c2 values from Goody and Yung (2nd edition)  MKS units
    ! These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4
    real(r8), parameter :: c1 = 3.741832d-16 ! W m-2
    real(r8), parameter :: c2 = 1.438786d-2  ! m K
    
    real(r8) :: x(12) = [-0.981560634246719D0,  -0.904117256370475D0, &
                         -0.769902674194305D0,  -0.587317954286617D0, &
                         -0.367831498998180D0,  -0.125233408511469D0, &
                          0.125233408511469D0,   0.367831498998180D0, &
                          0.587317954286617D0,   0.769902674194305D0, &
                          0.904117256370475D0,   0.981560634246719D0]

    real(r8) :: w(12) = [ 0.047175336386512D0,   0.106939325995318D0, &
                          0.160078328543346D0,   0.203167426723066D0, &
                          0.233492536538355D0,   0.249147045813403D0, &
                          0.249147045813403D0,   0.233492536538355D0, &
                          0.203167426723066D0,   0.160078328543346D0, &
                          0.106939325995318D0,   0.047175336386512D0]

    call gomars_v1_rad_final()

    allocate(wnoi     (nspecti))
    allocate(dwni     (nspecti))
    allocate(wavei    (nspecti))
    allocate(wnov     (nspectv))
    allocate(dwnv     (nspectv))
    allocate(wavev    (nspectv))
    allocate(solar    (nspectv))
    allocate(tauray   (nspectv))

    allocate(co2i(ntref,npint,nrefh2o,nspecti,ngauss))
    allocate(co2v(ntref,npint,nrefh2o,nspectv,ngauss))

    allocate(fzeroi(nspecti))
    allocate(fzerov(nspectv))

    allocate(wv    (nspectv))
    allocate(wi    (nspecti))

    allocate(planckir(nspecti,8501))

    allocate(pfgasref(npint))

    ! Set visible spectrum.
    do is = 1, nspectv
      wnov  (is) = 0.5_r8 * (bwnv(is+1) + bwnv(is))
      dwnv  (is) = bwnv(is+1) - bwnv(is)
      wavev (is) = 1.0e+4_r8 / wnov(is)
      tauray(is) = (8.7_r8 / g) * (1.527_r8 * (1.0 + 0.013 / wavev(is)**2)/ wavev(is)**4) / pray0
    end do

    ! Set infrared spectrum.
    do is = 1, nspecti
      wnoi  (is) = 0.5_r8 * (bwni(is+1) + bwni(is))
      dwni  (is) = bwni(is+1) - bwni(is)
      wavei (is) = 1.0e+4_r8 / wnoi(is)
    end do

    ! For each IR wavelength interval, compute the integral of B(T), the
    ! Planck function, divided by the wavelength interval, in cm-1.  The
    ! integration is in MKS units, the final answer is the same as the
    ! original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.
    do is = 1, nspecti
      a = 1.0e-2_r8 / bwni(is+1)
      b = 1.0e-2_r8 / bwni(is  )
      bpa = 0.5_r8 * (b + a)
      bma = 0.5_r8 * (b - a)
      do it = 500, 9000
        ans = 0
        do i = 1, 12
          y   = bma * x(i) + bpa
          ans = ans + w(i) * c1 / (y**5 * (exp(c2 / (y * dble(it) / 1.0d+1)) - 1))
        end do
        planckir(is,it-499) = ans * bma / (pi * dwni(is))
      end do
    end do

    ! Fill the (VISIBLE) arrays Qextv, Qscatv, WV, GV
    do is = 1, nspectv
      if (qscatv(is) >= qextv(is)) then
        qscatv(is) = 0.99999_r8 * qextv(is)
      end if
      wv(is) = qscatv(is) / qextv(is)
    end do

    ! Fill the (INFRARED) arrays Qexti, Qscati, WI, GI
    do is = 1, nspecti
      if (qscati(is) >= qexti(is)) then
        qscati(is) = 0.99999_r8 * qexti(is)
      end if
      wi(is) = qscati(is) / qexti(is)
    end do

    ! Interpolate CO2 k coefficients to the finer pressure grid.
    call laginterp(pgasref, pfgasref, co2i, co2v, fzeroi, fzerov)

    ptop = 10**pfgasref(1) * 100 ! Pa

  end subroutine gomars_v1_rad_init

  subroutine gomars_v1_rad_final()

    if (allocated(wnoi     )) deallocate(wnoi     )
    if (allocated(dwni     )) deallocate(dwni     )
    if (allocated(wavei    )) deallocate(wavei    )
    if (allocated(wnov     )) deallocate(wnov     )
    if (allocated(dwnv     )) deallocate(dwnv     )
    if (allocated(wavev    )) deallocate(wavev    )
    if (allocated(solar    )) deallocate(solar    )
    if (allocated(tauray   )) deallocate(tauray   )
    if (allocated(co2i     )) deallocate(co2i     )
    if (allocated(co2v     )) deallocate(co2v     )
    if (allocated(fzeroi   )) deallocate(fzeroi   )
    if (allocated(fzerov   )) deallocate(fzerov   )
    if (allocated(wv       )) deallocate(wv       )
    if (allocated(wi       )) deallocate(wi       )
    if (allocated(planckir )) deallocate(planckir )
    if (allocated(pfgasref )) deallocate(pfgasref )

  end subroutine gomars_v1_rad_final

  subroutine gomars_v1_rad_run(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k
    real(r8), dimension(2*nlev+3                 ) :: pmid ! p_lev + p_lev
    real(r8), dimension(2*nlev+3                 ) :: plev ! p     + p_lev
    real(r8), dimension(2*nlev+3                 ) :: tmid ! t_lev + t_lev
    real(r8), dimension(2*nlev+3                 ) :: tlev ! t     + t_lev
    real(r8), dimension(2*nlev+3                 ) :: qh2o
    real(r8), dimension(2*nlev+4,nspectv         ) :: qxvdst    , qxvcld
    real(r8), dimension(2*nlev+4,nspectv         ) :: qsvdst    , qsvcld
    real(r8), dimension(2*nlev+4,nspectv         ) :: gvdst     , gvcld
    real(r8), dimension(2*nlev+4,nspecti         ) :: qxidst    , qxicld
    real(r8), dimension(2*nlev+4,nspecti         ) :: qsidst    , qsicld
    real(r8), dimension(2*nlev+4,nspecti         ) :: gidst     , gicld
    real(r8), dimension(2*nlev+4                 ) :: qextrefdst, qextrefcld
    real(r8), dimension(2*nlev+4                 ) :: taurefdst , taurefcld
    real(r8), dimension(2*nlev+4,2               ) :: taudst    , taucld
    real(r8), dimension(nlayrad ,nspectv,ngauss  ) :: wbarv
    real(r8), dimension(nlayrad ,nspectv,ngauss  ) :: cosbv
    real(r8), dimension(nlayrad ,nspecti,ngauss  ) :: wbari
    real(r8), dimension(nlayrad ,nspecti,ngauss  ) :: cosbi
    real(r8), dimension(nlayrad ,nspectv,ngauss  ) :: dtauv
    real(r8), dimension(nlayrad ,nspectv,ngauss  ) :: tauv
    real(r8), dimension(2*nlev+3,nspectv,ngauss  ) :: taucumv
    real(r8), dimension(nlayrad ,nspecti,ngauss  ) :: dtaui
    real(r8), dimension(nlayrad ,nspecti,ngauss  ) :: taui
    real(r8), dimension(2*nlev+3,nspecti,ngauss  ) :: taucumi
    real(r8), dimension(2*nlev+3                 ) :: taucum
    ! Total gas opacity through the column for each spectral interval and each Gauss point
    real(r8), dimension(         nspectv,ngauss-1) :: taugsurf
    real(r8), dimension(nlayrad                  ) :: fluxdnv
    real(r8), dimension(nlayrad                  ) :: fluxupv
    real(r8), dimension(nlayrad                  ) :: fmnetv
    real(r8), dimension(nlayrad                  ) :: fluxdni
    real(r8), dimension(nlayrad                  ) :: fluxupi
    real(r8), dimension(nlayrad                  ) :: fmneti
    real(r8), dimension(2*nlev+3                 ) :: suntot
    real(r8), dimension(2*nlev+3                 ) :: irtot
    real(r8) nfluxtopv, diffvt
    real(r8) nfluxtopi, albi, nonlte, ptstrat

    qxvdst = 0; qxidst = 0
    qsvdst = 0; qsidst = 0
    gvdst  = 0; gidst  = 0
    qxvcld = 0; qxicld = 0
    qsvcld = 0; qsicld = 0
    gvcld  = 0; gicld  = 0
    qextrefcld = 1
    taurefcld  = 0

    associate (mesh         => state%mesh        , &
               lat          => state%mesh%lat    , & ! in
               cosz         => state%cosz        , & ! in
               alsp         => state%alsp        , & ! in
               als          => state%als         , & ! inout
               ps           => state%ps          , & ! in
               dp           => state%dp_dry      , & ! in
               pk           => state%pk          , & ! in
               p            => state%p           , & ! in
               p_lev        => state%p_lev       , & ! in
               t            => state%t           , & ! in
               t_lev        => state%t_lev       , & ! in
               tg           => state%tg          , & ! in
               tstrat       => state%tstrat      , & ! inout
               q            => state%q           , & ! in
               qsfc         => state%qsfc        , & ! in
               co2ice_sfc   => state%co2ice_sfc  , & ! in
               npcflag      => state%npcflag     , & ! in
               tausurf      => state%tausurf     , & ! out
               detau        => state%detau       , & ! out
               vsflx_sfc_dn => state%vsflx_sfc_dn, & ! out
               vsdif_sfc_dn => state%vsdif_sfc_dn, & ! out
               irflx_sfc_dn => state%irflx_sfc_dn, & ! out
               fluxsfc      => state%fluxsfc     , & ! out
               fuptopv      => state%fuptopv     , & ! out
               fdntopv      => state%fdntopv     , & ! out
               fupsfcv      => state%fupsfcv     , & ! out
               fdnsfcv      => state%fdnsfcv     , & ! out
               fuptopi      => state%fuptopi     , & ! out
               fupsfci      => state%fupsfci     , & ! out
               fdnsfci      => state%fdnsfci     , & ! out
               qrad         => state%qrad        )   ! out
    do i = 1, mesh%ncol
      ! ------------------------------------------------------------------------
      ! Update surface albedo.
      als(i) = alsp(i)
      if (co2ice_sfc(i) > 0) then
        als(i) = merge(alices, alicen, lat(i) < 0)
      else if (albfeed .and. qsfc(i,iMa_vap) > icethresh_kgm2 .and. .not. npcflag(i)) then
        als(i) = icealb
      end if
      ! ------------------------------------------------------------------------
      ! Fill the radiation variables.
      call fillpt(p(i,:), p_lev(i,:), t(i,:), t_lev(i,:), tg(i), tstrat(i), plev, tlev, pmid, tmid)
      qh2o(1:3) = 0
      if (active_water) then
        do k = 1, nlev
          qh2o(2*k+2) = m_co2 / m_h2o * q(i,k,iMa_vap)
          qh2o(2*k+3) = qh2o(2*k+2)
        end do
      else
        do k = 1, nlev
          qh2o(2*k+2) = 1.0e-7_r8
          qh2o(2*k+3) = qh2o(2*k+2)
        end do
      end if
      ! ------------------------------------------------------------------------
      ! Calculate optical properties of dust.
      if (active_dust) then
        call opt_dst(q(i,:,:), plev, qxvdst, qsvdst, gvdst, qxidst, qsidst, gidst, qextrefdst, taurefdst, taudst)
        taurefdst(1:3) = 0
        taucum   (1:3) = 0
        do k = 4, 2 * nlev + 3
          taucum(k) = taucum(k-1) + taurefdst(k)
        end do
      end if
      ! Fill special bottom radiation level to zero.
      taurefdst(2*nlev+4) = 0
      tausurf(i) = taucum(2*nlev+3)
      ! ------------------------------------------------------------------------
      ! Calculate optical properties of water ice cloud.
      if (cloudon) then
        call opt_cld(q(i,:,:), plev, qxvcld, qsvcld, gvcld, qxicld, qsicld, gicld, qextrefcld, taurefcld, taucld)
      else
        taurefcld = 0
      end if
      ! ------------------------------------------------------------------------
      ! Set up and solve for the visible fluxes.
      if (cosz(i) >= 1.0e-5_r8) then
        ! Sun is up.
        call optcv(plev, pmid, tmid, qh2o, qxvdst, qsvdst, gvdst, qxvcld, qsvcld, gvcld, qextrefcld, &
                   wbarv, cosbv, dtauv, tauv, taucumv, taugsurf, taurefdst, taurefcld)
        call sfluxv(dtauv, tauv, taucumv, taugsurf, cosz(i), als(i), wbarv, cosbv, &
                    fluxupv, fluxdnv, fmnetv, nfluxtopv, diffvt, detau(i,:,:))
        suntot(3) = fmnetv(1) - nfluxtopv
        do k = 2, nlayrad
          suntot(2*k+1) = fmnetv(k) - fmnetv(k-1)
        end do
      else
        ! Sun is down.
        do k = 1, nlayrad
          suntot(2*k+1) = 0
          fluxdnv(k)    = 0
        end do
        diffvt           = 0
        fluxupv(1)       = 0
        fluxdnv(1)       = 0
        fluxupv(nlayrad) = 0
        fluxdnv(nlayrad) = 0
      end if
      vsflx_sfc_dn(i) = fluxdnv(nlayrad)
      vsdif_sfc_dn(i) = diffvt
      fuptopv     (i) = fluxupv(1)
      fdntopv     (i) = fluxdnv(1)
      fupsfcv     (i) = fluxupv(nlayrad)
      fdnsfcv     (i) = fluxdnv(nlayrad)
      ! ------------------------------------------------------------------------
      ! Set up and solve for the infrared fluxes.
      ! Check for ground ice, and change infrared albedo if there is any ice.
      albi = 1 - egognd
      if (co2ice_sfc(i) > 0) then
        albi = 1 - merge(egoco2s, egoco2n, mesh%lat(i) < 0)
      end if
      ! Calculate the optical depth due to all sources in the infrared bands.
      call optci(plev, pmid, tmid, qh2o, qxidst, qsidst, gidst, qextrefdst, &
                 qxicld, qsicld, gicld, qextrefcld, wbari, cosbi, dtaui, taui, &
                 taucumi, taugsurf, taurefdst, taurefcld)
      ! Calculate the fluxes in the infrared bands.
      call sfluxi(plev, tlev, dtaui, taucumi, taugsurf, albi, wbari, cosbi, &
                  fluxupi, fluxdni, fmneti, nfluxtopi)
      irtot(3) = fmneti(1) - nfluxtopi
      do k = 2, nlayrad
        irtot(2*k+1) = fmneti(k) - fmneti(k-1)
      end do
      irflx_sfc_dn(i) = fluxdni(nlayrad)
      fluxsfc     (i) = (1 - als(i)) * fluxdnv(nlayrad)
      fuptopi     (i) = fluxupi(1)
      fupsfci     (i) = fluxupi(nlayrad)
      fdnsfci     (i) = fluxdni(nlayrad)
      ! Set the solar and infrared heating. Also add the non-LTE correction.
      do k = 1, nlev
        nonlte = 2.2e2_r8 * plev(2*k+2) / (1 + 2.2e2_r8 * plev(2*k+2))
        qrad(i,k) = (suntot(2*k+3) * nonlte + irtot(2*k+3)) / (cpd * dp(i,k) / g * pk(i,k))
      end do
      ! Update stratospheric temperature.
      nonlte = 2.2e2_r8 * pstrat / (1 + 2.2e2_r8 * pstrat)
      tstrat(i) = tstrat(i) + dt * (suntot(3) * nonlte + irtot(3)) / (cpd * ptrop / g)
    end do
    end associate

  end subroutine gomars_v1_rad_run

  subroutine update_toa_solar_flux(ls)

    real(r8), intent(in) :: ls

    integer is
    real(r8) rsdist

    rsdist = solar_dist(ls)**2
    
    ! Calculate solar flux at the current Mars distance.
    do is = 1, nspectv
      solar(is) = solar_1au(is) / rsdist
    end do

  end subroutine update_toa_solar_flux

  subroutine update_sfc_solar_flux(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, is, ig
    real(r8) factor

    associate (mesh         => state%mesh        , &
               cosz         => state%cosz        , & ! in
               detau        => state%detau       , & ! in
               vsdif_sfc_dn => state%vsdif_sfc_dn, & ! in
               solar_sfc_dn => state%solar_sfc_dn, & ! out
               vsflx_sfc_dn => state%vsflx_sfc_dn)   ! inout
    do i = 1, mesh%ncol
      solar_sfc_dn(i) = 0
      if (cosz(i) >= 1.0e-5_r8) then
        do is = 1, nspectv
          factor = cosz(i) * solar(is)
          do ig = 1, ngauss - 1
            if (detau(i,is,ig) <= 5) then
              solar_sfc_dn(i) = solar_sfc_dn(i) + factor * exp(-detau(i,is,ig) / cosz(i)) * gweight(ig) * (1 - fzerov(is))
            end if
          end do
          ig = ngauss
          if (detau(i,is,ig) <= 5) then
            solar_sfc_dn(i) = solar_sfc_dn(i) + factor * exp(-detau(i,is,ig) / cosz(i)) * fzerov(is)
          end if
        end do
      end if
      vsflx_sfc_dn(i) = vsdif_sfc_dn(i) + solar_sfc_dn(i)
    end do
    end associate

  end subroutine update_sfc_solar_flux

  subroutine fillpt(p, p_lev, t, t_lev, tg, tstrat, plev, tlev, pmid, tmid)

    ! Put the pressure and temperature arrays onto the radiation grid.

    real(r8), intent(in ) :: p    (nlev  )
    real(r8), intent(in ) :: p_lev(nlev+1)
    real(r8), intent(in ) :: t    (nlev  )
    real(r8), intent(in ) :: t_lev(nlev+1)
    real(r8), intent(in ) :: tg
    real(r8), intent(in ) :: tstrat
    real(r8), intent(out) :: plev(2*nlev+3)
    real(r8), intent(out) :: tlev(2*nlev+3)
    real(r8), intent(out) :: pmid(2*nlev+3)
    real(r8), intent(out) :: tmid(2*nlev+3)

    integer k

    plev(1:2) = ptrop * 0.5_r8
    do k = 1, nlev
      plev(2*k+2) = p(k)
    end do
    do k = 1, nlev + 1
      plev(2*k+1) = p_lev(k)
    end do

    tlev(1:3) = tstrat
    do k = 1, nlev
      tlev(2*k+2) = t(k)
    end do
    do k = 2, nlev
      tlev(2*k+1) = t_lev(k)
    end do
    tlev(2*nlev+3) = tg

    ! tmid_rad and pmid_rad are used by optci and optcv subroutines to get the
    ! index for CO2 K-coefficient interpolation.
    tmid(1) = tlev(2)
    tmid(2) = tlev(2)
    pmid(1) = plev(1)
    pmid(2) = plev(2)
    do k = 1, nlev
      tmid(2*k+1) = tlev(2*k+1)
      tmid(2*k+2) = tlev(2*k+1)
      pmid(2*k+1) = plev(2*k+1)
      pmid(2*k+2) = plev(2*k+1)
    end do
    tmid(2*nlev+3) = tlev(2*nlev+3)
    pmid(2*nlev+3) = plev(2*nlev+3)

  end subroutine fillpt

  subroutine opt_dst(q, plev, qxv, qsv, gv, qxi, qsi, gi, qextrefdst, taurefdst, taudst)

    real(r8), intent(in ) :: q         (nlev,ntracers)
    real(r8), intent(in ) :: plev      (2*nlev+3)
    real(r8), intent(out) :: qxv       (2*nlev+4,nspectv)
    real(r8), intent(out) :: qsv       (2*nlev+4,nspectv)
    real(r8), intent(out) :: gv        (2*nlev+4,nspectv)
    real(r8), intent(out) :: qxi       (2*nlev+4,nspecti)
    real(r8), intent(out) :: qsi       (2*nlev+4,nspecti)
    real(r8), intent(out) :: gi        (2*nlev+4,nspecti)
    real(r8), intent(out) :: qextrefdst(2*nlev+4)
    real(r8), intent(out) :: taurefdst (2*nlev+4)
    real(r8), intent(out) :: taudst    (2)

    integer k, l, is, ib
    real(r8) dev2, cst
    real(r8) Mo, No, Rs, Ao
    real(r8) surf(nbin_rt)

    ! Set default values.
    do l = 1, 2 * nlev + 4
      qextrefdst(l) = 1
      taurefdst (l) = 0
    end do
    do is = 1, nspectv
      do l = 1, 2 * nlev + 4
        qxv(l,is) = qextrefdst(l)
        qsv(l,is) = qextrefdst(l) * 0.99_r8
        gv (l,is) = 0
      end do
    end do
    do is = 1, nspecti
      do l = 1, 2 * nlev + 4
        qxi(l,is) = qextrefdst(l)
        qsi(l,is) = qextrefdst(l) * 0.99_r8
        gi (l,is) = 0
      end do
    end do

    dev2 = 1.0_r8 / (sqrt2 * dev_dst)
    cst = 0.75_r8 / (pi * rho_dst)

    taudst = 0

    do k = 1, nlev
      if (q(k,iMa_dst) > 1.0e-8_r8 .and. q(k,iNb_dst) > 1) then
        ! Calculate the cross-section mean radius (Rs) of the log-normal distribution.
        Mo = q(k,iMa_dst)
        No = q(k,iNb_dst) + 1
        Rs = (Mo / No * cst)**athird * exp(-0.5_r8 * dev_dst**2)
        ! Calculate the total cross sectional area (Ao) of water ice particles.
        Ao = No * pi * Rs**2
        ! Define the cross-section weighted distribution, i.e. surface/size bin. Change Rs to Reff.
        Rs = 1.0_r8 / min(max(Rs * exp(1.5_r8 * dev_dst**2), 1.0e-7_r8), 50.0e-6_r8)
        do ib = 1, nbin_rt
          surf(ib) = 0.5_r8 * (erf(log(radb_rt(ib+1) * Rs) * dev2) - erf(log(radb_rt(ib) * Rs) * dev2))
        end do
        ! Calculate the averaged values of the optical properties for the whole distribution.
        do l = 2 * k + 2, 2 * k + 3
          do is = 1, nspectv
            qxv(l,is) = 0
            qsv(l,is) = 0
            do ib = 1, nbin_rt
              qxv(l,is) = qxv(l,is) + qextv_dst (ib,is) * surf(ib)
              qsv(l,is) = qsv(l,is) + qscatv_dst(ib,is) * surf(ib)
              gv (l,is) = gv (l,is) + gv_dst    (ib,is) * surf(ib)
            end do
            qsv(l,is) = min(qsv(l,is), 0.99999_r8 * qxv(l,is))
          end do
          do is = 1, nspecti
            qxi(l,is) = 0
            qsi(l,is) = 0
            do ib = 1, nbin_rt
              qxi(l,is) = qxi(l,is) + qexti_dst (ib,is) * surf(ib)
              qsi(l,is) = qsi(l,is) + qscati_dst(ib,is) * surf(ib)
              gi (l,is) = gi (l,is) + gi_dst    (ib,is) * surf(ib)
            end do
            qsi(l,is) = min(qsi(l,is), 0.99999_r8 * qxi(l,is))
          end do
          qextrefdst(l) = qxv(l,nrefv)
          taurefdst (l) = Ao * qextrefdst(l) * (plev(l) - plev(l-1)) / g
          ! For diagnostics: dust opacity at reference wavelengths (vis and ir).
          taudst(1) = taudst(1) + taurefdst(l)
          taudst(2) = taudst(2) + taurefdst(l) / qextrefdst(l) * (qxi(l,4) - qsi(l,4))
        end do
      end if
    end do

  end subroutine opt_dst

  subroutine opt_cld(q, plev, qxv, qsv, gv, qxi, qsi, gi, qextrefcld, taurefcld, taucld)

    real(r8), intent(in ) :: q         (nlev,ntracers)
    real(r8), intent(in ) :: plev      (2*nlev+3)
    real(r8), intent(out) :: qxv       (2*nlev+4,nspectv)
    real(r8), intent(out) :: qsv       (2*nlev+4,nspectv)
    real(r8), intent(out) :: gv        (2*nlev+4,nspectv)
    real(r8), intent(out) :: qxi       (2*nlev+4,nspecti)
    real(r8), intent(out) :: qsi       (2*nlev+4,nspecti)
    real(r8), intent(out) :: gi        (2*nlev+4,nspecti)
    real(r8), intent(out) :: qextrefcld(2*nlev+4)
    real(r8), intent(out) :: taurefcld (2*nlev+4)
    real(r8), intent(out) :: taucld    (2)

    integer n, k, l, i, is, ib
    integer irap
    real(r8) dev2, cst
    real(r8) mantletocore
    real(r8) Mo, No, Rs, Ao
    real(r8) surf(nbin_rt)

    n = 2 * nlev + 3

    do l = 1, n + 1
      qextrefcld(l) = 1
      taurefcld (l) = 0
    end do

    do is = 1, nspectv
      do l = 1, n + 1
        qxv(l,is) = qextrefcld(l)
        qsv(l,is) = qextrefcld(l) * 0.99_r8
        gv (l,is) = 0
      end do
    end do

    do is = 1, nspecti
      do l = 1, n + 1
        qxi(l,is) = qextrefcld(l)
        qsi(l,is) = qextrefcld(l) * 0.99_r8
        gi (l,is) = 0
      end do
    end do

    dev2 = 1.0_r8 / (sqrt2 * dev_ice)

    taucld = 0

    do k = 1, nlev
      if (q(k,iMa_cld) + q(k,iMa_cor) > 1.0e-7 .and. q(k,iNb_cld) > 1) then
        ! Determine the ratio of dust core over that of ice mantle.
        mantletocore = (q(k,iMa_cor) / rho_dst / (q(k,iMa_cld) / rho_ice + q(k,iMa_cor) / rho_dst))**athird
        ! Find the index to which corresponds the optical properties of the core to mantle ratio.
        ! Those properties were determined off-line using the Toon and Ackerman coated spheres code.
        irap = nratio
        do i = 1, nratio
          if (mantletocore < cor_ratio(i) .and. i /= 1) then
            irap = i - 1
            exit
          else if (mantletocore == cor_ratio(i)) then
            irap = i
            exit
          else if (mantletocore < cor_ratio(i) .and. i == 1) then
            irap = i
            exit
          end if
        end do
        ! Calculate the cross-section mean radius (Rs) of the log-normal distribution.
        Mo = q(k,iMa_cld) + q(k,iMa_cor)
        No = q(k,iNb_cld)
        cst = 0.75_r8 / (pi * (q(k,iMa_cld) / rho_ice + q(k,iMa_cor) / rho_dst) / Mo)
        Rs = (Mo / No * cst)**athird * exp(-0.5_r8 * dev_ice**2)
        ! Calculate the total cross sectional area (Ao) of water ice particles.
        Ao = No * pi * Rs**2
        ! Define the cross-section weighted distribution, i.e. surface/size bin. Change Rs to Reff.
        Rs = Rs * exp(1.5_r8 * dev_ice**2)
        Rs = 1.0_r8 / min(max(Rs, 1.0e-7_r8), 100.0e-6_r8)
        do ib = 1, nbin_rt
          surf(ib) = 0.5_r8 * (erf(log(radb_rt(ib+1) * Rs) * dev2) - erf(log(radb_rt(ib) * Rs) * dev2))
        end do
        ! Calculate the averaged values of the optical properties for the whole distribution.
        do l = 2 * k + 2, 2 * k + 3
          do is = 1, nspectv
            qxv(l,is) = 0
            qsv(l,is) = 0
            do ib = 1, nbin_rt
              qxv(l,is) = qxv(l,is) + surf(ib) * qextv_cld (irap,ib,is)
              qsv(l,is) = qsv(l,is) + surf(ib) * qscatv_cld(irap,ib,is)
              gv (l,is) = gv (l,is) + surf(ib) * gv_cld    (irap,ib,is)
            end do
            qsv(l,is) = min(qsv(l,is), 0.99999_r8 * qxv(l,is))
          end do
          do is = 1, nspecti
            qxi(l,is) = 0
            qsi(l,is) = 0
            do ib = 1, nbin_rt
              qxi(l,is) = qxi(l,is) + surf(ib) * qexti_cld (irap,ib,is)
              qsi(l,is) = qsi(l,is) + surf(ib) * qscati_cld(irap,ib,is)
              gi (l,is) = gi (l,is) + surf(ib) * gi_cld    (irap,ib,is)
            end do
            qsi(l,is) = min(qsi(l,is), 0.99999_r8 * qxi(l,is))
          end do
          qextrefcld(l) = qxv(l,nrefv)
          taurefcld (l) = Ao * qextrefcld(l) * (plev(l) - plev(l-1)) / g
          ! For diagnostics: cloud opacity at reference wavelengths (vis and ir).
          taucld(1) = taucld(1) + taurefcld(l)
          taucld(2) = taucld(2) + taurefcld(l) / qextrefcld(l) * (qxi(l,4) - qsi(l,4))
        end do
      end if
    end do

  end subroutine opt_cld

  subroutine optcv(plev, pmid, tmid, qh2o, qxvdst, qsvdst, gvdst, qxvcld, qsvcld, gvcld, qextrefcld, &
                   wbarv, cosbv, dtauv, tauv, taucumv, taugsurf, taurefdst, taurefcld)

    real(r8), intent(in   ) :: plev      (2*nlev+3)
    real(r8), intent(in   ) :: pmid      (2*nlev+3)
    real(r8), intent(in   ) :: tmid      (2*nlev+3)
    real(r8), intent(in   ) :: qh2o      (2*nlev+3)
    real(r8), intent(in   ) :: qxvdst    (2*nlev+4,nspectv)
    real(r8), intent(in   ) :: qsvdst    (2*nlev+4,nspectv)
    real(r8), intent(in   ) :: gvdst     (2*nlev+4,nspectv)
    real(r8), intent(in   ) :: qxvcld    (2*nlev+4,nspectv)
    real(r8), intent(in   ) :: qsvcld    (2*nlev+4,nspectv)
    real(r8), intent(in   ) :: gvcld     (2*nlev+4,nspectv)
    real(r8), intent(in   ) :: qextrefcld(2*nlev+4)
    real(r8), intent(  out) :: wbarv     (nlayrad ,nspectv,ngauss  )
    real(r8), intent(  out) :: cosbv     (nlayrad ,nspectv,ngauss  )
    real(r8), intent(  out) :: dtauv     (nlayrad ,nspectv,ngauss  )
    real(r8), intent(  out) :: tauv      (nlevrad ,nspectv,ngauss  )
    real(r8), intent(  out) :: taucumv   (2*nlev+3,nspectv,ngauss  )
    real(r8), intent(  out) :: taugsurf  (         nspectv,ngauss-1)
    real(r8), intent(inout) :: taurefdst (2*nlev+4)
    real(r8), intent(inout) :: taurefcld (2*nlev+4)

    integer k, l, is, ig
    integer , dimension(  2*nlev+3               ) :: idx_t  
    integer , dimension(  2*nlev+3               ) :: idx_p  
    integer , dimension(  2*nlev+3               ) :: idx_h2o
    real(r8), dimension(  2*nlev+3               ) :: wratio
    real(r8), dimension(4,2*nlev+3               ) :: lcoef 
    real(r8), dimension(  2*nlev+4               ) :: taurefdst_save
    real(r8), dimension(  2*nlev+4               ) :: taurefcld_save
    real(r8), dimension(  2*nlev+3,nspectv       ) :: tray
    real(r8), dimension(  2*nlev+3,nspectv       ) :: tdst
    real(r8), dimension(  2*nlev+3,nspectv       ) :: tcld
    real(r8), dimension(  2*nlev+4,nspectv,ngauss) :: dtaukv
    real(r8) taugas
    real(r8) trayaer ! Tau Rayleigh scattering plus aerosol opacity
    real(r8) kcoef(4)

    taurefdst_save = taurefdst
    taurefcld_save = taurefcld

    ! Determine the total gas opacity throughout the column.
    taugsurf(:,1:ngauss-1) = 0

    do k = 2, 2 * nlev + 3
      call tpindex(pmid(k), tmid(k), qh2o(k), lcoef(:,k), idx_t(k), idx_p(k), idx_h2o(k), wratio(k))
      taurefdst(k) = taurefdst(k) / qxvdst(k,nrefv)
      taurefcld(k) = taurefcld(k) / qextrefcld(k)
      do is = 1, nspectv
        tray(k,is) = tauray(is) * (plev(k) - plev(k-1))
        tdst(k,is) = taurefdst(k) * qxvdst(k,is)
        tcld(k,is) = taurefcld(k) * qxvcld(k,is)
      end do
    end do

    do k = 2, 2 * nlev + 3
      do is = 1, nspectv
        trayaer = tray(k,is) + tdst(k,is) + tcld(k,is)
        do ig = 1, ngauss - 1
          ! Interpolate between water mixing ratios.
          ! wratio is zero if the requested water amount is equal to, or outside
          ! the range of the reference values.
          kcoef(1) = co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                     co2v(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig))
          kcoef(2) = co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                     co2v(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig))
          kcoef(3) = co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                     co2v(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig))
          kcoef(4) = co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                     co2v(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig))
          ! Interpolate the CO2 k-coefficients to the requested temperature and pressure.
          taugas = (lcoef(1,k) * kcoef(1) + lcoef(2,k) * kcoef(2) + &
                    lcoef(3,k) * kcoef(3) + lcoef(4,k) * kcoef(4)) * cmk * (plev(k) - plev(k-1))
          taugsurf(is,ig) = taugas + taugsurf(is,ig)
          dtaukv(k,is,ig) = taugas + trayaer
        end do
        ! Now fill in the "clear" part of the spectrum (ig = ngauss), which holds
        ! continuum opacity only.
        dtaukv(k,is,ngauss) = trayaer
      end do
    end do

    do is = 1, nspectv
      ! First, the special "clear" channel
      ig = ngauss
      do l = 1, nlayrad - 1
        k = 2 * l + 1
        dtauv(l,is,ig) = dtaukv(k,is,ig) + dtaukv(k+1,is,ig)
        cosbv(l,is,ig) = (gvdst(k  ,is) * taurefdst(k  ) * qsvdst(k  ,is)  + &
                          gvdst(k+1,is) * taurefdst(k+1) * qsvdst(k+1,is)  + &
                          gvcld(k  ,is) * taurefcld(k  ) * qsvcld(k  ,is)  + &
                          gvcld(k+1,is) * taurefcld(k+1) * qsvcld(k+1,is)) / &
                         (tray(k,is) + tray(k+1,is) + &
                          taurefdst(k  ) * qsvdst(k  ,is) + &
                          taurefdst(k+1) * qsvdst(k+1,is) + &
                          taurefcld(k  ) * qsvcld(k  ,is) + &
                          taurefcld(k+1) * qsvcld(k+1,is))
        wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefdst(k+1) * qsvdst(k+1,is) + &
                          taurefcld(k) * qsvcld(k,is) + taurefcld(k+1) * qsvcld(k+1,is) + &
                          (tray(k,is) + tray(k+1,is)) * 0.9999_r8) / dtauv(l,is,ig)
      end do
      ! Special bottom layer
      l = nlayrad
      k = 2 * l + 1
      dtauv(l,is,ig) = dtaukv(k,is,ig)
      cosbv(l,is,ig) = (gvdst(k,is) * taurefdst(k) * qsvdst(k,is)  + &
                        gvcld(k,is) * taurefcld(k) * qsvcld(k,is)) / &
                      (tray(k,is) + taurefdst(k) * qsvdst(k,is) + taurefcld(k) * qsvcld(k,is))
      wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefcld(k) * qsvcld(k,is) + &
                        tray(k,is) * 0.9999_r8) / dtauv(l,is,ig)
      ! Now the other Gauss points, if needed
      do ig = 1, ngauss - 1
        do l = 1, nlayrad - 1
          k = 2 * l + 1
          dtauv(l,is,ig) = dtaukv(k,is,ig) + dtaukv(k+1,is,ig)
          cosbv(l,is,ig) = cosbv(l,is,ngauss)
          wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefdst(k+1) * qsvdst(k+1,is) + &
                            taurefcld(k) * qsvcld(k,is) + taurefcld(k+1) * qsvcld(k+1,is) + &
                            (tray(k,is) + tray(k+1,is)) * 0.9999_r8) / dtauv(l,is,ig)
        end do
        ! Special bottom layer
        l = nlayrad
        k = 2 * l + 1
        dtauv(l,is,ig) = dtaukv(k,is,ig)
        cosbv(l,is,ig) = cosbv(l,is,ngauss)
        wbarv(l,is,ig) = (taurefdst(k) * qsvdst(k,is) + taurefcld(k) * qsvcld(k,is) + &
                          tray(k,is) * 0.9999_r8) / dtauv(l,is,ig)
      end do
    end do

    ! Calculate the total extinction optical depths.
    do is = 1, nspectv
      ig = ngauss
      tauv(1,is,ig) = 0
      do l = 1, nlayrad
        tauv(l+1,is,ig) = tauv(l,is,ig) + dtauv(l,is,ig)
      end do
      taucumv(1,is,ig) = 0
      do k = 2, 2 * nlev + 3
        taucumv(k,is,ig) = taucumv(k-1,is,ig) + dtaukv(k,is,ig)
      end do
      do ig = 1, ngauss - 1
        tauv(1,is,ig) = 0
        do l = 1, nlayrad
          tauv(l+1,is,ig) = tauv(l,is,ig) + dtauv(l,is,ig)
        end do
        taucumv(1,is,ig) = 0
        do k = 2, 2 * nlev + 3
          taucumv(k,is,ig) = taucumv(k-1,is,ig) + dtaukv(k,is,ig)
        end do
      end do
    end do

    ! Restore the original values of tauref and taurefcld.
    taurefdst = taurefdst_save
    taurefcld = taurefcld_save

  end subroutine optcv

  subroutine sfluxv(dtauv, tauv, taucumv, taugsurf, cosz, als, wbarv, cosbv, &
                    fluxupv, fluxdnv, fmnetv, nfluxtopv, diffvt, detau)

    real(r8), intent(in ) :: dtauv   (nlayrad ,nspectv,ngauss  )
    real(r8), intent(in ) :: tauv    (nlevrad ,nspectv,ngauss  )
    real(r8), intent(in ) :: taucumv (2*nlev+3,nspectv,ngauss  )
    real(r8), intent(in ) :: taugsurf(         nspectv,ngauss-1)
    real(r8), intent(in ) :: cosz
    real(r8), intent(in ) :: als
    real(r8), intent(in ) :: wbarv   (nlayrad ,nspectv,ngauss  )
    real(r8), intent(in ) :: cosbv   (nlayrad ,nspectv,ngauss  )
    real(r8), intent(out) :: fluxupv (nlayrad)
    real(r8), intent(out) :: fluxdnv (nlayrad)
    real(r8), intent(out) :: fmnetv  (nlayrad)
    real(r8), intent(out) :: nfluxtopv
    real(r8), intent(out) :: diffvt
    real(r8), intent(out) :: detau   (         nspectv,ngauss  )

    integer n, k, l, is, ig
    real(r8) fzero
    real(r8) btop
    real(r8) bsfc
    real(r8) fluxup
    real(r8) fluxdn
    real(r8) fmupv(nlayrad)
    real(r8) fmdnv(nlayrad)
    real(r8) diffv

    n = 2 * nlev + 3

    fluxupv   = 0
    fluxdnv   = 0
    fmnetv    = 0
    nfluxtopv = 0
    diffvt    = 0

    ! Calculate the net flux in each spectral interval.
    do is = 1, nspectv
      fzero = fzerov(is)
      if (fzerov(is) < 0.99_r8) then
        do ig = 1, ngauss - 1
          call getdetau(      &
            dtauv  (1,is,ig), &
            tauv   (1,is,ig), &
            taucumv(1,is,ig), &
            wbarv  (1,is,ig), &
            cosbv  (1,is,ig), &
            detau  (  is,ig))
          if (taugsurf(is,ig) < tlimits) then
            fzero = fzero + (1 - fzerov(is)) * gweight(ig)
          else
            ! Set up the top and bottom boundary conditions.
            btop = 0
            bsfc = als * cosz * solar(is) * exp(-min(detau(is,ig) / cosz, maxexp))
            ! We can now solve for the coefficients of the two-stream problem.
            call gfluxv(        &
              solar  (  is   ), &
              dtauv  (1,is,ig), &
              tauv   (1,is,ig), &
              taucumv(1,is,ig), &
              wbarv  (1,is,ig), &
              cosbv  (1,is,ig), &
              cosz            , &
              als             , &
              btop            , &
              bsfc            , &
              fmupv           , &
              fmdnv           , &
              diffv           , &
              fluxup          , &
              fluxdn          , &
              detau    (is,ig))
            ! Calculate the cumulative visible net flux.
            nfluxtopv = nfluxtopv + (fluxup - fluxdn) * gweight(ig) * (1 - fzerov(is))
            do l = 1, nlayrad
              fmnetv (l) = fmnetv (l) + (fmupv(l) - fmdnv(l)) * gweight(ig) * (1 - fzerov(is))
              fluxupv(l) = fluxupv(l) + fmupv(l) * gweight(ig) * (1 - fzerov(is))
              fluxdnv(l) = fluxdnv(l) + fmdnv(l) * gweight(ig) * (1 - fzerov(is))
            end do
            ! The diffuse component of the downward solar flux.
            diffvt = diffvt + diffv * gweight(ig) * (1 - fzerov(is))
          end if
        end do
      end if
      ! Special 17th Gauss point
      ig = ngauss
      call getdetau(      &
        dtauv  (1,is,ig), &
        tauv   (1,is,ig), &
        taucumv(1,is,ig), &
        wbarv  (1,is,ig), &
        cosbv  (1,is,ig), &
        detau  (  is,ig))
      btop = 0
      bsfc = als * cosz * solar(is) * exp(-min(detau(is,ig) / cosz, maxexp))
      call gfluxv(        &
        solar  (  is   ), &
        dtauv  (1,is,ig), &
        tauv   (1,is,ig), &
        taucumv(1,is,ig), &
        wbarv  (1,is,ig), &
        cosbv  (1,is,ig), &
        cosz            , &
        als             , &
        btop            , &
        bsfc            , &
        fmupv           , &
        fmdnv           , &
        diffv           , &
        fluxup          , &
        fluxdn          , &
        detau    (is,ig))
      nfluxtopv = nfluxtopv + (fluxup - fluxdn) * fzero
      do l = 1, nlayrad
        fmnetv (l) = fmnetv (l) + (fmupv(l) - fmdnv(l)) * fzero
        fluxupv(l) = fluxupv(l) + fmupv(l) * fzero
        fluxdnv(l) = fluxdnv(l) + fmdnv(l) * fzero
      end do
      diffvt = diffvt + diffv * fzero
    end do

  end subroutine sfluxv

  subroutine getdetau(dtdel, tdel, taucumin, wdel, cdel, detau)

    real(r8), intent(in ) :: dtdel   (nlayrad)
    real(r8), intent(in ) :: tdel    (nlevrad)
    real(r8), intent(in ) :: taucumin(2*nlev+3)
    real(r8), intent(in ) :: wdel    (nlayrad)
    real(r8), intent(in ) :: cdel    (nlayrad)
    real(r8), intent(out) :: detau

    integer k, l
    real(r8) factor
    real(r8) tau   (nlevrad)
    real(r8) taucum(2*nlev+3)

    ! Delta-Eddington scaling
    factor    = 1 - wdel(1) * cdel(1)**2
    tau   (1) = tdel(1) * factor
    taucum(1) = 0
    taucum(2) = taucumin(2) * factor
    taucum(3) = taucum(2) + (taucumin(3) - taucumin(2)) * factor

    do l = 1, nlayrad - 1
      factor      = 1 - wdel(l) * cdel(l)**2
      tau   (l+1) = tau(l) + dtdel(l) * factor
      k           = 2 * l + 2
      taucum(k  ) = tau(l+1)
      taucum(k+1) = taucum(k) + (taucumin(k+1) - taucumin(k)) * factor
    end do

    l         = nlayrad
    factor    = 1 - wdel(l) * cdel(l)**2
    tau(l+1)  = tau(l) + dtdel(l) * factor
    k         = 2 * l + 1
    taucum(k) = tau(l+1)
    detau     = taucum(k)

  end subroutine getdetau

  subroutine gfluxv(sol, dtdel, tdel, taucumin, wdel, cdel, cosz, als, &
                    btop, bsfc, fmidp, fmidm, diffv, fluxup, fluxdn, detau)

    real(r8), intent(in ) :: sol
    real(r8), intent(in ) :: dtdel    (nlayrad)
    real(r8), intent(in ) :: tdel     (nlevrad)
    real(r8), intent(in ) :: taucumin (2*nlev+3)
    real(r8), intent(in ) :: wdel     (nlayrad)
    real(r8), intent(in ) :: cdel     (nlayrad)
    real(r8), intent(in ) :: cosz
    real(r8), intent(in ) :: als
    real(r8), intent(in ) :: btop
    real(r8), intent(in ) :: bsfc
    real(r8), intent(out) :: fmidp    (nlayrad)
    real(r8), intent(out) :: fmidm    (nlayrad)
    real(r8), intent(out) :: diffv
    real(r8), intent(out) :: fluxup
    real(r8), intent(out) :: fluxdn
    real(r8), intent(out) :: detau

    integer , parameter :: nlp = 101 ! Must be larger than 2*nlev+3
    integer k, l
    real(r8) factor
    real(r8) tau   (nlevrad)
    real(r8) w0    (nlayrad)
    real(r8) cosbar(nlayrad)
    real(r8) dtau  (nlayrad)
    real(r8) taucum(2*nlev+3)
    real(r8) alpha (nlp)
    real(r8) lambda(nlp)
    real(r8) gamma (nlp)
    real(r8) g1    (nlp)
    real(r8) g2    (nlp)
    real(r8) g3    (nlp)
    real(r8) g4
    real(r8) denom
    real(r8) am
    real(r8) ap
    real(r8) cpm1  (nlp)
    real(r8) cmm1  (nlp)
    real(r8) cp    (nlp)
    real(r8) cm    (nlp)
    real(r8) ep
    real(r8) e1    (nlp)
    real(r8) e2    (nlp)
    real(r8) e3    (nlp)
    real(r8) e4    (nlp)
    real(r8) x1    (nlp)
    real(r8) x2    (nlp)
    real(r8) taumid
    real(r8) cpmid
    real(r8) cmmid

    ! Delta-Eddington scaling
    factor    = 1 - wdel(1) * cdel(1)**2
    tau   (1) = tdel(1) * factor
    taucum(1) = 0
    taucum(2) = taucumin(2) * factor
    taucum(3) = taucum(2) + (taucumin(3) - taucumin(2)) * factor

    do l = 1, nlayrad - 1
      factor      = 1 - wdel(l) * cdel(l)**2
      w0    (l  ) = wdel(l) * (1 - cdel(l)**2) / factor
      cosbar(l  ) = cdel(l) / (1 + cdel(l))
      dtau  (l  ) = dtdel(l) * factor
      tau   (l+1) = tau(l) + dtau(l)
      k           = 2 * l + 2
      taucum(k  ) = tau(l+1)
      taucum(k+1) = taucum(k) + (taucumin(k+1) - taucumin(k)) * factor
    end do

    l           = nlayrad
    factor      = 1 - wdel(l) * cdel(l)**2
    w0    (l  ) = wdel(l) * (1 - cdel(l)**2) / factor
    cosbar(l  ) = cdel(l) / (1 + cdel(l))
    dtau  (l  ) = dtdel(l) * factor
    tau   (l+1) = tau(l) + dtau(l)
    k           = 2 * l + 1
    taucum(k  ) = tau(l+1)
    detau       = taucum(k)

    do l = 1, nlayrad
      alpha (l) = sqrt((1 - w0(l)) / (1 - w0(l) * cosbar(l)))
      g1    (l) = (sqrt3 * 0.5_r8) * (2 - w0(l) * (1 + cosbar(l)))
      g2    (l) = (sqrt3 * 0.5_r8 * w0(l)) * (1 - cosbar(l))
      g3    (l) = 0.5_r8 * (1 - sqrt3 * cosbar(l) * cosz)
      lambda(l) = sqrt(g1(l)**2 - g2(l)**2)
      gamma (l) = (g1(l) - lambda(l)) / g2(l)
    end do

    do l = 1, nlayrad
      g4    = 1 - g3(l)
      denom = lambda(l)**2 - 1.0_r8 / cosz**2
      ! There is a potential problem if w0=0 and ubarv=cosz, then denom will
      ! vanish. This only happens physically when the scattering goes to zero.
      ! Prevent this with an if statement.
      if (denom == 0) denom = 1.0e-10_r8

      am = sol * w0(l) * (g4    * (g1(l) + 1.0_r8 / cosz) + g2(l) * g3(l)) / denom
      ap = sol * w0(l) * (g3(l) * (g1(l) - 1.0_r8 / cosz) + g2(l) * g4   ) / denom

      factor  = exp(-min(tau(l) / cosz, maxexp))
      cpm1(l) = ap * factor
      cmm1(l) = am * factor

      factor = exp(-min(tau(l+1) / cosz, maxexp))
      cp(l)  = ap * factor
      cm(l)  = am * factor
    end do

    ! Calculate the exponential terms for the tridiaonal rotated layered method.
    do l = 1, nlayrad
      ep = exp(min(taumax, lambda(l) * dtau(l))) ! Clipped exponential
      e1(l) = ep + gamma(l) / ep
      e2(l) = ep - gamma(l) / ep
      e3(l) = gamma(l) * ep + 1.0_r8 / ep
      e4(l) = gamma(l) * ep - 1.0_r8 / ep
    end do

    call dsolver(nlayrad, gamma, cp, cm, cpm1, cmm1, e1, e2, e3, e4, btop, bsfc, als, x1, x2)

    do l = 1, nlayrad - 1
      ep = exp(min(taumax, lambda(l) * (taucum(2*l+1) - taucum(2*l))))
      g4 = 1 - g3(l)
      denom = lambda(l)**2 - 1.0_r8 / cosz**2
      if (denom == 0) denom = 1.0e-10_r8
      am = sol * w0(l) * (g4    * (g1(l) + 1.0_r8 / cosz) + g2(l) * g3(l)) / denom
      ap = sol * w0(l) * (g3(l) * (g1(l) - 1.0_r8 / cosz) + g2(l) * g4   ) / denom

      taumid = taucum(2*l+1)
      cpmid = ap * exp(-min(taumid / cosz, maxexp))
      cmmid = am * exp(-min(taumid / cosz, maxexp))

      fmidp(l) = x1(l) * ep + gamma(l) * x2(l) / ep + cpmid
      fmidm(l) = x1(l) * ep * gamma(l) + x2(l) / ep + cmmid

      ! Add the direct flux to the downwelling term.
      fmidm(l) = fmidm(l) + cosz * sol * exp(-min(taumid / cosz, maxexp))
    end do

    ! Flux at the top layer
    ep = 1
    g4 = 1 - g3(1)
    denom = lambda(1)**2 - 1.0_r8 / cosz**2
    if (denom == 0) denom = 1.0e-10_r8
    am = sol * w0(1) * (g4    * (g1(1) + 1.0_r8 / cosz) + g2(1) * g3(1)) / denom
    ap = sol * w0(1) * (g3(1) * (g1(1) - 1.0_r8 / cosz) + g2(1) * g4   ) / denom

    cpmid = ap
    cmmid = am

    fluxup = x1(1) * ep + gamma(1) * x2(1) / ep + cpmid
    fluxdn = x1(1) * ep * gamma(1) + x2(1) / ep + cmmid

    ! Add the direct flux to the downwelling term.
    fluxdn = fluxdn + cosz * sol * exp(-min(taucum(1) / cosz, maxexp))

    ! This is for the "special" bottom layer, where we take dtau instead of dtau/2.
    l = nlayrad
    ep = exp(min(taumax, lambda(l) * (taucum(2*nlev+3) - taucum(2*nlev+2))))
    g4 = 1 - g3(l)
    denom = lambda(l)**2 - 1.0_r8 / cosz**2
    if (denom == 0) denom = 1.0e-10_r8
    am = sol * w0(l) * (g4    * (g1(l) + 1.0_r8 / cosz) + g2(l) * g3(l)) / denom
    ap = sol * w0(l) * (g3(l) * (g1(l) - 1.0_r8 / cosz) + g2(l) * g4   ) / denom

    taumid = min(taucum(2*nlev+3), taumax)
    cpmid = ap * exp(-min(taumid / cosz, maxexp))
    cmmid = am * exp(-min(taumid / cosz, maxexp))

    fmidp(l) = x1(l) * ep + gamma(l) * x2(l) / ep + cpmid
    fmidm(l) = x1(l) * ep * gamma(l) + x2(l) / ep + cmmid

    diffv = fmidm(l)

    fmidm(l) = fmidm(l) + cosz * sol * exp(-min(taumid / cosz, maxexp))

  end subroutine gfluxv

  subroutine optci(plev, pmid, tmid, qh2o, qxidst, qsidst, gidst, qextrefdst, qxicld, qsicld, gicld, qextrefcld, &
                   wbari, cosbi, dtaui, taui, taucumi, taugsurf, taurefdst, taurefcld)

    real(r8), intent(in   ) :: plev      (2*nlev+3)
    real(r8), intent(in   ) :: pmid      (2*nlev+3)
    real(r8), intent(in   ) :: tmid      (2*nlev+3)
    real(r8), intent(in   ) :: qh2o      (2*nlev+3)
    real(r8), intent(in   ) :: qxidst    (2*nlev+4,nspecti)
    real(r8), intent(in   ) :: qsidst    (2*nlev+4,nspecti)
    real(r8), intent(in   ) :: gidst     (2*nlev+4,nspecti)
    real(r8), intent(in   ) :: qextrefdst(2*nlev+4)
    real(r8), intent(in   ) :: qxicld    (2*nlev+4,nspecti)
    real(r8), intent(in   ) :: qsicld    (2*nlev+4,nspecti)
    real(r8), intent(in   ) :: gicld     (2*nlev+4,nspecti)
    real(r8), intent(in   ) :: qextrefcld(2*nlev+4)
    real(r8), intent(  out) :: wbari     (nlayrad ,nspecti,ngauss)
    real(r8), intent(  out) :: cosbi     (nlayrad ,nspecti,ngauss)
    real(r8), intent(  out) :: dtaui     (nlayrad ,nspecti,ngauss)
    real(r8), intent(  out) :: taui      (nlevrad ,nspecti,ngauss)
    real(r8), intent(  out) :: taucumi   (2*nlev+3,nspecti,ngauss)
    real(r8), intent(  out) :: taugsurf  (         nspecti,ngauss-1)
    real(r8), intent(inout) :: taurefdst (2*nlev+4)
    real(r8), intent(inout) :: taurefcld (2*nlev+4)

    integer k, l, is, ig
    integer idx_t          (2*nlev+3)
    integer idx_p          (2*nlev+3)
    integer idx_h2o        (2*nlev+3)
    real(r8) dtauki        (2*nlev+4,nspecti,ngauss)
    real(r8) taudstk       (2*nlev+4,nspecti)
    real(r8) taucldk       (2*nlev+4,nspecti)
    real(r8) taurefdst_save(2*nlev+4)
    real(r8) taurefcld_save(2*nlev+4)
    real(r8) wratio        (2*nlev+3)
    real(r8) lcoef       (4,2*nlev+3)
    real(r8) taudst        (2*nlev+3,nspecti)
    real(r8) taucld        (2*nlev+3,nspecti)
    real(r8) ans
    real(r8) kcoef(4)
    real(r8) tauac

    do is = 1, nspecti
      do ig = 1, ngauss
        dtauki(2*nlev+4,is,ig) = 0
      end do
      taudstk(2*nlev+4,is) = 0
      taucldk(2*nlev+4,is) = 0
    end do

    taurefdst_save = taurefdst
    taurefcld_save = taurefcld

    do ig = 1, ngauss - 1
      do is = 1, nspecti
        taugsurf(is,ig) = 0
      end do
    end do

    do k = 2, 2 * nlev + 3
      call tpindex(pmid(k), tmid(k), qh2o(k), lcoef(:,k), idx_t(k), idx_p(k), idx_h2o(k), wratio(k))
      taurefdst(k) = taurefdst(k) / qextrefdst(k)
      taurefcld(k) = taurefcld(k) / qextrefcld(k)
      do is = 1, nspecti
        taudst(k,is) = taurefdst(k) * qxidst(k,is)
        taucld(k,is) = taurefcld(k) * qxicld(k,is)
      end do
    end do

    do k = 2, 2 * nlev + 3
      do is = 1, nspecti
        do ig = 1, ngauss - 1
          kcoef(1) = co2i(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2i(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                     co2i(idx_t(k)  ,idx_p(k)  ,idx_h2o(k)  ,is,ig))
          kcoef(2) = co2i(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2i(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                     co2i(idx_t(k)  ,idx_p(k)+1,idx_h2o(k)  ,is,ig))
          kcoef(3) = co2i(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2i(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)+1,is,ig) - &
                     co2i(idx_t(k)+1,idx_p(k)+1,idx_h2o(k)  ,is,ig))
          kcoef(4) = co2i(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig) + wratio(k) * &
                    (co2i(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)+1,is,ig) - &
                     co2i(idx_t(k)+1,idx_p(k)  ,idx_h2o(k)  ,is,ig))
          ! Interpolate the CO2 k-coefficients to the requested temperature and pressure.
          ans = (lcoef(1,k) * kcoef(1) + lcoef(2,k) * kcoef(2) + &
                 lcoef(3,k) * kcoef(3) + lcoef(4,k) * kcoef(4)) * cmk * (plev(k) - plev(k-1))
          taugsurf(is,ig) = taugsurf(is,ig) + ans
          dtauki(k,is,ig) = taudst(k,is) + taucld(k,is) + ans
        end do
        ! Now fill in the "clear" part of the spectrum (ig = ngauss), which holds
        ! continuum opacity only.
        dtauki(k,is,ngauss) = taudst(k,is) + taucld(k,is)
      end do
    end do

    do is = 1, nspecti
      do k = 2, 2 * nlev + 4
        taudstk(k,is) = taurefdst(k) * qsidst(k,is)
        taucldk(k,is) = taurefcld(k) * qsicld(k,is)
      end do
    end do

    do is = 1, nspecti
      ig = ngauss
      do l = 1, nlayrad
        k = 2 * l + 1
        dtaui(l,is,ig) = dtauki(k,is,ig) + dtauki(k+1,is,ig) + 1.0d-50
        if (dtaui(l,is,ig) > 1.0e-9_r8) then
          wbari(l,is,ig) = (taudstk(k,is) + taudstk(k+1,is) + taucldk(k,is) + taucldk(k+1,is)) / dtaui(l,is,ig)
        else
          wbari(l,is,ig) = 0
          dtaui(l,is,ig) = 1.0e-9_r8
        end if
        tauac = taudstk(k,is) + taudstk(k+1,is) + taucldk(k,is) + taucldk(k+1,is)
        if (tauac > 0) then
          cosbi(l,is,ig) = (gidst(k,is) * taudstk(k,is) + gidst(k+1,is) * taudstk(k+1,is)  + &
                            gicld(k,is) * taucldk(k,is) + gicld(k+1,is) * taucldk(k+1,is)) / &
                           (taudstk(k,is) + taudstk(k+1,is) + taucldk(k,is) + taucldk(k+1,is))
        else
          cosbi(l,is,ig) = 0
        end if
      end do
      do ig = 1, ngauss - 1
        do l = 1, nlayrad
          k = 2 * l + 1
          dtaui(l,is,ig) = dtauki(k,is,ig) + dtauki(k+1,is,ig) + 1.0d-50
          if (dtaui(l,is,ig) > 1.0e-9_r8) then
            wbari(l,is,ig) = (taudstk(k,is) + taudstk(k+1,is) + taucldk(k,is) + taucldk(k+1,is)) / dtaui(l,is,ig)
          else
            wbari(l,is,ig) = 0
            dtaui(l,is,ig) = 1.0e-9_r8
          end if
          cosbi(l,is,ig) = cosbi(l,is,ngauss)
        end do
      end do
    end do

    ! Total extinction optical depths
    do is = 1, nspecti
      ig = ngauss
      taui(1,is,ig) = 0
      do l = 1, nlayrad
        taui(l+1,is,ig) = taui(l,is,ig) + dtaui(l,is,ig)
      end do
      taucumi(1,is,ig) = 0
      do k = 2, 2 * nlev + 3
        taucumi(k,is,ig) = taucumi(k-1,is,ig) + dtauki(k,is,ig)
      end do
      do ig = 1, ngauss - 1
        taui(1,is,ig) = 0
        do l = 1, nlayrad
          taui(l+1,is,ig) = taui(l,is,ig) + dtaui(l,is,ig)
        end do
        taucumi(1,is,ig) = 0
        do k = 2, 2 * nlev + 3
          taucumi(k,is,ig) = taucumi(k-1,is,ig) + dtauki(k,is,ig)
        end do
      end do
    end do

    taurefdst = taurefdst_save
    taurefcld = taurefcld_save

  end subroutine optci

  subroutine sfluxi(plev, tlev, dtaui, taucumi, taugsurf, albi, wbari, cosbi, fluxupi, fluxdni, fmneti, nfluxtopi)

    real(r8), intent(in ) :: plev     (2*nlev+3)
    real(r8), intent(in ) :: tlev     (2*nlev+3)
    real(r8), intent(in ) :: dtaui    (nlayrad ,nspecti,ngauss  )
    real(r8), intent(in ) :: taucumi  (2*nlev+3,nspecti,ngauss  )
    real(r8), intent(in ) :: taugsurf (         nspecti,ngauss-1)
    real(r8), intent(in ) :: albi
    real(r8), intent(in ) :: wbari    (nlayrad ,nspecti,ngauss  )
    real(r8), intent(in ) :: cosbi    (nlayrad ,nspecti,ngauss  )
    real(r8), intent(out) :: fluxupi  (nlayrad)
    real(r8), intent(out) :: fluxdni  (nlayrad)
    real(r8), intent(out) :: fmneti   (nlayrad)
    real(r8), intent(out) :: nfluxtopi

    integer k, l, is, ig, nts, ntt
    real(r8) ttop
    real(r8) tsfc
    real(r8) btop
    real(r8) bsfc
    real(r8) fzero
    real(r8) ftopup
    real(r8) fmupi(nlevrad)
    real(r8) fmdni(nlevrad)

    ! Zero fluxes.
    nfluxtopi = 0
    fluxupi   = 0
    fluxdni   = 0
    fmneti    = 0

    ttop = tlev(2)
    tsfc = tlev(2*nlev+3)

    ntt = ttop * 10 - 499
    nts = tsfc * 10 - 499

    do is = 1, nspecti
      bsfc  = (1 - albi) * planckir(is,nts)
      fzero = fzeroi(is)
      if (fzero < 0.99) then
        do ig = 1, ngauss - 1
          if (taugsurf(is,ig) < tlimiti) then
            fzero = fzero + (1 - fzeroi(is)) * gweight(ig)
          else
            btop = (1 - exp(-dtaui(1,is,ig) * plev(2) / (plev(4) - plev(2)) / ubari)) * planckir(is,ntt)
            call gfluxi(        &
              is              , &
              tlev            , &
              dtaui  (:,is,ig), &
              taucumi(:,is,ig), &
              wbari  (:,is,ig), &
              cosbi  (:,is,ig), &
              albi            , &
              btop            , &
              bsfc            , &
              ftopup          , &
              fmupi           , &
              fmdni           )
            nfluxtopi = nfluxtopi + ftopup * dwni(is) * gweight(ig) * (1 - fzeroi(is))
            do l = 1, nlevrad - 1
              fluxupi(l) = fluxupi(l) + fmupi(l) * dwni(is) * gweight(ig) * (1 - fzeroi(is))
              fluxdni(l) = fluxdni(l) + fmdni(l) * dwni(is) * gweight(ig) * (1 - fzeroi(is))
              fmneti (l) = fmneti (l) + (fmupi(l) - fmdni(l)) * dwni(is) * gweight(ig) * (1 - fzeroi(is))
            end do
          end if
        end do
      end if
      ! Special 17th Gauss-point
      ig   = ngauss
      btop = (1 - exp(-dtaui(1,is,ig) * plev(2) / (plev(4) - plev(2)) / ubari)) * planckir(is,ntt)
      call gfluxi(        &
        is              , &
        tlev            , &
        dtaui  (:,is,ig), &
        taucumi(:,is,ig), &
        wbari  (:,is,ig), &
        cosbi  (:,is,ig), &
        albi            , &
        btop            , &
        bsfc            , &
        ftopup          , &
        fmupi           , &
        fmdni           )
      nfluxtopi = nfluxtopi + ftopup * dwni(is) * fzero
      do l = 1, nlevrad - 1
        fluxupi(l) = fluxupi(l) + fmupi(l) * dwni(is) * fzero
        fluxdni(l) = fluxdni(l) + fmdni(l) * dwni(is) * fzero
        fmneti (l) = fmneti (l) + (fmupi(l) - fmdni(l)) * dwni(is) * fzero
      end do
    end do

  end subroutine sfluxi

  subroutine gfluxi(is, tlev, dtau, taucum, wbar, cosb, albi, btop, bsfc, ftopup, fmup, fmdn)

    integer , intent(in ) :: is
    real(r8), intent(in ) :: tlev  (2*nlev+3)
    real(r8), intent(in ) :: dtau  (nlayrad )
    real(r8), intent(in ) :: taucum(2*nlev+3)
    real(r8), intent(in ) :: wbar  (nlayrad )
    real(r8), intent(in ) :: cosb  (nlayrad )
    real(r8), intent(in ) :: albi
    real(r8), intent(in ) :: btop
    real(r8), intent(in ) :: bsfc
    real(r8), intent(out) :: ftopup
    real(r8), intent(out) :: fmup  (nlayrad)
    real(r8), intent(out) :: fmdn  (nlayrad)

    integer, parameter :: nlp = 101 ! Must be larger than 2*nlev+3
    integer l, nt, nt2
    real(r8) term, ep, em, dtauk, cpmid, cmmid, fluxup, fluxdn
    real(r8) alpha (nlp)
    real(r8) lambda(nlp)
    real(r8) gamma (nlp)
    real(r8) b0    (nlp)
    real(r8) b1    (nlp)
    real(r8) cp    (nlp)
    real(r8) cm    (nlp)
    real(r8) cpm1  (nlp)
    real(r8) cmm1  (nlp)
    real(r8) e1    (nlp)
    real(r8) e2    (nlp)
    real(r8) e3    (nlp)
    real(r8) e4    (nlp)
    real(r8) x1    (nlp)
    real(r8) x2    (nlp)

    do l = 1, nlayrad - 1
      alpha (l) = sqrt((1 - wbar(l)) / (1 - wbar(l) * cosb(l)))
      lambda(l) = alpha(l) * (1 - wbar(l) * cosb(l)) / ubari
      nt2       = tlev(2*l+2) * 10 - 499
      nt        = tlev(2*l  ) * 10 - 499
      b1    (l) = (planckir(is,nt2) - planckir(is,nt)) / dtau(l)
      b0    (l) = planckir(is,nt)
    end do

    l         = nlayrad
    alpha (l) = sqrt((1 - wbar(l)) / (1 - wbar(l) * cosb(l)))
    lambda(l) = alpha(l) * (1 - wbar(l) * cosb(l)) / ubari
    nt        = tlev(2*l+1) * 10 - 499
    nt2       = tlev(2*l  ) * 10 - 499
    b1    (l) = (planckir(is,nt) - planckir(is,nt2)) / dtau(l)
    b0    (l) = planckir(is,nt2)

    do l = 1, nlayrad
      gamma(l) = (1 - alpha(l)) / (1 + alpha(l))
      term    = ubari / (1 - wbar(l) * cosb(l))
      cp  (l) = b0(l) + b1(l) * dtau(l) + b1(l) * term
      cm  (l) = b0(l) + b1(l) * dtau(l) - b1(l) * term
      cpm1(l) = b0(l) + b1(l) * term
      cmm1(l) = b0(l) - b1(l) * term
    end do

    do l = 1, nlayrad
      ep    = exp(min(lambda(l) * dtau(l), maxexp))
      em    = 1.0_r8 / ep
      e1(l) = ep + gamma(l) * em
      e2(l) = ep - gamma(l) * em
      e3(l) = gamma(l) * ep + em
      e4(l) = gamma(l) * ep - em
    end do

    call dsolver(nlayrad, gamma, cp, cm, cpm1, cmm1, e1, e2, e3, e4, btop, bsfc, albi, x1, x2)

    do l = 1, nlayrad - 1
      dtauk   = taucum(2*l+1) - taucum(2*l)
      ep      = exp(min(lambda(l) * dtauk, maxexp))
      em      = 1.0_r8 / ep
      term    = ubari / (1 - wbar(l) * cosb(l))
      cpmid   = b0(l) + b1(l) * dtauk + b1(l) * term
      cmmid   = b0(l) + b1(l) * dtauk - b1(l) * term
      fmup(l) = x1(l) * ep + gamma(l) * x2(l) * em + cpmid
      fmdn(l) = x1(l) * ep * gamma(l) + x2(l) * em + cmmid
      ! For flux, integrate over the hemisphere treating intensity constant.
      fmup(l) = fmup(l) * pi
      fmdn(l) = fmdn(l) * pi
    end do

    l       = nlayrad
    ep      = exp(min(lambda(l) * dtau(l), taumax))
    em      = 1.0_r8 / ep
    term    = ubari / (1 - wbar(l) * cosb(l))
    cpmid   = b0(l) + b1(l) * dtau(l) + b1(l) * term
    cmmid   = b0(l) + b1(l) * dtau(l) - b1(l) * term
    fmup(l) = x1(l) * ep + gamma(l) * x2(l) * em + cpmid
    fmdn(l) = x1(l) * ep * gamma(l) + x2(l) * em + cmmid
    fmup(l) = fmup(l) * pi
    fmdn(l) = fmdn(l) * pi

    ep     = 1
    em     = 1
    term   = ubari / (1 - wbar(1) * cosb(1))
    cpmid  = b0(1) + b1(1) * term
    cmmid  = b0(1) - b1(1) * term
    fluxup = x1(1) * ep + gamma(1) * x2(1) * em + cpmid
    fluxdn = x1(1) * ep * gamma(1) + x2(1) * em + cmmid

    ftopup = (fluxup - fluxdn) * pi

  end subroutine gfluxi

  subroutine dsolver(nl, gamma, cp, cm, cpm1, cmm1, e1, e2, e3, e4, btop, bsfc, als, x1, x2)

    integer, parameter :: nmax = 201

    integer , intent(in ) :: nl
    real(r8), intent(in ) :: gamma(nl)
    real(r8), intent(in ) :: cp   (nl)
    real(r8), intent(in ) :: cm   (nl)
    real(r8), intent(in ) :: cpm1 (nl)
    real(r8), intent(in ) :: cmm1 (nl)
    real(r8), intent(in ) :: e1   (nl)
    real(r8), intent(in ) :: e2   (nl)
    real(r8), intent(in ) :: e3   (nl)
    real(r8), intent(in ) :: e4   (nl)
    real(r8), intent(in ) :: btop
    real(r8), intent(in ) :: bsfc
    real(r8), intent(in ) :: als
    real(r8), intent(out) :: x1   (nl)
    real(r8), intent(out) :: x2   (nl)

    integer i, l

    real(r8) a(nmax), b(nmax), c(nmax), d(nmax), x(nmax)

    a(1) = 0
    b(1) = gamma(1) + 1
    c(1) = gamma(1) - 1
    d(1) = btop - cmm1(1)

    l = 0
    do i = 2, 2 * nl - 2, 2
      l    = l + 1
      a(i) = (e1(l) + e3(l)) * (gamma(l+1) - 1)
      b(i) = (e2(l) + e4(l)) * (gamma(l+1) - 1)
      c(i) = 2 * (1 - gamma(l+1)**2)
      d(i) = (gamma(l+1) - 1) * (cpm1(l+1) - cp(l)) + (1 - gamma(l+1)) * (cm(l) - cmm1(l+1))
    end do

    l = 0
    do i = 3, 2 * nl - 1, 2
      l    = l + 1
      a(i) = 2 * (1 - gamma(l)**2)
      b(i) = (e1(l) - e3(l)) * (gamma(l+1) + 1)
      c(i) = (e1(l) + e3(l)) * (gamma(l+1) - 1)
      d(i) = e3(l) * (cpm1(l+1) - cp(l)) + e1(l) * (cm(l) - cmm1(l+1))
    end do

    a(2*nl) = e1(nl) - als * e3(nl)
    b(2*nl) = e2(nl) - als * e4(nl)
    c(2*nl) = 0
    d(2*nl) = bsfc - cp(nl) + als * cm(nl)

    call tridiag_thomas(a(:2*nl), b(:2*nl), c(:2*nl), d(:2*nl), x(:2*nl))

    do l = 1, nl
      i = 2 * l
      x1(l) = x(i-1) + x(i)
      x2(l) = x(i-1) - x(i)

      if (x2(l) /= 0) then
        if (abs(x2(l) / (x(i-1) + 1.0e-20_r8)) < 1.0e-30) x2(l) = 0
      end if
    end do

  end subroutine dsolver

  subroutine tpindex(p, t, qh2o, coef, idx_t, idx_p, idx_h2o, wratio)

    ! Interpolate the CO2 K-coefficients to the current P, T values.

    real(r8), intent(in ) :: p          ! FIXME: Check units.
    real(r8), intent(in ) :: t
    real(r8), intent(in ) :: qh2o
    real(r8), intent(out) :: coef(4)
    integer , intent(out) :: idx_t      ! Temperature-grid index
    integer , intent(out) :: idx_p      ! Pressure-grid index
    integer , intent(out) :: idx_h2o    ! Water abundance index
    real(r8), intent(out) :: wratio     ! Water abundance ratio

    integer i
    real(r8) ct, cp, plog

    ! Get the upper and lower Temperature-grid indicies that bound the
    ! requested temperature.  If the requested temperature is outside
    ! the T-grid, set up to extrapolate from the appropriate end.
    if (t <= tgasref(1)) then
      idx_t = 1
    else
      idx_t = 0
      do i = 1, ntref - 1
        if (t > tgasref(i) .and. t <= tgasref(i+1)) then
          idx_t = i
          exit
        end if
      end do
      if (idx_t == 0) then
        idx_t = ntref - 1
      end if
    end if
    ct = (t - tgasref(idx_t)) / (tgasref(idx_t+1) - tgasref(idx_t))

    ! Get the upper and lower Pressure-grid indicies that bound the
    ! requested pressure.  If the requested pressure is outside
    ! the P-grid, set up to extrapolate from the appropiate end.
    plog = log10(p / 100) ! hPa
    idx_p = 0
    do i = 2, npint - 1
      if (plog <= pfgasref(i)) then
        idx_p = i - 1
        exit
      end if
    end do
    if (idx_p == 0) then
      idx_p = npint - 1
    end if
    cp = (plog - pfgasref(idx_p)) / (pfgasref(idx_p+1) - pfgasref(idx_p))

    ! Fill the interpolation coefficients.
    coef(1) = (1 - cp) * (1 - ct)
    coef(2) = cp * (1 - ct)
    coef(3) = cp * ct
    coef(4) = (1 - cp) * ct

    ! Get the indices for water abundance. There are 10 sets of k-coefficients
    ! with differing amounts of water vs CO2.
    if (qh2o <= wrefh2o(1)) then
      idx_h2o = 1
      wratio = 0
    else if (qh2o >= wrefh2o(nrefh2o)) then
      idx_h2o = nrefh2o
      wratio = 0
    else
      do i = 2, nrefh2o
        if (qh2o >= wrefh2o(i-1) .and. qh2o < wrefh2o(i)) then
          idx_h2o = i - 1
          wratio = (qh2o - wrefh2o(i-1)) / (wrefh2o(i) - wrefh2o(i-1))
          exit
        end if
      end do
    end if

  end subroutine tpindex

  subroutine lagrange(x, xi, yi, ans)

    !  Lagrange interpolation - Polynomial interpolation at point x 
    !  xi(1) <= x <= xi(4).  Yi(n) is the functional value at XI(n).
    
    real(8) x, xi(4), yi(4), ans
    real(8) fm1, fm2, fm3, fm4
    
    fm1 = x - xi(1)
    fm2 = x - xi(2)
    fm3 = x - xi(3)
    fm4 = x - xi(4)
    
    !  Get the "answer" at the requested X
    ans = fm2 * fm3 * fm4 * yi(1) / ((xi(1) - xi(2)) * (xi(1) - xi(3)) * (xi(1) - xi(4))) + &
          fm1 * fm3 * fm4 * yi(2) / ((xi(2) - xi(1)) * (xi(2) - xi(3)) * (xi(2) - xi(4))) + &
          fm1 * fm2 * fm4 * yi(3) / ((xi(3) - xi(1)) * (xi(3) - xi(2)) * (xi(3) - xi(4))) + &
          fm1 * fm2 * fm3 * yi(4) / ((xi(4) - xi(1)) * (xi(4) - xi(2)) * (xi(4) - xi(3))) 
    
  end subroutine lagrange

  subroutine laginterp(pgref, pint, co2i, co2v, fzeroi, fzerov)

    real(r8), intent(in ) :: pgref (      npref                       )
    real(r8), intent(out) :: pint  (      npint                       )
    real(r8), intent(out) :: co2i  (ntref,npint,nrefh2o,nspecti,ngauss)
    real(r8), intent(out) :: co2v  (ntref,npint,nrefh2o,nspectv,ngauss)
    real(r8), intent(out) :: fzeroi(                    nspecti       )
    real(r8), intent(out) :: fzerov(                    nspectv       )
    
    real(r8) co2i8(ntref,npref,nrefh2o,nspecti,ngauss)
    real(r8) co2v8(ntref,npref,nrefh2o,nspectv,ngauss)
    real(r8) x, xi(4), yi(4), ans
    real(r8) pref(npref), p
    integer n, nt, np, nh, ng, nw, m, i

    pint = [                                  &
      -6.0d0, -5.8d0, -5.6d0, -5.4d0, -5.2d0, &
      -5.0d0, -4.8d0, -4.6d0, -4.4d0, -4.2d0, &
      -4.0d0, -3.8d0, -3.6d0, -3.4d0, -3.2d0, &
      -3.0d0, -2.8d0, -2.6d0, -2.4d0, -2.2d0, &
      -2.0d0, -1.8d0, -1.6d0, -1.4d0, -1.2d0, &
      -1.0d0, -0.8d0, -0.6d0, -0.4d0, -0.2d0, &
       0.0d0,  0.2d0,  0.4d0,  0.6d0,  0.8d0, &
       1.0d0,  1.2d0,  1.4d0,  1.6d0,  1.8d0, &
       2.0d0,  2.2d0,  2.4d0,  2.6d0,  2.8d0, &
       3.0d0,  3.2d0,  3.4d0,  3.6d0,  3.8d0, &
       4.0d0                                  &
    ]

    ! Take log of the reference pressures.
    do n = 1, npref
      pref(n) = log10(pgref(n))
    end do

    ! Get CO2 k coefficients.
    open(20, file=trim(data_root)//'/old/CO2H2O_V_12_95_INTEL', form='unformatted', status='old')
    read(20) co2v8
    read(20) fzerov
    close(20)

    open(20, file=trim(data_root)//'/old/CO2H2O_IR_12_95_INTEL', form='unformatted', status='old')
    read(20) co2i8
    read(20) fzeroi
    close(20)

    ! Take Log10 of the values - we interpolate the log10 of the values,
    ! not the values themselves.   Smallest value is 1.0E-200.
    do nt = 1, ntref
      do np = 1, npref
        do nh = 1, nrefh2o
          do ng = 1, ngauss
            do nw = 1, nspectv
              if (co2v8(nt,np,nh,nw,ng) > 1.0d-200) then
                co2v8(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
              else
                co2v8(nt,np,nh,nw,ng) = -200.0
              end if
            end do
            do nw = 1, nspecti
              if (co2i8(nt,np,nh,nw,ng) > 1.0d-200) then
                co2i8(nt,np,nh,nw,ng) = log10(co2i8(nt,np,nh,nw,ng))
              else
                co2i8(nt,np,nh,nw,ng) = -200.0
              end if
            end do
          end do
        end do
      end do
    end do
    !  Interpolate the values:  first the IR
    do nt = 1, ntref
      do nh = 1, nrefh2o
        do nw = 1, nspecti
          do ng = 1, ngauss
            ! First, the initial interval (P=1e-6 to 1e-5)
            n = 1 
            do m = 1, 5
              x     = pint(m)
              xi(1) = pref(n)
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2i8(nt,n  ,nh,nw,ng)
              yi(2) = co2i8(nt,n+1,nh,nw,ng)
              yi(3) = co2i8(nt,n+2,nh,nw,ng)
              yi(4) = co2i8(nt,n+3,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2i(nt,m,nh,nw,ng) = 10.0**ans
            end do 
            do n = 2, npref - 2
              do m = 1, 5
                i     = (n - 1) * 5 + m
                x     = pint(i)
                xi(1) = pref(n-1)
                xi(2) = pref(n  )
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2i8(nt,n-1,nh,nw,ng)
                yi(2) = co2i8(nt,n  ,nh,nw,ng)
                yi(3) = co2i8(nt,n+1,nh,nw,ng)
                yi(4) = co2i8(nt,n+2,nh,nw,ng)
                call lagrange(x, xi, yi, ans)
                co2i(nt,i,nh,nw,ng) = 10.0**ans
              end do 
            end do
            !  Now, get the last interval (P=1e+3 to 1e+4)
            n = npref - 1
            do m = 1, 5
              i     = (n - 1) * 5 + m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n)
              xi(4) = pref(n+1)
              yi(1) = co2i8(nt,n-2,nh,nw,ng)
              yi(2) = co2i8(nt,n-1,nh,nw,ng)
              yi(3) = co2i8(nt,n,nh,nw,ng)
              yi(4) = co2i8(nt,n+1,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2i(nt,i,nh,nw,ng) = 10.0**ans
            end do  
            !  Fill the last pressure point
            co2i(nt,npint,nh,nw,ng) = 10.0**co2i8(nt,npref,nh,nw,ng)
          end do
        end do
      end do
    end do
    ! Interpolate the values:  now the visible
    do nt = 1, ntref
      do nh = 1, nrefh2o
        do nw = 1, nspectv
          do ng = 1, ngauss
            ! First, the initial interval (P=1e-6 to 1e-5)
            n = 1 
            do m = 1, 5
              x     = pint(m  )
              xi(1) = pref(n  )
              xi(2) = pref(n+1)
              xi(3) = pref(n+2)
              xi(4) = pref(n+3)
              yi(1) = co2v8(nt,n  ,nh,nw,ng)
              yi(2) = co2v8(nt,n+1,nh,nw,ng)
              yi(3) = co2v8(nt,n+2,nh,nw,ng)
              yi(4) = co2v8(nt,n+3,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2v(nt,m,nh,nw,ng) = 10.0**ans
            end do 
            do n = 2, npref - 2
              do m = 1, 5
                i     = (n - 1) * 5 + m
                x     = pint(i  )
                xi(1) = pref(n-1)
                xi(2) = pref(n  )
                xi(3) = pref(n+1)
                xi(4) = pref(n+2)
                yi(1) = co2v8(nt,n-1,nh,nw,ng)
                yi(2) = co2v8(nt,n  ,nh,nw,ng)
                yi(3) = co2v8(nt,n+1,nh,nw,ng)
                yi(4) = co2v8(nt,n+2,nh,nw,ng)
                call lagrange(x, xi, yi, ans)
                co2v(nt,i,nh,nw,ng) = 10.0**ans
              end do 
            end do
            !  Now, get the last interval (P=1e+3 to 1e+4)
            n = npref - 1
            do m = 1, 5
              i     = (n - 1) * 5 + m
              x     = pint(i)
              xi(1) = pref(n-2)
              xi(2) = pref(n-1)
              xi(3) = pref(n  )
              xi(4) = pref(n+1)
              yi(1) = co2v8(nt,n-2,nh,nw,ng)
              yi(2) = co2v8(nt,n-1,nh,nw,ng)
              yi(3) = co2v8(nt,n  ,nh,nw,ng)
              yi(4) = co2v8(nt,n+1,nh,nw,ng)
              call lagrange(x, xi, yi, ans)
              co2v(nt,i,nh,nw,ng) = 10.0**ans
            end do  
            !  Fill the last pressure point
            co2v(nt,npint,nh,nw,ng) = 10.0**co2v8(nt,npref,nh,nw,ng)
          end do
        end do
      end do
    end do

  end subroutine laginterp

end module gomars_v1_rad_mod
