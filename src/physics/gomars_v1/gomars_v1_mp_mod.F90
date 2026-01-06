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

module gomars_v1_mp_mod

  use flogger
  use math_mod
  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_tracers_mod
  use gomars_v1_types_mod

  implicit none

  private

  public gomars_v1_mp_init
  public gomars_v1_mp_final
  public gomars_v1_mp_run
  public rho_aer, std_aer
  public rad_rt, radb_rt
  public qextv_dst, qscatv_dst, gv_dst
  public qexti_dst, qscati_dst, gi_dst
  public qextv_cld, qscatv_cld, gv_cld
  public qexti_cld, qscati_cld, gi_cld
  public cor_ratio

  real(r8), parameter :: rmin  = 0.1e-6_r8
  real(r8), parameter :: rmax  = 10.0e-6_r8
  real(r8), parameter :: rbmin = 0.0001e-6_r8
  real(r8), parameter :: rbmax = 1.0e-2_r8
  real(r8), parameter :: cor_ratio(nratio) = [   &
    0.10_r8, 0.20_r8, 0.25_r8, 0.30_r8, 0.35_r8, &
    0.40_r8, 0.45_r8, 0.50_r8, 0.55_r8, 0.60_r8, &
    0.65_r8, 0.70_r8, 0.80_r8, 0.90_r8, 0.99_r8  &
  ]

  real(r8), allocatable, dimension(:    ) :: rho_aer
  real(r8), allocatable, dimension(:    ) :: std_aer
  real(r8), allocatable, dimension(:    ) :: rad_aer
  real(r8), allocatable, dimension(:    ) :: radb_aer
  real(r8), allocatable, dimension(:    ) :: rad_rt       ! Water ice cloud particle size bins
  real(r8), allocatable, dimension(:    ) :: radb_rt      ! Water ice cloud particle size bin boundaries
  real(r8), allocatable, dimension(:,:,:) :: qextv_cld
  real(r8), allocatable, dimension(:,:,:) :: qscatv_cld
  real(r8), allocatable, dimension(:,:,:) :: gv_cld
  real(r8), allocatable, dimension(:,:,:) :: qexti_cld
  real(r8), allocatable, dimension(:,:,:) :: qscati_cld
  real(r8), allocatable, dimension(:,:,:) :: gi_cld
  real(r8), allocatable, dimension(:,:  ) :: qextv_dst
  real(r8), allocatable, dimension(:,:  ) :: qscatv_dst
  real(r8), allocatable, dimension(:,:  ) :: gv_dst
  real(r8), allocatable, dimension(:,:  ) :: qexti_dst
  real(r8), allocatable, dimension(:,:  ) :: qscati_dst
  real(r8), allocatable, dimension(:,:  ) :: gi_dst

contains

  subroutine gomars_v1_mp_init()

    integer i, j, l
    real(r8) vrat
    real(r8) factor

    call gomars_v1_mp_final()

    allocate(rho_aer(ntracers))
    allocate(std_aer(ntracers))

    rho_aer(iMa_dst) = rho_dst
    rho_aer(iNb_dst) = rho_dst
    rho_aer(iMa_cld) = rho_ice
    rho_aer(iNb_cld) = rho_ice
    rho_aer(iMa_cor) = rho_ice

    std_aer(iMa_dst) = dev_dst
    std_aer(iNb_dst) = dev_dst
    std_aer(iMa_cld) = dev_ice
    std_aer(iNb_cld) = dev_ice
    std_aer(iMa_cor) = dev_ice

    ! Initialize the aerosol size bins for microphysics.
    allocate(rad_aer (nbin_aer  ))
    allocate(radb_aer(nbin_aer+1))

    rad_aer (1) = rmin
    radb_aer(1) = rbmin

    vrat = log(rmax / rmin) / (nbin_aer - 1.0_r8) * 3
    vrat = exp(vrat)

    do i = 1, nbin_aer - 1
      rad_aer(i+1) = rad_aer(i) * vrat**athird
    end do

    do i = 1, nbin_aer
      radb_aer(i+1) = ((2 * vrat) / (vrat + 1))**athird * rad_aer(i)
    end do

    ! Initialize the water ice cloud radiation properties.
    allocate(rad_rt  (nbin_rt   ))
    allocate(radb_rt (nbin_rt+1 ))

    rad_rt (1        ) = 1.0e-7_r8
    rad_rt (nbin_rt  ) = 50.0e-6_r8
    radb_rt(1        ) = rbmin
    radb_rt(nbin_rt+1) = rbmax

    vrat = log(rad_rt(nbin_rt) / rad_rt(1)) / (nbin_rt - 1.0_r8) * 3
    vrat = exp(vrat)

    do i = 1, nbin_rt - 1
      rad_rt (i+1) = rad_rt(i) * vrat**athird
      radb_rt(i+1) = ((2 * vrat) / (vrat + 1))**athird * rad_rt(i)
    end do

    allocate(qextv_cld (nratio,nbin_rt,nspectv))
    allocate(qscatv_cld(nratio,nbin_rt,nspectv))
    allocate(gv_cld    (nratio,nbin_rt,nspectv))
    allocate(qexti_cld (nratio,nbin_rt,nspecti))
    allocate(qscati_cld(nratio,nbin_rt,nspecti))
    allocate(gi_cld    (nratio,nbin_rt,nspecti))
    allocate(qextv_dst (       nbin_rt,nspectv))
    allocate(qscatv_dst(       nbin_rt,nspectv))
    allocate(gv_dst    (       nbin_rt,nspectv))
    allocate(qexti_dst (       nbin_rt,nspecti))
    allocate(qscati_dst(       nbin_rt,nspecti))
    allocate(gi_dst    (       nbin_rt,nspecti))

    open(60, file=trim(data_root)//'/old/waterCoated_vis_JD_12bands.dat')
    open(61, file=trim(data_root)//'/old/waterCoated_ir_JD_12bands.dat')
    open(62, file=trim(data_root)//'/old/Dust_vis_wolff2010_JD_12bands.dat')
    open(63, file=trim(data_root)//'/old/Dust_ir_wolff2010_JD_12bands.dat')

    do j = 1, nratio
      do i = 1, nbin_rt
        read(60, '(7(e12.7,x))') (qextv_cld (j,i,l), l=1,nspectv)
        read(60, '(7(e12.7,x))') (qscatv_cld(j,i,l), l=1,nspectv)
        read(60, '(7(e12.7,x))') (gv_cld    (j,i,l), l=1,nspectv)
        read(61, '(5(e12.7,x))') (qexti_cld (j,i,l), l=1,nspecti)
        read(61, '(5(e12.7,x))') (qscati_cld(j,i,l), l=1,nspecti)
        read(61, '(5(e12.7,x))') (gi_cld    (j,i,l), l=1,nspecti)
      end do
    end do

    do i = 1, nbin_rt
      read(62, '(7(e11.5,x))') (qextv_dst (i,l), l=1,nspectv)
      read(62, '(7(e11.5,x))') (qscatv_dst(i,l), l=1,nspectv)
      read(62, '(7(e11.5,x))') (gv_dst    (i,l), l=1,nspectv)
      read(63, '(5(e11.5,x))') (qexti_dst (i,l), l=1,nspecti)
      read(63, '(5(e11.5,x))') (qscati_dst(i,l), l=1,nspecti)
      read(63, '(5(e11.5,x))') (gi_dst    (i,l), l=1,nspecti)
      ! factor = qextv_dst(i,6) / (vistoir * qexti_dst(i,4))
      factor = 1
      do l = 1, nspecti
        qexti_dst (i,l) = qexti_dst (i,l) * factor 
        qscati_dst(i,l) = qscati_dst(i,l) * factor 
      end do
    end do

    close(60)
    close(61)
    close(62)
    close(63)

  end subroutine gomars_v1_mp_init

  subroutine gomars_v1_mp_final()

    if (allocated(rho_aer   )) deallocate(rho_aer   )
    if (allocated(std_aer   )) deallocate(std_aer   )
    if (allocated(rad_aer   )) deallocate(rad_aer   )
    if (allocated(radb_aer  )) deallocate(radb_aer  )
    if (allocated(rad_rt    )) deallocate(rad_rt    )
    if (allocated(radb_rt   )) deallocate(radb_rt   )
    if (allocated(qextv_cld )) deallocate(qextv_cld )
    if (allocated(qscatv_cld)) deallocate(qscatv_cld)
    if (allocated(gv_cld    )) deallocate(gv_cld    )
    if (allocated(qexti_cld )) deallocate(qexti_cld )
    if (allocated(qscati_cld)) deallocate(qscati_cld)
    if (allocated(gi_cld    )) deallocate(gi_cld    )
    if (allocated(qextv_dst )) deallocate(qextv_dst )
    if (allocated(qscatv_dst)) deallocate(qscatv_dst)
    if (allocated(gv_dst    )) deallocate(gv_dst    )
    if (allocated(qexti_dst )) deallocate(qexti_dst )
    if (allocated(qscati_dst)) deallocate(qscati_dst)
    if (allocated(gi_dst    )) deallocate(gi_dst    )

  end subroutine gomars_v1_mp_final

  subroutine gomars_v1_mp_run(state, dt)

    type(gomars_v1_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    integer i, k, m
    real(r8) mo, no
    real(r8) dens(nlev,naer)
    real(r8) rcor(nlev,naer)

    associate (mesh        => state%mesh       , &
               p           => state%p          , & ! in
               p_lev       => state%p_lev      , & ! in
               dz          => state%dz         , & ! in
               dz_lev      => state%dz_lev     , & ! in
               t           => state%t          , & ! in
               t_lev       => state%t_lev      , & ! in
               kh          => state%kh         , & ! in
               q           => state%q          , & ! inout
               rho         => state%rhod       , & ! out
               rho_lev     => state%rhod_lev   , & ! out
               deposit     => state%deposit    , & ! out
               qsfc        => state%qsfc       , & ! inout
               qflx_sfc_dn => state%qflx_sfc_dn)   ! inout
    do i = 1, mesh%ncol
      ! Update air density since temperature may be changed.
      do k = 1, mesh%nlev
        rho    (i,k) = dry_air_density(t    (i,k), p    (i,k))
        rho_lev(i,k) = dry_air_density(t_lev(i,k), p_lev(i,k))
      end do
      rho_lev(i,nlev+1) = dry_air_density(t(i,nlev), p(i,nlev))
      do k = 1, mesh%nlev
        do m = 1, naer
          q(i,k,m) = max(q(i,k,m), 0.0_r8)
          select case (m)
          case (1)
            mo = q(i,k,iMa_dst)
            no = q(i,k,iNb_dst) + 1.0d-50
            dens(k,m) = rho_aer(m)
            rcor(k,m) = (mo / no * 0.75_r8 / pi / dens(k,m))**athird * exp(3 * std_aer(m)**2)
            if (mo < 1.0e-20_r8) rcor(k,m) = 1.0e-8_r8
          case (2)
            dens(k,m) = rho_aer(iMa_dst)
            rcor(k,m) = rcor(k,iMa_dst) * exp(-3 * std_aer(m)**2)
          case (3)
            mo = q(i,k,iMa_cld) + q(i,k,iMa_cor) + 1.0d-50
            no = q(i,k,iNb_cld) + 1.0d-50
            dens(k,m) = q(i,k,iMa_cld) / mo * rho_aer(iMa_cld) + &
                        q(i,k,iMa_cor) / mo * rho_aer(iMa_dst)
            dens(k,m) = min(max(dens(k,m), rho_ice), rho_dst)
            rcor(k,m) = (mo / no * 0.75_r8 / pi / dens(k,m))**athird * exp(3 * std_aer(m)**2)
            if (mo < 1.0e-20_r8) rcor(k,m) = 1.0e-8_r8
          case (4)
            dens(k,m) = dens(k,iMa_cld)
            rcor(k,m) = rcor(k,iMa_cld) * exp(-3 * std_aer(m)**2)
          case (5)
            dens(k,m) = dens(k,iMa_cld)
            rcor(k,m) = rcor(k,iMa_cld)
          end select
        end do
      end do
      ! Now compute the microphysical processes.
      ! Always check the order of sedimentation and nucleation condensation.
      ! PBL is called before microphysics, so do sedimentation first.
      call sedim(dt, t_lev(i,:), rho(i,:), rho_lev(i,:), dz(i,:), dz_lev(i,:), &
                 kh(i,:), q(i,:,:), rcor(:,:), dens(:,:), deposit(i,:), qflx_sfc_dn(i,:))

      ! Update tracers on the surface.
      qsfc(i,:) = qsfc(i,:) + deposit(i,:)

      if (qsfc(i,iMa_dst) > 1) then
        do k = 1, nlev
          print *, k, rcor(k,iMa_dst), q(i,k,iMa_dst)
        end do
        print *, mesh%gid(i), qsfc(i,iMa_dst)
        stop 999
      end if
    end do
    end associate

  end subroutine gomars_v1_mp_run

  subroutine sedim(dt, t_lev, rho, rho_lev, dz, dz_lev, kh, q, rcor, dens, deposit, qflx_sfc_dn)

    ! Sedimentation is based on the standard Stokes-Cunningham relationships for
    ! particle fall velocity with the slip correction for the thin Martian
    ! atmosphere.
    !
    !   vf = 2 * g * r^2 * ρ_p / (9 * μ) * (1 + cc * K_n)
    ! 
    ! where K_n is Knudsen number
    !
    !   K_n = λ / r
    !
    ! and λ is the mean free path of gas molecules.
    !
    !   cc = 1.246 + 0.42 * exp(-0.87 / K_n)

    real(r8), intent(in   ) :: dt
    real(r8), intent(in   ) :: t_lev       (nlev+1)
    real(r8), intent(in   ) :: rho         (nlev)
    real(r8), intent(in   ) :: rho_lev     (nlev+1)
    real(r8), intent(in   ) :: dz          (nlev)
    real(r8), intent(in   ) :: dz_lev      (nlev+1)
    real(r8), intent(in   ) :: kh          (nlev+1)        ! Vertical eddy coefficient from PBL (m2 s-1)
    real(r8), intent(inout) :: q           (nlev,ntracers)
    real(r8), intent(inout) :: rcor        (nlev,naer    ) ! Particle radius (m)
    real(r8), intent(inout) :: dens        (nlev,naer    ) ! Particle density (kg m-3)
    real(r8), intent(inout) :: deposit     (     ntracers)
    real(r8), intent(inout) :: qflx_sfc_dn (     ntracers) ! Tracer mass downward flux at the surface (kg m-2 s-1)

    integer k, m
    real(r8) q_old(nlev)
    real(r8) wt       ! Terminal velocity (m s-1)
    real(r8) dv       ! Gas molecule mean diffusion velocity (m s-1)
    real(r8) mfp      ! Gas molecule mean free path (m)
    real(r8) kn       ! Knudsen number (<<1: continuum regime, ~1: transition regime >>1: free molecular regime)
    real(r8) cc       ! Cunningham slip correction factor
    real(r8) w
    real(r8) lnrho    ! log(rho(k-1) / rho(k))
    real(r8) theta
    real(r8) sigma
    real(r8) rap
    real(r8) cfl
    real(r8) exp2t
    real(r8) ft(nlev+1)
    real(r8) fs(nlev+1)
    real(r8) as(nlev)
    real(r8) bs(nlev)
    real(r8) cs(nlev)
    real(r8) ds(nlev)
    logical solve_implicit

    real(r8) m1, m2

    deposit = 0
    theta = 0

    do m = 1, naer
      q_old = rho * q(:,m)
      q(:,m) = q_old
      do k = 2, nlev
        ! Calculate fall velocity on half levels.
        ! - Thermal velocity of CO2 molecules
        wt = sqrt(scale_co2 * t_lev(k))
        ! - CO2 gas diffusion velocity
        dv = 1.59e-6_r8 * t_lev(k)**1.5_r8 / (t_lev(k) + 244.4_r8)
        ! - CO2 gas mean free path
        mfp = 2 * dv / (rho_lev(k) * wt)
        ! - Knudsen number
        kn = mfp / rcor(k,m)
        ! - Cunningham slip correction factor (Alvarez et al., 2024 gave 1.168 + 0.552 * exp(-0.990 / kn)
        cc = 1.246_r8 + 0.42_r8 * exp(-0.87_r8 / kn)
        ! - Fall velocity
        w  = 2 * g * rcor(k,m)**2 * dens(k,m) / (9 * dv) * (1 + cc * kn)
        w  = -w * exp(-std_aer(m)**2)
        ! - Courant number
        cfl = w * dt / dz_lev(k)
        ! Correct the fall velocity accounting for mixing.
        lnrho = log(rho(k-1) / rho(k))
        if (kh(k) /= 0) then
          theta = 0.5_r8 * (w * dz_lev(k) / kh(k) + lnrho)
          if (theta /= 0) then
            sigma = 1.0_r8 / tanh(theta) - 1.0_r8 / theta
          else
            sigma = 1
          end if
        else
          sigma = 1
        end if
        if (q_old(k) == 0) then
          rap = 10
          if (q_old(k-1) == 0) then
            rap = 1
          end if
        else
          rap = min(max(q_old(k-1) / q_old(k), 0.1_r8), 10.0_r8)
        end if
        if (.not. (rap > 0.9 .and. rap < 1.1 .or. abs(w * dt) > dz(k))) then
          if (w < 0) then
            w = w * exp(sigma * log((exp(-cfl * log(rap)) - 1) / (cfl * (1 - rap))))
          else
            w = 0
          end if
        end if
        ! Calculate fluxes.
        if (kh(k) /= 0) then
          if (theta /= 0) then ! Has turbulent mixing and sedimentation or density gradient.
            exp2t = exp(2 * theta)
            ft(k) = (w + kh(k) / dz_lev(k) * lnrho) / (exp2t - 1)
            fs(k) = ft(k) * exp2t
          else ! Only has turbulent mixing and no density gradient.
            ft(k) = kh(k) / dz_lev(k)
            fs(k) = ft(k)
          end if
        else
          if (w < 0) then ! Only has sedimentation.
            ft(k) = -w
            fs(k) = 0
          else ! Neither has turbulent mixing nor sedimentation.
            ! FIXME: w should never be positive.
            ft(k) = 0
            fs(k) = w
          end if
        end if
      end do
      ! Boundary conditions for the fluxes.
      ft(1) = 0
      fs(1) = 0
      ft(nlev+1) = -w
      fs(nlev+1) = 0
      ! Calculate the coefficients of the continuity equation.
      solve_implicit = .false.
      do k = 1, nlev
        cs(k) = ft(k+1) + fs(k) - dz(k) / dt
        if (cs(k) > 0) then
          solve_implicit = .true.
          exit
        end if
        as(k) = -dz(k) / dt
        bs(k) = -ft(k)
        ds(k) = -fs(k+1)
      end do
      if (solve_implicit) then
        do k = 1, nlev
          as(k) = ft(k)
          bs(k) = -(ft(k+1) + fs(k) + dz(k) / dt)
          cs(k) = fs(k+1)
          ds(k) = -dz(k) / dt * q_old(k)
        end do
        call tridiag_thomas(as, bs, cs, ds, q(:,m))
      else
        q(1,m) = (cs(1) * q_old(1) + ds(1) * q_old(2)) / as(1)
        do k = 2, nlev - 1
          q(k,m) = (bs(k) * q_old(k-1) + cs(k) * q_old(k) + ds(k) * q_old(k+1)) / as(k)
        end do
      end if
      ! Calculate the tracer mass falling on the ground.
      if (m == iMa_cld) then
        deposit    (iMa_vap) = deposit    (iMa_vap) + q(nlev,m) * ft(nlev+1) * dt
        qflx_sfc_dn(iMa_cld) = qflx_sfc_dn(iMa_cld) + q(nlev,m) * ft(nlev+1)
      else if (m == iMa_dst) then
        deposit    (iMa_dst) = deposit    (iMa_dst) + q(nlev,m) * ft(nlev+1) * dt
        qflx_sfc_dn(iMa_dst) = qflx_sfc_dn(iMa_dst) + q(nlev,m) * ft(nlev+1)
      else if (m == iMa_cor) then
        deposit    (iMa_cor) = deposit    (iMa_cor) + q(nlev,m) * ft(nlev+1) * dt
        qflx_sfc_dn(iMa_cor) = qflx_sfc_dn(iMa_cor) + q(nlev,m) * ft(nlev+1)
      end if
      if (.not. solve_implicit) then
        q(nlev,m) = (bs(nlev) * q_old(nlev-1) + cs(nlev) * q_old(nlev)) / as(nlev)
      end if
      q(:,m) = q(:,m) / rho
    end do

  end subroutine sedim

  subroutine nucleacond(dt, p, t, rho, q)

    real(r8), intent(in   ) :: dt
    real(r8), intent(in   ) :: p           (nlev)
    real(r8), intent(in   ) :: t           (nlev)
    real(r8), intent(in   ) :: rho         (nlev)
    real(r8), intent(inout) :: q           (nlev,ntracers)

    integer i, k
    real(r8) qsat
    real(r8) ph2o
    real(r8) sat_ratio
    real(r8) Mo
    real(r8) No
    real(r8) Rn
    real(r8) Rm
    real(r8) dev2
    real(r8) n_aer(nbin_aer)
    real(r8) m_aer(nbin_aer)
    real(r8) rate (nbin_aer)

    do k = 1, nlev
      qsat = water_vapor_saturation_mixing_ratio_mars(t(k), p(k))
      ph2o = q(k,iMa_vap) * (m_co2 / m_h2o) * p(k)
      sat_ratio = q(k,iMa_vap) / qsat

      Mo = q(k,iMa_dst)
      No = q(k,iNb_dst) + 1.0d-50
      Rn = (Mo / No * 0.75_r8 / pi / rho_aer(iMa_dst))**athird * exp(-1.5_r8 * std_aer(iMa_dst)**2)
      Rn = min(max(Rn, rad_aer(1)), rad_aer(nbin_aer))
      Rm = Rn * exp(3 * std_aer(iMa_dst)**2)
      Rn = 1.0_r8 / Rn
      Rm = 1.0_r8 / Rm
      dev2 = 1.0_r8 / (sqrt2 * std_aer(iMa_dst))
      do i = 1, nbin_aer
        m_aer(i) = 0.5_r8 * Mo * (erf(log(radb_aer(i+1) * Rm) * dev2) - erf(log(radb_aer(i) * Rm) * dev2))
        n_aer(i) = 0.5_r8 * No * (erf(log(radb_aer(i+1) * Rn) * dev2) - erf(log(radb_aer(i) * Rn) * dev2))
      end do

      ! Calculate the nucleation rates.
    end do

  end subroutine nucleacond

  subroutine nuclea(ph2o, t, sat_ratio, n_aer, rate)

    real(r8), intent(in ) :: ph2o
    real(r8), intent(in ) :: t
    real(r8), intent(in ) :: sat_ratio
    real(r8), intent(in ) :: n_aer(nbin_aer)
    real(r8), intent(out) :: rate (nbin_aer)

  end subroutine nuclea

  pure real(r8) function sig(t) result(res)

    real(r8), intent(in) :: t

    res = (141 - 0.15_r8 * t) * 1.0e-3_r8

  end function sig

end module gomars_v1_mp_mod
