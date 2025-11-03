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

module gomars_v1_lsm_mod

  use flogger
  use formula_mod
  use math_mod
  use gomars_v1_const_mod
  use gomars_v1_objects_mod
  use gomars_v1_types_mod
  use gomars_v1_namelist_mod
  use gomars_v1_tracers_mod

  implicit none

  private

  public gomars_v1_lsm_init
  public gomars_v1_lsm_final
  public gomars_v1_lsm_run
  public sthick, sthick_lev, sdepth, sdepth_lev

  real(r8), public, parameter :: gidn  = 0.0545_r8
  real(r8), public, parameter :: gids  = 0.0805_r8
  real(r8), public, parameter :: factl = 0.25_r8
  real(r8), public, parameter :: factm = 1.2_r8
  real(r8), public, parameter :: skind = 0.06_r8
  real(r8), parameter :: t1 = 50
  real(r8), parameter :: t2 = 350
  integer , parameter :: max_iter = 30

  real(r8) sdepth    (nsoil  )
  real(r8) sdepth_lev(nsoil+1)
  real(r8) sthick    (nsoil  )
  real(r8) sthick_lev(nsoil+1)

contains

  subroutine gomars_v1_lsm_init()

    integer i, k, iblk

    call gomars_v1_lsm_final()

    ! Thickness between full levels
    sthick(1) = factl * skind
    do k = 2, nsoil
      sthick(k) = sthick(k-1) * factm
    end do

    ! Depth of half levels
    sdepth_lev(1) = 0
    do k = 2, nsoil + 1
      sdepth_lev(k) = sdepth_lev(k-1) + sthick(k-1)
    end do

    ! Depth of full levels
    do k = 1, nsoil
      sdepth(k) = 0.5_r8 * (sdepth_lev(k) + sdepth_lev(k+1))
    end do

    ! Thickness between half levels
    sthick_lev(1) = 0.5_r8 * sthick(1)
    do k = 2, nsoil
      sthick_lev(k) = sdepth(k) - sdepth(k-1)
    end do
    sthick_lev(nsoil+1) = 0.5_r8 * sthick(nsoil)

    ! Initialize soil model.
    do iblk = 1, size(objects)
      associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state)
      do i = 1, mesh%ncol
        ! Default values
        do k = 1, nsoil
          state%rhosoil(i,k) = 1481.39_r8
          state%cpsoil (i,k) = 840.0_r8
        end do
        if (state%grss(i)) then
          ! Southern polar cap of ground ice
          do k = 1, nsoil
            if (sdepth(k) > gids) then
              state%zin    (i,k) = 2236.995_r8
              state%rhosoil(i,k) = 1781.99_r8
              state%cpsoil (i,k) = 1404.09_r8
            end if
          end do
        else if (state%grsn(i)) then
          ! Northern polar cap of ground ice
          do k = 1, nsoil
            if (sdepth(k) > gidn) then
              state%zin    (i,k) = 1100.00_r8
              state%rhosoil(i,k) = 1781.99_r8
              state%cpsoil (i,k) = 1404.09_r8
            end if
          end do
        end if
        ! Half levels
        state%scond(i,1) = state%zin(i,1)**2 / (state%rhosoil(i,1) * state%cpsoil(i,1))
        do k = 2, nsoil
          state%scond(i,k) = 0.5_r8 * (                                          &
            state%zin(i,k-1)**2 / (state%rhosoil(i,k-1) * state%cpsoil(i,k-1)) + &
            state%zin(i,k  )**2 / (state%rhosoil(i,k  ) * state%cpsoil(i,k  )))
        end do

        ! Soil temperature for cold run varying from 170K to 200K.
        do k = 1, nsoil
          state%stemp(i,k) = state%zavgtg(i)
        end do

        ! Initialize some parameters.
        state%h2osub_sfc(i) = 0
        state%h2oice_sfc(i) = state%qsfc(i,iMa_vap)
      end do
      end associate
    end do

  end subroutine gomars_v1_lsm_init

  subroutine gomars_v1_lsm_final()

  end subroutine gomars_v1_lsm_final

  subroutine gomars_v1_lsm_run(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k, iter
    logical done
    real(r8) dmgdt, tsat, emg15, emgout, downir, astar
    real(r8) tl, th, fl, fh, f, df, dx_old, dx, tmp
    real(r8) fcdn, tgp, tinp, wflux, qflx, qsat
    real(r8) flux(nsoil+1), a(nsoil), b(nsoil), c(nsoil), d(nsoil)

    associate (mesh         => state%mesh        , &
               lat          => state%mesh%lat    , & ! in
               alsp         => state%alsp        , & ! in
               als          => state%als         , & ! out
               qsfc         => state%qsfc        , & ! in
               npcflag      => state%npcflag     , & ! in
               rhosoil      => state%rhosoil     , & ! in
               cpsoil       => state%cpsoil      , & ! in
               scond        => state%scond       , & ! in
               stemp        => state%stemp       , & ! in
               zin          => state%zin         , & ! in
               rhouch       => state%rhouch      , & ! in
               ps           => state%ps          , & ! inout
               tg           => state%tg          , & ! out
               t            => state%t           , & ! in
               q            => state%q           , & ! in
               ht_sfc       => state%ht_sfc      , & ! in
               irflx_sfc_dn => state%irflx_sfc_dn, & ! in
               vsflx_sfc_dn => state%vsflx_sfc_dn, & ! in
               co2ice_sfc   => state%co2ice_sfc  , & ! in
               h2osub_sfc   => state%h2osub_sfc  , & ! in
               h2oice_sfc   => state%h2oice_sfc    & ! in
              )
    do i = 1, mesh%ncol
      dmgdt = 0
      ! Set surface albedo.
      als(i) = alsp(i)
      if (co2ice_sfc(i) > 0) then
        als(i) = merge(alices, alicen, lat(i) < 0)
      else if (albfeed .and. qsfc(i,iMa_vap) > icethresh_kgm2 .and. .not. npcflag(i)) then
        als(i) = icealb
      end if
      tsat = dewpoint_temperature_mars(ps(i))
      ! ------------------------------------------------------------------------
      ! No CO2 ice on the ground
      if (co2ice_sfc(i) <= 0) then
        co2ice_sfc(i) = 0
        ! Emissivities for bare ground
        emg15   = eg15gnd
        emgout  = egognd
        downir  = emg15 * irflx_sfc_dn(i)
        done = .false.
        astar = (1 - als(i)) * vsflx_sfc_dn(i)
        call funcd(astar, downir, rhouch(i), scond(i,1), stemp(i,1), ps(i), t1, &
                   t(i,nlev), q(i,nlev,iMa_vap), h2oice_sfc(i), npcflag(i), fl, df)
        call funcd(astar, downir, rhouch(i), scond(i,1), stemp(i,1), ps(i), t2, &
                   t(i,nlev), q(i,nlev,iMa_vap), h2oice_sfc(i), npcflag(i), fh, df)
        if (fl == 0) then
          tg = t1
          done = .true.
        else if (fh == 0) then
          tg = t2
          done = .true.
        else if (fl < 0) then
          tl = t1
          th = t2
        else
          tl = t2
          th = t1
        end if
        if (.not. done) then
          tg(i) = 0.5_r8 * (t1 + t2)
          dx_old = abs(t2 - t1)
          dx = dx_old
          call funcd(astar, downir, rhouch(i), scond(i,1), stemp(i,1), ps(i), tg(i), &
                     t(i,nlev), q(i,nlev,iMa_vap), h2oice_sfc(i), npcflag(i), f, df)
          do iter = 1, max_iter
            if (((tg(i) - th) * df - f) * ((tg(i) - tl) * df - f) >= 0 .or. abs(2 * f) > abs(dx_old * df)) then
              dx_old = dx
              dx = 0.5 * (th - tl)
              tg(i) = tl + dx
              if (tl == tg(i)) then
                done = .true.
                exit
              end if
            else
              dx_old = dx
              dx = f / df
              tmp = tg(i)
              tg(i) = tg(i) - dx
              if (tmp == tg(i)) then
                done = .true.
                exit
              end if
            end if
            if (abs(dx) < 1.0e-3) then
              done = .true.
              exit
            end if
            call funcd(astar, downir, rhouch(i), scond(i,1), stemp(i,1), ps(i), tg(i), &
                       t(i,nlev), q(i,nlev,iMa_vap), h2oice_sfc(i), npcflag(i), f, df)
            if (f < 0) then
              tl = tg(i)
            else
              th = tg(i)
            end if
          end do
        end if
        if (done) then
          qflx = 0
          if (h2oice_sfc(i) > 0 .and. .not. npcflag(i)) then
            qsat = water_vapor_saturation_mixing_ratio_mars(tg(i), ps(i))
            qflx = -rhouch(i) * (q(i,nlev,iMa_vap) - qsat) / cpd
            if (qflx * dt >= h2oice_sfc(i)) then
              qflx = h2oice_sfc(i) / dt
              h2oice_sfc(i) = 0
            else
              h2oice_sfc(i) = h2oice_sfc(i) - qflx * dt
            end if
          else if (npcflag(i)) then
            qsat = water_vapor_saturation_mixing_ratio_mars(tg(i), ps(i))
            qflx = -rhouch(i) * (q(i,nlev,iMa_vap) - qsat) / cpd
            h2oice_sfc(i) = h2oice_sfc(i) - qflx * dt
          end if
        else
          call log_error('Subroutine newtg did not converge!', __FILE__, __LINE__)
        end if
        ! NOTE: qflx > 0 is sublimation, qflx < 0 is condensation.
        h2osub_sfc(i) = h2osub_sfc(i) + qflx * dt

        if (tg(i) < tsat) then
          tg(i)  = tsat
          emg15  = eg15gnd
          emgout = egognd

          fcdn = -scond(i,1) * (stemp(i,1) - tsat) / sthick_lev(1)
          tgp = dt * ((1 - als(i)) * vsflx_sfc_dn(i) + ht_sfc(i) - emg15 * (stbo * tsat**4) - fcdn) / xlhtc

          ! Check if there is any CO2 ice accumulation.
          if (tgp < 0) then
            dmgdt = -tgp / dt
            co2ice_sfc(i) = -tgp
          else
            dmgdt = 0
            ! This term represents the last amounts of ice evaporating resulting in an
            ! increase in Tg. It still depends on TINP until I can figure out what to do with it.
            tinp  = sqrdy / zin(i,1)
            tg(i) = tsat + tgp * xlhtc * tinp
            co2ice_sfc(i) = 0
          end if
        end if
      ! ------------------------------------------------------------------------
      ! No CO2 ice on the ground
      else
        tg(i) = tsat

        ! Only modify if we have water condensation.
        qsat = water_vapor_saturation_mixing_ratio_mars(tg(i), ps(i))
        ! See Eq. (1) in Haberle et al. (2019), which used a bulk transfer approach,
        ! but note rhouch contains cpd, so here divides cpd.
        wflux = -rhouch(i) * (q(i,nlev,iMa_vap) - qsat) / cpd
        if (wflux < 0) then
          h2oice_sfc(i) = h2oice_sfc(i) - wflux * dt
          h2osub_sfc(i) = h2osub_sfc(i) + wflux * dt
        end if

        if (.not. npcflag(i) .and. h2osub_sfc(i) > qsfc(i,iMa_vap)) then
          h2osub_sfc(i) = qsfc(i,iMa_vap)
        end if

        if (lat(i) < 0) then
          emg15  = eg15co2s
          emgout = egoco2s
        else
          emg15  = eg15co2n
          emgout = egoco2n
        end if

        ! New soil scheme: surface boundary condition with ice on the ground.
        fcdn = -2 * scond(i,1) * (stemp(i,1) - tsat) / sthick(1)
        tgp  = -co2ice_sfc(i) + dt * ((1 - als(i)) * vsflx_sfc_dn(i) + ht_sfc(i) - emg15 * (stbo * tsat**4) - fcdn) / xlhtc

        ! Check if there is still CO2 ice left.
        if (tgp < 0) then
          dmgdt = -(co2ice_sfc(i) + tgp) / dt
          co2ice_sfc(i) = -tgp
        else
          dmgdt = -co2ice_sfc(i) / dt
          tinp  = sqrdy / zin(i,1)
          tg(i) = tsat + tgp * xlhtc * tinp
          co2ice_sfc(i) = 0
        end if
      end if
      ! Setup the tridiagonal matrix.
      a(1) = 0
      do k = 2, nsoil
        a(k) = -dt * scond(i,k  ) / (rhosoil(i,k) * cpsoil(i,k) * sthick(k) * sthick_lev(k  ))
      end do
      do k = 1, nsoil - 1
        c(k) = -dt * scond(i,k+1) / (rhosoil(i,k) * cpsoil(i,k) * sthick(k) * sthick_lev(k+1))
      end do
      c(nsoil) = 0
      do k = 1, nsoil
        b(k) = 1 - a(k) - c(k)
      end do
      d = stemp(i,:)
      ! Append the surface boundary condition.
      b(1) = b(1) + dt * scond(i,1) / (rhosoil(i,1) * cpsoil(i,1) * sthick(1) * sthick_lev(1))
      d(1) = d(1) + dt * scond(i,1) / (rhosoil(i,1) * cpsoil(i,1) * sthick(1) * sthick_lev(1)) * tg(i)
      ! Call tridiagonal solver to update soil temperature.
      call tridiag_thomas(a, b, c, d, stemp(i,:))
      ! Calculate the CO2 condensation temperature at the surface.
      tsat = dewpoint_temperature_mars(ps(i))
      if (tg(i) < tsat .or. co2ice_sfc(i) > 0) then
        tg(i) = tsat
      end if
    end do
    end associate

  end subroutine gomars_v1_lsm_run

  subroutine funcd(astar, downir, rhouch, scond, stemp, ps, tg, tbot, qbot, h2oice_sfc, npcflag, f, df)

    real(r8), intent(in ) :: astar
    real(r8), intent(in ) :: downir
    real(r8), intent(in ) :: rhouch
    real(r8), intent(in ) :: scond
    real(r8), intent(in ) :: stemp
    real(r8), intent(in ) :: ps
    real(r8), intent(in ) :: tg
    real(r8), intent(in ) :: tbot
    real(r8), intent(in ) :: qbot
    real(r8), intent(in ) :: h2oice_sfc
    logical , intent(in ) :: npcflag
    real(r8), intent(out) :: f
    real(r8), intent(out) :: df

    real(r8) qg

    f = astar + downir + rhouch * tbot - rhouch * tg + 2 * scond * (stemp - tg) / sthick(1) - stbo * tg**4
    df = -rhouch - 2 * scond / sthick(1) - 4 * stbo * tg**3

    if (latent_heat) then
      if (h2oice_sfc > 0 .or. npcflag) then
        qg = water_vapor_saturation_mixing_ratio_mars(tg, ps)
        f  = f + rhouch * (2.8e6_r8 / cpd) * (qbot - qg)
        df = df - 6146.1_r8 * rhouch * (2.8e6_r8 / cpd) * qg / tg**2
      end if
    end if

  end subroutine funcd

end module gomars_v1_lsm_mod