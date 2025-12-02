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

module gomars_v1_co2cyc_mod

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_types_mod
  use gomars_v1_namelist_mod, only: co2scav
  use gomars_v1_tracers_mod

  implicit none

  private

  public gomars_v1_co2cyc_run

contains

  subroutine gomars_v1_co2cyc_run(state, tend)

    type(gomars_v1_state_type), intent(inout) :: state
    type(gomars_v1_tend_type ), intent(inout) :: tend

    integer i, k, m
    integer k_scavup, k_scavdn
    real(r8) tsat, tgp, dm(0:nlev)

    associate (mesh        => state%mesh       , &
               zin         => state%zin        , & ! in
               tstrat      => state%tstrat     , & ! inout
               ps          => state%ps         , & ! in
               p           => state%p          , & ! in
               dp          => state%dp_dry     , & ! in
               t           => state%t          , & ! inout
               tg          => state%tg         , & ! inout
               q           => state%q          , & ! inout
               qsfc        => state%qsfc       , & ! inout
               qflx_sfc_dn => state%qflx_sfc_dn, & ! inout
               co2ice_sfc  => state%co2ice_sfc , & ! in
               dpsdt       => tend %dpsdt      )   ! inout
    do i = 1, mesh%ncol
      dm = 0
      ! ------------------------------------------------------------------------
      ! Calculate stratospheric condensation.
      if (tstrat(i) < tsat_strat) then
        dm     (  0) = cpd * (tsat_strat - tstrat(i)) * ptrop / g / xlhtc ! [kg m-2]
        tstrat (i  ) = tsat_strat
        dpsdt  (i  ) = dpsdt(i) + dm(0) / dt * g
      end if
      ! ------------------------------------------------------------------------
      ! Calculate tropospheric condensation.
      do k = 1, mesh%nlev
        tsat = dewpoint_temperature_mars(p(i,k))
        if (t(i,k) < tsat) then
          dm   (  k) = cpd * (tsat - t(i,k)) * dp(i,k) / g / xlhtc
          t    (i,k) = tsat
          dpsdt(i  ) = dpsdt(i) + dm(k) / dt * g
        end if
      end do
      ! ------------------------------------------------------------------------
      if (dpsdt(i) > 0) then
        ! CO2 frost point at this surface pressure.
        tsat = dewpoint_temperature_mars(ps(i))
        if (co2ice_sfc(i) > 0) then
          ! Case 1: CO2 ice already on the ground.
          ! Add condensation to existing CO2 ice mass.
          co2ice_sfc(i) = co2ice_sfc(i) + dpsdt(i) * dt / g
          tg(i) = tsat
        else
          ! Case 2: No CO2 ice on the ground; Ground is warmer.
          ! Ground temperature drops when CO2 ice sublimes on warmer surface.
          tgp  = tg(i) - sqrdy / zin(i,1) * dpsdt(i) * xlhtc * dt / g
          ! Calculate how much CO2 ice sublimes on hitting the ground.
          if (tgp >= tsat) then
            ! Case 2A: All CO2 ice sublimes; No net condensation.
            dpsdt     (i) = 0
            tg        (i) = tgp
            co2ice_sfc(i) = 0
          else
            ! Case 2B: Ground cooled to CO2 frost point and some CO2 ice remains.
            ! Calculate how much CO2 ice remains.
            dpsdt     (i) = dpsdt(i) * (tsat - tgp) / (tg(i) - tgp)
            tg        (i) = tsat
            co2ice_sfc(i) = dpsdt(i) * dt / g
          end if
        end if
      end if
      ! ------------------------------------------------------------------------
      ! Aerosol scavenging by CO2 snow fall.
      if (co2scav) then
        k_scavup = 0
        k_scavdn = 0
        ! Determine the portion of atmosphere affected by CO2 condensation.
        do k = 1, mesh%nlev
          if (dm(k) > 0) then
            k_scavup = k
            exit
          end if
        end do
        do k = mesh%nlev, 1, -1
          if (dm(k) > 0) then
            k_scavdn = k
            exit
          end if
        end do
        if (k_scavup /= 0 .and. k_scavdn /= 0) then
          do k = k_scavup, k_scavdn
            if (k_scavdn == mesh%nlev) then
              ! If condensation occurs down to the surface, put all aerosols on the surface (in fact, only the cloud mass matters).
              qsfc       (i,iMa_vap) = qsfc       (i,iMa_vap) + scaveff * q(i,k,iMa_cld) * dp(i,k) / g
              qsfc       (i,iMa_dst) = qsfc       (i,iMa_dst) + scaveff * q(i,k,iMa_dst) * dp(i,k) / g
              qsfc       (i,iMa_cor) = qsfc       (i,iMa_cor) + scaveff * q(i,k,iMa_cor) * dp(i,k) / g
              qflx_sfc_dn(i,iMa_dst) = qflx_sfc_dn(i,iMa_dst) + scaveff * q(i,k,iMa_dst) * dp(i,k) / g / dt
              qflx_sfc_dn(i,iMa_cor) = qflx_sfc_dn(i,iMa_cor) + scaveff * q(i,k,iMa_cor) * dp(i,k) / g / dt
              qflx_sfc_dn(i,iMa_cld) = qflx_sfc_dn(i,iMa_cld) + scaveff * q(i,k,iMa_cld) * dp(i,k) / g / dt
            else
              ! If condensation occurs in a restricted portion, put aerosols in the highest layer unaffected by CO2 condensation.
              do m = 1, ntracers
                q(i,k_scavdn+1,m) = q(i,k_scavdn+1,m) + scaveff * q(i,k,m) * dp(i,k) / dp(i,k_scavdn+1)
              end do
            end if
            do m = 1, naer
              q(i,k,m) = q(i,k,m) * (1 - scaveff)
            end do
            dm(k) = 0
          end do
        end if
      end if
    end do
    end associate

  end subroutine gomars_v1_co2cyc_run

end module gomars_v1_co2cyc_mod