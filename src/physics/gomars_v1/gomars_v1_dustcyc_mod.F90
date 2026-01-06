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

module gomars_v1_dustcyc_mod

  use formula_mod
  use gomars_v1_const_mod
  use gomars_v1_namelist_mod
  use gomars_v1_types_mod
  use gomars_v1_tracers_mod

  implicit none

  private

  public gomars_v1_dustcyc_run

  real(r8), parameter :: cpp       = 744.5_r8
  real(r8), parameter :: reff_lift = 2.0e-6_r8

contains

  subroutine gomars_v1_dustcyc_run(state, dt)

    type(gomars_v1_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    if (use_ddl) call ddl_lift(state, dt)
    if (use_wsl) call wsl_lift(state, dt)

  end subroutine gomars_v1_dustcyc_run

  subroutine ddl_lift(state, dt)

    type(gomars_v1_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    integer i, k
    real(r8) x, b, psx, pconx, rm, mp, madd, nadd

    x = rd / cpp + 1

    associate (mesh   => state%mesh      , &
               ps     => state%ps        , & ! in
               pcon   => state%pcon      , & ! in
               p      => state%p         , & ! in
               dp     => state%dp_dry    , & ! in
               ht_pbl => state%ht_pbl    , & ! in
               dstflx => state%dstflx_ddl, & ! out
               qsfc   => state%qsfc      , & ! inout
               q      => state%q         )   ! inout
    do i = 1, mesh%ncol
      dstflx(i) = 0
      if (ht_pbl(i) > 0) then
        psx   = ps  (i)**x
        pconx = pcon(i)**x
        b = (psx - pconx) / ((ps(i) - pcon(i)) * x * psx / ps(i))
        dstflx(i) = max(alpha_d * (1 - b) * ht_pbl(i) * dt, 0.0_r8)
        if (dstflx(i) > 0) then
          rm = reff_lift * exp(-0.5_r8 * dev_dst**2)
          mp = 4.0_r8 / 3.0_r8 * pi * rm**3 * rho_dst
          madd = dstflx(i) / dp(i,nlev) * g
          nadd = madd / mp
          ! Find index of PBL top.
          do k = nlev, 1, -1
            if (p(i,k) >= pcon(i)) then
              q(i,k,iMa_dst) = q(i,k,iMa_dst) + madd
              q(i,k,iNb_dst) = q(i,k,iNb_dst) + nadd
            else
              exit
            end if
          end do
          qsfc(i,iMa_dst) = qsfc(i,iMa_dst) - dstflx(i)
        end if
      end if
    end do
    end associate

  end subroutine ddl_lift

  subroutine wsl_lift(state, dt)

    type(gomars_v1_state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    integer i
    real(r8) rhos, tau, hflx, rm, mp, madd, nadd

    associate (mesh       => state%mesh      , &
               ps         => state%ps        , & ! in
               tg         => state%tg        , & ! in
               dp         => state%dp_dry    , & ! in
               taux       => state%taux      , & ! in
               tauy       => state%tauy      , & ! in
               co2ice_sfc => state%co2ice_sfc, & ! in
               dstflx     => state%dstflx_wsl, & ! out
               qsfc       => state%qsfc      , & ! inout
               q          => state%q         )   ! inout
    do i = 1, mesh%ncol
      rhos = dry_air_density(tg(i), ps(i))
      dstflx(i) = 0
      tau = sqrt(taux(i)**2 + tauy(i)**2)
      if (tau <= tau_thresh .or. co2ice_sfc(i) > 0) then
        dstflx(i) = 0
      else
        hflx = 2.61_r8 * (tau**1.5_r8 / g / sqrt(rhos)) * &
          (1 - sqrt(tau_thresh / tau)) * (1 + sqrt(tau_thresh / tau))**2
        dstflx(i) = max(alpha_n * hflx * dt, 0.0_r8)
        if (dstflx(i) > 0) then
          rm = reff_lift * exp(-0.5_r8 * dev_dst**2)
          mp = 4.0_r8 / 3.0_r8 * pi * rm**3 * rho_dst
          madd = dstflx(i) / dp(i,nlev) * g
          nadd = madd / mp
          q(i,nlev,iMa_dst) = q(i,nlev,iMa_dst) + madd
          q(i,nlev,iNb_dst) = q(i,nlev,iNb_dst) + nadd
          qsfc(i,iMa_dst) = qsfc(i,iMa_dst) - dstflx(i)
        end if
      end if
    end do
    end associate

  end subroutine wsl_lift

end module gomars_v1_dustcyc_mod
