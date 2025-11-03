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

module gomars_v1_cnvadj_mod

  use gomars_v1_const_mod
  use gomars_v1_types_mod

  implicit none

  private

  public gomars_v1_cnvadj_run

contains

  subroutine gomars_v1_cnvadj_run(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i, k, k2, k3, m
    real(r8) wgt(nlev), pcexp, wgtcon, dpcon
    real(r8) qcon(ntracers)

    associate (mesh   => state%mesh  , &
               ps     => state%ps    , & ! in
               p      => state%p     , & ! in
               p_lev  => state%p_lev , & ! in
               dp     => state%dp_dry, & ! in
               pk     => state%pk    , & ! in
               pt     => state%pt    , & ! inout
               pt_lev => state%pt_lev, & ! in
               t      => state%t     , & ! inout
               q      => state%q     , & ! inout
               pcon   => state%pcon  , & ! out
               ptcon  => state%ptcon )   ! out
    do i = 1, mesh%ncol
      ! Weighting factor for convection. wgt * pt is proportional to the head
      ! energy of the layer.
      do k = 1, mesh%nlev
        wgt(k) = dp(i,k) * pk(i,k)
      end do

      loop_k1: do k = 1, mesh%nlev - 1
        if (pt(i,k) > pt(i,k+1)) then
          ! --------------------------------------------------------------------
          ! No instability found. No adjustments are needed.
          if (k + 1 /= nlev) cycle loop_k1
          ! --------------------------------------------------------------------
          ! Assume a little convection near the ground.
          if (pt_lev(i,nlev) > pt(i,nlev)) then
            ptcon(i) = 0.5_r8 * (pt_lev(i,nlev+1) + pt(i,nlev))
            pcon (i) = p(i,nlev)
          ! --------------------------------------------------------------------
          ! Shallow conective layer presents.
          else
            ptcon(i) = pt(i,nlev)
            pcexp    = 2 * (ptcon(i) - pt_lev(i,nlev)) / (pt(i,nlev-1) - pt_lev(i,nlev))
            pcon (i) = p_lev(i,nlev) * (p(i,nlev-1) / p_lev(i,nlev))**pcexp
            pcon (i) = max(pcon(i), p_lev(i,nlev-1))
          end if
          cycle loop_k1
        else
          ! Instability found. Do adjustments required to make pt non-decreasing
          ! in altitude.
          wgtcon   = wgt(  k+1)
          dpcon    = dp (i,k+1)
          qcon      = q  (i,k+1,:)
          ptcon(i)  = pt (i,k+1)
          ! Search upwards.
          loop_k2: do k2 = k, 1, -1
            if (pt(i,k2) <= ptcon(i)) then
              ! ----------------------------------------------------------------
              ! Found new instability.
              if (pt(i,k2) /= ptcon(i)) then
                ptcon(i) = (pt(i,k2) * wgt(k2) + ptcon(i) * wgtcon) / (wgt(k2) + wgtcon)
              end if
              ! Mix tracers.
              do m = 1, ntracers
                if (q(i,k2,m) /= qcon(m)) then
                  qcon(m) = (q(i,k2,m) * dp(i,k2) + qcon(m) * dpcon) / (dp(i,k2) + dpcon)
                end if
              end do
              ! Adjusts the layers from k2 to k.
              loop_k3: do k3 = k2, k + 1
                pt    (i,k3  ) = ptcon(i)
                pt_lev(i,k3+1) = ptcon(i)
                q     (i,k3,:) = qcon
              end do loop_k3
              pcon(i) = p_lev(i,k2)
              wgtcon  = wgtcon + wgt(k2)
              dpcon   = dpcon + dp(i,k2)
            else
              ! ----------------------------------------------------------------
              ! Found no more instabilities. Only need to check for pcon adjustment.
              if (pt_lev(i,k2+1) >= ptcon(i)) cycle loop_k1
              if (k + 1 /= nlev) cycle loop_k1
              pcexp   = 2 * (ptcon(i) - pt_lev(i,k2+1)) / (pt(i,k2) - pt_lev(i,k2+1))
              pcon(i) = p_lev(i,k2+1) * (p(i,k2) / p_lev(i,k2+1))**pcexp
              pcon(i) = max(pcon(i), p_lev(i,k2))
            end if
          end do loop_k2
        end if
      end do loop_k1
      ! ------------------------------------------------------------------------
      ! Update temperature.
      do k = 1, mesh%nlev
        t(i,k) = pt(i,k) * (p(i,k) / ps(i))**rd_o_cpd
      end do
    end do
    end associate

  end subroutine gomars_v1_cnvadj_run

end module gomars_v1_cnvadj_mod