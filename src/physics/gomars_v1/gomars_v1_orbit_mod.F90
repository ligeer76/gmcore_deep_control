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

module gomars_v1_orbit_mod

  use datetime
  use gomars_v1_const_mod
  use gomars_v1_types_mod

  implicit none

  private

  public gomars_v1_orbit_init
  public solar_dist
  public solar_decl_angle
  public update_solar_decl_angle
  public gomars_v1_orbit_cosz
  public gomars_v1_orbit_cosz_avg
  public decl_angle
  public cos_decl
  public sin_decl
  public time_of_day

  ! Convert Mkm to AU
  real(r8), parameter :: au        = 1.0_r8 / 149.597927_r8
  ! Aphelion Sun-Mars distance (Mkm)
  real(r8), parameter :: raphe     = mars_raphe
  ! Perihelion Sun-Mars distance (Mkm)
  real(r8), parameter :: rperi     = mars_rperi
  ! Obliquity (rad)
  real(r8), parameter :: obliq     = 25.2193_r8 * rad
  real(r8), parameter :: sin_obliq = sin(obliq)
  ! Eccentricity (rad)
  real(r8), parameter :: eccen     = (raphe - rperi) / (raphe + rperi)
  ! real(r8), parameter :: eccen     = 0.093379
  ! Semimajor axis (AU)
  real(r8), parameter :: semia     = rperi / (1 - eccen) * au
  ! real(r8), parameter :: semia     = 1.52369
  ! Semi-laus rectum
  real(r8), parameter :: pelip     = 0.5_r8 * (raphe + rperi) * (1 - eccen**2) * au
  ! Sol day per Martian year
  real(r8), parameter :: year_sol  = mars_sol_per_mars_year
  ! Sol day at perihelion (in datetime library)
  ! real(r8), parameter :: peri_sol  = 485
  ! Difference of solar longiutde between Ls~0 and perihelion (rad)
  real(r8) peri_dls

  ! Solar declination angle (rad), updated periodically
  real(r8) decl_angle
  real(r8) cos_decl
  real(r8) sin_decl
  real(r8) time_of_day

contains

  subroutine gomars_v1_orbit_init()

    type(datetime_type) time
    real(r8) peri_ls

    call time%init(my=0, sols=peri_sol, planet='mars')
    peri_dls = pi2 - time%solar_longitude()

  end subroutine gomars_v1_orbit_init

  ! Caluclate the distance between Sun and Mars in AU units.
  real(r8) function solar_dist(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    res = semia * (1 - eccen * cos(eccen_anomaly('mars', ls) + peri_dls))

  end function solar_dist

  pure real(r8) function solar_decl_angle(ls) result(res)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    res = asin(sin(ls) * sin_obliq)

  end function solar_decl_angle

  subroutine update_solar_decl_angle(ls)

    real(r8), intent(in) :: ls ! Solar longitude (rad)

    decl_angle = solar_decl_angle(ls)
    cos_decl = cos(decl_angle)
    sin_decl = sin(decl_angle)

  end subroutine update_solar_decl_angle

  subroutine gomars_v1_orbit_cosz(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i

    associate (mesh    => state%mesh        , &
               lon     => state%mesh%lon    , & ! in
               cos_lat => state%mesh%cos_lat, & ! in
               sin_lat => state%mesh%sin_lat, & ! in
               cosz    => state%cosz        )   ! out
    do i = 1, mesh%ncol
      cosz(i) = sin_lat(i) * sin_decl + cos_lat(i) * cos_decl * cos(pi2 * time_of_day + lon(i) - pi)
      if (cosz(i) < 1.0e-5_r8) cosz(i) = 0
    end do
    end associate

  end subroutine gomars_v1_orbit_cosz

  subroutine gomars_v1_orbit_cosz_avg(state)

    type(gomars_v1_state_type), intent(inout) :: state

    integer i
    real(r8) lon0, adt, c1, c2, c3, c4
    real(r8) t1, cosz1
    real(r8) t2, cosz2
    real(r8) t0, tr1, trr

    associate (mesh    => state%mesh        , &
               lon     => state%mesh%lon    , & ! in
               cos_lat => state%mesh%cos_lat, & ! in
               sin_lat => state%mesh%sin_lat, & ! in
               cosz    => state%cosz        )   ! out
    do i = 1, mesh%ncol
      lon0 = lon(i) - pi
      c1  = sin_decl * sin_lat(i)
      c2  = cos_decl * cos_lat(i)
      adt = pi2 * (dt / earth_day_seconds)
      t1  = pi2 * time_of_day
      t2  = t1 + adt
      c3  = t1 + lon0
      c4  = t2 + lon0
      cosz1 = c1 + c2 * cos(c3)
      cosz2 = c1 + c2 * cos(c4)

      if (cosz1 >= 0 .and. cosz2 >= 0) then
        ! Sun above the horizon for the entire time period.
        cosz(i) = c1 + c2 * (sin(c4) - sin(c3)) / adt
      else if (cosz1 <= 0 .and. cosz2 <= 0) then
        ! Sun below the horizon for the entire time period.
        cosz(i) = 0
      else
        ! Sun rises or sets during the time period.
        tr1 = acos(-c1 / c2)
        trr = tr1 - lon0
        if (trr > t2) then
          tr1 = pi2 - tr1
          trr = tr1 - lon0
          if (trr < t1) trr = trr + pi2
          if (trr > t2) trr = trr - pi2
        else if (trr < t1) then
          trr = trr + pi2
          if (.not. (trr >= t1 .and. trr <= t2)) then
            trr = trr + pi2
            if (.not. (trr >= t1 .and. trr <= t2)) then
              tr1 = pi2 - tr1
              trr = tr1 - lon0
              if (trr < t1) trr = trr + pi2
              if (trr > t2) trr = trr - pi2
            end if
          end if
        end if
        if (cosz1 < 0 .and. cosz2 > 0) then
          ! Sun rises after time t1.
          if (trr < t1 .or. trr > t2) trr = t2
          cosz(i) = (c1 * (t2 - trr) + c2 * (sin(t2 + lon0) - sin(trr + lon0))) / adt
        else
          ! Sun sets before time t2.
          if (trr < t1 .or. trr > t2) trr = t1
          cosz(i) = (c1 * (trr - t1) + c2 * (sin(trr + lon0) - sin(t1 + lon0))) / adt
        end if
      end if
    end do
    end associate

  end subroutine gomars_v1_orbit_cosz_avg

end module gomars_v1_orbit_mod
