!=======================================================================
!  Tropical Kelvin / Equatorial Rossby test case IC
!  - Matches the "associate(...)" + work arrays + halo fill + A2C wind
!    conversion style of baroclinic_wave_test_set_ic
!
!  pertt:
!    0 -> Kelvin-like: symmetric, v=0, add (u',T') in-phase
!    1 -> Equatorial Rossby-like: antisymmetric streamfunction ψ -> (u',v')
!
!  Notes:
!  - Background is taken from evaluate_z_temperature(...) (same as your
!    baroclinic wave), but background winds are set to 0 to avoid jet.
!  - Vertical taper uses the same cubic polynomial taper you already use.
!  - Uses p = mg(i,j,k) from your hybrid coordinate state.
!=======================================================================

MODULE tropical_wave_test_mod

  IMPLICIT NONE

  private

  public tropical_wave_test_init
  public tropical_wave_test_set_ic

  REAL(8), PARAMETER ::               &
  a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
  Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
  g     = 9.80616d0,             & ! Gravity (m s^2)
  cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
  Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
  Rvap  = 461.5d0,               & ! Ideal gas constnat for water vapor
  Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
  pi    = 3.14159265358979d0,    & ! pi
  p0    = 100000.0d0,            & ! surface pressure (Pa)
  kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
  omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
  deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

  REAL(8), PARAMETER ::               &
  T0E        = 310.d0     ,      & ! temperature at equatorial surface (K)
  T0P        = 240.d0     ,      & ! temperature at polar surface (K)
  B          = 2.d0       ,      & ! jet half-width parameter
  K          = 3.d0       ,      & ! jet width parameter
  lapse      = 0.005d0             ! lapse rate parameter

  ! REAL(8), PARAMETER :: pi = 3.1415926535897932384626433832795d0
  ! REAL(8), PARAMETER :: p0       = 100000.0d0    ! Surface pressure (Pa)

  ! --- Tunable parameters (ideally read from namelist_mod) ---
  INTEGER :: pertt = 0                  ! 0 Kelvin, 1 Rossby-like
  REAL(8) :: pert_uamp   = 0.5d0        ! m/s   Kelvin u' amplitude
  REAL(8) :: pert_tamp   = 0.5d0        ! K     Kelvin T' amplitude
  REAL(8) :: pert_z      = 20000.d0     ! m     vertical taper top
  REAL(8) :: pert_lat0   = 10.d0*pi/180.d0  ! rad equatorial trapping width
  INTEGER :: pert_m      = 2            ! zonal wavenumber
  REAL(8) :: pert_lon0   = pi           ! rad phase center
  REAL(8) :: pert_psiamp = 2.0d7        ! m^2/s streamfunction amplitude scale (Rossby)

  REAL(8) :: dxepsilon = 1.0d-6         ! rad for numerical derivatives
  
  integer :: moist = 0
  integer :: deep  = 0
  ! integer :: pertt = 0

  namelist /tropical_wave_control/ moist
  namelist /tropical_wave_control/ deep
  namelist /tropical_wave_control/ pertt

CONTAINS

  subroutine tropical_wave_test_init(namelist_path)

    use namelist_mod
    use tracer_mod

    character(*), intent(in) :: namelist_path

    integer ignore

    open(11, file=namelist_path, status='old')
    read(11, nml=tropical_wave_control, iostat=ignore)
    close(11)

    if (moist == 1) then
      call tracer_add('moist', dt_adv, 'qv', 'Water vapor', 'kg kg-1')
    end if

    ptop = 226.0d0 ! Recommended pressure at the model top

  end subroutine tropical_wave_test_init
  !---------------------------------------------------------------------
  !  Public IC setter: matches your baroclinic_wave_test_set_ic format
  !---------------------------------------------------------------------
  SUBROUTINE tropical_wave_test_set_ic(block)

    USE namelist_mod
    USE block_mod
    USE tracer_mod
    USE operators_mod
    USE formula_mod
    USE latlon_parallel_mod
    USE latlon_operators_mod

    TYPE(block_type), INTENT(INOUT), TARGET :: block

    INTEGER :: i, j, k
    REAL(8) :: thetav, rho, z, qv

    ASSOCIATE (mesh   => block%mesh            , &
               gzs    => block%static%gzs      , & ! out
               mgs    => block%dstate(1)%mgs   , & ! out
               u      => block%aux%u           , & ! work array (A-grid)
               v      => block%aux%v           , & ! work array (A-grid)
               u_lon  => block%dstate(1)%u_lon , & ! out (C-grid)
               v_lat  => block%dstate(1)%v_lat , & ! out (C-grid)
               t      => block%aux%t           , & ! work array
               pt     => block%dstate(1)%pt    , & ! out
               mg     => block%dstate(1)%mg    , & ! work array (pressure at full levels)
               ph_lev => block%dstate(1)%ph_lev, & ! work array (geopotential at half/lev)
               q      => tracers(block%id)%q   )   ! out (tracers)

      !---------------------------------------------------------------
      ! Surface fields
      !---------------------------------------------------------------
      gzs%d = 0.d0
      mgs%d = p0

      CALL calc_mg(block, block%dstate(1))
      CALL calc_ph(block, block%dstate(1))

      z = 0.d0
      DO k = mesh%full_kds, mesh%full_kde
        DO j = mesh%full_jds, mesh%full_jde
          DO i = mesh%full_ids, mesh%full_ide

            CALL tropical_wave_test( &
              deep   = deep             , &
              moist  = 0                , & ! recommend dry
              pertt  = pertt            , &
              X      = scale_X          , &
              lon    = mesh%full_lon(i) , &
              lat    = mesh%full_lat(j) , &
              p      = mg%d(i,j,k)      , &
              z      = z                , &
              zcoords= 0                , &
              u      = u%d(i,j,k)       , &
              v      = v%d(i,j,k)       , &
              t      = t%d(i,j,k)       , &
              thetav = thetav           , &
              phis   = gzs%d(i,j)       , &
              ps     = mgs%d(i,j)       , &
              rho    = rho              , &
              q      = qv               )

            ! No moisture in this case, but keep interface consistent
            IF (moist == 1) THEN
              q%d(i,j,k,1) = qv
            END IF

          END DO
        END DO
      END DO

      !---------------------------------------------------------------
      ! Moist branch (kept for consistency; for this test, moist=0)
      !---------------------------------------------------------------
      IF (moist == 1) THEN
        DO k = mesh%full_kds, mesh%full_kde
          DO j = mesh%full_jds, mesh%full_jde
            DO i = mesh%full_ids, mesh%full_ide
              mgs%d(i,j) = mgs%d(i,j) - (ph_lev%d(i,j,k+1) - ph_lev%d(i,j,k)) * q%d(i,j,k,1)
            END DO
          END DO
        END DO
        CALL fill_halo(mgs)

        DO k = mesh%full_kds, mesh%full_kde
          DO j = mesh%full_jds, mesh%full_jde
            DO i = mesh%full_ids, mesh%full_ide
              q%d(i,j,k,1) = q%d(i,j,k,1) / (1.d0 - q%d(i,j,k,1))
            END DO
          END DO
        END DO
        CALL fill_halo(q, 1)
      END IF

      !---------------------------------------------------------------
      ! Potential temperature (modified, consistent with your code)
      !---------------------------------------------------------------
      CALL calc_mg(block, block%dstate(1))
      DO k = mesh%full_kds, mesh%full_kde
        DO j = mesh%full_jds, mesh%full_jde
          DO i = mesh%full_ids, mesh%full_ide
            IF (moist == 1) THEN
              qv = q%d(i,j,k,1)
            ELSE
              qv = 0.d0
            END IF
            pt%d(i,j,k) = modified_potential_temperature(t%d(i,j,k), mg%d(i,j,k), qv)
          END DO
        END DO
      END DO
      CALL fill_halo(pt)

      !---------------------------------------------------------------
      ! Convert A-grid wind to C-grid wind
      !---------------------------------------------------------------
      CALL fill_halo(u)
      CALL fill_halo(v)
      CALL wind_a2c_operator(u, v, u_lon, v_lat)
      CALL fill_halo(u_lon)
      CALL fill_halo(v_lat)

      init_hydrostatic_gz = .true.

    END ASSOCIATE

  END SUBROUTINE tropical_wave_test_set_ic


  !---------------------------------------------------------------------
  ! Pointwise test case (BIND(C) just like your baroclinic_wave_test)
  !---------------------------------------------------------------------
  SUBROUTINE tropical_wave_test(deep,moist,pertt,X,lon,lat,p,z,zcoords, &
                                u,v,t,thetav,phis,ps,rho,q) &
    BIND(c, name = "tropical_wave_test")

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: deep, moist, pertt
    REAL(8), INTENT(IN) :: X, lon, lat
    REAL(8), INTENT(INOUT) :: p, z
    INTEGER, INTENT(IN) :: zcoords

    REAL(8), INTENT(OUT) :: u, v, t, thetav, phis, ps, rho, q

    REAL(8) :: t_bg, uu, vv, tt, qv

    !------------------------------------------------
    ! Background: reuse your balanced p-z-T relation
    !------------------------------------------------
    IF (zcoords .EQ. 1) THEN
      ! (Not used by current IC loop; keep consistent)
      CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t_bg)
    ELSE
      CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t_bg)
    END IF

    !------------------------------------------------
    ! Base fields: quiet background (no jet)
    !------------------------------------------------
    ps   = p0
    phis = 0.d0
    u    = 0.d0
    v    = 0.d0
    t    = t_bg

    ! Dry by default
    qv = 0.d0
    q  = qv

    !------------------------------------------------
    ! Add perturbations
    !------------------------------------------------
    IF (pertt .EQ. 0) THEN
      CALL kelvin_perturbation(lon, lat, z, uu, vv, tt)
      u = u + uu
      v = v + vv
      t = t + tt

    ELSEIF (pertt .EQ. 1) THEN
      CALL equatorial_rossby_perturbation(lon, lat, z, uu, vv)
      u = u + uu
      v = v + vv
      ! Optional: add a small temperature perturbation (commented)
      ! t = t + 0.25d0*pert_tamp * antisym_lat(lat) * cos(dble(pert_m)*(lon-pert_lon0)) * vertical_taper(z)
    END IF

    !------------------------------------------------
    ! Density / virtual potential temperature
    !------------------------------------------------
    rho = p / (Rd * t)
    thetav = t * (1.d0 + 0.61d0*q) * (p0 / p)**(Rd / cp)

  END SUBROUTINE tropical_wave_test


  !---------------------------------------------------------------------
  ! Kelvin: symmetric G(lat), v=0, (u',T') in-phase
  !---------------------------------------------------------------------
  SUBROUTINE kelvin_perturbation(lon, lat, z, du, dv, dT)
    IMPLICIT NONE
    REAL(8), INTENT(IN)  :: lon, lat, z
    REAL(8), INTENT(OUT) :: du, dv, dT

    REAL(8) :: Glat, Wlon, S

    Glat = EXP( - (lat/pert_lat0)**2 )
    Wlon = COS( DBLE(pert_m) * (lon - pert_lon0) )
    S    = vertical_taper(z)

    du = pert_uamp * Glat * Wlon * S
    dv = 0.d0
    dT = pert_tamp * Glat * Wlon * S
  END SUBROUTINE kelvin_perturbation


  !---------------------------------------------------------------------
  ! Equatorial Rossby-like: antisymmetric streamfunction ψ -> u,v
  ! (uses the same central-difference pattern as your baroclinic pertt=1)
  !---------------------------------------------------------------------
  SUBROUTINE equatorial_rossby_perturbation(lon, lat, z, du, dv)
    IMPLICIT NONE
    REAL(8), INTENT(IN)  :: lon, lat, z
    REAL(8), INTENT(OUT) :: du, dv

    REAL(8) :: psi_p, psi_m, psi_e, psi_w

    psi_p = equatorial_streamfunction(lon, lat + dxepsilon, z)
    psi_m = equatorial_streamfunction(lon, lat - dxepsilon, z)
    psi_e = equatorial_streamfunction(lon + dxepsilon, lat, z)
    psi_w = equatorial_streamfunction(lon - dxepsilon, lat, z)

    du = - (psi_p - psi_m) / (2.d0 * dxepsilon)
    dv =   (psi_e - psi_w) / (2.d0 * dxepsilon * COS(lat))
  END SUBROUTINE equatorial_rossby_perturbation


  !---------------------------------------------------------------------
  ! Antisymmetric latitude envelope
  !---------------------------------------------------------------------
  REAL(8) FUNCTION antisym_lat(lat)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lat
    antisym_lat = (lat/pert_lat0) * EXP( - (lat/pert_lat0)**2 )
  END FUNCTION antisym_lat


  !---------------------------------------------------------------------
  ! Streamfunction ψ for equatorial Rossby-like perturbation
  !---------------------------------------------------------------------
  REAL(8) FUNCTION equatorial_streamfunction(lon, lat, z)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: lon, lat, z

    REAL(8) :: Hlat, Wlon, S

    Hlat = antisym_lat(lat)
    Wlon = COS( DBLE(pert_m) * (lon - pert_lon0) )
    S    = vertical_taper(z)

    equatorial_streamfunction = - pert_psiamp * Hlat * Wlon * S
  END FUNCTION equatorial_streamfunction


  !---------------------------------------------------------------------
  ! Vertical taper: same cubic polynomial taper as your existing code
  !  - 1 at z=0
  !  - smoothly to 0 at z=pert_z
  !  - 0 above pert_z
  !---------------------------------------------------------------------
  REAL(8) FUNCTION vertical_taper(z)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z
    REAL(8) :: zz
    IF (z < pert_z) THEN
      zz = z / pert_z
      vertical_taper = 1.d0 - 3.d0*zz**2 + 2.d0*zz**3
    ELSE
      vertical_taper = 0.d0
    END IF
  END FUNCTION vertical_taper

  !-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    
    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    REAL(8), INTENT(OUT) :: &
                z,          & ! Altitude (m)
                t             ! Temperature (K)

    INTEGER :: ix

    REAL(8) :: z0, z1, z2
    REAL(8) :: p0, p1, p2

    z0 = 0.d0
    z1 = 10000.d0

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z0, p0, t)
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z1, p1, t)

    DO ix = 1, 100
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)

      CALL evaluate_pressure_temperature(deep, X, lon, lat, z2, p2, t)

      IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
        EXIT
      END IF

      z0 = z1
      p0 = p1

      z1 = z2
      p1 = p2
    END DO

    z = z2

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p0, t)

  END SUBROUTINE evaluate_z_temperature

  SUBROUTINE evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    REAL(8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    REAL(8), INTENT(OUT) :: &
                p,          & ! Pressure (Pa)
                t             ! Temperature (K)

    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constA, constB, constC, constH, scaledZ
    REAL(8) :: tau1, tau2, inttau1, inttau2
    REAL(8) :: rratio, inttermT

    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
            + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**K &
             - K / (K + 2.d0) * (rratio * cos(lat))**(K + 2.d0)

    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    t = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

  END SUBROUTINE evaluate_pressure_temperature

END MODULE tropical_wave_test_mod
