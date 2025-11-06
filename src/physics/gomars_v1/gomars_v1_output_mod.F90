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

module gomars_v1_output_mod

  use fiona
  use gomars_v1_const_mod
  use gomars_v1_objects_mod

  implicit none

  private

  public gomars_v1_add_output
  public gomars_v1_output

contains

  subroutine gomars_v1_add_output(tag, dtype)

    character(*), intent(in) :: tag
    character(*), intent(in) :: dtype

    character(4) :: dims_2d  (3) = ['lon ', 'lat ', 'time']
    character(4) :: dims_3d  (4) = ['lon ', 'lat ', 'lev ', 'time']
    character(4) :: soil_dims(4) = ['lon ', 'lat ', 'soil', 'time']

    call fiona_add_dim(tag, 'soil', size=nsoil)

    call fiona_add_var(tag, 'alsp'        , long_name='Surface albedo'                               , units=''       , dim_names=dims_2d  (1:2), dtype=dtype)
    call fiona_add_var(tag, 'zin'         , long_name='Surface thermal inertia'                      , units=''       , dim_names=dims_2d  (1:2), dtype=dtype)
    call fiona_add_var(tag, 'npcflag'     , long_name='Flag of northern polar cap of water ice'      , units=''       , dim_names=dims_2d  (1:2), dtype=dtype)
    call fiona_add_var(tag, 'gnd_ice'     , long_name='Flag of GRS ground ice'                       , units=''       , dim_names=dims_2d  (1:2), dtype=dtype)
    call fiona_add_var(tag, 'tg'          , long_name='Ground temperature'                           , units='K'      , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'tstrat'      , long_name='Stratospheric temperature'                    , units='K'      , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'solar_sfc_dn', long_name='Downward solar flux at the surface'           , units='W m-2'  , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'vsflx_sfc_dn', long_name='Downward visible flux at the surface'         , units='W m-2'  , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'vsdif_sfc_dn', long_name='Downward diffused visible flux at the surface', units='W m-2'  , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'cosz'        , long_name='Cosine of solar zenith angle'                 , units=''       , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'qrad'        , long_name='Heat rate due to radiation'                   , units='K s-1'  , dim_names=dims_3d  (1:4), dtype=dtype)
    call fiona_add_var(tag, 'shr2'        , long_name='Square of wind shear'                         , units='s-2'    , dim_names=dims_3d  (1:4), dtype=dtype)
    call fiona_add_var(tag, 'ustar'       , long_name='u* in similarity theory'                      , units='m s-1'  , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'tstar'       , long_name='Temperature scale'                            , units='K'      , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'co2ice_sfc'  , long_name='Ground CO2 ice'                               , units='kg m-2' , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'pcon'        , long_name='Pressure of PBL'                              , units='Pa'     , dim_names=dims_2d  (1:3), dtype=dtype)
    call fiona_add_var(tag, 'ptcon'       , long_name='Potential temperature of PBL'                 , units='K'      , dim_names=dims_2d  (1:3), dtype=dtype)

    ! Soil variables
    call fiona_add_var(tag, 'scond'       , long_name='Soil thermal conductivity'                    , units=''       , dim_names=soil_dims(1:3), dtype=dtype)
    call fiona_add_var(tag, 'stemp'       , long_name='Soil temperature'                             , units='K'      , dim_names=soil_dims     , dtype=dtype)

    ! Tendencies
    call fiona_add_var(tag, 'dudt_phys'   , long_name='Physics tendency of zonal wind'               , units='m s-2'  , dim_names=dims_3d  (1:4), dtype=dtype)
    call fiona_add_var(tag, 'dvdt_phys'   , long_name='Physics tendency of meridional wind'          , units='m s-2'  , dim_names=dims_3d  (1:4), dtype=dtype)
    call fiona_add_var(tag, 'dtdt_phys'   , long_name='Physics tendency of temperature'              , units='K s-1'  , dim_names=dims_3d  (1:4), dtype=dtype)
    call fiona_add_var(tag, 'dpsdt_phys'  , long_name='Physics tendency of surface pressure'         , units='Pa s-1' , dim_names=dims_2d  (1:3), dtype=dtype)

  end subroutine gomars_v1_add_output

  subroutine gomars_v1_output(tag, iblk)

    character(*), intent(in) :: tag
    integer, intent(in) :: iblk
    integer soil_start(3), soil_count(3)
    real(8), allocatable :: tmp(:)

    associate (mesh => objects(iblk)%mesh, state => objects(iblk)%state, tend => objects(iblk)%tend)
    soil_start = [mesh%cell_start_2d(1), mesh%cell_start_2d(2), 1    ]
    soil_count = [mesh%cell_count_2d(1), mesh%cell_count_2d(2), nsoil]
    call fiona_output(tag, 'alsp'        , reshape(state%alsp        , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'zin'         , reshape(state%zin         , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    allocate(tmp(mesh%ncol))
    where (state%npcflag)
      tmp = 1
    else where
      tmp = 0
    end where
    call fiona_output(tag, 'npcflag'     , reshape(tmp               , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    where (state%grss .or. state%grsn)
      tmp = 1
    else where
      tmp = 0
    end where
    call fiona_output(tag, 'gnd_ice'     , reshape(tmp               , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    deallocate(tmp)
    call fiona_output(tag, 'tg'          , reshape(state%tg          , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'tstrat'      , reshape(state%tstrat      , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'solar_sfc_dn', reshape(state%solar_sfc_dn, mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'vsflx_sfc_dn', reshape(state%vsflx_sfc_dn, mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'vsdif_sfc_dn', reshape(state%vsdif_sfc_dn, mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'cosz'        , reshape(state%cosz        , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'qrad'        , reshape(state%qrad        , mesh%cell_count_3d), start=mesh%cell_start_3d, count=mesh%cell_count_3d)
    call fiona_output(tag, 'shr2'        , reshape(state%shr2        , mesh%cell_count_3d), start=mesh%cell_start_3d, count=mesh%cell_count_3d)
    call fiona_output(tag, 'ustar'       , reshape(state%ustar       , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'tstar'       , reshape(state%tstar       , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'co2ice_sfc'  , reshape(state%co2ice_sfc  , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'pcon'        , reshape(state%pcon        , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    call fiona_output(tag, 'ptcon'       , reshape(state%ptcon       , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)

    ! Soil variables
    call fiona_output(tag, 'scond'       , reshape(state%scond       , soil_count        ), start=soil_start        , count=soil_count        )
    call fiona_output(tag, 'stemp'       , reshape(state%stemp       , soil_count        ), start=soil_start        , count=soil_count        )

    ! Tendencies
    call fiona_output(tag, 'dudt_phys'   , reshape(tend%dudt         , mesh%cell_count_3d), start=mesh%cell_start_3d, count=mesh%cell_count_3d)
    call fiona_output(tag, 'dvdt_phys'   , reshape(tend%dvdt         , mesh%cell_count_3d), start=mesh%cell_start_3d, count=mesh%cell_count_3d)
    call fiona_output(tag, 'dtdt_phys'   , reshape(tend%dtdt         , mesh%cell_count_3d), start=mesh%cell_start_3d, count=mesh%cell_count_3d)
    call fiona_output(tag, 'dpsdt_phys'  , reshape(tend%dpsdt        , mesh%cell_count_2d), start=mesh%cell_start_2d, count=mesh%cell_count_2d)
    end associate

  end subroutine gomars_v1_output

end module gomars_v1_output_mod