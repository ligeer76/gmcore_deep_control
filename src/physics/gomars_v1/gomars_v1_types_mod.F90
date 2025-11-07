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

module gomars_v1_types_mod

  use physics_types_mod
  use tracer_mod
  use gomars_v1_const_mod

  implicit none

  private

  public gomars_v1_state_type
  public gomars_v1_tend_type
  public physics_use_wet_tracers

  type, extends(physics_state_type) :: gomars_v1_state_type
    ! Surface pressure at time level n (Pa)
    real(r8), allocatable, dimension(:      ) :: ps_old
    ! [restart] Pressure of planetary boundary layer top (Pa)
    real(r8), allocatable, dimension(:      ) :: pcon
    ! Potential temperature of planetary boundary layer top (K)
    real(r8), allocatable, dimension(:      ) :: ptcon
    ! Stratospheric temperature (K)
    real(r8), allocatable, dimension(:      ) :: tstrat
    ! Surface CO2 ice (?)
    real(r8), allocatable, dimension(:      ) :: co2ice_sfc
    !
    real(r8), allocatable, dimension(:      ) :: latheat
    ! Tracer mass on the surface (kg m-2)
    real(r8), allocatable, dimension(:,    :) :: qsfc
    !
    real(r8), allocatable, dimension(  :    ) :: fluxupv
    !
    real(r8), allocatable, dimension(  :    ) :: fluxdnv
    !
    real(r8), allocatable, dimension(  :    ) :: fmnetv
    !
    real(r8), allocatable, dimension(  :    ) :: fluxupi
    !
    real(r8), allocatable, dimension(  :    ) :: fluxdni
    !
    real(r8), allocatable, dimension(  :    ) :: fmneti
    !
    real(r8), allocatable, dimension(:      ) :: fuptopv
    !
    real(r8), allocatable, dimension(:      ) :: fdntopv
    !
    real(r8), allocatable, dimension(:      ) :: fupsfcv
    !
    real(r8), allocatable, dimension(:      ) :: fdnsfcv
    !
    real(r8), allocatable, dimension(:      ) :: fuptopi
    !
    real(r8), allocatable, dimension(:      ) :: fupsfci
    !
    real(r8), allocatable, dimension(:      ) :: fdnsfci
    ! [restart] Downward diffuse visible (solar) flux at the surface (W m-2)
    real(r8), allocatable, dimension(:      ) :: vsdif_sfc_dn
    ! Downward visible flux at the surface (W m-2)
    real(r8), allocatable, dimension(:      ) :: vsflx_sfc_dn
    ! Downward IR flux at the surface (W m-2)
    real(r8), allocatable, dimension(:      ) :: irflx_sfc_dn
    !
    real(r8), allocatable, dimension(:      ) :: tausurf
    ! Delta-Eddington optical depth (???)
    real(r8), allocatable, dimension(:,  :,:) :: detau
    ! Downward solar flux at the surface (W m-2)
    real(r8), allocatable, dimension(    :  ) :: solar_sfc_dn   ! solar
    ! Total absorption of solar energy by the atmosphere (?)
    real(r8), allocatable, dimension(:      ) :: ssun
    !
    real(r8), allocatable, dimension(:      ) :: fluxsfc
    ! Heat rate on full levels (including TOA) due to radiation (K s-1)
    real(r8), allocatable, dimension(:,:    ) :: qrad
    ! Heat rate at the surface due to PBL (positive upward) (W m-2)
    real(r8), allocatable, dimension(:      ) :: ht_pbl
    ! Exchange of heat between surface and air (positive downward) (W m-2)
    real(r8), allocatable, dimension(:      ) :: ht_sfc         ! fa
    ! Surface albedo from data
    real(r8), allocatable, dimension(:      ) :: alsp
    ! Surface albedo considering surface CO2 ice
    real(r8), allocatable, dimension(:      ) :: als
    ! Flag of northern GRS ground ice
    logical , allocatable, dimension(:      ) :: grsn
    ! Flag of southern GRS ground ice
    logical , allocatable, dimension(:      ) :: grss
    ! Flag of northern polar cap of water ice
    logical , allocatable, dimension(:      ) :: npcflag
    ! cpd * rho * ustar * cdh
    real(r8), allocatable, dimension(:      ) :: rhouch
    ! Squared wind shear (s-1)
    real(r8), allocatable, dimension(:,:    ) :: shr2
    ! Gradient Richardson number
    real(r8), allocatable, dimension(:,:    ) :: ri
    ! Eddy mixing coefficient for momentum (m2 s-1)
    real(r8), allocatable, dimension(:,:    ) :: km
    ! Eddy mixing coefficient for heat (m2 s-1)
    real(r8), allocatable, dimension(:,:    ) :: kh
    ! Thermal inertia (J m-2 K-1 s-1/2)
    real(r8), allocatable, dimension(:,:    ) :: zin
    real(r8), allocatable, dimension(:,:    ) :: rhosoil
    real(r8), allocatable, dimension(:,:    ) :: cpsoil
    !
    real(r8), allocatable, dimension(:,:    ) :: scond
    !
    real(r8), allocatable, dimension(:,:    ) :: stemp
    !
    real(r8), allocatable, dimension(:      ) :: zavgtg
    ! Water ice upward sublimation amount per unit area at the surface (kg m-2)
    real(r8), allocatable, dimension(:      ) :: h2osub_sfc     ! subflux
    ! [restart] Water ice at the surface (???)
    real(r8), allocatable, dimension(:      ) :: h2oice_sfc     ! gndice
    ! Wind stress dust lifting flux (kg m-2 s-1)
    real(r8), allocatable, dimension(:      ) :: dstflx_wsl
    ! Dust devil lifting flux (kg m-2 s-1)
    real(r8), allocatable, dimension(:      ) :: dstflx_ddl
    ! Square of Mars distance from sun
    real(r8) :: rsdist
    ! Volume mixing ratio of H2O on full and half levels (1)
    real(r8), allocatable, dimension(    :  ) :: qh2o_rad
    !
    real(r8), allocatable, dimension(:,    :) :: deposit
  contains
    procedure :: init  => gomars_v1_state_init
    procedure :: clear => gomars_v1_state_clear
    final gomars_v1_state_final
  end type gomars_v1_state_type

  type, extends(physics_tend_type) :: gomars_v1_tend_type
  contains
    procedure :: init  => gomars_v1_tend_init
    procedure :: clear => gomars_v1_tend_clear
    procedure :: reset => gomars_v1_tend_reset
    final gomars_v1_tend_final
  end type gomars_v1_tend_type

contains

  subroutine gomars_v1_state_init(this, mesh)

    class(gomars_v1_state_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%ps_old        (mesh%ncol                   )); this%ps_old        = 0
    allocate(this%pcon          (mesh%ncol                   )); this%pcon          = 0
    allocate(this%ptcon         (mesh%ncol                   )); this%ptcon         = 0
    allocate(this%tstrat        (mesh%ncol                   )); this%tstrat        = 0
    allocate(this%co2ice_sfc    (mesh%ncol                   )); this%co2ice_sfc    = 0
    allocate(this%latheat       (mesh%ncol                   )); this%latheat       = 0
    allocate(this%qsfc          (mesh%ncol,ntracers          )); this%qsfc          = 0
    allocate(this%fluxupv       (nlayrad                     )); this%fluxupv       = 0
    allocate(this%fluxdnv       (nlayrad                     )); this%fluxdnv       = 0
    allocate(this%fmnetv        (nlayrad                     )); this%fmnetv        = 0
    allocate(this%fluxupi       (nlayrad                     )); this%fluxupi       = 0
    allocate(this%fluxdni       (nlayrad                     )); this%fluxdni       = 0
    allocate(this%fmneti        (nlayrad                     )); this%fmneti        = 0
    allocate(this%fuptopv       (mesh%ncol                   )); this%fuptopv       = 0
    allocate(this%fdntopv       (mesh%ncol                   )); this%fdntopv       = 0
    allocate(this%fupsfcv       (mesh%ncol                   )); this%fupsfcv       = 0
    allocate(this%fdnsfcv       (mesh%ncol                   )); this%fdnsfcv       = 0
    allocate(this%fuptopi       (mesh%ncol                   )); this%fuptopi       = 0
    allocate(this%fupsfci       (mesh%ncol                   )); this%fupsfci       = 0
    allocate(this%fdnsfci       (mesh%ncol                   )); this%fdnsfci       = 0
    allocate(this%vsdif_sfc_dn  (mesh%ncol                   )); this%vsdif_sfc_dn  = 0
    allocate(this%vsflx_sfc_dn  (mesh%ncol                   )); this%vsflx_sfc_dn  = 0
    allocate(this%irflx_sfc_dn  (mesh%ncol                   )); this%irflx_sfc_dn  = 0
    allocate(this%tausurf       (mesh%ncol                   )); this%tausurf       = 0
    allocate(this%detau         (mesh%ncol,nspectv,ngauss    )); this%detau         = 0
    allocate(this%solar_sfc_dn  (mesh%ncol                   )); this%solar_sfc_dn  = 0
    allocate(this%ssun          (mesh%ncol                   )); this%ssun          = 0
    allocate(this%fluxsfc       (mesh%ncol                   )); this%fluxsfc       = 0
    allocate(this%qrad          (mesh%ncol,mesh%nlev         )); this%qrad          = 0
    allocate(this%ht_pbl        (mesh%ncol                   )); this%ht_pbl        = 0
    allocate(this%ht_sfc        (mesh%ncol                   )); this%ht_sfc        = 0
    allocate(this%alsp          (mesh%ncol                   )); this%alsp          = 0
    allocate(this%als           (mesh%ncol                   )); this%als           = 0
    allocate(this%grsn          (mesh%ncol                   )); this%grsn          = .false.
    allocate(this%grss          (mesh%ncol                   )); this%grss          = .false.
    allocate(this%npcflag       (mesh%ncol                   )); this%npcflag       = .false.
    allocate(this%rhouch        (mesh%ncol                   )); this%rhouch        = 0
    allocate(this%shr2          (mesh%ncol,mesh%nlev+1       )); this%shr2          = 0
    allocate(this%ri            (mesh%ncol,mesh%nlev+1       )); this%ri            = 0
    allocate(this%km            (mesh%ncol,mesh%nlev+1       )); this%km            = 0
    allocate(this%kh            (mesh%ncol,mesh%nlev+1       )); this%kh            = 0
    allocate(this%zin           (mesh%ncol,nsoil             )); this%zin           = 0
    allocate(this%rhosoil       (mesh%ncol,nsoil             )); this%rhosoil       = 0
    allocate(this%cpsoil        (mesh%ncol,nsoil             )); this%cpsoil        = 0
    allocate(this%scond         (mesh%ncol,nsoil+1           )); this%scond         = 0
    allocate(this%stemp         (mesh%ncol,nsoil             )); this%stemp         = 0
    allocate(this%zavgtg        (mesh%ncol                   )); this%zavgtg        = 0
    allocate(this%h2osub_sfc    (mesh%ncol                   )); this%h2osub_sfc    = 0
    allocate(this%h2oice_sfc    (mesh%ncol                   )); this%h2oice_sfc    = 0
    allocate(this%dstflx_wsl    (mesh%ncol                   )); this%dstflx_wsl    = 0
    allocate(this%dstflx_ddl    (mesh%ncol                   )); this%dstflx_ddl    = 0
    allocate(this%qh2o_rad      (2*mesh%nlev+3               )); this%qh2o_rad      = 0
    allocate(this%deposit       (mesh%ncol,ntracers          )); this%deposit       = 0

    call this%physics_state_init(mesh)

  end subroutine gomars_v1_state_init

  subroutine gomars_v1_state_clear(this)

    class(gomars_v1_state_type), intent(inout) :: this

    if (allocated(this%ps_old       )) deallocate(this%ps_old       )
    if (allocated(this%pcon         )) deallocate(this%pcon         )
    if (allocated(this%ptcon        )) deallocate(this%ptcon        )
    if (allocated(this%tstrat       )) deallocate(this%tstrat       )
    if (allocated(this%co2ice_sfc   )) deallocate(this%co2ice_sfc   )
    if (allocated(this%latheat      )) deallocate(this%latheat      )
    if (allocated(this%qsfc         )) deallocate(this%qsfc         )
    if (allocated(this%fluxupv      )) deallocate(this%fluxupv      )
    if (allocated(this%fluxdnv      )) deallocate(this%fluxdnv      )
    if (allocated(this%fmnetv       )) deallocate(this%fmnetv       )
    if (allocated(this%fluxupi      )) deallocate(this%fluxupi      )
    if (allocated(this%fluxdni      )) deallocate(this%fluxdni      )
    if (allocated(this%fmneti       )) deallocate(this%fmneti       )
    if (allocated(this%fuptopv      )) deallocate(this%fuptopv      )
    if (allocated(this%fdntopv      )) deallocate(this%fdntopv      )
    if (allocated(this%fupsfcv      )) deallocate(this%fupsfcv      )
    if (allocated(this%fdnsfcv      )) deallocate(this%fdnsfcv      )
    if (allocated(this%fuptopi      )) deallocate(this%fuptopi      )
    if (allocated(this%fupsfci      )) deallocate(this%fupsfci      )
    if (allocated(this%fdnsfci      )) deallocate(this%fdnsfci      )
    if (allocated(this%vsdif_sfc_dn )) deallocate(this%vsdif_sfc_dn )
    if (allocated(this%vsflx_sfc_dn )) deallocate(this%vsflx_sfc_dn )
    if (allocated(this%irflx_sfc_dn )) deallocate(this%irflx_sfc_dn )
    if (allocated(this%tausurf      )) deallocate(this%tausurf      )
    if (allocated(this%detau        )) deallocate(this%detau        )
    if (allocated(this%solar_sfc_dn )) deallocate(this%solar_sfc_dn )
    if (allocated(this%ssun         )) deallocate(this%ssun         )
    if (allocated(this%fluxsfc      )) deallocate(this%fluxsfc      )
    if (allocated(this%qrad         )) deallocate(this%qrad         )
    if (allocated(this%ht_pbl       )) deallocate(this%ht_pbl       )
    if (allocated(this%ht_sfc       )) deallocate(this%ht_sfc       )
    if (allocated(this%alsp         )) deallocate(this%alsp         )
    if (allocated(this%als          )) deallocate(this%als          )
    if (allocated(this%grsn         )) deallocate(this%grsn         )
    if (allocated(this%grss         )) deallocate(this%grss         )
    if (allocated(this%npcflag      )) deallocate(this%npcflag      )
    if (allocated(this%rhouch       )) deallocate(this%rhouch       )
    if (allocated(this%shr2         )) deallocate(this%shr2         )
    if (allocated(this%ri           )) deallocate(this%ri           )
    if (allocated(this%km           )) deallocate(this%km           )
    if (allocated(this%kh           )) deallocate(this%kh           )
    if (allocated(this%zin          )) deallocate(this%zin          )
    if (allocated(this%rhosoil      )) deallocate(this%rhosoil      )
    if (allocated(this%cpsoil       )) deallocate(this%cpsoil       )
    if (allocated(this%scond        )) deallocate(this%scond        )
    if (allocated(this%stemp        )) deallocate(this%stemp        )
    if (allocated(this%zavgtg       )) deallocate(this%zavgtg       )
    if (allocated(this%h2osub_sfc   )) deallocate(this%h2osub_sfc   )
    if (allocated(this%h2oice_sfc   )) deallocate(this%h2oice_sfc   )
    if (allocated(this%dstflx_wsl   )) deallocate(this%dstflx_wsl   )
    if (allocated(this%dstflx_ddl   )) deallocate(this%dstflx_ddl   )
    if (allocated(this%qh2o_rad     )) deallocate(this%qh2o_rad     )
    if (allocated(this%deposit      )) deallocate(this%deposit      )

    call this%physics_state_clear()

  end subroutine gomars_v1_state_clear

  subroutine gomars_v1_state_final(this)

    type(gomars_v1_state_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v1_state_final

  subroutine gomars_v1_tend_init(this, mesh)

    class(gomars_v1_tend_type), intent(inout) :: this
    type(physics_mesh_type), intent(in), target :: mesh

    call this%clear()

    call this%physics_tend_init(mesh)

  end subroutine gomars_v1_tend_init

  subroutine gomars_v1_tend_clear(this)

    class(gomars_v1_tend_type), intent(inout) :: this

    call this%physics_tend_clear()

  end subroutine gomars_v1_tend_clear

  subroutine gomars_v1_tend_reset(this)

    class(gomars_v1_tend_type), intent(inout) :: this

    call this%physics_tend_reset()

  end subroutine gomars_v1_tend_reset

  subroutine gomars_v1_tend_final(this)

    type(gomars_v1_tend_type), intent(inout) :: this

    call this%clear()

  end subroutine gomars_v1_tend_final

end module gomars_v1_types_mod