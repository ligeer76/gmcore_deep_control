! ==============================================================================
! This file is part of GMCORE since 2019.
!
! GMCORE is a dynamical core for atmospheric model.
!
! GMCORE is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY. You may contact authors for helping or cooperation.
! ==============================================================================
! Description:
!
!   The tracer advection can be seperated into different batches. Each batch
!   can have different time step size. The wind and mass flux are accumulated
!   along model integration, and averaged to middle time level of advection time
!   step cycle.
!
!   The batch type allocates necessary arrays, and provides wind accumulation
!   subroutines.
!
! Note:
!
!   - It needs to verify the wind and mass flux accumulation manners:
!     Averaging wind and mass flux on n + 1/2 time level, or n time level.
!
! Authors:
!
!   - Li Dong (Institute of Atmospheric Physics, Chinese Academy of Sciences)
! ==============================================================================

module adv_batch_mod

  use flogger
  use container
  use const_mod
  use namelist_mod
  use math_mod
  use time_mod
  use perf_mod
  use latlon_mesh_mod
  use latlon_field_types_mod
  use latlon_operators_mod
  use latlon_parallel_mod
  use vert_coord_mod

  implicit none

  private

  public adv_batch_type
  public adv_fill_vhalo
  public adv_batch_use_deep_xy
  public adv_batch_area
  public adv_batch_face_length_x
  public adv_batch_face_length_y
  public adv_batch_edge_dx
  public adv_batch_edge_dy

  ! Different tracers can be combined into one batch, and advected in different
  ! time steps.
  type adv_batch_type
    character(strlen_scheme) :: scheme_h = 'N/A'
    character(strlen_scheme) :: scheme_v = 'N/A'
    character(strlen_loc   ) :: loc  = 'cell'
    character(strlen_name  ) :: name = ''
    logical  :: initialized = .false.
    logical  :: dynamic     = .false.
    logical  :: passive     = .true.
    logical  :: semilag     = .false.
    integer  :: ntracers    = 1
    integer  :: nstep       = 0     ! Number of dynamic steps for one adv step
    integer  :: step        = 0     ! Step counter
    real(r8) :: dt                  ! Advection time step size in seconds
    integer , allocatable :: idx(:) ! Global index of tracers in this batch
    type(latlon_mesh_type), pointer :: mesh => null()
    type(array_type) fields
    type(latlon_field3d_type) m     ! Dry-air weight at n time level
    type(latlon_field3d_type) mfx   ! Mass flux in x direction
    type(latlon_field3d_type) mfy   ! Mass flux in y direction
    type(latlon_field3d_type) mfz   ! Mass flux in z direction
    type(latlon_field3d_type) u
    type(latlon_field3d_type) v
    type(latlon_field3d_type) w     ! Only for advection tests
    type(latlon_field3d_type) qmfx
    type(latlon_field3d_type) qmfy
    type(latlon_field3d_type) qmfz
    type(latlon_field3d_type) rdp     ! Radius at batch center
    type(latlon_field3d_type) rdp_x   ! Radius at x-face
    type(latlon_field3d_type) rdp_y   ! Radius at y-face
    type(latlon_field3d_type) rdp_z   ! Radius at z-face
    ! Semi-Lagrangian variables
    type(latlon_field3d_type) u_frac
    type(latlon_field3d_type) w_frac
    type(latlon_field3d_type) mfx_frac
    type(latlon_field3d_type) mfz_frac
    type(latlon_field3d_type) cflx
    type(latlon_field3d_type) cfly
    type(latlon_field3d_type) cflz
    ! FFSL variables
    type(latlon_field3d_type) divx
    type(latlon_field3d_type) divy
    type(latlon_field3d_type) qmfx0 ! Inner tracer mass flux in x direction
    type(latlon_field3d_type) qmfy0 ! Inner tracer mass flux in y direction
    type(latlon_field3d_type) qx
    type(latlon_field3d_type) qy
    ! Background batch
    type(adv_batch_type), pointer :: bg => null()
  contains
    procedure :: init       => adv_batch_init
    procedure :: clear      => adv_batch_clear
    procedure :: copy_m_old => adv_batch_copy_m_old
    procedure :: set_wind   => adv_batch_set_wind
    procedure :: accum_wind => adv_batch_accum_wind
    procedure :: calc_cflxy_mass   => adv_batch_calc_cflxy_mass
    procedure :: calc_cflxy_tracer => adv_batch_calc_cflxy_tracer
    procedure :: calc_cflz_mass    => adv_batch_calc_cflz_mass
    procedure :: calc_cflz_tracer  => adv_batch_calc_cflz_tracer
    procedure, private :: prepare_semilag_h => adv_batch_prepare_semilag_h
    procedure, private :: prepare_semilag_v => adv_batch_prepare_semilag_v
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, filter_mesh, filter_halo, mesh, halo, scheme, batch_loc, batch_name, dt, dynamic, passive, semilag, idx, bg)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: filter_mesh
    type(latlon_halo_type), intent(in) :: filter_halo(:)
    type(latlon_mesh_type), intent(in), target :: mesh
    type(latlon_halo_type), intent(in) :: halo(:)
    character(*), intent(in) :: scheme
    character(*), intent(in) :: batch_loc
    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    logical, intent(in) :: dynamic
    logical, intent(in) :: passive
    logical, intent(in) :: semilag
    integer, intent(in), optional :: idx(:)
    class(adv_batch_type), intent(in), target, optional :: bg

    character(10) loc0, locx, locy, locz

    call this%clear()

    this%mesh => filter_mesh
    this%fields = array(30)

    this%loc      = batch_loc
    this%name     = batch_name
    this%dt       = dt
    this%dynamic  = dynamic
    this%passive  = passive
    this%semilag  = semilag
    this%nstep    = dt / dt_dyn
    this%step     = 0

    if (present(bg)) this%bg => bg

    ! Discriminate horizontal and vertical schemes.
    if (count_string(scheme, ':') == 1) then
      this%scheme_h = split_string(scheme, ':', 1)
      this%scheme_v = split_string(scheme, ':', 2)
    else
      this%scheme_h = scheme
      this%scheme_v = scheme
    end if

    select case (this%loc)
    case ('cell')
      loc0 = 'cell'
      locx = 'lon'
      locy = 'lat'
      locz = 'lev'
    case ('lev')
      loc0 = 'lev'
      locx = 'lev_lon'
      locy = 'lev_lat'
      locz = 'cell'
    case ('vtx')
      loc0 = 'vtx'
      locx = 'lat'
      locy = 'lon'
      locz = 'none'
    case default
      call log_error('Invalid grid location ' // trim(this%loc) // '!', __FILE__, __LINE__, pid=proc%id_model)
    end select

    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_m'                             , &
      long_name         ='Dry-air weight'                                    , &
      units             ='Pa'                                                , &
      loc               =loc0                                                , &
      mesh              =filter_mesh                                         , &
      halo              =filter_halo                                         , &
      halo_cross_pole   =.true.                                              , &
      output            =merge('h0', '  ', advection)                        , &
      restart           =.true.                                              , &
      field             =this%m                                              )
    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_u'                             , &
      long_name         ='U wind component'                                  , &
      units             ='m s-1'                                             , &
      loc               =locx                                                , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =merge('h0', '  ', advection)                        , &
      restart           =.false.                                             , &
      field             =this%u                                              )
    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_v'                             , &
      long_name         ='V wind component'                                  , &
      units             ='m s-1'                                             , &
      loc               =locy                                                , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            =merge('h0', '  ', advection)                        , &
      restart           =.false.                                             , &
      field             =this%v                                              )
    if (advection .and. nlev > 1 .and. locz /= 'none') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_w'                             , &
        long_name       ='Vertical coordinate velocity'                      , &
        units           ='s-1'                                               , &
        loc             =locz                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          =merge('h0', '  ', advection)                        , &
        restart         =.false.                                             , &
        field           =this%w                                              )
    end if
    if (locz /= 'none') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_mfz'                           , &
        long_name       ='Vertical mass flux'                                , &
        units           ='Pa s-1'                                            , &
        loc             =locz                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%mfz                                            )
    end if
    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_mfx'                           , &
      long_name         ='Mass flux in x direction'                          , &
      units             ='Pa m s-1'                                          , &
      loc               =locx                                                , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%mfx                                            )
    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_mfy'                           , &
      long_name         ='Mass flux in y direction'                          , &
      units             ='Pa m s-1'                                          , &
      loc               =locy                                                , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%mfy                                            )
    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_qmfx'                          , &
      long_name         ='Tracer mass flux in x direction'                   , &
      units             ='Pa kg kg-1 m s-1'                                  , &
      loc               =locx                                                , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%qmfx                                           )
    call append_field(this%fields                                            , &
      name              =trim(this%name) // '_qmfy'                          , &
      long_name         ='Tracer mass flux in y direction'                   , &
      units             ='Pa kg kg-1 m s-1'                                  , &
      loc               =locy                                                , &
      mesh              =mesh                                                , &
      halo              =halo                                                , &
      output            ='h1'                                                , &
      restart           =.false.                                             , &
      field             =this%qmfy                                           )
    if (locz /= 'none') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_qmfz'                          , &
        long_name       ='Tracer mass flux in z direction'                   , &
        units           ='Pa kg kg-1 m s-1'                                  , &
        loc             =locz                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%qmfz                                           )
    end if
    if (.not. this%passive) then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_u_frac'                        , &
        long_name       ='Fractional U wind component'                       , &
        units           ='m s-1'                                             , &
        loc             =locx                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%u_frac                                         )
    end if
    if (advection .and. nlev > 1 .and. locz /= 'none') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_w_frac'                        , &
        long_name       ='Fractional vertical coordinate velocity'           , &
        units           ='s-1'                                               , &
        loc             =locz                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%w_frac                                         )
    end if
    if (.not. this%passive .or. this%semilag) then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_mfx_frac'                      , &
        long_name       ='Fractional mass flux in x direction'               , &
        units           ='Pa m-1 s-1'                                        , &
        loc             =locx                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%mfx_frac                                       )
    end if
    if ((.not. this%passive .or. this%semilag) .and. locz /= 'none') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_mfz_frac'                      , &
        long_name       ='Fractional vertical mass flux'                     , &
        units           ='Pa m-2 s-1'                                        , &
        loc             =locz                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%mfz_frac                                       )
    end if
    if (.not. this%passive .or. this%semilag) then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_cflx'                          , &
        long_name       ='CFL number in x direction'                         , &
        units           =''                                                  , &
        loc             =locx                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%cflx                                           )
    end if
    if (.not. this%passive .or. this%semilag) then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_cfly'                          , &
        long_name       ='CFL number in y direction'                         , &
        units           =''                                                  , &
        loc             =locy                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%cfly                                           )
    end if
    if ((.not. this%passive .or. this%semilag) .and. locz /= 'none') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_cflz'                          , &
        long_name       ='CFL number in z direction'                         , &
        units           =''                                                  , &
        loc             =locz                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%cflz                                           )
    end if
    select case (this%scheme_h)
    case ('ffsl')
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_divx'                          , &
        long_name       ='Mass flux divergence in x direction'               , &
        units           ='Pa s-1'                                            , &
        loc             =loc0                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%divx                                           )
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_divy'                          , &
        long_name       ='Mass flux divergence in y direction'               , &
        units           ='Pa s-1'                                            , &
        loc             =loc0                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        output          ='h1'                                                , &
        restart         =.false.                                             , &
        field           =this%divy                                           )
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_qmfx0'                         , &
        long_name       ='Inner tracer mass flux in x direction'             , &
        units           ='Pa kg kg-1 m s-1'                                  , &
        loc             =locx                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        restart         =.false.                                             , &
        field           =this%qmfx0                                          )
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_qmfy0'                         , &
        long_name       ='Inner tracer mass flux in y direction'             , &
        units           ='Pa kg kg-1 m s-1'                                  , &
        loc             =locy                                                , &
        mesh            =mesh                                                , &
        halo            =halo                                                , &
        restart         =.false.                                             , &
        field           =this%qmfy0                                          )
    end select
    if (.not. this%passive .or. this%scheme_h == 'ffsl') then
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_qx'                            , &
        long_name       ='Tracer dry mixing ratio after advection in x direction', &
        units           ='kg kg-1'                                           , &
        loc             =loc0                                                , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        restart         =.false.                                             , &
        halo_cross_pole =.true.                                              , &
        field           =this%qx                                             )
      call append_field(this%fields                                          , &
        name            =trim(this%name) // '_qy'                            , &
        long_name       ='Tracer dry mixing ratio after advection in y direction', &
        units           ='kg kg-1'                                           , &
        loc             =loc0                                                , &
        mesh            =filter_mesh                                         , &
        halo            =filter_halo                                         , &
        restart         =.false.                                             , &
        field           =this%qy                                             )
    end if

    if (present(idx)) then
      this%ntracers = size(idx)
      allocate(this%idx(this%ntracers))
      this%idx = idx
    end if

    call time_add_alert(batch_name, seconds=dt/time_scale)

    this%initialized = .true.

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    class(*), pointer :: field
    integer i

    do i = 1, this%fields%size
      field => this%fields%value_at(i)
      select type (field)
      type is (latlon_field2d_type)
        call field%clear()
      type is (latlon_field3d_type)
        call field%clear()
      type is (latlon_field4d_type)
        call field%clear()
      end select
    end do
    call this%fields%clear()

    if (allocated(this%idx)) deallocate(this%idx)

    this%loc         = 'cell'
    this%name        = ''
    this%dt          = 0
    this%initialized = .false.
    this%dynamic     = .false.
    this%ntracers    = 0
    this%nstep       = 0
    this%step        = 0
    this%bg          => null()
    call this%rdp  %clear()
    call this%rdp_x%clear()
    call this%rdp_y%clear()
    call this%rdp_z%clear()

  end subroutine adv_batch_clear

  subroutine adv_batch_copy_m_old(this, m)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: m

    call this%m%copy(m)
    call adv_fill_vhalo(this%m, no_negvals=.true.)
    if (.not. this%passive .or. this%semilag) call fill_halo(this%m, async=.true.)

  end subroutine adv_batch_copy_m_old

  subroutine adv_batch_set_wind(this, u, v, w, mfx, mfy, mfz, m, dt, rdp, rdp_x, rdp_y, rdp_z)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: u
    type(latlon_field3d_type), intent(in) :: v
    type(latlon_field3d_type), intent(in), optional :: w
    type(latlon_field3d_type), intent(in), optional :: mfx
    type(latlon_field3d_type), intent(in), optional :: mfy
    type(latlon_field3d_type), intent(in), optional :: mfz
    type(latlon_field3d_type), intent(in), optional :: m
    real(r8), intent(in), optional :: dt
    type(latlon_field3d_type), intent(in), optional :: rdp
    type(latlon_field3d_type), intent(in), optional :: rdp_x
    type(latlon_field3d_type), intent(in), optional :: rdp_y
    type(latlon_field3d_type), intent(in), optional :: rdp_z

    call perf_start('adv_batch_set_wind')

    call this%u%link(u)
    call this%v%link(v)
    if (present(w  )) call this%w  %link(w  )
    if (present(mfx)) call this%mfx%link(mfx)
    if (present(mfy)) call this%mfy%link(mfy)
    if (present(mfz)) call this%mfz%link(mfz)
    if (present(m  )) call this%copy_m_old(m)
    if (present(rdp  )) call this%rdp  %link(rdp  )
    if (present(rdp_x)) call this%rdp_x%link(rdp_x)
    if (present(rdp_y)) call this%rdp_y%link(rdp_y)
    if (present(rdp_z)) call this%rdp_z%link(rdp_z)

    if (this%semilag) then
      call this%prepare_semilag_h(dt)
      call this%prepare_semilag_v(dt)
    end if

    call perf_stop('adv_batch_set_wind')

  end subroutine adv_batch_set_wind

  subroutine adv_batch_accum_wind(this, u, v, mfx, mfy, mfz, dt, rdp, rdp_x, rdp_y, rdp_z)

    ! FIXME: We do not need to accumulate u and v.

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: u
    type(latlon_field3d_type), intent(in) :: v
    type(latlon_field3d_type), intent(in) :: mfx
    type(latlon_field3d_type), intent(in) :: mfy
    type(latlon_field3d_type), intent(in), optional :: mfz
    real(r8), intent(in), optional :: dt
    type(latlon_field3d_type), intent(in), optional :: rdp
    type(latlon_field3d_type), intent(in), optional :: rdp_x
    type(latlon_field3d_type), intent(in), optional :: rdp_y
    type(latlon_field3d_type), intent(in), optional :: rdp_z

    if (present(rdp  )) call this%rdp  %link(rdp  )
    if (present(rdp_x)) call this%rdp_x%link(rdp_x)
    if (present(rdp_y)) call this%rdp_y%link(rdp_y)
    if (present(rdp_z)) call this%rdp_z%link(rdp_z)

    if (this%step == 0) then
      ! Reset step.
      this%u  %d = 0
      this%v  %d = 0
      this%mfx%d = 0
      this%mfy%d = 0
      if (present(mfz)) this%mfz%d = 0
      this%step = 1
    end if
    if (this%step == this%nstep) then
      ! This is the end step.
      this%u%d = (this%u%d + u%d) / this%nstep
      this%v%d = (this%v%d + v%d) / this%nstep
      call this%mfx%add(mfx, with_halo=.true.); this%mfx%d = this%mfx%d / this%nstep
      call this%mfy%add(mfy, with_halo=.true.); this%mfy%d = this%mfy%d / this%nstep
      if (present(mfz)) this%mfz%d = (this%mfz%d + mfz%d) / this%nstep
    else
      ! Accumulating.
      this%u%d = this%u%d + u%d
      this%v%d = this%v%d + v%d
      call this%mfx%add(mfx, with_halo=.true.)
      call this%mfy%add(mfy, with_halo=.true.)
      if (present(mfz)) this%mfz%d = this%mfz%d + mfz%d
    end if
    this%step = this%step + 1
    if (this%step > this%nstep) then
      this%step = 0
      if (this%semilag) then
        call this%prepare_semilag_h(dt)
        call this%prepare_semilag_v(dt)
      end if
    end if

  end subroutine adv_batch_accum_wind

  logical function adv_batch_use_deep_xy(this)

    class(adv_batch_type), intent(in) :: this

    adv_batch_use_deep_xy = deepwater .and. use_mesh_change .and. associated(this%rdp%d) .and. &
      associated(this%rdp_x%d) .and. associated(this%rdp_y%d)

  end function adv_batch_use_deep_xy

  real(r8) function adv_batch_area(this, i, j, k)

    class(adv_batch_type), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(r8) scale

    scale = 1.0_r8
    if (adv_batch_use_deep_xy(this)) scale = this%rdp%d(i,j,k) / radius

    select case (this%loc)
    case ('cell', 'lev')
      adv_batch_area = this%mesh%area_cell(j) * scale**2
    case ('vtx')
      adv_batch_area = this%mesh%area_vtx(j) * scale**2
    case default
      adv_batch_area = 0.0_r8
    end select

  end function adv_batch_area

  real(r8) function adv_batch_face_length_x(this, i, j, k)

    class(adv_batch_type), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(r8) scale

    scale = 1.0_r8
    if (adv_batch_use_deep_xy(this)) scale = this%rdp_x%d(i,j,k) / radius

    select case (this%loc)
    case ('cell', 'lev')
      adv_batch_face_length_x = this%mesh%le_lon(j) * scale
    case ('vtx')
      adv_batch_face_length_x = this%mesh%de_lat(j) * scale
    case default
      adv_batch_face_length_x = 0.0_r8
    end select

  end function adv_batch_face_length_x

  real(r8) function adv_batch_face_length_y(this, i, j, k)

    class(adv_batch_type), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(r8) scale

    scale = 1.0_r8
    if (adv_batch_use_deep_xy(this)) scale = this%rdp_y%d(i,j,k) / radius

    select case (this%loc)
    case ('cell', 'lev')
      adv_batch_face_length_y = this%mesh%le_lat(j) * scale
    case ('vtx')
      adv_batch_face_length_y = this%mesh%de_lon(j) * scale
    case default
      adv_batch_face_length_y = 0.0_r8
    end select

  end function adv_batch_face_length_y

  real(r8) function adv_batch_edge_dx(this, i, j, k)

    class(adv_batch_type), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(r8) scale

    scale = 1.0_r8
    if (adv_batch_use_deep_xy(this)) scale = this%rdp_x%d(i,j,k) / radius

    select case (this%loc)
    case ('cell', 'lev')
      adv_batch_edge_dx = this%mesh%de_lon(j) * scale
    case ('vtx')
      adv_batch_edge_dx = this%mesh%de_lat(j) * scale
    case default
      adv_batch_edge_dx = 0.0_r8
    end select

  end function adv_batch_edge_dx

  real(r8) function adv_batch_edge_dy(this, i, j, k)

    class(adv_batch_type), intent(in) :: this
    integer, intent(in) :: i, j, k
    real(r8) scale

    scale = 1.0_r8
    if (adv_batch_use_deep_xy(this)) scale = this%rdp_y%d(i,j,k) / radius

    select case (this%loc)
    case ('cell', 'lev')
      adv_batch_edge_dy = this%mesh%de_lat(j) * scale
    case ('vtx')
      adv_batch_edge_dy = this%mesh%de_lon(j) * scale
    case default
      adv_batch_edge_dy = 0.0_r8
    end select

  end function adv_batch_edge_dy

  subroutine adv_batch_calc_cflxy_mass(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    integer i, j, k
    real(r8) dx, dy

    associate (mesh   => this%mesh  , &
               u      => this%u     , & ! in
               v      => this%v     , & ! in
               cflx   => this%cflx  , & ! out
               cfly   => this%cfly  , & ! out
               u_frac => this%u_frac)   ! out
    do k = mesh%full_kds, mesh%full_kde
      do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
        do i = mesh%half_ids - 1, mesh%half_ide
          dx = adv_batch_edge_dx(this, i, j, k)
          cflx%d(i,j,k) = u%d(i,j,k) * dt / dx
          u_frac%d(i,j,k) = u%d(i,j,k) - int(cflx%d(i,j,k)) * dx / dt
        end do
      end do
      do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
        do i = mesh%full_ids, mesh%full_ide
          dy = adv_batch_edge_dy(this, i, j, k)
          cfly%d(i,j,k) = v%d(i,j,k) * dt / dy
        end do
      end do
    end do
    end associate

  end subroutine adv_batch_calc_cflxy_mass

  subroutine adv_batch_calc_cflxy_tracer(this, mx, my, mfx, mfy, cflx, cfly, mfx_frac, dt)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(inout) :: mx
    type(latlon_field3d_type), intent(inout) :: my
    type(latlon_field3d_type), intent(in) :: mfx
    type(latlon_field3d_type), intent(in) :: mfy
    type(latlon_field3d_type), intent(inout) :: cflx
    type(latlon_field3d_type), intent(inout) :: cfly
    type(latlon_field3d_type), intent(inout) :: mfx_frac
    real(r8), intent(in) :: dt

    real(r8) mc, dm, face_len
    integer ks, ke, i, j, k, l

    call perf_start('adv_batch_calc_cflxy_tracer')

    associate (mesh => this%mesh)
    select case (this%loc)
    case ('cell', 'lev')
      ks = merge(mesh%full_kds, mesh%half_kds, this%loc == 'cell')
      ke = merge(mesh%full_kde, mesh%half_kde, this%loc == 'cell')
      call wait_halo(mx)
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole
          do i = mesh%half_ids - 1, mesh%half_ide
            face_len = adv_batch_face_length_x(this, i, j, k)
            dm = mfx%d(i,j,k) * face_len * dt
            mfx_frac%d(i,j,k) = mfx%d(i,j,k)
            if (dm >= 0) then
              do l = i, mesh%full_ims, -1
                mc = mx%d(l,j,k) * adv_batch_area(this, l, j, k)
                if (dm < mc) exit
                dm = dm - mc
                mfx_frac%d(i,j,k) = mfx_frac%d(i,j,k) - mc / face_len / dt
              end do
              cflx%d(i,j,k) = i - l + dm / mc
            else
              do l = i + 1, mesh%full_ime
                mc = mx%d(l,j,k) * adv_batch_area(this, l, j, k)
                if (-dm < mc) exit
                dm = dm + mc
                mfx_frac%d(i,j,k) = mfx_frac%d(i,j,k) + mc / face_len / dt
              end do
              cflx%d(i,j,k) = i + 1 - l + dm / mc
            end if
          end do
        end do
      end do
      call wait_halo(my)
      do k = ks, ke
        do j = mesh%half_jds - merge(0, 1, mesh%has_south_pole()), mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide
            face_len = adv_batch_face_length_y(this, i, j, k)
            dm = mfy%d(i,j,k) * face_len * dt
            if (dm >= 0) then
              do l = j, mesh%full_jms, -1
                mc = my%d(i,l,k) * adv_batch_area(this, i, l, k)
                if (dm < mc) exit
                dm = dm - mc
              end do
              cfly%d(i,j,k) = j - l + dm / mc
            else
              do l = j + 1, mesh%full_jme
                mc = my%d(i,l,k) * adv_batch_area(this, i, l, k)
                if (-dm < mc) exit
                dm = dm + mc
              end do
              cfly%d(i,j,k) = j + 1 - l + dm / mc
            end if
            ! Clip CFL number that are out of range. This should be very rare and in the polar region.
            cfly%d(i,j,k) = min(max(cfly%d(i,j,k), -1.0_r8), 1.0_r8)
          end do
        end do
      end do
    case ('vtx')
      ks = mesh%full_kds
      ke = mesh%full_kde
      call wait_halo(mx)
      do k = ks, ke
        do j = mesh%half_jds, mesh%half_jde
          do i = mesh%full_ids, mesh%full_ide + 1
            face_len = adv_batch_face_length_x(this, i, j, k)
            dm = mfx%d(i,j,k) * face_len * dt
            if (dm >= 0) then
              do l = i, mesh%half_ims, -1
                mc = mx%d(l,j,k) * adv_batch_area(this, l, j, k)
                if (dm < mc) exit
                dm = dm - mc
              end do
              cflx%d(i,j,k) = i - l + dm / mc
            else
              do l = i + 1, mesh%half_ime
                mc = mx%d(l,j,k) * adv_batch_area(this, l, j, k)
                if (-dm < mc) exit
                dm = dm + mc
              end do
              cflx%d(i,j,k) = i + 1 - l + dm / mc
            end if
          end do
        end do
      end do
      call wait_halo(my)
      do k = ks, ke
        do j = mesh%full_jds_no_pole, mesh%full_jde_no_pole + merge(0, 1, mesh%has_north_pole())
          do i = mesh%half_ids, mesh%half_ide
            face_len = adv_batch_face_length_y(this, i, j, k)
            dm = mfy%d(i,j,k) * face_len * dt
            if (dm >= 0) then
              do l = j, mesh%half_jms, -1
                mc = my%d(i,l,k) * adv_batch_area(this, i, l, k)
                if (dm < mc) exit
                dm = dm - mc
              end do
              cfly%d(i,j,k) = j - l + dm / mc
            else
              do l = j + 1, mesh%half_jme
                mc = my%d(i,l,k) * adv_batch_area(this, i, l, k)
                if (-dm < mc) exit
                dm = dm + mc
              end do
              cfly%d(i,j,k) = j + 1 - l + dm / mc
            end if
            ! Clip CFL number that are out of range. This should be very rare and in the polar region.
            cfly%d(i,j,k) = min(max(cfly%d(i,j,k), -1.0_r8), 1.0_r8)
          end do
        end do
      end do
    end select
    end associate

    call perf_stop('adv_batch_calc_cflxy_tracer')

  end subroutine adv_batch_calc_cflxy_tracer

  subroutine adv_batch_calc_cflz_mass(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: dt

    real(r8) dz
    integer i, j, k, l

    if (.not. associated(this%w%d)) return

    associate (mesh   => this%mesh   , &
               w      => this%w      , & ! in
               cflz   => this%cflz   , & ! out
               w_frac => this%w_frac )   ! out
    do k = mesh%half_kds + 1, mesh%half_kde - 1
      do j = mesh%full_jds, mesh%full_jde
        do i = mesh%full_ids, mesh%full_ide
          dz = w%d(i,j,k) * dt
          w_frac%d(i,j,k) = w%d(i,j,k)
          if (dz >= 0) then
            do l = k - 1, mesh%full_kms, -1
              if (dz < mesh%full_dlev(l)) exit
              dz = dz - mesh%full_dlev(l)
              w_frac%d(i,j,k) = w_frac%d(i,j,k) - mesh%full_dlev(l) / dt
            end do
            cflz%d(i,j,k) = k - 1 - l + dz / mesh%full_dlev(l)
          else
            do l = k, mesh%full_kme
              if (dz > -mesh%full_dlev(l)) exit
              dz = dz + mesh%full_dlev(l)
              w_frac%d(i,j,k) = w_frac%d(i,j,k) + mesh%full_dlev(l) / dt
            end do
            cflz%d(i,j,k) = k - l + dz / mesh%full_dlev(l)
          end if
        end do
      end do
    end do
    end associate

  end subroutine adv_batch_calc_cflz_mass

  subroutine adv_batch_calc_cflz_tracer(this, m, mfz, cflz, mfz_frac, dt)

    class(adv_batch_type), intent(inout) :: this
    type(latlon_field3d_type), intent(in) :: m
    type(latlon_field3d_type), intent(in) :: mfz
    type(latlon_field3d_type), intent(inout) :: cflz
    type(latlon_field3d_type), intent(inout) :: mfz_frac
    real(r8), intent(in) :: dt

    real(r8) mc, dm
    integer i, j, k, l

    associate (mesh => this%mesh)
    select case (this%loc)
    case ('cell')
      do k = mesh%half_kds + 1, mesh%half_kde - 1
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            ! NOTE: Here we ignore the horizontal area of the cell since it cancels out in shallow-atmosphere approximation.
            dm = mfz%d(i,j,k) * dt ! 𝜹π/𝜹η dη/dt * dt
            mfz_frac%d(i,j,k) = mfz%d(i,j,k)
            if (dm >= 0) then
              do l = k - 1, mesh%full_kms, -1
                mc = m%d(i,j,l) ! 𝜹π/𝜹η 𝜹η = 𝜹π
                if (dm < mc) exit
                dm = dm - mc
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) - mc / dt
              end do
              cflz%d(i,j,k) = k - 1 - l + dm / mc
            else
              do l = k, mesh%full_kme
                mc = m%d(i,j,l)
                if (-dm < mc) exit
                dm = dm + mc
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) + mc / dt
              end do
              cflz%d(i,j,k) = k - l + dm / mc
            end if
          end do
        end do
      end do
    case ('lev')
      do k = mesh%full_kds, mesh%full_kde
        do j = mesh%full_jds, mesh%full_jde
          do i = mesh%full_ids, mesh%full_ide
            dm = mfz%d(i,j,k) * dt
            mfz_frac%d(i,j,k) = mfz%d(i,j,k)
            if (dm >= 0) then
              do l = k, mesh%half_kms, -1
                mc = m%d(i,j,l)
                if (dm < mc) exit
                dm = dm - mc
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) - mc / dt
              end do
              cflz%d(i,j,k) = k - l + dm / mc
            else
              do l = k + 1, mesh%half_kme
                mc = m%d(i,j,l)
                if (-dm < mc) exit
                dm = dm + mc
                mfz_frac%d(i,j,k) = mfz_frac%d(i,j,k) + mc / dt
              end do
              cflz%d(i,j,k) = k + 1 - l + dm / mc
            end if
          end do
        end do
      end do
    end select
    end associate

  end subroutine adv_batch_calc_cflz_tracer

  subroutine adv_batch_prepare_semilag_h(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt

    dt_opt = this%dt; if (present(dt)) dt_opt = dt

    associate (m        => this%m       , & ! in
               mfx      => this%mfx     , & ! in
               mfy      => this%mfy     , & ! in
               u        => this%u       , & ! in
               v        => this%v       , & ! in
               mfx_frac => this%mfx_frac, & ! out
               cflx     => this%cflx    , & ! out
               cfly     => this%cfly    , & ! out
               divx     => this%divx    , & ! out
               divy     => this%divy    )   ! out
    ! Calculate horizontal CFL number and divergence along each axis.
    if (this%passive) then
      call this%calc_cflxy_tracer(m, m, mfx, mfy, cflx, cfly, mfx_frac, dt_opt)
    else
      call this%calc_cflxy_mass(dt_opt)
    end if
    if (this%scheme_h == 'ffsl') then
      if (adv_batch_use_deep_xy(this)) then
        call divx_operator_deep(u, this%rdp_x, divx)
        call divy_operator_deep(v, this%rdp_y, divy)
      else
        call divx_operator(u, divx)
        call divy_operator(v, divy)
      end if
    end if
    end associate

  end subroutine adv_batch_prepare_semilag_h

  subroutine adv_batch_prepare_semilag_v(this, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in), optional :: dt

    real(r8) dt_opt

    dt_opt = this%dt; if (present(dt)) dt_opt = dt

    associate (m        => this%m       , & ! in
               mfz      => this%mfz     , & ! in
               mfz_frac => this%mfz_frac, & ! out
               cflz     => this%cflz    )   ! out
    if (this%passive) then
      call this%calc_cflz_tracer(m, mfz, cflz, mfz_frac, dt_opt)
    else
      call this%calc_cflz_mass(dt_opt)
    end if
    end associate

  end subroutine adv_batch_prepare_semilag_v

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

  subroutine adv_fill_vhalo(f, no_negvals)

    type(latlon_field3d_type), intent(inout) :: f
    logical, intent(in) :: no_negvals

    integer kds, kde, kms, kme, i, j, k

    select case (f%loc)
    case ('cell')
      kds = f%mesh%full_kds
      kde = f%mesh%full_kde
      kms = f%mesh%full_kms
      kme = f%mesh%full_kme
    case ('lev')
      kds = f%mesh%half_kds
      kde = f%mesh%half_kde
      kms = f%mesh%half_kms
      kme = f%mesh%half_kme
    case default
      stop 'Unhandled branch in adv_fill_vhalo!'
    end select

    ! Set upper and lower boundary conditions.
    do k = kds - 1, kms, -1
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          ! f%d(i,j,k) = f%d(i,j,kds)
          ! f%d(i,j,k) = 2 * f%d(i,j,k+1) - f%d(i,j,k+2)
          f%d(i,j,k) = 3 * f%d(i,j,k+1) - 3 * f%d(i,j,k+2) + f%d(i,j,k+3)
          ! f%d(i,j,k) = 4 * f%d(i,j,k+1) - 6 * f%d(i,j,k+2) + 4 * f%d(i,j,k+3) - f%d(i,j,k+4)
        end do
      end do
    end do
    do k = kde + 1, kme
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          ! f%d(i,j,k) = f%d(i,j,kde)
          ! f%d(i,j,k) = 2 * f%d(i,j,k-1) - f%d(i,j,k-2)
          f%d(i,j,k) = 3 * f%d(i,j,k-1) - 3 * f%d(i,j,k-2) + f%d(i,j,k-3)
          ! f%d(i,j,k) = 4 * f%d(i,j,k-1) - 6 * f%d(i,j,k-2) + 4 * f%d(i,j,k-3) - f%d(i,j,k-4)
        end do
      end do
    end do
    if (no_negvals) then
      do j = f%mesh%full_jds, f%mesh%full_jde
        do i = f%mesh%full_ids, f%mesh%full_ide
          if (any(f%d(i,j,kms:kds-1) < 0)) f%d(i,j,kms:kds-1) = f%d(i,j,kds)
          if (any(f%d(i,j,kde+1:kme) < 0)) f%d(i,j,kde+1:kme) = f%d(i,j,kde)
        end do
      end do
    end if

  end subroutine adv_fill_vhalo

end module adv_batch_mod
