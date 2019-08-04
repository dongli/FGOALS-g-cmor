! ------------------------------------------------------------------------------
! Description:
!
!   This module convert the output of FGOALS-g into CMOR compliant netCDF files.
!   Two JSON files will be read. One is the standard CMOR table file, e.g.
!   CMIP6_Amon.json, the other is the customized variable mapping between model
!   and table, e.g. gamil_vars.Amon.json.
!
! Authors:
!
!   - Li Dong <dongli@lasg.iap.ac.cn>
! ------------------------------------------------------------------------------

module cmor_fgoals_g_mod

  use cmor_users_functions
  use json_module
  use flogger
  use fiona
  use datetime
  use namelist_mod, &
    start_time_str => start_time, end_time_str => end_time
  use gamil_reader_mod, &
    gamil_lon  => lon , gamil_lon_bnds      => lon_bnds, &
    gamil_lat  => lat , gamil_lat_bnds      => lat_bnds, &
    gamil_lev  => lev , gamil_lev_bnds      => lev_bnds, &
    gamil_ptop => ptop, gamil_time_interval => time_interval

  implicit none

  private

  public cmor_fgoals_g_init
  public cmor_fgoals_g_run
  public cmor_fgoals_g_final

  character(256), parameter :: table_root = '../tables/cmip6/Tables'

  ! Axis indices
  integer, parameter :: time_axis_idx      = 1
  integer, parameter :: time1_axis_idx     = 2
  integer, parameter :: time2_axis_idx     = 3
  integer, parameter :: lon_axis_idx       = 4
  integer, parameter :: lat_axis_idx       = 5
  integer, parameter :: lev_axis_idx       = 6
  integer, parameter :: ilev_axis_idx      = 7
  integer, parameter :: plev19_axis_idx    = 8
  integer, parameter :: plev8_axis_idx     = 9
  integer, parameter :: height2m_axis_idx  = 10
  integer, parameter :: height10m_axis_idx = 11

  type var_info_type
    integer var_id
    character(30) table_var_name
    character(30) model_var_name
    character(30) :: units    = ''
    character(30) :: vinterp  = ''
    character(4 ) :: positive = ''
    character(10), allocatable :: dims(:)
  contains
    final :: var_info_final
  end type var_info_type

  type axis_bundle_type
    integer, dimension(3) :: axis_ids_2d           ! lon, lat, time
    integer, dimension(4) :: axis_ids_3d_plev19    ! lon, lat, plev19, time
    integer, dimension(4) :: axis_ids_3d_plev8     ! lon, lat, plev8 , time
    integer, dimension(4) :: axis_ids_3d_full      ! lon, lat, lev   , time
    integer, dimension(4) :: axis_ids_3d_half      ! lon, lat, ilev  , time
    integer, dimension(4) :: axis_ids_2d_height2m  ! lon, lat, time  , height2m
    integer, dimension(4) :: axis_ids_2d_height10m ! lon, lat, time  , height10m
  end type axis_bundle_type

  type model_info_type
    character(50) :: time_units = 'N/A'
    integer table_id
    integer axis_ids(11) ! See axis indices for ordering.
    type(axis_bundle_type) axes_time  ! Axes with time axis
    type(axis_bundle_type) axes_time1 ! Axes with time1 axis
    type(axis_bundle_type) axes_time2 ! Axes with time2 axis
    integer zfactor_id
    ! Variables
    integer num_var
    type(var_info_type), allocatable :: var_info(:)
  contains
    procedure :: init => model_info_init
    procedure :: create_gamil_cmor_objects => model_info_create_gamil_cmor_objects
    procedure :: clear => model_info_clear
    final model_info_final
  end type model_info_type

  type(datetime_type) start_time
  type(datetime_type) end_time
  type(timedelta_type) time_period
  integer num_time
  type(model_info_type) gamil

  real(8), parameter :: cmor_plev19(19) = [      &
    100000.0d0, 92500.0d0, 85000.0d0, 70000.0d0, &
    60000.0d0, 50000.0d0, 40000.0d0, 30000.0d0,  &
    25000.0d0, 20000.0d0, 15000.0d0, 10000.0d0,  &
    7000.0d0, 5000.0d0, 3000.0d0, 2000.0d0,      &
    1000.0d0, 500.0d0, 100.0d0                   &
  ]

  real(8), parameter :: cmor_plev8(8) = [        &
    100000.0d0, 85000.0d0, 70000.0d0, 50000.0d0, &
    25000.0d0, 10000.0d0, 5000.0d0, 1000.0d0     &
  ]

contains

  subroutine cmor_fgoals_g_init()

    integer ierr, i

    call io_init()

    ierr = cmor_setup(inpath=table_root, netcdf_file_action=cmor_replace) ! , exit_control=cmor_exit_on_warning)
    call handle_cmor_error(ierr, __FILE__, __LINE__)

    start_time = create_datetime(start_time_str, '%Y', calendar=datetime_noleap_calendar)
    end_time   = create_datetime(end_time_str  , '%Y', calendar=datetime_noleap_calendar)

    time_period = end_time - start_time

    ! Get the first model data to inquire dimension information.
    call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // start_time%format('%Y-%m') // '.nc')
    call gamil_reader_get_grids()
    call gamil_reader_get_att('time', 'units', gamil%time_units) ! NOTE: Here we assume all the time units are the same.
    call gamil_reader_close()

  end subroutine cmor_fgoals_g_init

  subroutine cmor_fgoals_g_run()

    integer i

    do i = 1, size(frequencies)
      if (frequencies(i) /= '') then
        call gamil%init(trim(table_root) // '/CMIP6_' // trim(frequencies(i)) // '.json', &
                        '../src/gamil_vars.' // trim(frequencies(i)) // '.json')
        call gamil_write(frequencies(i), gamil_hist_tags(i), gamil_steps_per_file(i), gamil_time_formats(i))
      end if
    end do

  end subroutine cmor_fgoals_g_run

  subroutine cmor_fgoals_g_final()

    integer ierr

    ierr = cmor_close()
    call handle_cmor_error(ierr, __FILE__, __LINE__)

    call gamil_reader_final()

  end subroutine cmor_fgoals_g_final

  subroutine handle_cmor_error(ierr, file, line)

    integer, intent(in) :: ierr
    character(*), intent(in) :: file
    integer, intent(in) :: line

    character(8) line_s

    write(line_s, "(I0)") line

    if (ierr /= 0) then
      call log_error('CMOR action failed!', file, line)
      stop 1
    end if

  end subroutine handle_cmor_error

  subroutine gamil_write(frequency, hist_tag, steps_per_file, time_format)

    character(*), intent(in) :: frequency
    character(*), intent(in) :: hist_tag
    integer     , intent(in) :: steps_per_file
    character(*), intent(in) :: time_format

    type(datetime_type) time
    type(timedelta_type) dt
    integer ierr, ivar, itime, last_year
    integer time_step
    real(8) time_axis_value(1), time_axis_bnds(2,1)
    real(8), allocatable, dimension(:,:  ) :: ps
    real(8), allocatable, dimension(:,:  ) :: array_2d
    real(8), allocatable, dimension(:,:,:) :: array_3d_full
    real(8), allocatable, dimension(:,:,:) :: array_3d_half
    real(8), allocatable, dimension(:,:,:) :: array_3d_plev19
    real(8), allocatable, dimension(:,:,:) :: array_3d_plev8

    allocate(ps             (size(gamil_lon),size(gamil_lat)))
    allocate(array_2d       (size(gamil_lon),size(gamil_lat)))
    allocate(array_3d_full  (size(gamil_lon),size(gamil_lat),size(gamil_lev)))
    allocate(array_3d_half  (size(gamil_lon),size(gamil_lat),size(gamil_lev_bnds)))
    allocate(array_3d_plev19(size(gamil_lon),size(gamil_lat),size(cmor_plev19)))
    allocate(array_3d_plev8 (size(gamil_lon),size(gamil_lat),size(cmor_plev8 )))

    call gamil%create_gamil_cmor_objects(frequency)

    select case (frequency)
    case ('Amon')
      dt = timedelta(months=1)
      num_time = time_period%total_seconds() / 86400 / 365 * 12
    case ('day')
      dt = timedelta(days=1)
      num_time = time_period%total_seconds() / 86400
    case ('6hrLev', '6hrPlev', '6hrPlevPt')
      dt = timedelta(hours=6)
      num_time = time_period%total_seconds() / 86400 * 4
    case ('3hr')
      dt = timedelta(hours=3)
      num_time = time_period%total_seconds() / 86400 * 8
    end select
    do ivar = 1, gamil%num_var
      if (gamil%var_info(ivar)%model_var_name == 'XXX') cycle ! Skip the incomplete variable.
      call log_notice('Convert variable ' // trim(gamil%var_info(ivar)%table_var_name) // ' ...')
      time = start_time
      time_step = 1
      last_year = time%year
      do itime = 1, num_time
        call log_print(time%isoformat())
        ! Close previous file.
        if (time%year /= last_year) then
          ierr = cmor_close(gamil%var_info(ivar)%var_id, preserve=1)
          last_year = time%year
        end if
        if (time_step == 1) then
          call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.' // trim(hist_tag) // '.' // time%format(time_format) // '.nc')
        end if
        call gamil_reader_get_var('time'     , time_axis_value(1), time_step=time_step)
        call gamil_reader_get_var('time_bnds', time_axis_bnds(:,1), time_step=time_step)
        ! FIXME: Reset time_axis_value?
        time_axis_value(1) = (time_axis_bnds(2,1) + time_axis_bnds(1,1)) * 0.5d0
        select case (size(gamil%var_info(ivar)%dims))
        case (3) ! 2D variable
          call gamil_reader_get_var(gamil%var_info(ivar)%model_var_name, array_2d, time_step=time_step)
          ierr = cmor_write(                    &
            var_id=gamil%var_info(ivar)%var_id, &
            data=array_2d                     , &
            ntimes_passed=1                   , &
            time_vals=time_axis_value         , &
            time_bnds=time_axis_bnds)
        case (4) ! 3D variable
          if (any(gamil%var_info(ivar)%dims == 'alevel')) then
            call gamil_reader_get_var('PS', ps, time_step=time_step)
            call gamil_reader_get_var(gamil%var_info(ivar)%model_var_name, array_3d_full, time_step=time_step)
            ierr = cmor_write(                     &
              var_id=gamil%var_info(ivar)%var_id , &
              data=array_3d_full                 , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds)
            ierr = cmor_write(                     &
              var_id=gamil%zfactor_id            , &
              data=ps                            , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds           , &
              store_with=gamil%var_info(ivar)%var_id)
          else if (any(gamil%var_info(ivar)%dims == 'alevhalf')) then
            call gamil_reader_get_var('PS', ps, time_step=time_step)
            call gamil_reader_get_var(gamil%var_info(ivar)%model_var_name, array_3d_half, time_step=time_step)
            ierr = cmor_write(                     &
              var_id=gamil%var_info(ivar)%var_id , &
              data=array_3d_half                 , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds)
            ierr = cmor_write(                     &
              var_id=gamil%zfactor_id            , &
              data=ps                            , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds           , &
              store_with=gamil%var_info(ivar)%var_id)
          else if (any(gamil%var_info(ivar)%dims == 'plev19')) then
            call gamil_reader_get_var(             &
              gamil%var_info(ivar)%model_var_name, &
              array_3d_plev19                    , &
              time_step=time_step                , &
              plev=cmor_plev19                   , &
              use_log_linear=gamil%var_info(ivar)%vinterp == 'log_linear')
            ierr = cmor_write(                     &
              var_id=gamil%var_info(ivar)%var_id , &
              data=array_3d_plev19               , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds)
          else if (any(gamil%var_info(ivar)%dims == 'plev8')) then
            call gamil_reader_get_var(             &
              gamil%var_info(ivar)%model_var_name, &
              array_3d_plev8                     , &
              time_step=time_step                , &
              plev=cmor_plev8                    , &
              use_log_linear=gamil%var_info(ivar)%vinterp == 'log_linear')
            ierr = cmor_write(                     &
              var_id=gamil%var_info(ivar)%var_id , &
              data=array_3d_plev8                , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds)
          else if (any(gamil%var_info(ivar)%dims == 'height2m')) then
            call gamil_reader_get_var(gamil%var_info(ivar)%model_var_name, array_2d, time_step=time_step)
            ierr = cmor_write(                     &
              var_id=gamil%var_info(ivar)%var_id , &
              data=array_2d                      , &
              ntimes_passed=1                    , &
              time_vals=time_axis_value          , &
              time_bnds=time_axis_bnds)
          end if
        end select
        time_step = time_step + 1
        if (time_step > steps_per_file) then
          time_step = 1
          call gamil_reader_close()
        end if
        time = time + dt
      end do
      ierr = cmor_close(gamil%var_info(ivar)%var_id)
    end do

    deallocate(ps             )
    deallocate(array_2d       )
    deallocate(array_3d_full  )
    deallocate(array_3d_half  )
    deallocate(array_3d_plev19)
    deallocate(array_3d_plev8 )

  end subroutine gamil_write

  subroutine model_info_init(this, table_json_file_path, model_json_file_path)

    class(model_info_type), intent(inout) :: this
    character(*)          , intent(in   ) :: table_json_file_path
    character(*)          , intent(in   ) :: model_json_file_path

    type(json_core) json
    type(json_file) f_table
    type(json_value), pointer :: table_var_entry
    type(json_value), pointer :: table_var
    type(json_file) f_model
    type(json_value), pointer :: model_var_entry
    type(json_value), pointer :: model_var
    logical found
    integer ivar, idim, num_dim
    character(kind=json_CK, len=:), allocatable :: str

    call this%clear()

    ! NOTE: Here we open TWO json files.
    ! Load table json file to inquire information.
    call f_table%load_file(filename=table_json_file_path)
    call f_table%get('variable_entry', table_var_entry)
    if (.not. associated(table_var_entry)) then
      call log_error('Failed to parse ' // trim(table_json_file_path) // ' to get variable_entry!')
    end if
    ! Load model json file to inquire variable mapping between CMOR and model.
    call f_model%load_file(filename=model_json_file_path)
    call f_model%get('.', model_var_entry)
    if (.not. associated(model_var_entry)) then
      call log_error('Failed to parse ' // trim(model_json_file_path) // ' to get model_var_entry!')
    end if
    call json%info(model_var_entry, n_children=this%num_var)
    allocate(this%var_info(this%num_var))
    do ivar = 1, this%num_var
      call json%get_child(model_var_entry, ivar, model_var)
      call json%info(model_var, name=str)
      this%var_info(ivar)%table_var_name = str
      ! Find out the variable in CMOR table.
      call json%get(table_var_entry, this%var_info(ivar)%table_var_name, table_var)
      ! Get the dimenions from CMOR table.
      call json%get(table_var, 'dimensions', str)
      num_dim = string_count(str, ' ') + 1
      if (allocated(this%var_info(ivar)%dims)) deallocate(this%var_info(ivar)%dims)
      allocate(this%var_info(ivar)%dims(num_dim))
      do idim = 1, num_dim
        this%var_info(ivar)%dims(idim) = string_split(str, idim, ' ')
      end do
      call json%get(model_var, 'var_name', str)
      this%var_info(ivar)%model_var_name = str
      call json%get(model_var, 'units', str)
      this%var_info(ivar)%units = str
      call json%get(model_var, 'vinterp', str, found)
      if (found) this%var_info(ivar)%vinterp = str
      call json%get(model_var, 'positive', str, found)
      if (found) this%var_info(ivar)%positive = str
      if (this%var_info(ivar)%model_var_name /= 'XXX') then
        write(*, '("* ", A8, " -> ", A8)', advance='no') this%var_info(ivar)%model_var_name, this%var_info(ivar)%table_var_name
        write(*, '(A10)', advance='no') this%var_info(ivar)%units
        select case (this%var_info(ivar)%positive)
        case ('up')
          write(*, '(A)', advance='no') '↑'
        case ('down')
          write(*, '(A)', advance='no') '↓'
        case default
          write(*, '(A)', advance='no') ' '
        end select
        write(*, *) this%var_info(ivar)%dims
      else
        write(*, '("* ", A8, " -> ", A8, " ", A)') this%var_info(ivar)%model_var_name, this%var_info(ivar)%table_var_name, &
          'CHECK MODEL DATA PLEASE!'
      end if
    end do

  end subroutine model_info_init

  subroutine model_info_create_gamil_cmor_objects(this, frequency)

    class(model_info_type), intent(inout) :: this
    character(*)          , intent(in   ) :: frequency

    integer ierr, i
    integer, dimension(3) :: local_axis_ids_2d
    integer, dimension(4) :: local_axis_ids_3d_full
    integer, dimension(4) :: local_axis_ids_3d_half
    integer, dimension(4) :: local_axis_ids_3d_plev19
    integer, dimension(4) :: local_axis_ids_3d_plev8
    integer, dimension(4) :: local_axis_ids_2d_height2m
    integer, dimension(4) :: local_axis_ids_2d_height10m
    type(json_value), pointer :: var

    ierr = cmor_dataset_json('../src/CMIP6_' // trim(experiment_id) // '.json')
    call handle_cmor_error(ierr, __FILE__, __LINE__)
    call log_notice('Create dataset for ' // trim(experiment_id) // '.')

    this%table_id = cmor_load_table(trim(table_root) // '/CMIP6_' // trim(frequency) // '.json')
    call log_notice('Load table CMIP6_' // trim(frequency) // '.json.')

    call cmor_set_table(this%table_id)

    ! Time axis
    this%axis_ids(time_axis_idx) = cmor_axis(   &
      table_entry='time'                      , &
      units=this%time_units                   , &
      interval=gamil_time_interval)

    ! Time1 axis
    this%axis_ids(time1_axis_idx) = cmor_axis(  &
      table_entry='time1'                     , &
      units=this%time_units                   , &
      interval=gamil_time_interval)

    ! Time2 axis
    this%axis_ids(time2_axis_idx) = cmor_axis(  &
      table_entry='time2'                     , &
      units=this%time_units                   , &
      interval=gamil_time_interval)

    ! Longitude axis
    this%axis_ids(lon_axis_idx) = cmor_axis(    &
      table_entry='longitude'                 , &
      length=size(gamil_lon)                  , &
      units='degrees_east'                    , &
      coord_vals=gamil_lon                    , &
      cell_bounds=gamil_lon_bnds)

    ! Latitude axis
    this%axis_ids(lat_axis_idx) = cmor_axis(    &
      table_entry='latitude'                  , &
      length=size(gamil_lat)                  , &
      units='degrees_north'                   , &
      coord_vals=gamil_lat                    , &
      cell_bounds=gamil_lat_bnds)

    ! Full sigma level axis
    this%axis_ids(lev_axis_idx) = cmor_axis(    &
      table_entry='standard_sigma'            , &
      units='1'                               , &
      length=size(gamil_lev)                  , &
      coord_vals=gamil_lev                    , &
      cell_bounds=gamil_lev_bnds)

    ! Half sigma level axis
    this%axis_ids(ilev_axis_idx) = cmor_axis(   &
      table_entry='standard_sigma_half'       , &
      units='1'                               , &
      length=size(gamil_lev_bnds)             , &
      coord_vals=gamil_lev_bnds)

    ! Standard 19 pressure level axis
    this%axis_ids(plev19_axis_idx) = cmor_axis( &
      table_entry='plev19'                    , &
      length=size(cmor_plev19)                , &
      units='Pa'                              , &
      coord_vals=cmor_plev19)

    ! Standard 8 pressure level axis
    this%axis_ids(plev8_axis_idx) = cmor_axis(  &
      table_entry='plev8'                     , &
      length=size(cmor_plev8)                 , &
      units='Pa'                              , &
      coord_vals=cmor_plev8)

    ! Height 2m level axis
    this%axis_ids(height2m_axis_idx) = cmor_axis( &
      table_entry='height2m'                    , &
      length=1                                  , &
      units='m'                                 , &
      coord_vals=[2])

    ! Height 10m level axis
    this%axis_ids(height10m_axis_idx) = cmor_axis( &
      table_entry='height10m'                    , &
      length=1                                   , &
      units ='m'                                 , &
      coord_vals=[10])

    this%axes_time%axis_ids_2d = [          &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(time_axis_idx)          &
    ]
    this%axes_time%axis_ids_3d_plev19 =   [ &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(plev19_axis_idx)      , &
      this%axis_ids(time_axis_idx)          &
    ]
    this%axes_time%axis_ids_3d_plev8 = [    &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(plev8_axis_idx)       , &
      this%axis_ids(time_axis_idx)          &
    ]
    this%axes_time%axis_ids_3d_full = [     &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(lev_axis_idx)         , &
      this%axis_ids(time_axis_idx)          &
    ]
    this%axes_time%axis_ids_3d_half = [     &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(ilev_axis_idx)        , &
      this%axis_ids(time_axis_idx)          &
    ]
    this%axes_time%axis_ids_2d_height2m = [ &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(time_axis_idx)        , &
      this%axis_ids(height2m_axis_idx)      &
    ]
    this%axes_time%axis_ids_2d_height10m = [&
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(time_axis_idx)        , &
      this%axis_ids(height10m_axis_idx)     &
    ]
    this%axes_time1%axis_ids_2d = [         &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(time1_axis_idx)         &
    ]
    this%axes_time1%axis_ids_3d_plev19 = [  &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(plev19_axis_idx)      , &
      this%axis_ids(time1_axis_idx)         &
    ]
    this%axes_time1%axis_ids_3d_plev8 =   [ &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(plev8_axis_idx)       , &
      this%axis_ids(time1_axis_idx)         &
    ]
    this%axes_time1%axis_ids_3d_full = [    &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(lev_axis_idx)         , &
      this%axis_ids(time1_axis_idx)         &
    ]
    this%axes_time1%axis_ids_3d_half = [    &
      this%axis_ids(lon_axis_idx)         , &
      this%axis_ids(lat_axis_idx)         , &
      this%axis_ids(ilev_axis_idx)        , &
      this%axis_ids(time1_axis_idx)         &
    ]
    this%axes_time1%axis_ids_2d_height2m = [ &
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(time1_axis_idx)        , &
      this%axis_ids(height2m_axis_idx)       &
    ]
    this%axes_time1%axis_ids_2d_height10m = [&
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(time1_axis_idx)        , &
      this%axis_ids(height10m_axis_idx)      &
    ]
    this%axes_time2%axis_ids_2d = [          &
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(time2_axis_idx)          &
    ]
    this%axes_time2%axis_ids_3d_plev19 =    [&
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(plev19_axis_idx)       , &
      this%axis_ids(time2_axis_idx)          &
    ]
    this%axes_time2%axis_ids_3d_plev8 =    [ &
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(plev8_axis_idx)        , &
      this%axis_ids(time2_axis_idx)          &
    ]
    this%axes_time2%axis_ids_3d_full = [     &
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(lev_axis_idx)          , &
      this%axis_ids(time2_axis_idx)          &
    ]
    this%axes_time2%axis_ids_3d_half = [     &
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(ilev_axis_idx)         , &
      this%axis_ids(time2_axis_idx)          &
    ]
    this%axes_time2%axis_ids_2d_height2m = [ &
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(time2_axis_idx)        , &
      this%axis_ids(height2m_axis_idx)       &
    ]
    this%axes_time2%axis_ids_2d_height10m = [&
      this%axis_ids(lon_axis_idx)          , &
      this%axis_ids(lat_axis_idx)          , &
      this%axis_ids(time2_axis_idx)        , &
      this%axis_ids(height10m_axis_idx)      &
    ]
    ierr = cmor_zfactor(                     &
      zaxis_id=this%axis_ids(lev_axis_idx) , &
      zfactor_name='ptop'                  , &
      units='Pa'                           , &
      zfactor_values=gamil_ptop)

    this%zfactor_id = cmor_zfactor(          &
      zaxis_id=this%axis_ids(lev_axis_idx) , &
      zfactor_name='ps'                    , &
      units='Pa'                           , &
      axis_ids=this%axes_time%axis_ids_2d)

    do i = 1, this%num_var
      if (this%var_info(i)%model_var_name == 'XXX') cycle
      if (any(this%var_info(i)%dims == 'time1')) then
        local_axis_ids_2d           = this%axes_time1%axis_ids_2d
        local_axis_ids_3d_full      = this%axes_time1%axis_ids_3d_full
        local_axis_ids_3d_half      = this%axes_time1%axis_ids_3d_half
        local_axis_ids_3d_plev19    = this%axes_time1%axis_ids_3d_plev19
        local_axis_ids_3d_plev8     = this%axes_time1%axis_ids_3d_plev8
        local_axis_ids_2d_height2m  = this%axes_time1%axis_ids_2d_height2m
        local_axis_ids_2d_height10m = this%axes_time1%axis_ids_2d_height10m
      else if (any(this%var_info(i)%dims == 'time2')) then
        local_axis_ids_2d           = this%axes_time2%axis_ids_2d
        local_axis_ids_3d_full      = this%axes_time2%axis_ids_3d_full
        local_axis_ids_3d_half      = this%axes_time2%axis_ids_3d_half
        local_axis_ids_3d_plev19    = this%axes_time2%axis_ids_3d_plev19
        local_axis_ids_3d_plev8     = this%axes_time2%axis_ids_3d_plev8
        local_axis_ids_2d_height2m  = this%axes_time2%axis_ids_2d_height2m
        local_axis_ids_2d_height10m = this%axes_time2%axis_ids_2d_height10m
      else
        local_axis_ids_2d           = this%axes_time%axis_ids_2d
        local_axis_ids_3d_full      = this%axes_time%axis_ids_3d_full
        local_axis_ids_3d_half      = this%axes_time%axis_ids_3d_half
        local_axis_ids_3d_plev19    = this%axes_time%axis_ids_3d_plev19
        local_axis_ids_3d_plev8     = this%axes_time%axis_ids_3d_plev8
        local_axis_ids_2d_height2m  = this%axes_time%axis_ids_2d_height2m
        local_axis_ids_2d_height10m = this%axes_time%axis_ids_2d_height10m
      end if
      select case (size(this%var_info(i)%dims))
      case (3) ! 2D variable
        if (any(['pr  ','prc ','prl ','prsn'] == this%var_info(i)%table_var_name)) then
          if (this%var_info(i)%units == 'm s-1') then
            call log_warning('Change ' // trim(this%var_info(i)%model_var_name) // ' units from m s-1 to kg m-2 s-1, 1000 will be muplied to the values.')
            this%var_info(i)%units = 'kg m-2 s-1'
          end if
        end if
        this%var_info(i)%var_id = cmor_variable(         &
          table_entry=this%var_info(i)%table_var_name  , &
          units=this%var_info(i)%units                 , &
          axis_ids=local_axis_ids_2d                   , &
          positive=this%var_info(i)%positive)
      case (4) ! 3D variable
        if (any(this%var_info(i)%dims == 'alevel')) then
          this%var_info(i)%var_id = cmor_variable(       &
            table_entry=this%var_info(i)%table_var_name, &
            units=this%var_info(i)%units               , &
            axis_ids=local_axis_ids_3d_full            , &
            positive=this%var_info(i)%positive)
        else if (any(this%var_info(i)%dims == 'alevhalf')) then
          this%var_info(i)%var_id = cmor_variable(       &
            table_entry=this%var_info(i)%table_var_name, &
            units=this%var_info(i)%units               , &
            axis_ids=local_axis_ids_3d_half            , &
            positive=this%var_info(i)%positive)
        else if (any(this%var_info(i)%dims == 'plev19')) then
          this%var_info(i)%var_id = cmor_variable(       &
            table_entry=this%var_info(i)%table_var_name, &
            units=this%var_info(i)%units               , &
            axis_ids=local_axis_ids_3d_plev19          , &
            positive=this%var_info(i)%positive)
        else if (any(this%var_info(i)%dims == 'plev8')) then
          this%var_info(i)%var_id = cmor_variable(       &
            table_entry=this%var_info(i)%table_var_name, &
            units=this%var_info(i)%units               , &
            axis_ids=local_axis_ids_3d_plev8           , &
            positive=this%var_info(i)%positive)
        else if (any(this%var_info(i)%dims == 'height2m')) then
          this%var_info(i)%var_id = cmor_variable(       &
            table_entry=this%var_info(i)%table_var_name, &
            units=this%var_info(i)%units               , &
            axis_ids=local_axis_ids_2d_height2m        , &
            positive=this%var_info(i)%positive)
        else if (any(this%var_info(i)%dims == 'height10m')) then
          this%var_info(i)%var_id = cmor_variable(       &
            table_entry=this%var_info(i)%table_var_name, &
            units=this%var_info(i)%units               , &
            axis_ids=local_axis_ids_2d_height10m       , &
            positive=this%var_info(i)%positive)
        end if
      case default
        call log_error('Internal error!')
      end select
    end do

  end subroutine model_info_create_gamil_cmor_objects

  subroutine model_info_clear(this)

    class(model_info_type), intent(inout) :: this

    if (allocated(this%var_info)) deallocate(this%var_info)

  end subroutine model_info_clear

  subroutine model_info_final(this)

    type(model_info_type), intent(inout) :: this

    call this%clear()

  end subroutine model_info_final

  subroutine var_info_final(this)

    type(var_info_type), intent(inout) :: this

    if (allocated(this%dims)) deallocate(this%dims)

  end subroutine var_info_final

end module cmor_fgoals_g_mod
