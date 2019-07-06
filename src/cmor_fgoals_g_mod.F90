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

  type model_info_type
    character(50) :: time_units = 'N/A'
    integer num_time
    integer table_id
    integer, allocatable, dimension(:) :: axis_ids
    integer, allocatable, dimension(:) :: axis_ids_2d
    integer, allocatable, dimension(:) :: axis_ids_3d_std
    integer, allocatable, dimension(:) :: var_ids
    integer zfactor_id
    ! Variables
    integer num_var
    character(30), allocatable, dimension(:) :: cmor_var_names
    character(30), allocatable, dimension(:) :: model_var_names
    character(30), allocatable, dimension(:) :: units
    logical      , allocatable, dimension(:) :: is_var_2d
    character(30), allocatable, dimension(:) :: vinterp
  contains
    procedure :: init => model_info_init
    procedure :: create_gamil_cmor_objects => model_info_create_gamil_cmor_objects
    final model_info_final
  end type model_info_type

  type(datetime_type) start_time
  type(datetime_type) end_time
  type(model_info_type) gamil

  real(8) :: cmor_plev(19) = [ &
    100000.0d0, 92500.0d0, 85000.0d0, 70000.0d0, &
    60000.0d0, 50000.0d0, 40000.0d0, 30000.0d0,  &
    25000.0d0, 20000.0d0, 15000.0d0, 10000.0d0,  &
    7000.0d0, 5000.0d0, 3000.0d0, 2000.0d0,      &
    1000.0d0, 500.0d0, 100.0d0                   &
  ]

contains

  subroutine cmor_fgoals_g_init()

    type(timedelta_type) dt
    integer ierr, i

    call io_init()

    ierr = cmor_setup(inpath=table_root, netcdf_file_action=cmor_replace) ! , exit_control=cmor_exit_on_warning)
    call handle_cmor_error(ierr, __FILE__, __LINE__)

    start_time = create_datetime(start_time_str, '%Y', calendar=datetime_noleap_calendar)
    end_time   = create_datetime(end_time_str  , '%Y', calendar=datetime_noleap_calendar)

    dt = end_time - start_time
    gamil%num_time = dt%total_seconds() / 86400 / 365 * 12

    ! Get the first model data to inquire dimension information.
    call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // start_time%format('%Y-%m') // '.nc')
    call gamil_reader_get_grids()
    call gamil_reader_get_att('time', 'units', gamil%time_units)
    call gamil_reader_close()

    call gamil%init('../src/gamil_vars.json')

  end subroutine cmor_fgoals_g_init

  subroutine cmor_fgoals_g_run()

    type(datetime_type) time
    type(timedelta_type) dt
    integer ierr, i, j, last_year
    real(8) time_axis_value(1), time_axis_bnds(2,1)
    real(8), allocatable :: array_2d(:,:)
    real(8), allocatable :: array_3d(:,:,:)

    allocate(array_2d(size(gamil_lon),size(gamil_lat)))
    allocate(array_3d(size(gamil_lon),size(gamil_lat),size(gamil_lev)))

    call gamil%create_gamil_cmor_objects()

    time = start_time
    last_year = time%year
    dt = timedelta(months=1)
    do j = 1, 60 ! gamil_num_time
      ! Close previous file.
      if (time%year /= last_year) then
        ! Redefine the CMOR objects.
        call gamil%create_gamil_cmor_objects()
        do i = 1, gamil%num_var
          ierr = cmor_close(gamil%var_ids(i))
          call handle_cmor_error(ierr, __FILE__, __LINE__)
        end do
        last_year = time%year
      end if
      call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // time%format('%Y-%m') // '.nc')
      call gamil_reader_get_var('time', time_axis_value(1))
      call gamil_reader_get_var('time_bnds', time_axis_bnds(:,1))
      do i = 1, gamil%num_var
        if (gamil%is_var_2d(i)) then
          call gamil_reader_get_var(gamil%model_var_names(i), array_2d)
          ierr = cmor_write(            &
            var_id=gamil%var_ids(i)   , &
            data=array_2d             , &
            ntimes_passed=1           , &
            time_vals=time_axis_value , &
            time_bnds=time_axis_bnds)
          if (gamil%model_var_names(i) == 'PS') then
            ierr = cmor_write(            &
              var_id=gamil%zfactor_id   , &
              data=array_2d             , &
              ntimes_passed=1           , &
              time_vals=time_axis_value , &
              time_bnds=time_axis_bnds  , &
              store_with=gamil%var_ids(i))
          end if
        else
          call gamil_reader_get_var(gamil%model_var_names(i), array_3d, plev=cmor_plev, use_log_linear=gamil%vinterp(i) == 'log_linear')
          ierr = cmor_write(            &
            var_id=gamil%var_ids(i)   , &
            data=array_3d             , &
            ntimes_passed=1           , &
            time_vals=time_axis_value , &
            time_bnds=time_axis_bnds)
        end if
        call handle_cmor_error(ierr, __FILE__, __LINE__)
      end do
      call gamil_reader_close()
      time = time + dt
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

  subroutine model_info_init(this, json_file_path)

    class(model_info_type), intent(inout) :: this
    character(*)          , intent(in   ) :: json_file_path

    type(json_core) json
    type(json_file) f
    type(json_value), pointer :: root, child
    logical found
    integer i
    character(kind=json_CK, len=:), allocatable :: str

    call f%load_file(filename=json_file_path)
    call f%get('.', root)
    call json%info(root, n_children=this%num_var)
    allocate(this%var_ids        (this%num_var))
    allocate(this%cmor_var_names (this%num_var))
    allocate(this%model_var_names(this%num_var))
    allocate(this%units          (this%num_var))
    allocate(this%is_var_2d      (this%num_var))
    allocate(this%vinterp        (this%num_var))
    do i = 1, this%num_var
      call json%get_child(root, i, child)
      call json%info(child, name=str)
      this%cmor_var_names(i) = str
      call json%get(child, 'var_name', str)
      this%model_var_names(i) = str
      call json%get(child, 'units', str)
      this%units(i) = str
      call json%get(child, 'is_var_2d', this%is_var_2d(i), found)
      call json%get(child, 'vinterp', str, found)
      if (found) this%vinterp(i) = str
      write(*, '("- ", A8, " -> ", A8)') this%model_var_names(i), this%cmor_var_names(i)
    end do

  end subroutine model_info_init

  subroutine model_info_create_gamil_cmor_objects(this)

    class(model_info_type), intent(inout) :: this

    integer ierr, i
    type(json_value), pointer :: var

    ierr = cmor_dataset_json('../src/CMIP6_' // trim(experiment_id) // '.json')
    call handle_cmor_error(ierr, __FILE__, __LINE__)
    call log_notice('Create dataset for ' // trim(experiment_id) // '.')

    this%table_id = cmor_load_table(trim(table_root) // '/CMIP6_Amon.json')
    call log_notice('Load table CMIP6_Amon.json.')

    call cmor_set_table(this%table_id)

    if (.not. allocated(this%axis_ids       )) allocate(this%axis_ids       (4)) ! time, lon, lat, lev
    if (.not. allocated(this%axis_ids_2d    )) allocate(this%axis_ids_2d    (3)) ! lon, lat, time
    if (.not. allocated(this%axis_ids_3d_std)) allocate(this%axis_ids_3d_std(4)) ! lon, lat, lev, time

    this%axis_ids(1) = cmor_axis(       &
      table_entry='time'              , &
      length=this%num_time            , &
      units=this%time_units           , &
      interval=gamil_time_interval)

    this%axis_ids(2) = cmor_axis(       &
      table_entry='longitude'         , &
      length=size(gamil_lon)          , &
      units='degrees_east'            , &
      coord_vals=gamil_lon            , &
      cell_bounds=gamil_lon_bnds)

    this%axis_ids(3) = cmor_axis(       &
      table_entry='latitude'          , &
      length=size(gamil_lat)          , &
      units='degrees_north'           , &
      coord_vals=gamil_lat            , &
      cell_bounds=gamil_lat_bnds)

    this%axis_ids(4) = cmor_axis(       &
      table_entry='standard_sigma'    , &
      units='1'                       , &
      length=size(gamil_lev)          , &
      coord_vals=gamil_lev            , &
      cell_bounds=gamil_lev_bnds)

    this%axis_ids(5) = cmor_axis(       &
      table_entry='plev19'            , &
      length=size(cmor_plev)          , &
      units='Pa'                      , &
      coord_vals=cmor_plev)

    this%axis_ids_2d  = [               &
      this%axis_ids(2)                , &
      this%axis_ids(3)                , &
      this%axis_ids(1)                  &
    ]
    this%axis_ids_3d_std = [            &
      this%axis_ids(2)                , &
      this%axis_ids(3)                , &
      this%axis_ids(5)                , &
      this%axis_ids(1)                  &
    ]

    ierr = cmor_zfactor(                &
      zaxis_id=this%axis_ids(4)       , &
      zfactor_name='ptop'             , &
      units='Pa'                      , &
      zfactor_values=gamil_ptop)

    this%zfactor_id = cmor_zfactor( &
      zaxis_id=this%axis_ids(4)       , &
      zfactor_name='ps'               , &
      units='Pa'                      , &
      axis_ids=this%axis_ids_2d)

    do i = 1, size(this%cmor_var_names)
      if (this%is_var_2d(i)) then
        this%var_ids(i) = cmor_variable(        &
          table_entry=this%cmor_var_names(i)  , &
          units=this%units(i)                 , &
          axis_ids=this%axis_ids_2d)
      else
        this%var_ids(i) = cmor_variable(  &
          table_entry=this%cmor_var_names(i)  , &
          units=this%units(i)                 , &
          axis_ids=this%axis_ids_3d_std)
      end if
    end do

  end subroutine model_info_create_gamil_cmor_objects

  subroutine model_info_final(this)

    type(model_info_type), intent(inout) :: this

    if (allocated(this%cmor_var_names )) deallocate(this%cmor_var_names )
    if (allocated(this%model_var_names)) deallocate(this%model_var_names)
    if (allocated(this%units          )) deallocate(this%units          )
    if (allocated(this%is_var_2d      )) deallocate(this%is_var_2d      )
    if (allocated(this%vinterp        )) deallocate(this%vinterp        )
    if (allocated(this%axis_ids       )) deallocate(this%axis_ids       )
    if (allocated(this%axis_ids_2d    )) deallocate(this%axis_ids_2d    )
    if (allocated(this%axis_ids_3d_std)) deallocate(this%axis_ids_3d_std)

  end subroutine model_info_final

end module cmor_fgoals_g_mod
