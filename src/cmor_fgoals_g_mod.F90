module cmor_fgoals_g_mod

  use cmor_users_functions
  use flogger
  use fiona
  use datetime
  use namelist_mod, start_time_str => start_time, end_time_str => end_time
  use gamil_reader_mod

  implicit none

  private

  public cmor_fgoals_g_init
  public cmor_fgoals_g_run
  public cmor_fgoals_g_final

  character(256), parameter :: table_root = '../tables/cmip6/Tables'

  type(datetime_type) start_time
  type(datetime_type) end_time
  character(50) :: time_units = 'N/A'

  character(30), parameter :: gamil_time_interval = '10 minutes'
  real(8), parameter :: gamil_ptop = 219.4d0 ! Pa
  integer gamil_table_id
  integer gamil_lon_axis_id
  integer gamil_lat_axis_id
  integer gamil_lev_axis_id
  integer gamil_plev_axis_id
  integer gamil_time_axis_id
  integer gamil_ps_zfactor_id
  integer gamil_2d_axis_ids(3)
  integer gamil_3dp_axis_ids(4)

  character(8) :: gamil_vars(5,2)
  integer gamil_var_id(2)

  real(8), allocatable :: gamil_lon(:)
  real(8), allocatable :: gamil_lon_bnds(:)
  real(8), allocatable :: gamil_lat(:)
  real(8), allocatable :: gamil_lat_bnds(:)
  real(8), allocatable :: gamil_lev(:)
  real(8), allocatable :: gamil_lev_bnds(:)
  real(8) :: cmor_plev(19) = [ &
    100000.0d0, 92500.0d0, 85000.0d0, 70000.0d0, &
    60000.0d0, 50000.0d0, 40000.0d0, 30000.0d0,  &
    25000.0d0, 20000.0d0, 15000.0d0, 10000.0d0,  &
    7000.0d0, 5000.0d0, 3000.0d0, 2000.0d0,      &
    1000.0d0, 500.0d0, 100.0d0                   &
  ]
  integer gamil_num_time

contains

  subroutine cmor_fgoals_g_init()

    type(timedelta_type) dt
    integer ierr

    call io_init()

    ierr = cmor_setup(inpath=table_root, netcdf_file_action=cmor_replace) !, exit_control=cmor_exit_on_warning)
    call handle_cmor_error(ierr, __FILE__, __LINE__)

    start_time = create_datetime(start_time_str, '%Y', calendar=datetime_noleap_calendar)
    end_time   = create_datetime(end_time_str  , '%Y', calendar=datetime_noleap_calendar)

    dt = end_time - start_time
    gamil_num_time = dt%total_seconds() / 86400 / 365 * 12

    ! Get the first model data to inquire dimension information.
    call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // start_time%format('%Y-%m') // '.nc')
    call gamil_reader_get_grids(gamil_lon, gamil_lon_bnds, gamil_lat, gamil_lat_bnds, gamil_lev, gamil_lev_bnds)
    call gamil_reader_get_att('time', 'units', time_units)
    call gamil_reader_close()

    gamil_vars(:,1) = ['PS      ', 'ps      ', 'Pa      ', '2D      ', '        ']
    gamil_vars(:,2) = ['T       ', 'ta      ', 'K       ', '3D      ', 'LOG     ']

  end subroutine cmor_fgoals_g_init

  subroutine cmor_fgoals_g_run()

    type(datetime_type) time
    type(timedelta_type) dt
    integer ierr, i, j

    real(8) time_axis_value(1), time_axis_bnds(2,1)
    real(8), allocatable :: array_2d(:,:)
    real(8), allocatable :: array_3d(:,:,:)

    ierr = cmor_dataset_json('../src/CMIP6_' // trim(experiment_id) // '.json')
    call handle_cmor_error(ierr, __FILE__, __LINE__)
    call log_notice('Create dataset for ' // trim(experiment_id) // '.')

    gamil_table_id = cmor_load_table(trim(table_root) // '/CMIP6_Amon.json')
    call log_notice('Load table CMIP6_Amon.json.')

    call cmor_set_table(gamil_table_id)

    gamil_time_axis_id = cmor_axis(     &
      table_entry='time'              , &
      length=gamil_num_time           , &
      units=time_units                , &
      interval=gamil_time_interval)

    gamil_lon_axis_id = cmor_axis(      &
      table_entry='longitude'         , &
      length=size(gamil_lon)          , &
      units='degrees_east'            , &
      coord_vals=gamil_lon            , &
      cell_bounds=gamil_lon_bnds)

    gamil_lat_axis_id = cmor_axis(      &
      table_entry='latitude'          , &
      length=size(gamil_lat)          , &
      units='degrees_north'           , &
      coord_vals=gamil_lat            , &
      cell_bounds=gamil_lat_bnds)

    gamil_lev_axis_id = cmor_axis(      &
      table_entry='standard_sigma'    , &
      units='1'                       , &
      length=size(gamil_lev)          , &
      coord_vals=gamil_lev            , &
      cell_bounds=gamil_lev_bnds)

    gamil_plev_axis_id = cmor_axis(     &
      table_entry='plev19'            , &
      length=size(cmor_plev)          , &
      units='Pa'                      , &
      coord_vals=cmor_plev)

    gamil_2d_axis_ids  = [gamil_lon_axis_id,gamil_lat_axis_id,gamil_time_axis_id]
    gamil_3dp_axis_ids = [gamil_lon_axis_id,gamil_lat_axis_id,gamil_plev_axis_id,gamil_time_axis_id]

    ierr = cmor_zfactor(                &
      zaxis_id=gamil_lev_axis_id      , &
      zfactor_name='ptop'             , &
      units='Pa'                      , &
      zfactor_values=gamil_ptop)

    gamil_ps_zfactor_id = cmor_zfactor( &
      zaxis_id=gamil_lev_axis_id      , &
      zfactor_name='ps'               , &
      units='Pa'                      , &
      axis_ids=gamil_2d_axis_ids)

    do i = 1, size(gamil_vars, 2)
      write(*, *) '* ', trim(gamil_vars(1,i)), ' -> ', trim(gamil_vars(2,i))
      select case (gamil_vars(4,i))
      case ('2D')
        gamil_var_id(i) = cmor_variable(  &
          table_entry=gamil_vars(2,i)   , &
          units=gamil_vars(3,i)         , &
          axis_ids=gamil_2d_axis_ids)
      case ('3D')
        gamil_var_id(i) = cmor_variable(  &
          table_entry=gamil_vars(2,i)   , &
          units=gamil_vars(3,i)         , &
          axis_ids=gamil_3dp_axis_ids)
      end select
    end do

    allocate(array_2d(size(gamil_lon),size(gamil_lat)))
    allocate(array_3d(size(gamil_lon),size(gamil_lat),size(gamil_lev)))

    time = start_time
    dt = timedelta(months=1)
    do j = 1, 10 ! gamil_num_time
      call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // time%format('%Y-%m') // '.nc')
      call gamil_reader_get_var('time', time_axis_value(1))
      call gamil_reader_get_var('time_bnds', time_axis_bnds(:,1))
      do i = 1, size(gamil_vars, 2)
        select case (gamil_vars(4,i))
        case ('2D')
          call gamil_reader_get_var(gamil_vars(1,i), array_2d)
          ierr = cmor_write(            &
            var_id=gamil_var_id(i)    , &
            data=array_2d             , &
            ntimes_passed=1           , &
            time_vals=time_axis_value , &
            time_bnds=time_axis_bnds)
          if (gamil_vars(1,i) == 'PS') then
            ierr = cmor_write(            &
              var_id=gamil_ps_zfactor_id, &
              data=array_2d             , &
              ntimes_passed=1           , &
              time_vals=time_axis_value , &
              time_bnds=time_axis_bnds  , &
              store_with=gamil_var_id(i))
          end if
        case ('3D')
          call gamil_reader_get_var(gamil_vars(1,i), array_3d, plev=cmor_plev, use_log_linear=gamil_vars(5,i) == 'LOG')
          ierr = cmor_write(            &
            var_id=gamil_var_id(i)    , &
            data=array_3d             , &
            ntimes_passed=1           , &
            time_vals=time_axis_value , &
            time_bnds=time_axis_bnds)
        end select
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

end module cmor_fgoals_g_mod
