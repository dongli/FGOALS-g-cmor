module cmor_fgoals_g_mod

  use cmor_users_functions
  use netcdf

  implicit none

  private

  public cmor_fgoals_g_init
  public cmor_fgoals_g_run
  public cmor_fgoals_g_final

  character(256), parameter :: table_root = __FILE__ // '../table/cmip6/Tables'

  character(30) :: gamil_time_interval = '10 minutes'
  integer gamil_lon_axis_id
  integer gamil_lat_axis_id
  integer gamil_plev_axis_id
  integer gamil_time_axis_id

  integer gamil_num_lon
  real(8), allocatable :: gamil_lon(:)
  real(8), allocatable :: gamil_lon_bnds(:)
  integer gamil_num_lat
  real(8), allocatable :: gamil_lat(:)
  real(8), allocatable :: gamil_lat_bnds(:)
  integer, parameter :: cmor_num_plev = 19
  real(8) :: cmor_plev(cmor_num_plev) = [100000., 92500., 85000., 70000., &
              60000., 50000., 40000., 30000., 25000., 20000., &
              15000., 10000., 7000., 5000., 3000., 2000., 1000., 500., 100.]
  integer gamil_num_time

contains

  subroutine cmor_fgoals_g_init()

    integer ierr

    ierr = cmor_setup(inpath=table_root, netcdf_file_action=cmor_replace, exit_control=cmor_exit_on_warning)
    call handle_cmor_error(ierr, __FILE__, __LINE__)

    ierr = cmor_dataset_json(__FILE__ // '../src/CMIP6_piControl.json')
    call handle_cmor_error(ierr, __FILE__, __LINE__)

    gamil_lon_axis_id = cmor_axis(   &
      table='CMIP6_Amon.json', &
      table_entry='longitude', &
      length=gamil_num_lon   , &
      units='degrees_east'   , &
      coord_vals=gamil_lon   , &
      cell_bounds=gamil_lon_bnds)

    gamil_lat_axis_id = cmor_axis(   &
      table='CMIP6_Amon.json', &
      table_entry='latitude' , &
      length=gamil_num_lat   , &
      units='degrees_north'  , &
      coord_vals=gamil_lat   , &
      cell_bounds=gamil_lat_bnds)

    gamil_plev_axis_id = cmor_axis(  &
      table='CMIP6_Amon.json', &
      table_entry='plev19'   , &
      length=cmor_num_plev   , &
      units='Pa'             , &
      coord_vals=cmor_plev)

    gamil_time_axis_id = cmor_axis(  &
      table='CMIP6_Amon.json', &
      table_entry='time'     , &
      length=gamil_num_time  , &
      units=''               , &
      interval=gamil_time_interval)

  end subroutine cmor_fgoals_g_init

  subroutine cmor_fgoals_g_run(experiment_id, experiment_path)

    character(*), intent(in) :: experiment_id
    character(*), intent(in) :: experiment_path



  end subroutine cmor_fgoals_g_run

  subroutine cmor_fgoals_g_final()

  end subroutine cmor_fgoals_g_final

  subroutine handle_cmor_error(ierr, file, line)

    integer, intent(in) :: ierr
    character(*), intent(in) :: file
    integer, intent(in) :: line

    character(8) line_s

    write(line_s, "(I0)") line

    if (ierr /= 0) then
      write(*, *) '[Error]: ' // trim(file) // ':' // trim(line_s) // ': CMOR action failed!'
      stop 1
    end if

  end subroutine handle_cmor_error

end module cmor_fgoals_g_mod
