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
    integer table_id
    integer, allocatable, dimension(:) :: axis_ids
    integer, allocatable, dimension(:) :: axis_ids_2d
    integer, allocatable, dimension(:) :: axis_ids_3d_std
    integer, allocatable, dimension(:) :: axis_ids_3d_full
    integer, allocatable, dimension(:) :: axis_ids_3d_half
    integer, allocatable, dimension(:) :: var_ids
    integer zfactor_id
    ! Variables
    integer num_var
    character(30), allocatable, dimension(:) :: cmor_var_names
    character(30), allocatable, dimension(:) :: model_var_names
    character(30), allocatable, dimension(:) :: units
    logical      , allocatable, dimension(:) :: is_var_2d
    logical      , allocatable, dimension(:) :: is_model_full_levels
    logical      , allocatable, dimension(:) :: is_model_half_levels
    character(30), allocatable, dimension(:) :: vinterp
    character(4 ), allocatable, dimension(:) :: positive
  contains
    procedure :: init => model_info_init
    procedure :: create_gamil_cmor_objects => model_info_create_gamil_cmor_objects
    final model_info_final
  end type model_info_type

  type(datetime_type) start_time
  type(datetime_type) end_time
  integer num_time
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
    num_time = dt%total_seconds() / 86400 / 365 * 12

    ! Get the first model data to inquire dimension information.
    call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // start_time%format('%Y-%m') // '.nc')
    call gamil_reader_get_grids()
    call gamil_reader_get_att('time', 'units', gamil%time_units)
    call gamil_reader_close()

    call gamil%init('../src/gamil_vars.' // trim(frequency) // '.json')

  end subroutine cmor_fgoals_g_init

  subroutine cmor_fgoals_g_run()

    call gamil_write()

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

  subroutine gamil_write()

    type(datetime_type) time
    type(timedelta_type) dt
    integer ierr, ivar, itime, last_year
    real(8) time_axis_value(1), time_axis_bnds(2,1)
    real(8), allocatable, dimension(:,:  ) :: ps
    real(8), allocatable, dimension(:,:  ) :: array_2d
    real(8), allocatable, dimension(:,:,:) :: array_3d_std
    real(8), allocatable, dimension(:,:,:) :: array_3d_full
    real(8), allocatable, dimension(:,:,:) :: array_3d_half

    allocate(ps           (size(gamil_lon),size(gamil_lat)))
    allocate(array_2d     (size(gamil_lon),size(gamil_lat)))
    allocate(array_3d_std (size(gamil_lon),size(gamil_lat),size(cmor_plev)))
    allocate(array_3d_full(size(gamil_lon),size(gamil_lat),size(gamil_lev)))
    allocate(array_3d_half(size(gamil_lon),size(gamil_lat),size(gamil_lev_bnds)))

    call gamil%create_gamil_cmor_objects()

    dt = timedelta(months=1)
    do ivar = 1, gamil%num_var
      if (gamil%model_var_names(ivar) == 'XXX') cycle
      call log_notice('Convert variable ' // trim(gamil%cmor_var_names(ivar)) // ' ...')
      time = start_time
      last_year = time%year
      do itime = 1, num_time
        call log_print(time%isoformat())
        ! Close previous file.
        if (time%year /= last_year) then
          ierr = cmor_close(gamil%var_ids(ivar), preserve=1)
          last_year = time%year
        end if
        call gamil_reader_open(trim(experiment_path) // '/' // trim(case_id) // '.gamil.h0.' // time%format('%Y-%m') // '.nc')
        call gamil_reader_get_var('time'     , time_axis_value(1))
        call gamil_reader_get_var('time_bnds', time_axis_bnds(:,1))
        time_axis_value(1) = (time_axis_bnds(2,1) + time_axis_bnds(1,1)) * 0.5d0
        if (gamil%is_var_2d(ivar)) then
          call gamil_reader_get_var(gamil%model_var_names(ivar), array_2d)
          ierr = cmor_write(                  &
            var_id=gamil%var_ids(ivar)      , &
            data=array_2d                   , &
            ntimes_passed=1                 , &
            time_vals=time_axis_value       , &
            time_bnds=time_axis_bnds)
        else
          if (gamil%is_model_full_levels(ivar)) then
            call gamil_reader_get_var('PS', ps)
            call gamil_reader_get_var(gamil%model_var_names(ivar), array_3d_full)
            ierr = cmor_write(                &
              var_id=gamil%var_ids(ivar)    , &
              data=array_3d_full            , &
              ntimes_passed=1               , &
              time_vals=time_axis_value     , &
              time_bnds=time_axis_bnds)
            ierr = cmor_write(                &
              var_id=gamil%zfactor_id       , &
              data=ps                       , &
              ntimes_passed=1               , &
              time_vals=time_axis_value     , &
              time_bnds=time_axis_bnds      , &
              store_with=gamil%var_ids(ivar))
          else if (gamil%is_model_half_levels(ivar)) then
            call gamil_reader_get_var('PS', ps)
            call gamil_reader_get_var(gamil%model_var_names(ivar), array_3d_half)
            ierr = cmor_write(                &
              var_id=gamil%var_ids(ivar)    , &
              data=array_3d_half            , &
              ntimes_passed=1               , &
              time_vals=time_axis_value     , &
              time_bnds=time_axis_bnds)
            ierr = cmor_write(                &
              var_id=gamil%zfactor_id       , &
              data=ps                       , &
              ntimes_passed=1               , &
              time_vals=time_axis_value     , &
              time_bnds=time_axis_bnds      , &
              store_with=gamil%var_ids(ivar))
          else
            call gamil_reader_get_var(gamil%model_var_names(ivar), array_3d_std, &
                                      plev=cmor_plev, use_log_linear=gamil%vinterp(ivar) == 'log_linear')
            ierr = cmor_write(                &
              var_id=gamil%var_ids(ivar)    , &
              data=array_3d_std             , &
              ntimes_passed=1               , &
              time_vals=time_axis_value     , &
              time_bnds=time_axis_bnds)
          end if
        end if
        call gamil_reader_close()
        time = time + dt
      end do
      ierr = cmor_close(gamil%var_ids(ivar))
    end do

    deallocate(ps           )
    deallocate(array_2d     )
    deallocate(array_3d_std )
    deallocate(array_3d_full)
    deallocate(array_3d_half)

  end subroutine gamil_write

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
    allocate(this%var_ids             (this%num_var))
    allocate(this%cmor_var_names      (this%num_var))
    allocate(this%model_var_names     (this%num_var))
    allocate(this%units               (this%num_var))
    allocate(this%is_var_2d           (this%num_var))
    allocate(this%is_model_full_levels(this%num_var))
    allocate(this%is_model_half_levels(this%num_var))
    allocate(this%vinterp             (this%num_var))
    allocate(this%positive            (this%num_var)); this%positive = ''
    do i = 1, this%num_var
      call json%get_child(root, i, child)
      call json%info(child, name=str)
      this%cmor_var_names(i) = str
      call json%get(child, 'var_name', str)
      this%model_var_names(i) = str
      call json%get(child, 'units', str)
      this%units(i) = str
      call json%get(child, 'is_var_2d', this%is_var_2d(i), found)
      call json%get(child, 'is_model_full_levels', this%is_model_full_levels(i), found)
      call json%get(child, 'is_model_half_levels', this%is_model_half_levels(i), found)
      call json%get(child, 'vinterp', str, found)
      if (found) this%vinterp(i) = str
      call json%get(child, 'positive', str, found)
      if (found) this%positive(i) = str
      if (this%model_var_names(i) /= 'XXX') then
        write(*, '("* ", A8, " -> ", A8)', advance='no') this%model_var_names(i), this%cmor_var_names(i)
        write(*, '(A10)', advance='no') this%units(i)
        if (this%positive(i) == 'up') then
          write(*, '(A)', advance='no') '↑'
        else if (this%positive(i) == 'down') then
          write(*, '(A)', advance='no') '↓'
        end if
        write(*, *)
      else
        write(*, '("* ", A8, " -> ", A8, " ", A)') this%model_var_names(i), this%cmor_var_names(i), &
          'CHECK MODEL DATA PLEASE!'
      end if
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

    if (.not. allocated(this%axis_ids        )) allocate(this%axis_ids        (6)) ! time, lon, lat, lev, ilev, plev19
    if (.not. allocated(this%axis_ids_2d     )) allocate(this%axis_ids_2d     (3)) ! lon, lat, time
    if (.not. allocated(this%axis_ids_3d_std )) allocate(this%axis_ids_3d_std (4)) ! lon, lat, plev19, time
    if (.not. allocated(this%axis_ids_3d_full)) allocate(this%axis_ids_3d_full(4)) ! lon, lat, lev, time
    if (.not. allocated(this%axis_ids_3d_half)) allocate(this%axis_ids_3d_half(4)) ! lon, lat, ilev, time

    ! Time axis
    this%axis_ids(1) = cmor_axis(       &
      table_entry='time'              , &
      units=this%time_units           , &
      interval=gamil_time_interval)

    ! Longitude axis
    this%axis_ids(2) = cmor_axis(       &
      table_entry='longitude'         , &
      length=size(gamil_lon)          , &
      units='degrees_east'            , &
      coord_vals=gamil_lon            , &
      cell_bounds=gamil_lon_bnds)

    ! Latitude axis
    this%axis_ids(3) = cmor_axis(       &
      table_entry='latitude'          , &
      length=size(gamil_lat)          , &
      units='degrees_north'           , &
      coord_vals=gamil_lat            , &
      cell_bounds=gamil_lat_bnds)

    ! Full sigma level axis
    this%axis_ids(4) = cmor_axis(         &
      table_entry='standard_sigma'      , &
      units='1'                         , &
      length=size(gamil_lev)            , &
      coord_vals=gamil_lev              , &
      cell_bounds=gamil_lev_bnds)

    ! Half sigma level axis
    this%axis_ids(5) = cmor_axis(         &
      table_entry='standard_sigma_half' , &
      units='1'                         , &
      length=size(gamil_lev_bnds)       , &
      coord_vals=gamil_lev_bnds)

    ! Standard pressure level axis
    this%axis_ids(6) = cmor_axis(         &
      table_entry='plev19'              , &
      length=size(cmor_plev)            , &
      units='Pa'                        , &
      coord_vals=cmor_plev)

    this%axis_ids_2d  = [                 &
      this%axis_ids(2)                  , &
      this%axis_ids(3)                  , &
      this%axis_ids(1)                    &
    ]
    this%axis_ids_3d_std = [              &
      this%axis_ids(2)                  , &
      this%axis_ids(3)                  , &
      this%axis_ids(6)                  , &
      this%axis_ids(1)                    &
    ]
    this%axis_ids_3d_full = [             &
      this%axis_ids(2)                  , &
      this%axis_ids(3)                  , &
      this%axis_ids(4)                  , &
      this%axis_ids(1)                    &
    ]
    this%axis_ids_3d_half = [             &
      this%axis_ids(2)                  , &
      this%axis_ids(3)                  , &
      this%axis_ids(5)                  , &
      this%axis_ids(1)                    &
    ]

    ierr = cmor_zfactor(                  &
      zaxis_id=this%axis_ids(4)         , &
      zfactor_name='ptop'               , &
      units='Pa'                        , &
      zfactor_values=gamil_ptop)

    this%zfactor_id = cmor_zfactor(       &
      zaxis_id=this%axis_ids(4)         , &
      zfactor_name='ps'                 , &
      units='Pa'                        , &
      axis_ids=this%axis_ids_2d)

    do i = 1, size(this%cmor_var_names)
      if (this%model_var_names(i) == 'XXX') cycle
      if (this%is_var_2d(i)) then
        if (this%positive(i) /= '') then
          this%var_ids(i) = cmor_variable(        &
            table_entry=this%cmor_var_names(i)  , &
            units=this%units(i)                 , &
            axis_ids=this%axis_ids_2d           , &
            positive=this%positive(i))
        else
          this%var_ids(i) = cmor_variable(        &
            table_entry=this%cmor_var_names(i)  , &
            units=this%units(i)                 , &
            axis_ids=this%axis_ids_2d)
        end if
      else
        if (this%is_model_full_levels(i)) then
          this%var_ids(i) = cmor_variable(  &
            table_entry=this%cmor_var_names(i)  , &
            units=this%units(i)                 , &
            axis_ids=this%axis_ids_3d_full      , &
            positive=this%positive(i))
        else if (this%is_model_half_levels(i)) then
          this%var_ids(i) = cmor_variable(  &
            table_entry=this%cmor_var_names(i)  , &
            units=this%units(i)                 , &
            axis_ids=this%axis_ids_3d_half      , &
            positive=this%positive(i))
        else
          this%var_ids(i) = cmor_variable(  &
            table_entry=this%cmor_var_names(i)  , &
            units=this%units(i)                 , &
            axis_ids=this%axis_ids_3d_std      , &
            positive=this%positive(i))
        end if
      end if
    end do

  end subroutine model_info_create_gamil_cmor_objects

  subroutine model_info_final(this)

    type(model_info_type), intent(inout) :: this

    if (allocated(this%var_ids             )) deallocate(this%var_ids             )
    if (allocated(this%cmor_var_names      )) deallocate(this%cmor_var_names      )
    if (allocated(this%model_var_names     )) deallocate(this%model_var_names     )
    if (allocated(this%units               )) deallocate(this%units               )
    if (allocated(this%is_var_2d           )) deallocate(this%is_var_2d           )
    if (allocated(this%is_model_full_levels)) deallocate(this%is_model_full_levels)
    if (allocated(this%is_model_half_levels)) deallocate(this%is_model_half_levels)
    if (allocated(this%vinterp             )) deallocate(this%vinterp             )
    if (allocated(this%positive            )) deallocate(this%positive            )
    if (allocated(this%axis_ids            )) deallocate(this%axis_ids            )
    if (allocated(this%axis_ids_2d         )) deallocate(this%axis_ids_2d         )
    if (allocated(this%axis_ids_3d_std     )) deallocate(this%axis_ids_3d_std     )
    if (allocated(this%axis_ids_3d_full    )) deallocate(this%axis_ids_3d_full    )
    if (allocated(this%axis_ids_3d_half    )) deallocate(this%axis_ids_3d_half    )

  end subroutine model_info_final

end module cmor_fgoals_g_mod
