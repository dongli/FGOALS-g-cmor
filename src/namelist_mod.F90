module namelist_mod

  use flogger

  implicit none

  character(30) experiment_id
  character(256) experiment_path
  character(10) :: realization_index = '1'
  character(10) :: initialization_index = '1'
  character(10) :: physics_index = '1'
  character(10) :: forcing_index = '1'
  character(20) :: branch_time_in_parent = '0.0'
  character(30) case_id
  character(30) :: frequencies(10) = ''
  character(3 ) :: gamil_hist_tags(10) = ''
  integer       :: gamil_steps_per_file(10) = 1
  character(30) :: gamil_time_formats(10) = '%Y-%m'
  character(30) :: start_time(10) = ''
  character(30) :: end_time(10) = ''
  character(10) :: selected_vars(100) = ''
  integer       :: output_interval_in_years(10) = 1

  namelist /cmor_fgoals_g/ &
    experiment_id        , &
    experiment_path      , &
    realization_index    , &
    initialization_index , &
    physics_index        , &
    forcing_index        , &
    branch_time_in_parent, &
    case_id              , &
    frequencies          , &
    gamil_hist_tags      , &
    gamil_steps_per_file , &
    gamil_time_formats   , &
    start_time           , &
    end_time             , &
    selected_vars        , &
    output_interval_in_years

contains

  subroutine namelist_parse(file_path)

    character(*), intent(in) :: file_path

    logical flag
    integer i

    inquire(file=file_path, exist=flag)
    if (.not. flag) then
      call log_error('Namelist file ' // trim(file_path) // ' does not exist!')
    end if
    open(10, file=file_path, status='old')
    read(10, nml=cmor_fgoals_g)
    close(10)

    do i = 1, size(frequencies)
      if (frequencies(i) == '') cycle
      if (.not. any(['Amon   ','day    ', '6hrLev ', '6hrPlev', '3hr    '] == frequencies(i))) then
        call log_error('Invalid namelist option frequency "' // trim(frequencies(i)) // '"! Choose Amon, 3hr, 6hrLev and day.')
      end if
    end do

    do i = 2, size(start_time)
      if (start_time(i) == '') start_time(i) = start_time(i - 1)
      if (end_time  (i) == '') end_time  (i) = end_time  (i - 1)
    end do

  end subroutine namelist_parse

end module namelist_mod
