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
  character(30) :: frequencies(30) = ''
  character(3 ) :: gamil_hist_tags(30) = ''
  integer       :: gamil_steps_per_file(30) = 1
  character(30) :: gamil_time_formats(30) = '%Y-%m'
  character(30) start_time
  character(30) end_time

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
    end_time

contains

  subroutine namelist_parse(file_path)

    character(*), intent(in) :: file_path

    integer i

    open(10, file=file_path, status='old')
    read(10, nml=cmor_fgoals_g)
    close(10)

    do i = 1, size(frequencies)
      if (frequencies(i) == '') cycle
      if (.not. any(['Amon  ','day   ', '6hrLev', '3hr   '] == frequencies(i))) then
        call log_error('Invalid namelist option frequency "' // trim(frequencies(i)) // '"! Choose Amon, 3hr, 6hrLev and day.')
      end if
    end do

  end subroutine namelist_parse

end module namelist_mod
