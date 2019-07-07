module namelist_mod

  use flogger

  implicit none

  character(30) experiment_id
  character(256) experiment_path
  character(30) case_id
  character(30) frequency
  character(30) start_time
  character(30) end_time

  namelist /cmor_fgoals_g/  &
    experiment_id         , &
    experiment_path       , &
    case_id               , &
    frequency             , &
    start_time            , &
    end_time

contains

  subroutine namelist_parse(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=cmor_fgoals_g)
    close(10)

    if (frequency /= 'Amon' .and. frequency /= '3hr' .and. frequency /= '6hr' .and. frequency /= 'day') then
      call log_error('Invalid namelist option frequency "' // trim(frequency) // '"! Choose Amon, 3hr, 6hr and day.')
    end if

  end subroutine namelist_parse

end module namelist_mod
