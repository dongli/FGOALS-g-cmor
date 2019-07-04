module namelist_mod

  implicit none

  character(30) experiment_id
  character(256) experiment_path
  character(30) case_id
  character(30) start_time
  character(30) end_time

  namelist /cmor_fgoals_g/ &
    experiment_id,         &
    experiment_path,       &
    case_id,               &
    start_time,            &
    end_time

end module namelist_mod
