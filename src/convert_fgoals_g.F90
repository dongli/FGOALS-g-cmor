program convert_fgoals_g

  use cmor_fgoals_g_mod

  implicit none

  character(30) experiment_id
  character(256) experiment_path

  call get_command_argument(1, experiment_id)
  call get_command_argument(2, experiment_path)

  call cmor_fgoals_g_init()

  call cmor_fgoals_g_run(experiment_id, experiment_path)

  call cmor_fgoals_g_final()

end program convert_fgoals_g
