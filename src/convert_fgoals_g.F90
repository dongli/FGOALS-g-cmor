program convert_fgoals_g

  use namelist_mod
  use cmor_fgoals_g_mod

  implicit none

  character(256) namelist_file_path

  call get_command_argument(1, namelist_file_path)

  call namelist_parse(namelist_file_path)

  call cmor_fgoals_g_init()

  call cmor_fgoals_g_run()

  call cmor_fgoals_g_final()

end program convert_fgoals_g
