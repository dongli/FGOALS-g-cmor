module cmor_fgoals_g_mod

  use cmor_user_functions
  use netcdf

  implicit none

  private

  public cmor_fgoals_g_init
  public cmor_fgoals_g_run
  public cmor_fgoals_g_final

  character(256), parameter :: table_root = __FILE__ // '../table/cmip6/Tables'

contains

  subroutine cmor_fgoals_g_init()

    integer ierr

    ierr = cmor_setup(inpath=table_root, netcdf_file_action='replace')

    ierr = cmor_dataset_json(

  end subroutine cmor_fgoals_g_init

  subroutine cmor_fgoals_g_run(experiment_id, experiment_path)

    character(*), intent(in) :: experiment_id
    character(*), intent(in) :: experiment_path



  end subroutine cmor_fgoals_g_run

  subroutine cmor_fgoals_g_final()

  end subroutine cmor_fgoals_g_final

end module cmor_fgoals_g_mod
