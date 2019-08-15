#!/bin/bash
#SBATCH --job-name cmor_pi_r1
#SBATCH --partition yuyq
#SBATCH --ntasks 1
#SBATCH --nodes 1
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr

source /public1/wab/dongli/starman/setup/bashrc
starman load gcc netcdf udunits uuid cmor json-fortran -c gcc_8.2.0
./convert_fgoals_g.exe ./namelist.piControl.r1i1p1f1
