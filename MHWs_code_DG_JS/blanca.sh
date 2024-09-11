#!/bin/bash
# see https://curc.readthedocs.io/en/latest/running-jobs/job-resources.html#slurm-resource-flags for an explanation of these
#SBATCH --nodes=1
#SBATCH --time=168:00:00 #168
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --mem=258G #258
#SBATCH --ntasks=8 #8
#SBATCH --job-name= ECCO_daily_zlev0_gap0_2_1992_1997 # ECCO_daily_zlev0_gap0_2_1992_1997 ||| OISST_daily_1992_2017  ||| ECCO_daily_ohc_k0_k5_gap0_2_1992_1997
#SBATCH --output= ECCO_daily_zlev0_gap0_2_1992_1997.out # ECCO_daily_zlev0_gap0_2_1992_1997 ||| OISST_daily_1992_2017 ||| ECCO_daily_ohc_k0_k5_gap0_2_1992_1997
#SBATCH --mail-user=jasa1084@colorado.edu
#SBATCH --mail-type=ALL

module purge

module load matlab

matlab -nodesktop -nodisplay -r 'clear; run_daily_cases_BLANCA_ECCO_daily_zlev0;'



