#!/bin/bash
# see https://curc.readthedocs.io/en/latest/running-jobs/job-resources.html#slurm-resource-flags for an explanation of these
#SBATCH --nodes=1
#SBATCH --time=168:00:00 #168
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --mem=258G #258
#SBATCH --ntasks=8 #8
#SBATCH --job-name=job_name
#SBATCH --output=job_name.out
#SBATCH --mail-user=jasa1084@colorado.edu
#SBATCH --mail-type=ALL

module purge

module load matlab

matlab -nodesktop -nodisplay -r 'clear; NASAPO_2020_MHWs_detect_and_info_main;'




