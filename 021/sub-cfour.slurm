#!/bin/bash
#SBATCH --account=rrg-tzeng
#SBATCH --mem=25G      # memory; default unit is megabytes.
#SBATCH --time=0-23:00           # time (DD-HH:MM).
#SBATCH --cpus-per-task=4
#SBATCH --output=out.%j
#SBATCH --error=err.%j

module load StdEnv/2023
module load intel/2023.2.1
module load openmpi/4.1.5
module load cfour/2.1
module load mrcc/20230828
export CFOUR_NUM_CORES=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Starting run at: `date`"

xclean
xcfour

echo "Program finished with exit code $? at: `date`"
