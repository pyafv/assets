#!/bin/bash -l
#SBATCH --job-name=test
#SBATCH --time=0:10:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=end,fail,time_limit
#SBATCH --mail-user=wwang138@jh.edu
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=48
##SBATCH -p express


### exit on error
set -e

### needed for shared partition
module load openmpi
module load slurm

conda deactivate
conda activate myenv

python --version

mpiexec -n "$SLURM_NTASKS" \
  --map-by "slot:PE=$SLURM_CPUS_PER_TASK" \
  --bind-to core \
  python test_forces_match.py

wait
date

