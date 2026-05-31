#!/bin/bash -l
#SBATCH --job-name=multiprocessing
#SBATCH --time=1:0:0
#SBATCH --account=bcamley1
##SBATCH --mail-type=end,fail,time_limit
#SBATCH --mail-user=wwang138@jh.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH -p shared


### exit on error
set -e

### needed for shared partition
module load openmpi
module load slurm

conda deactivate
conda activate myenv

python --version

python multiprocess.py

wait
date
