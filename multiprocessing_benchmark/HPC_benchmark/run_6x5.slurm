#!/bin/bash -l
#SBATCH --job-name=pyafv_30
#SBATCH --time=0:10:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=end,fail,time_limit
#SBATCH --mail-user=wwang138@jh.edu
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
##SBATCH -o %A_%a.out

###SBATCH -p parallel
#SBATCH -p shared
###SBATCH -p bigmem
###SBATCH --qos=qos_bigmem

### exit on error
set -e

### needed for shared partition
module load slurm

conda deactivate
conda activate myenv

python --version

python benchmark_parallel_6x5.py

wait
date
