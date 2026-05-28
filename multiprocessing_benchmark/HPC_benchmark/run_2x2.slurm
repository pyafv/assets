#!/bin/bash -l
#SBATCH --job-name=pyafv_4
#SBATCH --time=2:0:0
#SBATCH --account=bcamley1
#SBATCH --mail-type=end,fail,time_limit
#SBATCH --mail-user=wwang138@jh.edu
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
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

python benchmark_parallel_2x2.py

wait
date
