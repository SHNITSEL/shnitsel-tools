#!/bin/bash -l 
###############################################
##serial calculation of trajectories. 
##Super small systems only
###############################################
#SBATCH --job-name=traj
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV

module purge
module load python/3.9-anaconda
module load intel/2021.4.0
module load intelmpi/2021.6.0
module load mkl/2021.4.0
#module load hdf5/1.12.2-intel-impi
conda activate molcas

sh all_run_traj.sh

wait
