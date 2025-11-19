#!/bin/bash

#SBATCH -p compute
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --job-name=ORCA
#SBATCH --constraint=tellurium
##SBATCH -w tellurium1
##SBATCH --exclude=tellurium[1-8]


#!/usr/bin/env bash

echo "traj_00001"

. /user/mai/anaconda/miniconda3/etc/profile.d/conda.sh
conda activate sharc4.0
 module load intel mkl
. /user/mai/SHARC/REPOSITORY_switchable/sharc_main/bin/sharcvars.sh

#!/usr/bin/env bash

echo "traj_00007"



PRIMARY_DIR=/global/mai/SHARC_4.0_beta_tests/TUTORIAL/Tutorial/traj/Singlet_2/TRAJ_00007/

cd $PRIMARY_DIR

. $SHARC/sharcvars.sh
$SHARC/driver.py -i molcas input
