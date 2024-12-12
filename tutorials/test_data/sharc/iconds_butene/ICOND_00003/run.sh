#!/bin/bash

#$-N init_00003



PRIMARY_DIR=/mnt/irisgpfs/users/camuller/04_cyanine-project/01_butene/namd_3Sing/geom_butene_s0/init/ICOND_00003//

cd $PRIMARY_DIR

if [ -d ../ICOND_00000/SAVE ];
then
  if [ -d ./SAVE ];
  then
    rm -r ./SAVE
  fi
  cp -r ../ICOND_00000/SAVE ./
else
  echo "Should do a reference overlap calculation, but the reference data in ../ICOND_00000/ seems not OK."
  exit 1
fi


$SHARC/SHARC_MOLCAS.py QM.in >> QM.log 2>> QM.err
