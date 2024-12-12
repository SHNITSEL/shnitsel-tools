#!/bin/bash

#$-N init_00000



PRIMARY_DIR=/mnt/irisgpfs/users/camuller/04_cyanine-project/01_butene/namd_3Sing/geom_butene_s0/init/ICOND_00000//

cd $PRIMARY_DIR


$SHARC/SHARC_MOLCAS.py QM.in >> QM.log 2>> QM.err
