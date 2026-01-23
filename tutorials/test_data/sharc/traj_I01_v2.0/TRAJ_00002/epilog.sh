#/bin/bash

PRIMARY_DIR=/global/anna/Methylenimmonium/Trajectories/Singlet_2/TRAJ_00002/
COPY_DIR=$TMPDIR/Singlet_2/TRAJ_00002/

cp $COPY_DIR/output.* $COPY_DIR/restart.* $PRIMARY_DIR
rm -r $COPY_DIR
