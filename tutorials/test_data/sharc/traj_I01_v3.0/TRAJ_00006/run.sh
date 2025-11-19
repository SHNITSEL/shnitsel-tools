#$ -v USER_EPILOG=Singlet_2/TRAJ_00006//epilog.sh


PRIMARY_DIR=/user/anna/Bachelor/Methylenimmonium/Trajectories/Singlet_2/TRAJ_00006/
COPY_DIR=$TMPDIR/Singlet_2/TRAJ_00006/

mkdir -p $COPY_DIR
cp -r $PRIMARY_DIR/* $COPY_DIR
cd $COPY_DIR
echo $HOSTNAME > $PRIMARY_DIR/host_info
echo $(pwd) >> $PRIMARY_DIR/host_info
echo $(date) >> $PRIMARY_DIR/host_info

$SHARC/sharc.x input
err=$?

cp -r $COPY_DIR/output.* $COPY_DIR/restart.* $COPY_DIR/restart/ $PRIMARY_DIR

if [ $err == 0 ]; 
then
  rm -r $COPY_DIR
else
  echo "The calculation crashed at
date = $(date)
with error code $err. 
Please inspect the trajectory on
host = $HOSTNAME
in
dir  = $(pwd)
" > $PRIMARY_DIR/README
fi
