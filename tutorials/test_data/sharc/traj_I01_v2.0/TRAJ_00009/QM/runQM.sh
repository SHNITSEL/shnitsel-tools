cd QM
/user/anna/anaconda2/bin/python $SHARC/SHARC_COLUMBUS.py QM.in >> QM.log 2>> QM.err
err=$?

exit $err