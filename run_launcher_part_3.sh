# See README.md for description of steps to set up a run

# This script runs ddcouplefm or waqmerge to merge multidomain in delwaq input files

# set some environment variables, shared across run launchers
source run_launcher_part_0.sh
     
# add ddcouple libraries to library path
export LD_LIBRARY_PATH="$DDCOUPLE_PATH/lib:$LD_LIBRARY_PATH"

# change to run directory
cd $RUN_DIR

# execute ddcouple to splice output together
echo "Running ddcouplefm to splice together DWAQ output across "$NPROC" domains"
$DDCOUPLE_PATH/bin/$DDCOUPLE_NAME $RUN_NAME".mdu" $NPROC

# now that output is coupled, make the mass conservation correction
if [ "$OLDCODE" = true ]
then
	echo "Calling sfb_dfm_postprocessor.py to make mass conservation correction at tributary/POTW inflow sites"
	$PYTHON $SFB_DFM_PARENT_PATH/sfb_dfm/sfb_dfm_postprocessor.py
fi
