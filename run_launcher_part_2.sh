# See README.md for description of steps to set up a run

# This script  partitions and executes the DFM run


# USER INPUT: for debugging, we run the code in serial and output to the screen
DEBUG=false

# set some environment variables, shared across run launchers
source run_launcher_part_0.sh

# add to path environment variables
export PATH=$DFM_PATH/bin:$PATH 
export LD_LIBRARY_PATH=$DFM_PATH/lib:$LD_LIBRARY_PATH

# change to run directory
cd $RUN_DIR

if [ "$DEBUG" = true ]
then
	dflowfm --autostartstop $RUN_NAME.mdu > out.txt
else
	# Create partitions for parallel run
	echo "Partitioning into "$NPROC" subdomains, check "$RUN_DIR"/partition.txt for status"
	echo ""
	
	# this is the usual method
	#dflowfm --partition:ndomains=$NPROC:icgsolver=6 $RUN_NAME.mdu >partition.txt
	#dflowfm --partition:ndomains=$NPROC:icgsolver=6:genpolygon=1 $RUN_NAME.mdu >partition.txt
	dflowfm --partition:icgsolver=6 $RUN_NAME.mdu sfei_v23_straightened_part.pol >partition.txt

	# can't get it to work for straight grid, so supplied polygons
	echo ""
	
	# Execute parallel run
	echo "Executing DFM run, check "$RUN_DIR"/out.txt and "$RUN_DIR"/err.txt for status"
	echo ""
	mpiexec -n $NPROC dflowfm --autostartstop $RUN_NAME.mdu > out.txt 2> err.txt
	echo ""
fi
