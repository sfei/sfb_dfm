# See README.md for description of steps to set up a run

# This script  partitions and executes the DFM run

# User input: 
export RUN_NAME="wy2022_r52184" # name of the run (this will be name of *.mdu file and folder it's stored in)
SFB_DFM_PARENT_PATH=/boisevol1/hpcshared/open_bay/hydro/full_res/wy2022_r52184     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created
NPROC=16                                    # number of processors (16 is a good number)

## take this out for the t140737 run, put it back for old school run
DFMV=/opt/software/delft/dfm/r52184-opt/bin # path to DFM binaries

# assumes run directory is in same parent directory as sfb_dfm package
# (if someone changes this in sfb_dfm.py, this will not be true anymore, so you'll
# have to point to the new run directory)
RUN_DIR=$SFB_DFM_PARENT_PATH/runs/$RUN_NAME

# take this out for the t140737 run, put it back for old school run
# Add DFM to PATH environment variable and check it only points to one version
export PATH=$DFMV:$PATH
echo "Make sure PATH points to only one version of DFM:"
echo "     PATH="$PATH
echo ""

# change to run directory
cd $RUN_DIR

# add this for the t140737 run, take it out for the old school run
#export LD_LIBRARY_PATH=/opt/anaconda3/envs/dfm_t140737/lib

# Execute parallel run
echo "Executing DFM run, check "$RUN_DIR"/out.txt and "$RUN_DIR"/err.txt for status"
echo ""
mpiexec -n $NPROC dflowfm --autostartstop $RUN_NAME.mdu > out.txt 2> err.txt
echo ""

