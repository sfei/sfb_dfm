# See README.md for description of steps to set up a run

# This script partitions and executes the DFM run

# User input: 
RUN_NAME="wy2018"                             # name of the run (this will be name of *.mdu file and folder it's stored in)
NPROC=16                                    # number of processors (16 is a good number)
DFMV=/opt/software/delft/dfm/r52184-opt/bin # path to DFM binaries

# get the parent directory of the sfb_dfm package (assumes this script 
# is located inside the sfb_dfm package, but script can be run from anywhere)
SFB_DFM_PARENT_PATH=$(dirname $(dirname $(readlink -f "$0")))

# assumes run directory is in same parent directory as sfb_dfm package
# (if someone changes this in sfb_dfm.py, this will not be true anymore, so you'll
# have to point to the new run directory)
RUN_DIR=$SFB_DFM_PARENT_PATH/runs/$RUN_NAME

# Add DFM to PATH environment variable and check it only points to one version
export PATH=$DFMV:$PATH
echo "Make sure PATH points to only one version of DFM:"
echo "     PATH="$PATH
echo ""

# change to run directory
cd $RUN_DIR

# Create partitions for parallel run
echo "Partitioning into "$NPROC" subdomains, check "$RUN_DIR"/partition.txt for status"
echo ""
dflowfm --partition:ndomains=$NPROC:icgsolver=6 $RUN_NAME.mdu >partition.txt
echo ""

# Execute parallel run
echo "Executing DFM run, check "$RUN_DIR"/out.txt and "$RUN_DIR"/err.txt for status"
echo ""
mpiexec -n $NPROC dflowfm --autostartstop $RUN_NAME.mdu > out.txt 2> err.txt
echo ""
