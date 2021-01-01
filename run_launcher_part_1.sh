# See README.md for description of steps to set up a run

# This script adds stompy to the pythonpath environment variable and then runs 
# sfb_dfm.py, which creates most of the input files for an open bay DFM run

# User input: 
export RUN_NAME="wy2018"         # name of the run (this will be name of *.mdu file and folder it's stored in)
export RUN_START="2017-08-01"  # run start time in YYYY-MM-DD format
export RUN_STOP="2018-10-01"   # run end time in YYYY-MM-DD format
export MAKE_PLOTS="True"       # flag to make and save plots of boundary conditions

# get the parent directory of the sfb_dfm package (assumes this script 
# is located inside the sfb_dfm package, but script can be run from anywhere)
SFB_DFM_PARENT_PATH=$(dirname $(dirname $(readlink -f "$0")))

# assume stompy is inside the same parent directory as the sfb_dfm package
STOMPY_PATH=$SFB_DFM_PARENT_PATH/stompy
echo ""
echo "User must make sure stompy package is installed here:"
echo $STOMPY_PATH
echo ""

# add stompy to pythonpath
echo "Adding stompy to PYTHONPATH":
#export PYTHONPATH=$STOMPY_PATH:$PYTHONPATH 
export PYTHONPATH=$STOMPY_PATH    # this overrides paths to pre-existing stompy (safest option)
echo "PYTHONPATH="$PYTHONPATH
echo ""

# now run sfb_dfm.py to set up the run
echo "Running sfb_dfm.py"
echo ""
python3 $SFB_DFM_PARENT_PATH/sfb_dfm/sfb_dfm.py
echo ""
echo "If desired, add wind and meteo forcing before partitioning and executing"
echo "the run in run_launcher_part_2.sh"
echo ""
