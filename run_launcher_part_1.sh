# See README.md for description of steps to set up a run

# This script adds stompy to the pythonpath environment variable and then runs 
# sfb_dfm.py, which creates most of the input files for an open bay DFM run

# User input: 
export RUN_NAME="wy2004"         # name of the run (this will be name of *.mdu file and folder it's stored in)
export RUN_START="2003-08-01"  # run start time in YYYY-MM-DD format
export RUN_STOP="2004-10-01"   # run end time in YYYY-MM-DD format
export MAKE_PLOTS="True"       # flag to make and save plots of boundary conditions
SFB_DFM_PARENT_PATH=/chicagovol1/hpcshared/open_bay/hydro/full_res/wy2004     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created
PPATH=/opt/anaconda3/envs/delft_env/bin/python # path to python executable (/opt/anaconda3/bin/python3 on richmond, /opt/anaconda3/envs/delft_env/bin/python on all other servers)

# get the parent directory of the sfb_dfm package 
#SFB_DFM_PARENT_PATH=$(dirname $(dirname $(readlink -f "$0"))) # <<< you can use this command to get the path automatically, assuming this shell script is in the sfb_dfm directory

# echo SFB_DFM_PARENT_PATH
echo ""
echo "sfb_dfm package is located here:"
echo "SFB_DFM_PARENT_PATH="$SFB_DFM_PARENT_PATH
echo ""

# assume stompy is inside the same parent directory as the sfb_dfm package
STOMPY_PATH=$SFB_DFM_PARENT_PATH/stompy
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
$PPATH $SFB_DFM_PARENT_PATH/sfb_dfm/sfb_dfm.py
echo ""
