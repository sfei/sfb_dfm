# See README.md for description of steps to set up a run

# This script adds stompy to the pythonpath environment variable and then runs 
# sfb_dfm.py, which creates most of the input files for an open bay DFM run
# make sure you have activated an appropriate conda environment before running
# this script -- on chicago, boise, or fortcollins enter "conda activate delft_env"
# or on richmond, you don't need to do anything because the base environment is fine

# User input: 
export RUN_NAME="wy2006_alliek"         # name of the run (this will be name of *.mdu file and folder it's stored in)
export RUN_START="2005-08-01"  # run start time in YYYY-MM-DD format
export RUN_STOP="2006-10-01"   # run end time in YYYY-MM-DD format
export MAKE_PLOTS="True"       # flag to make and save plots of boundary conditions
export SFB_DFM_PARENT_PATH=/boisevol1/hpcshared/open_bay/hydro/full_res/wy2006     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created

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
export PYTHONPATH=$STOMPY_PATH:$PYTHONPATH   
echo "PYTHONPATH="$PYTHONPATH
echo ""

# now run sfb_dfm.py to set up the run
echo "Running sfb_dfm.py"
echo ""
python sfb_dfm.py
echo ""
