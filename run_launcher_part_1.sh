# See README.md for description of steps to set up a run

# This script runs sfb_dfm.py, which creates most of the input files for an open bay DFM run

# set some environment variables, shared across run launchers
source run_launcher_part_0.sh

# add to path environment variables
export PYTHONPATH=$STOMPY_PATH       

# now run sfb_dfm.py to set up the run 
# (PYTHON variable gives path to the python executable and is set in run_launcher_part_0.sh)
echo "Running sfb_dfm.py using this version of python: "$PYTHON
echo ""
$PYTHON sfb_dfm.py
echo ""
