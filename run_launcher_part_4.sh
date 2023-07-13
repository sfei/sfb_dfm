# See README.md for description of steps to set up a run

# This script only needs to be run for the old code (r52184)
# it runs a stompy based postprocessor to add sources back into the dwaq input files b/c they are missing

# set some environment variables, shared across run launchers
source run_launcher_part_0.sh

# add stompy path
export PYTHONPATH=$STOMPY_PATH       

# path to *.hyd file, should be based on run path and run name, if not enter whatever is correct
HYDRO_PATH=$RUN_DIR/DFM_DELWAQ_$RUN_NAME/$RUN_NAME.hyd

# change to run directory
cd $RUN_DIR

# now that output is coupled, make the mass conservation correction
echo "Calling sfb_dfm_postprocessor.py to make mass conservation correction at tributary/POTW inflow sites"
$PYTHON $SFB_DFM_PARENT_PATH/sfb_dfm/sfb_dfm_postprocessor.py

