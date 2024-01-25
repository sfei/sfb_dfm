# See README.md for description of steps to set up a run

# This script only needs to be run for the old code (r52184)
# it runs a stompy based postprocessor to add sources back into the dwaq input files b/c they are missing

# set some environment variables, shared across run launchers
source run_launcher_part_0.sh

# change to run directory
cd $RUN_DIR

# now that output is coupled, make the mass conservation correction
echo "Calling postprocess_correct_point_source_error.py to make mass conservation correction at tributary/POTW inflow sites"
$PYTHON $SFB_DFM_PARENT_PATH/sfb_dfm/postprocess_correct_point_source_error.py

# now call the script that bounds the salinity and temperature
echo "Calling postprocess_cap_salinity_and_temperature.py to place bounds on salinity and temperature in dwaq hydro files"
$PYTHON $SFB_DFM_PARENT_PATH/sfb_dfm/postprocess_cap_salinity_and_temperature.py


