# alliek october 2022
# script to create tidally filtered aggregated grid hydro inputs from full resolution hydro inputs for wy2022_bloom run

# I ran this script in the the delft_env anaconda environment on chicago,
# enter this at command line to activate: "conda activate delft_env"

RUN_NAME="wy2021"                             # name of the run (this will be name of *.mdu file and folder it's stored in)
SFB_DFM_PARENT_PATH=/boisevol1/hpcshared/open_bay/hydro/full_res/wy2021     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created


# user sets these variables
runname='wy2022_r52184'
runpath_fr='/boisevol1/hpcshared/open_bay/hydro/full_res/wy2021/runs/'
runpath_ag='/boisevol1/hpcshared/open_bay/hydro/agg673/'
shppath='/richmondvol1/hpcshared/Shapefiles/open_bay_agg_grid_hand_drawn/boxes-v2_quadruple_resolution_straight_SHRINKWRAP_SNAPPED_NODUPS.shp'

# use latest version of stompy ... aggregation script is a single *.py file in the stompy/model/delft directory
export PYTHONPATH="/opt/software/rusty/stompy/newest_commit/stompy"

input1=$runpath_fr$runname'/DFM_DELWAQ_'$runname'/'$runname'.hyd'
output1=$runpath_ag$runname'/'$runname'_agg/'$runname'_agg'
input2=$runpath_ag$runname'/'$runname'_agg/com-'$runname'_agg.hyd'
output2=$runpath_ag$runname'/'$runname'_agg_lp/'$runname'_agg_lp'

echo "Aggregate the hydro output..."
#python -m stompy.model.delft.waq_hydro_editor -i /path_to_hyd_file/file.hyd -a /path_to_grid_shapefile/grid.shp -o /path_to_output_folder/output_agg/output
python -m stompy.model.delft.waq_hydro_editor -i $input1 -a $shppath -o $output1

echo "Tidally filter aggregated output..."
#python -m stompy.model.delft.waq_hydro_editor -i output_agg/com-output.hyd -l -p -o output_lp/output 
#*Note: -p option can be included or deleted -- 
# inclusion: passes tau and salinity to the output without tidal filtering it 
# exclusion: tidally filters all the parameters
python -m stompy.model.delft.waq_hydro_editor -i $input2 -l -p -o $output2
