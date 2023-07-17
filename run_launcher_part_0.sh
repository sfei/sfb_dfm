# See README.md for description of steps to set up a run

################################################
##### user may want to change the following ####
################################################

# specify parameters for this run 
export RUN_NAME="wy2022_t140737_new_stations" # name of the run (this will be name of *.mdu file and folder it's stored in)
export RUN_START="2021-08-01"  # run start time in YYYY-MM-DD format
export RUN_STOP="2022-10-01"   # run end time in YYYY-MM-DD format
export MAKE_PLOTS="True"       # flag to make and save plots of boundary conditions
export SFB_DFM_PARENT_PATH=/chicagovol1/hpcshared/open_bay/hydro/full_res/wy2022_t140737     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created

# specify number of processors
export NPROC=16

# specify path to DFM 
#export DFM_PATH=/opt/software/delft/dfm/r52184-opt # path to DFM binaries and libraries (bin and lib folders should be inside here) 
export DFM_PATH=/opt/anaconda3/envs/dfm_t140737/

# specify exact location of python and ddcouplefm/waqmerge executables
export PYTHON=/opt/anaconda3/envs/delft_env/bin/python # path to python executable

# for the old code, point to ddcouplefm, and for the new code use waqmerge that comes from a newer distribution
#export DDCOUPLE_PATH=/opt/software/delft/ddcouplefm/1.02.01.50002/lnx64/ # path to ddcouplefm or waqmerge bin and lib directories
#export DDCOUPLE_NAME=ddcouplefm # this is ddcouplefm for older code and waqmerge for newer code
export DDCOUPLE_PATH=/opt/software/delft/dfm/t141798/
export DDCOUPLE_NAME=waqmerge

###############################################################################################################
##### the following don't need to be altered, provided user organizes folders in a certain way (see README.md)
###############################################################################################################

# these paths require user to organize folders in a certain way, and to clone stompy into the sfb_dfm parent path folder
export RUN_DIR=$SFB_DFM_PARENT_PATH/runs/$RUN_NAME   # put run in the run folder
export SFB_DFM_DIR=$SFB_DFM_PARENT_PATH/sfb_dfm # sfb_dfm should be found here
export STOMPY_PATH=$SFB_DFM_PARENT_PATH/stompy # stompy should be found here



