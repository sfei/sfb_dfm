# See README.md for description of steps to set up a run

################################################
##### user may want to change the following ####
################################################

# specify parameters for this run 
export RUN_NAME="wy2019" # name of the run (this will be name of *.mdu file and folder it's stored in)
export RUN_START="2018-08-01"  # run start time in YYYY-MM-DD format
export RUN_STOP="2019-10-01"   # run end time in YYYY-MM-DD format
export MAKE_PLOTS="True"       # flag to make and save plots of boundary conditions
export SFB_DFM_PARENT_PATH=/boisevol1/hpcshared/open_bay/hydro/full_res/wy2019/     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created

# specify number of processors
export NPROC=16

# bound the salinity and temperautre in the DWAQ hydro input files
export MINTEMP=0
export MAXTEMP=28
export MINSALT=0
export MAXSALT=40

# we are using an old anaconda environment for now, need to iron out some bugs
# so sfb_dfm.py works in newer anaconda environments
# this is the direct path to python executable:
export PYTHON=/opt/anaconda3/envs/delft_env/bin/python 

# specify if you want to use the old or new version of the DFM code
# default should be true as we are still delaying a transfer to the newer version
# for complicated reasons
OLDCODE=true

# specify path to DFM
if [ "$OLDCODE" = true ]
then
	export DFM_PATH=/opt/software/delft/dfm/r52184-opt/
else
	export DFM_PATH=/opt/anaconda3/envs/dfm_t140737/
fi

# here is where you point to the program that stitches together the model output across
# the different subdomains (16 of them if you use the usual 16 processors). depending if you
# ran the old or new code, you will need to use ddcouplefm (for the old code) or waqmerge 
# (for the new code)
if [ "$OLDCODE" = true ]
then
	export DDCOUPLE_PATH=/opt/software/delft/ddcouplefm/1.02.01.50002/lnx64/ 
	export DDCOUPLE_NAME=ddcouplefm 
else
	export DDCOUPLE_PATH=/opt/software/delft/dfm/t141798/
	export DDCOUPLE_NAME=waqmerge
fi

###############################################################################################################
##### the following don't need to be altered, provided user organizes folders in a certain way (see README.md)
###############################################################################################################

# these paths require user to organize folders in a certain way, and to clone stompy into the sfb_dfm parent path folder
export RUN_DIR=$SFB_DFM_PARENT_PATH/runs/$RUN_NAME   # put run in the run folder
export HYDRO_PATH=$RUN_DIR/DFM_DELWAQ_$RUN_NAME/$RUN_NAME.hyd
export SFB_DFM_DIR=$SFB_DFM_PARENT_PATH/sfb_dfm # sfb_dfm should be found here
export STOMPY_PATH=$SFB_DFM_PARENT_PATH/stompy # stompy should be found here
export PYTHONPATH=$STOMPY_PATH



