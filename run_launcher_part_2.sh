# See README.md for description of steps to set up a run

# This script 
# 	(1) partitions and executes the DFM run
# Then executes two postprocessing steps:
# 	(2) stitching together the 16 domains into one for the dwaq hydro input files, 
# 	(3) add back the missing source flows (bug in this version of the DFM model is that they are left out, violating mass conservation)

# User input: 
RUN_NAME="wy2004"                             # name of the run (this will be name of *.mdu file and folder it's stored in)
SFB_DFM_PARENT_PATH=/chicagovol1/hpcshared/open_bay/hydro/full_res/wy2004     # this is the directory where the sfb_dfm and stompy are located, and it is where the "runs" folder will be created
NPROC=16                                    # number of processors (16 is a good number)
DFMV=/opt/software/delft/dfm/r52184-opt/bin # path to DFM binaries
DDCOUPLEV=/opt/software/delft/ddcouplefm/1.02.01.50002/lnx64 # path to ddcouple, the executable from deltares that stitches dwaq output together
PPATH=/opt/anaconda3/envs/delft_env/bin/python # path to python executable (this should be /opt/anaconda3/bin/python3 on richmond, /opt/anaconda3/envs/delft_env/bin/python on all other servers)

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

# At this point the run may or may not have completed, but if it even partially completed, the 
# following code will perform some postprocessing tricks
cd "$SFB_DFM_PARENT_PATH/sfb_dfm"

# set path to ddcouple and its libraries
export PATH="$DDCOUPLEV/bin:$PATH"
export LD_LIBRARY_PATH="$DDCOUPLEV/lib:$LD_LIBRARY_PATH"

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

# path to *.hyd file, should be based on run path and run name, if not enter whatever is correct
export HYDRO_PATH=$RUN_DIR/DFM_DELWAQ_$RUN_NAME/$RUN_NAME.hyd

# change to run direcgtory and execute ddcouple to splice output together
echo "Running ddcouplefm to splice together DWAQ output across "$NPROC" domains"
cd $RUN_DIR
ddcouplefm $RUN_NAME $NPROC

# now that output is coupled, make the mass conservation correction
echo "Calling sfb_dfm_postprocessor.py to make mass conservation correction at tributary/POTW inflow sites"
$PPATH $SFB_DFM_PARENT_PATH/sfb_dfm/sfb_dfm_postprocessor.py

