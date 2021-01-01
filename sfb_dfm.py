#!/usr/bin/env python
"""
Script to configure and generate mdu ready for running dflow fm.


Uses sfei_v20 grid, a Southbay-enhancement of the SF Bay/Delta community
model grid, with some areas deepened (LSB bathy), trimmed (Coyote Creek)
and dredge (see dredge_grid.py in this directory).

2018-01-05: fix some corrupt bathymetry in Napa- see update_bathy_to_v20.py
"""

import os
import glob
import pdb
import io
import shutil
import subprocess
import numpy as np
import logging
import xarray as xr
import six
from pathlib import Path

import sys
import stompy.model.delft.io as dio
from stompy.model.delft import dfm_grid
from stompy.spatial import wkb2shp
from stompy.io.local import usgs_nwis,noaa_coops
from stompy import utils

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

log=logging.getLogger('sfb_dfm')

# Change working directory to script location
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

print('dname = ' + dname)

import sfb_dfm_utils 
#%% 
DAY=np.timedelta64(86400,'s') # useful for adjusting times

# get name of run, start time, stop time, and flag to plot boundary conditions (or not) 
# from environment variables set in run launcher script
run_name = os.getenv('RUN_NAME')
run_start = np.datetime64(os.getenv('RUN_START'))
run_stop = np.datetime64(os.getenv('RUN_STOP'))
make_plots = bool(os.getenv('MAKE_PLOTS'))
#run_name="wy2016" 
#run_start=np.datetime64('2015-08-01')
#run_stop=np.datetime64('2015-08-02')
#make_plots = False # True = make plots :) 

ALL_FLOWS_UNIT = False # for debug, set all volumetric flow rates to 1m3/s if True

run_folder = Path(os.path.dirname(dname)) # set run folder to parent directory of sfb_dfm package

## --------------------------------------------------

# Derived parameters used in multiple places

# base_dir=os.path.dirname(__file__) # for command line invocation
base_dir = Path(dname)              # right now the base dir is where this script sfb_dfm lives. 
runs_dir =  run_folder / 'runs'     # Go to the run folder, check out if a seperate folder for runs exists... (parent directory)
runs_dir.exists() or runs_dir.mkdir()
run_base_dir =  runs_dir / run_name  
run_base_dir.exists() or run_base_dir.mkdir()  # Make sure run directory exists (if not, make it..)
print('')
print('YOUR RUN WILL BE FOUND HERE = %s' % run_base_dir)
print('')
#%% 

abs_static_dir = base_dir / 'inputs-static' # real location of static directory
rel_static_dir = os.path.relpath(str(abs_static_dir), str(run_base_dir)) # static directory relative to the run directory
rel_static_dir = Path(rel_static_dir)       # os can read Pathlib objects, so keeping formatting consistent 

# reference date - can only be specified to day precision, so # truncate to day precision (rounds down)
ref_date = run_start.astype('datetime64[D]')
net_file = base_dir / 'sfei_v20_net.nc' 

# No longer using any new-style boundary conditions
old_bc_fn = run_base_dir / 'FlowFMold_bnd.ext' 
obs_shp_fn = abs_static_dir / 'observation-points.shp'

dredge_depth = -0.5 # m NAVD88, depth to enforce at inflows and discharges



# clear any stale bc files:
for fn in [old_bc_fn]:
    fn.exists() and fn.unlink() 

#%% 
## --------------------------------------------------------------------------------
# Edits to the template mdu:
mdu = dio.MDUFile('template.mdu')

# set run base path
mdu.base_path = str(run_base_dir)   # converting Path object to string for stompy.

if 1: # set dates
    # RefDate can only be specified to day precision
    mdu['time','RefDate'] = utils.to_datetime(ref_date).strftime('%Y%m%d')
    mdu['time','Tunit']   = 'M' # minutes.  kind of weird, but stick with what was used already
    mdu['time','TStart']  = 0
    mdu['time','TStop']   = int( (run_stop - run_start) / np.timedelta64(1,'m') )

mdu['geometry','LandBoundaryFile'] = rel_static_dir / "deltabay.ldb"

mdu['geometry','Kmx'] = 10 # 10 layers

# update location of the boundary conditions
# this has the source/sinks which cannot be written in the new style file
mdu['external forcing','ExtForceFile'] = old_bc_fn.name

#%%
# Load the grid now -- it's used for clarifying some inputs, but
# is also modified to deepen areas near inflows, before being written
# out near the end of the script
grid = dfm_grid.DFMGrid(str(net_file)) 
    
## split into relative and absolute directories (alliek Dec 2020)
rel_bc_dir = 'bc_files'
abs_bc_dir = run_base_dir / rel_bc_dir
abs_bc_dir.exists() or abs_bc_dir.mkdir()

# features which have manually set locations for this grid
adjusted_pli_fn = base_dir / 'nudged_features.pli' 

# this line worked in emma's repo that she left on hpc because, but since she cloned rusty's 
# sfb_dfm_repo he made updates to add_sfbay_freshwater, so changing to work with rusty's updated 
# sfb_dfm_utils (alliek dec 2020)
#sfb_dfm_utils.add_sfbay_freshwater(bc_dir,
#                                   run_start,run_stop,ref_date,
#                                   adjusted_pli_fn,
#                                   freshwater_dir=os.path.join(base_dir, 'sfbay_freshwater'),
#                                   grid=grid,
#                                   dredge_depth=dredge_depth,
#                                   old_bc_fn=old_bc_fn,
#                                   all_flows_unit=ALL_FLOWS_UNIT)          
sfb_dfm_utils.add_sfbay_freshwater(mdu,
                         str(rel_bc_dir), # added rel_bc_dir alliek dec 2020
                         str(adjusted_pli_fn),
                         freshwater_dir = str(base_dir / 'sfbay_freshwater'),
                         grid = grid,
                         dredge_depth = dredge_depth,
                         all_flows_unit = ALL_FLOWS_UNIT,
                         time_offset = None)

#%% 


# POTW inputs:
# The new-style boundary inputs file (FlowFM_bnd_new.ext) cannot represent
# sources and sinks, so these come in via the old-style file.

## split into relative and absolute directories (alliek Dec 2020)
rel_src_dir = 'source_files'
src_dir     = run_base_dir / rel_src_dir
src_dir.exists() or src_dir.mkdir() 

# path to potw repository
potw_dir  =  base_dir / 'sfbay_potw'

# this line worked in emma's hpc repo but rusty made updates since then, so changing to 
# work with rusty's updated sfb_dfm_utils
#sfb_dfm_utils.add_sfbay_potw(src_dir,
#                             run_start,run_stop,ref_date,
#                             potw_dir,
#                             adjusted_pli_fn,
#                             grid,dredge_depth,
#                             old_bc_fn,
#                             all_flows_unit=ALL_FLOWS_UNIT)
sfb_dfm_utils.add_sfbay_potw(mdu, 
                             str(rel_src_dir), # added rel_src_dir alliek dec 2020
                             str(potw_dir), 
                             str(adjusted_pli_fn), 
                             grid,
                             dredge_depth, 
                             all_flows_unit = ALL_FLOWS_UNIT, 
                             time_offset = None, 
                             write_salt = True, 
                             write_temp = True)

##
# not sure but I beleive the following logical indices help to deal with missing temperature
# data at jersey point after 1/11/2016. this is from changes emma made to delta_inflows.py in
# sfb_dfm_utils
temp_jersey = run_start>np.datetime64('2009-12-01')<run_stop<np.datetime64('2016-11-01')
temp_rio    = run_start>np.datetime64('2010-01-01')<run_stop<np.datetime64('2020-01-20')
# Delta boundary conditions
# saved over rusty's delta_inflow.py with emma's version but then changed to be more like rusty's 
# for better handling of boundary condition directory
sfb_dfm_utils.add_delta_inflow(mdu,
                               str(rel_bc_dir),
                               static_dir = str(abs_static_dir),
                               grid = grid,
                               dredge_depth = dredge_depth,
                               all_flows_unit = ALL_FLOWS_UNIT,
                               temp_jersey = temp_jersey,
                               temp_rio = temp_rio)
##


# This factor seems to be about right for Point Reyes tides
# to show up at SF with the right amplitude.  Without
# attenuation, in runs/w2013b tides at SF are 1.10x observed.
# The lag is a bit less clear, with SF tides at -2.5 minutes (leading),
# but SF currents at about -15 minutes (leading).  All of these
# are likely wrapped up in some friction calibration, for another
# day.
sfb_dfm_utils.add_ocean(str(run_base_dir),
                        str(rel_bc_dir),
                        run_start,
                        run_stop,
                        ref_date,
                        static_dir = str(abs_static_dir),
                        grid = grid,
                        factor = 0.901,
                        lag_seconds = 120,
                        old_bc_fn = str(old_bc_fn),
                        all_flows_unit = ALL_FLOWS_UNIT)

## 
if 1:            
    lines=["QUANTITY=frictioncoefficient",
           "FILENAME=%s/friction12e.xyz" % rel_static_dir,
           "FILETYPE=7",
           "METHOD=5",
           "OPERAND=O",
           ""]
    with open(str(old_bc_fn) ,'at') as fp:
        fp.write("\n".join(lines))

if 1:  # Copy grid file into run directory and update mdu
    mdu['geometry','NetFile'] = net_file.name
    dest = run_base_dir / mdu['geometry','NetFile']
    # write out the modified grid
    dfm_grid.write_dfm(grid, str(dest) , overwrite=True)


# update to work with rusty's changes to sfb_dfm_utils
#sfb_dfm_utils.add_initial_salinity_dyn(run_base_dir,
#                                       abs_static_dir,
#                                       mdu,
#                                       run_start)
sfb_dfm_utils.add_initial_salinity(str(run_base_dir),
                                   str(abs_static_dir),
                                   str(old_bc_fn),
                                   all_flows_unit = ALL_FLOWS_UNIT)


# WIND
#ludwig_ok=sfb_dfm_utils.add_erddap_ludwig_wind(run_base_dir,
#                                               run_start,run_stop,
#                                               old_bc_fn)
#if not ludwig_ok:
#    const_ok=sfb_dfm_utils.add_constant_wind(run_base_dir,mdu,[0,0],run_start,run_stop)
#    assert const_ok
#else:
#    assert ludwig_ok # or see lsb_dfm.py for constant field.

##

if 1: # fixed weir file is just referenced as static input
    mdu['geometry','FixedWeirFile'] = str(rel_static_dir / 'SBlevees_tdk.pli') 

if 1: 
    # evaporation was a bit out of control in south bay - try scaling back just
    # the evaporation some.  This is a punt!
    
    # update to work with rusty's updates to sfb_dfm_utils (added scale_precip)
    #sfb_dfm_utils.add_cimis_evap_precip(run_base_dir,mdu,scale_evap=0.5)
    # print('DEBUGGING THIS!!!')
    sfb_dfm_utils.add_cimis_evap_precip(str(run_base_dir), mdu, scale_precip=1.0, scale_evap=0.5)
    
if 1: # output locations
    mdu['output','CrsFile'] = str (rel_static_dir / "SB-observationcrosssection.pli")

##
if 1:
    # Observation points taken from shapefile for easier editing/comparisons in GIS
    obs_pnts = wkb2shp.shp2geom(str(obs_shp_fn))
    obs_fn   = 'observation_pnts.xyn'
    
    with open( str(run_base_dir/obs_fn),'wt') as fp:
        for idx,row in enumerate(obs_pnts):
            xy = np.array(row['geom'])
            fp.write("%12g %12g '%s'\n"%(xy[0], xy[1], row['name']))
    mdu['output','ObsFile'] = obs_fn

    if run_name.startswith('short'):
        mdu['output','MapInterval'] = 3600
    
##
mdu_fn = str(run_base_dir / (run_name + ".mdu")) 
mdu.write(mdu_fn)
print('Just printed out %s.' % mdu_fn)


if make_plots:
    from sfb_dfm_utils import plot_mdu # SW added function here 
    print('Now making plots...')
    plot_mdu.plot_MDU(mdu_fn, net_file)
##

## As of r52184, explicitly built with metis support, partitioning can be done automatically
## from here.
#
#cmd="%s/mpiexec -n %d %s/dflowfm --partition:ndomains=%d %s"%(dfm_bin_dir,nprocs,dfm_bin_dir,nprocs,
#                                                              mdu['geometry','NetFile'])
#pwd=os.getcwd()
#try:
#    os.chdir(run_base_dir)
#    res=subprocess.call(cmd,shell=True)
#finally:
#    os.chdir(pwd)
#
#
## similar, but for the mdu:
#cmd="%s/generate_parallel_mdu.sh %s %d 6"%(dfm_bin_dir,os.path.basename(mdu_fn),nprocs)
#try:
#    os.chdir(run_base_dir)
#    res=subprocess.call(cmd,shell=True)
#finally:
#    os.chdir(pwd)
