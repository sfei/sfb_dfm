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
from stompy.grid import unstructured_grid # added by alliek august 2024
from stompy.spatial import wkb2shp
from stompy.io.local import usgs_nwis,noaa_coops
from stompy import utils
import sfb_dfm_utils 

if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)
log=logging.getLogger('sfb_dfm')

#%% 
DAY=np.timedelta64(86400,'s') # useful for adjusting times

# get name of run, start time, stop time, and flag to plot boundary conditions (or not) 
# from environment variables set in run launcher script
run_folder = os.getenv('SFB_DFM_PARENT_PATH')
run_name = os.getenv('RUN_NAME')
run_start = np.datetime64(os.getenv('RUN_START'))
run_stop = np.datetime64(os.getenv('RUN_STOP'))
make_plots_str = os.getenv('MAKE_PLOTS')

# it's kind of involved to interpret True/False strings
if make_plots_str=='True':
    make_plots=True
elif make_plots_str=='False':
    make_plots=False
else:
    raise Exception('MAKE_PLOTS environment variable set in run_launcher_part_1.sh must be "True" or "False"')
    
#run_name="wy2016" 
#run_start=np.datetime64('2015-08-01')
#run_stop=np.datetime64('2015-08-02')
#make_plots = False # True = make plots :) 

ALL_FLOWS_UNIT = False # for debug, set all volumetric flow rates to 1m3/s if True


## --------------------------------------------------

# Derived parameters used in multiple places

# base_dir=os.path.dirname(__file__) # for command line invocation
base_dir = os.path.join(run_folder,'sfb_dfm')                   # right now the base dir is where this script sfb_dfm lives. 
runs_dir =  os.path.join(run_folder,'runs')     # Go to the run folder, check out if a seperate folder for runs exists... (parent directory)
os.path.exists(runs_dir) or os.makedirs(runs_dir)
run_base_dir =  os.path.join(runs_dir, run_name)  
os.path.exists(run_base_dir) or os.makedirs(run_base_dir)  # Make sure run directory exists (if not, make it..)
print('')
print('YOUR RUN WILL BE FOUND HERE = %s' % run_base_dir)
print('')
#%% 

abs_static_dir = os.path.join(base_dir, 'inputs-static') # real location of static directory
rel_static_dir = os.path.relpath(abs_static_dir, run_base_dir) # static directory relative to the run directory

# this is a somewhat awkward requirement due to sloppy coding by allie, helps the intial condition and ocean boudnary
# temperature routines find everything they need
abs_init_dir = os.path.join(base_dir,'sfb_dfm_utils')

# reference date - can only be specified to day precision, so # truncate to day precision (rounds down)
ref_date = run_start.astype('datetime64[D]')
net_file = os.path.join(base_dir, 'sfei_v23_straightened_net.nc')

# No longer using any new-style boundary conditions
old_bc_fn = os.path.join(run_base_dir ,'FlowFMold_bnd.ext')
obs_shp_fn = os.path.join(abs_static_dir ,'observation-points.shp')

# path to grid boundary shapefile
grid_boundary_fn = os.path.join(base_dir, 'derived', 'sfei_v23_straightened_net_OUTLINE.shp')

# path to cimis input file
cimis_fn = os.path.join(base_dir,'sfbay_cimis','union_city-hourly.nc')

dredge_depth = -0.5 # m NAVD88, depth to enforce at inflows and discharges



# clear any stale bc files:
for fn in [old_bc_fn]:
    os.path.exists(fn) and os.unlink(fn) 

#%% 
## --------------------------------------------------------------------------------
# Edits to the template mdu:
mdu = dio.MDUFile('template.mdu')

# set run base path
mdu.base_path = run_base_dir   

if 1: # set dates
    # RefDate can only be specified to day precision
    mdu['time','RefDate'] = utils.to_datetime(ref_date).strftime('%Y%m%d')
    mdu['time','Tunit']   = 'M' # minutes.  kind of weird, but stick with what was used already
    mdu['time','TStart']  = 0
    mdu['time','TStop']   = int( (run_stop - run_start) / np.timedelta64(1,'m') )

mdu['geometry','LandBoundaryFile'] = os.path.join(rel_static_dir, "deltabay.ldb")

mdu['geometry','Kmx'] = 10 # 10 layers

# update location of the boundary conditions
# this has the source/sinks which cannot be written in the new style file
mdu['external forcing','ExtForceFile'] = old_bc_fn

#%%
# Load the grid now -- it's used for clarifying some inputs, but
# is also modified to deepen areas near inflows, before being written
# out near the end of the script
#grid = dfm_grid.DFMGrid(str(net_file)) 
grid = unstructured_grid.UnstructuredGrid.read_dfm(str(net_file)) # alliek update to ugrid august 2024
    
## split into relative and absolute directories (alliek Dec 2020)
rel_bc_dir = 'bc_files'
abs_bc_dir = os.path.join(run_base_dir,rel_bc_dir)
os.path.exists(abs_bc_dir) or os.makedirs(abs_bc_dir)

# features which have manually set locations for this grid
adjusted_pli_fn = os.path.join(base_dir,'nudged_features.pli') 

# copy polygon file to guide parallelization to run directory
shutil.copyfile(os.path.join(abs_static_dir,'sfei_v23_straightened_part.pol'), 
                os.path.join(run_base_dir,'sfei_v23_straightened_part.pol'))


# this line worked in emma's repo that she left on hpc because, but since she cloned rusty's 
# sfb_dfm_repo he made updates to add_sfbay_freshwater, so changing to work with rusty's updated 
# sfb_dfm_utils (alliek dec 2020)     
sfb_dfm_utils.add_sfbay_freshwater(mdu,
                         str(rel_bc_dir), # added rel_bc_dir alliek dec 2020
                         str(adjusted_pli_fn),
                         freshwater_dir = os.path.join(base_dir,'sfbay_freshwater'),
                         grid = grid,
                         dredge_depth = dredge_depth,
                         temperature_fn=os.path.join(base_dir,'alameda_creek_climatology','Alameda_Creek_Climatological_Temperature.csv'),
                         all_flows_unit = ALL_FLOWS_UNIT,
                         time_offset = None)

#%% 


# POTW inputs:
# The new-style boundary inputs file (FlowFM_bnd_new.ext) cannot represent
# sources and sinks, so these come in via the old-style file.

## split into relative and absolute directories (alliek Dec 2020)
rel_src_dir = 'source_files'
src_dir     = os.path.join(run_base_dir, rel_src_dir)
os.path.exists(src_dir) or os.makedirs(src_dir) 

# path to potw repository
potw_dir  =  os.path.join(base_dir,'sfbay_potw')

# this line worked in emma's hpc repo but rusty made updates since then, so changing to 
# work with rusty's updated sfb_dfm_utils
sfb_dfm_utils.add_sfbay_potw(mdu, 
                             rel_src_dir, # added rel_src_dir alliek dec 2020
                             potw_dir, 
                             adjusted_pli_fn, 
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
#temp_jersey = run_start>np.datetime64('2009-12-01')<run_stop<np.datetime64('2016-11-01')
#temp_rio    = run_start>np.datetime64('2010-01-01')<run_stop<np.datetime64('2020-01-20')

# Delta boundary conditions
# saved over rusty's delta_inflow.py with emma's version but then changed to be more like rusty's 
# for better handling of boundary condition directory
sfb_dfm_utils.add_delta_inflow(mdu, 
                               base_dir = base_dir,
                               rel_bc_dir = rel_bc_dir,
                               static_dir = abs_static_dir,
                               grid = grid,
                               dredge_depth = dredge_depth,
                               all_flows_unit = ALL_FLOWS_UNIT,
                               temp_jersey = False,
                               temp_rio = False)
##


## prior to adding the ocean boundary, need to generate the temperature time series
#sfb_dfm_utils.make_ROMS_temp_tim(abs_init_dir,abs_bc_dir,ref_date,run_start,run_stop)

## also copy the ROMS pli file into the bc directory
#shutil.copyfile(os.path.join(abs_static_dir,'sea_temp_ROMS.pli'), os.path.join(abs_bc_dir,'sea_temp_ROMS.pli'))

# This factor seems to be about right for Point Reyes tides
# to show up at SF with the right amplitude.  Without
# attenuation, in runs/w2013b tides at SF are 1.10x observed.
# The lag is a bit less clear, with SF tides at -2.5 minutes (leading),
# but SF currents at about -15 minutes (leading).  All of these
# are likely wrapped up in some friction calibration, for another
# day.
sfb_dfm_utils.add_ocean(run_base_dir,
                        rel_bc_dir,
                        run_start,
                        run_stop,
                        ref_date,
                        base_dir,
                        static_dir = abs_static_dir,
                        grid = grid,
                        factor = 0.901,
                        lag_seconds = 120,
                        old_bc_fn = old_bc_fn,
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
    mdu['geometry','NetFile'] = net_file
    
    # alliek commented out august 2024:
    #dest = os.path.join(run_base_dir, mdu['geometry','NetFile'])
    ## write out the modified grid
    #dfm_grid.write_dfm(grid, str(dest) , overwrite=True)

# prior to adding initial salinity to the external forcing, we have to generate it!!!
# in july 2023 added secchidepth to this function ... note that secchidepth is actually 
# not an initial condition, but dfm doesn't allow it to vary in time, so it has the same 
# xyz format, and it's based on the same usgs cruise dataset as the initial salinity and temperature,
# so i added it here
sfb_dfm_utils.make_initial_sal_temp(run_start, run_stop, abs_init_dir, abs_bc_dir, grid)


# update to work with rusty's changes to sfb_dfm_utils
#sfb_dfm_utils.add_initial_salinity_dyn(run_base_dir,
#                                       abs_static_dir,
#                                       mdu,
#                                       run_start)
# altered by Allie to use updated initial salinity and also temperature!!! and secchidepth too now
sfb_dfm_utils.add_initial_salinity(run_base_dir,
                                   abs_static_dir, 
                                   abs_bc_dir, # added by allie
                                   old_bc_fn,
                                   all_flows_unit = ALL_FLOWS_UNIT)


if 1: # fixed weir file is just referenced as static input
    mdu['geometry','FixedWeirFile'] = os.path.join(rel_static_dir,'SBlevees_tdk.pli') 

if 1: 
    # evaporation was a bit out of control in south bay - try scaling back just
    # the evaporation some.  This is a punt!
    
    # update to work with rusty's updates to sfb_dfm_utils (added scale_precip)
    #sfb_dfm_utils.add_cimis_evap_precip(run_base_dir,mdu,scale_evap=0.5)
    # print('DEBUGGING THIS!!!')
    sfb_dfm_utils.add_cimis_evap_precip(cimis_fn, run_base_dir, mdu, scale_precip=1.0, scale_evap=0.5)
    
if 1: # output locations
    mdu['output','CrsFile'] = os.path.join(rel_static_dir, "SB-observationcrosssection.pli")

##
if 1:
    # Observation points taken from shapefile for easier editing/comparisons in GIS
    obs_pnts = wkb2shp.shp2geom(str(obs_shp_fn))
    obs_fn   = 'observation_pnts.xyn'
    
    with open( os.path.join(run_base_dir,obs_fn),'wt') as fp:
        for idx,row in enumerate(obs_pnts):
            xy = np.array([row['geom'].x,row['geom'].y])
            fp.write("%12g %12g '%s'\n"%(xy[0], xy[1], row['name']))
    mdu['output','ObsFile'] = obs_fn

    if run_name.startswith('short'):
        mdu['output','MapInterval'] = 3600


# initialze a string to print at the end of this script
user_instruction_string = ''

if 1: # if using temperature model
    mdu['physics','Temperature'] = 5 # this selects the "composite" temperature model
    mdu['physics','Soiltempthick'] = 0.1 # controls heat exchange with the bed, can be used to tune temperature model if desired 

    # copy the meteo grid file into the run folder
    shutil.copyfile(os.path.join(abs_static_dir,'meteo_coarse.grd'), os.path.join(run_base_dir,'meteo_coarse.grd'))

    # add lines to the external forcing file to reference meteo input that the user needs to upload!!!
    lines=["QUANTITY=humidity_airtemperature_cloudiness",
           "FILENAME=%s/hac.tem" % rel_bc_dir,
           "FILETYPE=6",
           "METHOD=3",
           "OPERAND=O",
           ""]
    with open(old_bc_fn ,'at') as fp:
        fp.write("\n".join(lines))

    # tell user to upload meteo file!!!
    user_instruction_string += ("\nATTENTION!\n" + 
    'User must manually upload the meterological forcing file, named hac.tem, and put it here:\n%s\n' % str(abs_bc_dir) + 
    'You can find scripts for generating hac.tem in the following location on SFEI''s Google Drive:\n' + 
    '\\1_Nutrient_Share\\2_Data_NUTRIENTS\\SFEI_Meteo\\Meteo4DFlow-SFB-UTC\\\n' + 
    '(https://drive.google.com/drive/folders/16ILqtvXCRjAupDagCSRMLPFE-0oWMBum)\n' +
    'Make sure to change the permissions of hac.tem, using chmod, once you have uploaded it, so DFM can read it!\n')


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
# updated wind approach
if 1:
    # add lines to the external forcing file to reference wind input that the user needs to upload!!!
    lines=["QUANTITY=windx",
           "FILENAME=%s/windx.amu" % rel_bc_dir,
           "FILETYPE=4",
           "METHOD=1",
           "OPERAND=O",
           "QUANTITY=windy",
           "FILENAME=%s/windy.amv" % rel_bc_dir,
           "FILETYPE=4",
           "METHOD=1",
           "OPERAND=O",
           ""]
    with open(old_bc_fn ,'at') as fp:
        fp.write("\n".join(lines))

    # tell user to upload wind files!!
    user_instruction_string += ("\nATTENTION!\n" + 
    'User must manually upload the two wind forcing files, named windx.amu and windy.amv, and put them here:\n%s\n' % str(abs_bc_dir) + 
    'You can find scripts for generating windx.amu and windy.amv in the following location on SFEI''s Google Drive:\n' + 
    '\\1_Nutrient_Share\\2_Data_NUTRIENTS\\SFEI_Wind\\Wind4DFlow-SFB-UTC\\\n' + 
    '(https://drive.google.com/drive/folders/1WzqvNo0I0yoWI3KHeDkvtnCrp_MqY5t4)\n' +
    'Make sure to change the permissions of windx.amu and windy.amv, using chmod, once you have uploaded them, so DFM can read them!\n')


# add zero momentum boundary condition for ocean boundary
if 1:

    lines=["QUANTITY=uxuyadvectionvelocitybnd",
           "FILENAME=%s/seauxuy.pli"  % rel_bc_dir,
           "FILETYPE=9",
           "METHOD=2",
           "OPERAND=O"]
    with open(old_bc_fn ,'at') as fp:
        fp.write("\n".join(lines)) 

    shutil.copyfile(os.path.join(abs_static_dir,'seauxuy.pli'),
        os.path.join(run_base_dir,rel_bc_dir,'seauxuy.pli'))

    for i in [1,2,3,4,5,6]:
        shutil.copyfile(os.path.join(abs_static_dir,'seauxuy.tim'),
            os.path.join(run_base_dir,rel_bc_dir,'seauxuy_%04d.tim' % i))


##
mdu_fn = os.path.join(run_base_dir,(run_name + ".mdu")) 
mdu.write(mdu_fn)
print('Just printed out %s.' % mdu_fn)


# print user instructions
print(user_instruction_string)
dum = input('Press any key to continue once you have uploaded the meteo and/or wind files...')


# in addition to making plots, this checks for NaN!!!
if make_plots:
    from sfb_dfm_utils import plot_mdu # SW added function here 
    print('Now making plots...')
    plot_mdu.plot_MDU(mdu_fn, grid_boundary_fn)
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
