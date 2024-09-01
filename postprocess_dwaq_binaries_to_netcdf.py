#!/usr/bin/env python
# coding: utf-8

# THIS CRASHES FOR BIG RUNS, NEED TO DEBUG oct 2022 alliek

# extract essentials from the DWAQ binaries and save as netcdf, analogous to mapfile
# code by allie king, september 2022
# borrowed from rusty's old residual circulation plotting script
# can use multiple processors to speed up the slow part, which is calculating velocity 
# vectors at cell centers from flows through cell faces


##############
# import stuff
##############

import sys, os
sys.path.append('/opt/software/rusty/stompy/newest_commit/stompy')
import numpy as np
import stompy.model.delft.waq_scenario as dwaq
import matplotlib.pyplot as plt
from stompy.model.stream_tracer import U_perot
from stompy.grid import unstructured_grid
import datetime 
import xarray as xr 
import multiprocessing
from joblib import Parallel, delayed
import matplotlib.colors as mcolors
from scipy.sparse import vstack


#############
# user input
#############

# run name
run_name = 'wy2013-v24'

# path to hyd file (output will be saved in the same directory)
hyd_file = '../runs/%s/DFM_DELWAQ_%s/%s.hyd' % (run_name,run_name,run_name)

# desired output time step in seconds (note open bay dwaq output is every 30 min, and delta is every hour)
dt_out = 1800

# number of processors for computing velocitis in parallel
nproc = 16

# becasue the code crashes if we try to save everything in one file, let's 
# break things into separate time chunks -- ntchunk specifies the number of time
# steps in each chunk 
ntchunk = 48*7*4 # 4 weeks 

# start and end time steps (zero index, set to None for max range)
itstart = None
itstop  = None

# include temperutre? don't try to include temperature if it was not modeled UNLESS it has the same number of 
# time steps as the other model output, otherwise you will get into trouble
include_temperature = False

# plot last time step to see if things look ok???
plot_last_time_step = True

############
# main
############

# get a note about how the dataset was generated
if '__file__' in locals():
    file_path = os.path.realpath(__file__)
else:
    file_path = 'dwaq_binaries_to_netcdf.py'

# get the directory of the binary files and the run base name
hyd_dir = os.path.dirname(hyd_file)
run_name = os.path.basename(hyd_file).split('.')[0]

# output path
output_path = os.path.join(hyd_dir,'%s_from_dwaq_hydro_timechunk%s.nc' % (run_name,'%06d'))

# create a test file in hydro directory to make sure we have permissions before we get into the thick of things
with open(os.path.join(hyd_dir, 'test.txt'),'wt+') as f:
    f.write('hello')
os.remove(os.path.join(hyd_dir, 'test.txt'))

# load hydro 
hydro = dwaq.HydroFiles(hyd_file)

# get info about start time, stop time, time step in seconds, number of time steps
t_start = np.datetime64(hydro.time0)
t_stop = t_start + np.timedelta64(hydro.t_secs[-1],'s')
time_sim = t_start + hydro.t_secs.astype('timedelta64[s]')
dt_sec = (time_sim[1] - time_sim[0])/np.timedelta64(1,'s')
n_time_sim = len(time_sim)

# start and end time
if itstart is None:
    itstart = 0
if itstop is None:
    itstop = n_time_sim

# if user specifies output time step greater than simulation output time step
nskip = int(np.round(dt_out/dt_sec))
if nskip<1:
    print('desired output time step of %f seconds is less ' +
          'than time step of %f seconds, will use %s seconds' % (dt_out,dt_sec,dt_sec))
    nskip = 1

# now get output time and numebr of output time steps
time_all = time_sim[itstart:itstop:nskip]
t_secs_all = hydro.t_secs[itstart:itstop:nskip]
n_time_all = len(time_all)

# compute the number of output files 
numchunks = int(np.ceil(n_time_all/ntchunk))

# do some grid handling stuff 
hydro.infer_2d_elements()
hydro.infer_2d_links()
g = hydro.grid()
unstructured_grid.cleanup_dfm_multidomains(g) # not sure this makes a difference

# not sure what this does but it returns a sparse matrix used to compute velocity
M = hydro.flowlink_to_edge(g) 

# these may be in hydro but I can't find them
n_layers = int(hydro.n_seg/hydro.n_2d_elements)
n_flow_per_layer = int(hydro.n_exch_x / n_layers)

# get the netcdf grid file
grid = hydro.get_geom()

# create xarray dataset with data
nLay = np.arange(n_layers)

# computing velocities is slow, so parallelize it
def one_timestep(it, print_progress=True):

    # time in seconds
    t_sec = t_secs_all[it]

    # read flow and volume from dwaq binaries
    flow    = hydro.flows(t_sec)
    volumes = hydro.volumes(t_sec)

    # don't worry about the vertical exchanges when compute horizontal velocities
    flow = flow[:hydro.n_exch_x]

    # initialize big matrix to store cell velocity by layer
    vel_cell = np.zeros((hydro.n_2d_elements,2,n_layers))

    # sort the flows according to corresponding link number and multiply by sign
    ind = np.argsort(hydro.exch_to_2d_link['link'])
    flow = flow[ind] * hydro.exch_to_2d_link['sgn'][ind]
    
    # reshape the volumes so all segments in layer 0 are followed by all segments in layer 1, etc
    volumes = volumes.reshape((n_layers, hydro.n_2d_elements))

    # reorder the flows matrix so, like the volume matrix, all links in layer 0 are followed by  
    # all links in layer 1, etc
    flow = np.reshape(flow, (n_flow_per_layer, n_layers)).transpose()
    
    # go layer by layer vertically
    for ilay in range(n_layers):
    
        # not sure what this does
        flow_edge = M.dot(flow[ilay,:])
    
        # compute velocity at cell center, inferred from velocites normal to cell edges
        vel_cell[:,:,ilay] = U_perot(g, flow_edge, volumes[ilay,:])

    # print progress
    if print_progress:  
        timer_out = (datetime.datetime.now()-timer_start).seconds/3600
        print('%0.2f hours since start %s, on time step %d of %d' % (timer_out, timer_start_string, it, n_time_all))

    return vel_cell

print('')
print('')
print('netcdf files will be saved here: %s' % output_path)
print('there will be %d netcdf files with %d time steps each' % (numchunks,ntchunk))


# start a timer
timer_start = datetime.datetime.now()
timer_start_string = timer_start.strftime('%Y-%h-%m %H:%M:%S')

print('beginning to loop through %d chunks at time %s' % (numchunks, timer_start_string))

# loop through the chunks
for ichunk in range(numchunks):

    # range of time steps in this chunk
    it1 = ichunk*ntchunk
    it2 = (ichunk+1)*ntchunk
    if it2>n_time_all:
        it2 = n_time_all
    ntime = it2-it1

    # time
    time = time_all[it1:it2]

    # initialize output
    volume = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))
    salinity = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))
    vertical_diffusivity = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))
    shear_stress = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))
    temperature = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))
    velocity_x = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))
    velocity_y = np.nan*np.ones((ntime, hydro.n_2d_elements, n_layers))

    # first extract everything except velocity 
    for it in range(it1,it2):

        # chunk time step
        itc = it-it1

        # time in seconds
        t_sec = t_secs_all[it]
    
        # read the dwaq binaries at this time step and save in output matrix
        volume[itc,:,:] = np.reshape(hydro.volumes(t_sec) , (n_layers, hydro.n_2d_elements)).transpose()
        salinity[itc,:,:] = np.reshape(hydro.seg_func(t_sec,label='salinity-file'), (n_layers, hydro.n_2d_elements)).transpose()
        vertical_diffusivity[itc,:,:] = np.reshape(hydro.seg_func(t_sec,label='vert-diffusion-file'), (n_layers, hydro.n_2d_elements)).transpose()
        shear_stress[itc,:,:] = np.reshape(hydro.seg_func(t_sec,label='shear-stresses-file'), (n_layers, hydro.n_2d_elements)).transpose()
        if include_temperature:
            temperature[itc,:,:] = np.reshape(hydro.seg_func(t_sec,label='temperature-file'), (n_layers, hydro.n_2d_elements)).transpose()
    
    # compute the water depth from the volumes
    water_depth = np.sum(volume,axis=2) / np.tile(g.cells_area(),(volume.shape[0],1))

    # compute velocity in parallel
    pout = np.array(Parallel(n_jobs=nproc)(delayed(one_timestep)(it) for it in range(it1,it2)))
    velocity_x            = pout[:,:,0,:]
    velocity_y            = pout[:,:,1,:]
    
    # create and save as an xarray dataset
    ds = xr.Dataset({
        'layernum': xr.DataArray(
                    data   = nLay,   
                    dims   = ['nLay'],
                    coords = {'nLay' : nLay},
                    attrs  = {
                        'description' : 'vertical layer index (0 is surface)'
                        }
                    ),
        'FlowElem_bl': xr.DataArray(
                    data   = grid.FlowElem_bl.values,   
                    dims   = ['nFlowElem'],
                    coords = {'nFlowElem' : grid.nFlowElem.values},
                    attrs  = {
                        'description' : 'bed elevation w.r.t. MLLW', 
                        'units'     : 'm'
                        }
                    ),
        'FlowElem_xcc': xr.DataArray(
                    data   = grid.FlowElem_xcc.values,   
                    dims   = ['nFlowElem'],
                    coords = {'nFlowElem' : grid.nFlowElem.values},
                    attrs  = {
                        'description' : 'x coordinate (UTM Zone 10)',
                        'units'     : 'm'
                        }
                    ),
        'FlowElem_ycc': xr.DataArray(
                    data   = grid.FlowElem_ycc.values,  
                    dims   = ['nFlowElem'],
                    coords = {'nFlowElem' : grid.nFlowElem.values},
                    attrs  = {
                        'description' : 'y coordinate (UTM Zone 10)',
                        'units'     : 'm'
                        }
                    ),
        'salinity': xr.DataArray(
                    data   = salinity,   
                    dims   = ['time','nFlowElem','nLay'],
                    coords = {'time': time,
                              'nFlowElem' : grid.nFlowElem.values,
                              'nLay' : nLay},
                    attrs  = {
                        'units'     : 'ppt ???'
                        },
                    ),
        'temperature': xr.DataArray(
                    data   = temperature,   
                    dims   = ['time','nFlowElem','nLay'],
                    coords = {'time': time,
                              'nFlowElem' : grid.nFlowElem.values,
                              'nLay' : nLay},
                    attrs  = {
                        'units'     : 'oC'
                        },
                    ),
        'vertical_diffusivity': xr.DataArray(
                    data   = vertical_diffusivity,   
                    dims   = ['time','nFlowElem','nLay'],
                    coords = {'time': time,
                              'nFlowElem' : grid.nFlowElem.values,
                              'nLay' : nLay},
                    attrs  = {
                        'units'     : 'm2/s ???'
                        },
                    ),
        'shear_stress': xr.DataArray(
                    data   = shear_stress,   
                    dims   = ['time','nFlowElem','nLay'],
                    coords = {'time': time,
                              'nFlowElem' : grid.nFlowElem.values,
                              'nLay' : nLay},
                    attrs  = {
                        'units'     : 'N/m2 ???'
                        },
                    ),
        'water_depth': xr.DataArray(
                    data   = water_depth,   
                    dims   = ['time','nFlowElem'],
                    coords = {'time': time,
                              'nFlowElem' : grid.nFlowElem.values},
                    attrs  = {
                        'units'     : 'm'
                        },
                    ),
        'velocity_x' : xr.DataArray(
                    data   = velocity_x,   
                    dims   = ['time','nFlowElem','nLay'],
                    coords = {'time': time,
                              'nFlowElem' : grid.nFlowElem.values,
                              'nLay' : nLay},
                    attrs  = {
                          'units'     : 'm/s'
                          },
            ),
        'velocity_y' : xr.DataArray(
            data   = velocity_y,   
            dims   = ['time','nFlowElem','nLay'],
            coords = {'time': time,
                      'nFlowElem' : grid.nFlowElem.values,
                      'nLay' : nLay},
            attrs  = {
                'units'     : 'm/s'
                },
            )
    
        },
    
            attrs = {'origin' : 'generated at %s by %s' % (datetime.datetime.now(), file_path)}
    
        ) 

    # print progress
    timer_out = (datetime.datetime.now()-timer_start).seconds/3600
    print('saving chunk %d of %d, %0.2f hours since start %s' % (ichunk+1, numchunks, timer_out, timer_start_string))
            
    # saving as netcdf file
    ds.to_netcdf(output_path % ichunk)

###########################################################################################
# for troubleshooting, may want to plot stuff ... let's try final time step of first chunk
###########################################################################################

if 1:

    # load dataset 
    data = xr.open_dataset(output_path % 0)

    # plot water depth
    fig1, ax1 = plt.subplots()
    h1 = g.plot_cells(values=data.water_depth.values[0,:], ax=ax1, cmap='jet')
    h1.set_clim([0,10])
    ax1.set_title('water depth (m)\nat time %s' % data.time.values[-1])
    fig1.colorbar(h1, ax=ax1, fraction=0.02,pad=0.01)
    ax1.axis('off')

    # plot top and bottom of a bunch of things
    def plot_top_bottom(varname,vmin,vmax):

        fig1, ax1 = plt.subplots(1,2)
        h1 = g.plot_cells(values=data[varname].values[0,:,0], ax=ax1[0], cmap='jet')
        h1.set_clim([vmin,vmax])
        ax1[0].set_title('top')
        h2 = g.plot_cells(values=data[varname].values[0,:,-1], ax=ax1[1], cmap='jet')
        h2.set_clim([vmin,vmax])
        ax1[1].set_title('bottom')
        fig1.suptitle('%s (%s)\nat time %s' % (varname,data[varname].units,data.time.values[-1]))
        fig1.colorbar(h1, ax=ax1.ravel().tolist(), fraction=0.02,pad=0.01)
        ax1[0].axis('off')
        ax1[1].axis('off')


    plot_top_bottom('salinity',0,33)
    plot_top_bottom('vertical_diffusivity',0,0.01)
    plot_top_bottom('shear_stress',0,0.8)
    plot_top_bottom('velocity_x',-0.4,0.4)
    plot_top_bottom('velocity_y',-0.4,0.4)
