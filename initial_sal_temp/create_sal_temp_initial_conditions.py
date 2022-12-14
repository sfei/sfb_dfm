import os, sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from stompy.spatial import interp_4d, proj_utils
from stompy.grid import unstructured_grid
import utm
import xarray as xr
from scipy.interpolate import griddata

## user inputs the start time of the simulation here
#run_start = np.datetime64('2003-08-01')

def make_initial_sal_temp(run_start, abs_init_dir, abs_bc_dir):

    ''' start time is a numpy datetime object correspondign to the start date and time of the model run
    and abs_init_dir is the absolute path to sfb_dfm/init_sal_temp, i.e. the folder where this script is found'''

    # cruise data filename -- prior to running this script, download it manually 
    # from here: https://sfbay.wr.usgs.gov/water-quality-database/
    cruise_data_fn = os.path.join(abs_init_dir,'wqdata.csv')
    
    # here's a file with the cruise station coordinates -- found this in /richmondvol1/hpcshared/inputs/stations/
    # hope they're correct!
    cruise_coords_fn = os.path.join(abs_init_dir,'station_coords.csv')
    
    # ROMS data filename -- this is derived from the CASCaDE model ocean temperature boundary condition, which is 
    # set to the climatalogical average from a 1980 - 2010 ROMS-based reanalysis, cited in Vroom et al 2017... 
    # Vroom says it is from the ROMS-based reanalysis but does not mention it is an average from 1980 - 2010... this
    # info is from personal communication with Noah Knowles. Rusty originally set up sfb_dfm to use temperatures from
    # Point Reyes -- these appear to be higher than the ROMS-based temperatures. Going forward, may want to evaluate 
    # various options, possibly a combo of ROMS and Point Reyes, maybe ask Noah how we can get original ROMS data 
    # and find an offset between Point Reyes and different points along the boundary????? If we can model tempearature
    # well enough without ocean temperature varying from year to year, don't have to worry about it, but if we are having
    # problems this is something we can troubleshoot
    roms_fn = os.path.join(abs_init_dir,'..','inputs-static','sea_temp_ROMS.nc')
    
    # path to the model grid (doesn't matter if it's the same exact one we use for the simulation, just need somthing to 
    # interpolate onto)
    grid_fn = os.path.join(abs_init_dir,'..','sfei_v20_net.nc')
    
    # time window for including cruise data (search on either side of the start time, so 15 days means a total 30 day window)
    time_window = np.timedelta64(15,'D')
    
    # set a constant ocean boundary salinity 
    ocean_salinity = 33.0
    
    # parameters for Rusty's extrapolation algorithm (not sure of units of time scale)
    alpha = 1e-5 # smoothness in space
    beta = 0.5   # smoothness in time
    time_scale = 2.0 # time scale for temporal smoothing 
    
    # specify a grid resolution (meters)
    dxy = 1000
    
    ############## read in cruise station coordinates #########################
    
    # read cruise station coordinates and create dictionaries to map station number to latitude and longitude
    coords = pd.read_csv(cruise_coords_fn,sep='\t')
    stn2xutm = {}
    stn2yutm = {}
    for ic in range(len(coords)):
        stn = coords.iloc[ic]['Station']
        lat = coords.iloc[ic]['Latitude (degrees)']
        lon = coords.iloc[ic]['Longitude (degrees)']
        x, y, _, _ = utm.from_latlon(lat,lon) 
        stn2xutm[stn] = x 
        stn2yutm[stn] = y 
    
    ############## munge the cruise data #########################
    
    # read the USGS cruise data
    cruise_data = pd.read_csv(cruise_data_fn)
    
    # add a numpy datetime64 time column to wqdata, note that wqdata time stamps are in local time (PT), 
    # need to convert to UTC because model is in UTC, and then need to remove timezone info because I'm not sure
    # the stompy time handling utilities can cope (they probably can but I'm trying to just copy Rusty's grid
    # extrapolation example to a tee without thinking about it) 
    datetime = []
    for it in range(len(cruise_data)):
        date = cruise_data.iloc[it]['Date']
        time = cruise_data.iloc[it]['Time']
        if len(time)<5:
            time = '0' + time
        datetime.append(pd.Timestamp('%s %s' % (date,time),tz='US/Pacific').tz_convert('UTC').replace(tzinfo=None))
    cruise_data['time'] = np.array(datetime).astype('datetime64[s]')
    
    # only use the cruise data within the prescribed time window around the start date
    ind = np.abs(cruise_data['time'].values - run_start) <= time_window
    if np.any(ind):
        cruise_data = cruise_data.loc[ind]
    else:
        raise Exception('no data found in %s within %s of %s; go to '  % (cruise_data_fn,time_window,run_start) + 
                        'https://sfbay.wr.usgs.gov/water-quality-database/, download the data in this time window, and' + 
                        'manually place it in the sfb_dfm/initial_sal_temp/ folder. sorry this part is not automated yet!')
    
    # add latitude and longitude to cruise data, and rename the temperature column
    cruise_data['x'] = cruise_data['Station'].map(stn2xutm)
    cruise_data['y'] = cruise_data['Station'].map(stn2yutm)
    cruise_data.rename(columns={'Temperature (Degrees Celsius)':'Temperature'}, inplace=True)
    
    ############## put the cruise data into a dataframe with correct format for the interpolator #########################
    
    # initialize a master dataframe that will hold depth-averaged cruise data along with the ROMS data
    data = pd.DataFrame(columns=['time','x','y','Temperature','Salinity'])
    
    # make the assumption that all depths at a given station on a given sampling date have the same time stamp
    time = np.unique(cruise_data['time'])
    for it in range(len(time)):
    
        ind = cruise_data['time'] == time[it]
        data_append = cruise_data.loc[ind][['x','y','Temperature','Salinity']].mean()
        data_append['time'] = time[it]
        data.loc[len(data)] = data_append
    
    ############## add the ROMS-based sea temperatures at the ocean boundary and assume constant salinity there #########################
    
    # now, bring in the ROMS sea temperatures at the ocean boundary
    roms = xr.open_dataset(roms_fn)
    
    # time in ROMS dataset is days since January 1, so compute time axis w.r.t. year containing the start date
    start_year = pd.Timestamp(run_start).year
    time_ROMS = np.datetime64('%d-01-01' % start_year) + roms.time.values*np.timedelta64(1,'D')
    
    # grab the ROMS temperatures at the start date
    ist = np.argmin(np.abs((time_ROMS - run_start)))
    temp_ROMS = roms.temperature.values[ist,:]
    
    # note number of roms stations
    nst = len(roms.station)
    
    # add each station
    for ist in range(nst):
    
        dict_append = {'time' : '%s 00:00:00' % run_start.astype('datetime64[D]'),
                       'x' : roms.xutm.values[ist],
                       'y' : roms.yutm.values[ist],
                       'Temperature' : roms.temperature.values[it,ist],
                       'Salinity' : ocean_salinity}
    
        data.loc[len(data)] = dict_append
    
    # fix the dtype, not sure why everythign is an object
    data = data.astype( dtype={'time' : 'datetime64[ns]',
                           'x' : 'float',
                           'y' : 'float',
                           'Temperature' : 'float',
                           'Salinity' : 'float'})
    
    # load the grid, get the cell center coordinates 
    grid = unstructured_grid.UnstructuredGrid.read_dfm(grid_fn,cleanup=True)
    xy = grid.cells_center()
    xx = xy[:,0]
    yy = xy[:,1]
    
    
    # come up with a nice regular grid, based on the extent of the model grid and the resplution 
    # specified by the user, only including points within one grid cell from the model grid
    xmax, ymax = xy.max(axis=0)
    xmin, ymin = xy.min(axis=0)
    xmax = np.ceil(xmax/dxy)*dxy
    ymax = np.ceil(ymax/dxy)*dxy
    xmin = np.floor(xmin/dxy)*dxy
    ymin = np.floor(ymin/dxy)*dxy
    xg, yg = np.meshgrid(np.arange(xmin,xmax+dxy,dxy), np.arange(ymin,ymax+dxy,dxy))
    xg = xg.flatten()
    yg = yg.flatten()
    ng = len(xg)
    include = np.zeros(ng, dtype=bool)
    for ig in range(ng):
        if np.min(np.sqrt((xg[ig]-xx)**2 + (yg[ig]-yy)**2)) <= dxy:
            include[ig] = True
    xg = xg[include]
    yg = yg[include]
    
    # make initial condition files for temperature and salinity
    for value_col in ['Temperature', 'Salinity']:
    
        # here's where we apply Rusty's widget to interpolate as the fish swims, with some kind of weighting in time
        # data2d maps to the grid cell centers
        time_int=interp_4d.interp_to_time_per_ll(data, run_start, lat_col='y', lon_col='x', value_col=value_col)
        time_int['weight']= (1+time_int['time_offset']/time_scale)**(-beta)
        smooth = interp_4d.weighted_grid_extrapolation(grid, time_int, alpha=alpha, value_col=value_col)
        assert np.all(np.isfinite(smooth)),"Error -- getting some non-finite values"
        
        # now we need to interpolate onto the regular grid xg, yg
        ongrid = griddata(xy, smooth, (xg, yg), method='linear')
    
        # linear method doesn't extrapolate, so use nearest neighbor method to fill in nan's
        ongrid_nn = griddata(xy, smooth, (xg, yg), method='nearest')
        ind = np.isnan(ongrid)
        ongrid[ind] = ongrid_nn[ind]
    
        # find max and min for colorbar
        cmin = np.min(ongrid)
        cmax = np.max(ongrid)
    
        # now plot the gridded interpolated field along with the data used to generate it 
        # (note not all of the cruise data is used, only data close enough in time)
        fig, ax = plt.subplots(figsize=(11.5, 8))
        sc = ax.scatter(xg, yg, s=10, c=ongrid, cmap='jet', vmin=cmin, vmax=cmax)
        cb = plt.colorbar(sc)
        cb.set_label(value_col)
        ax.axis('off')
    
        # add the input data and time offset -- note if have multiple measuremetns at 
        # same point in space, these will write on top of each other, oh well
        ax.scatter(time_int['x'], time_int['y'], s=200, c=time_int[value_col], cmap='jet', vmin=cmin, vmax=cmax)
        for i in range(len(time_int)):
            time_int1 = time_int.loc[i]
            ax.text(time_int1['x'], time_int1['y'], '%0.1f' % time_int1['time_offset'], ha='center', va='center')
    
        # save plot
        ax.set_title(('%s\n%s initial condition (small dots) shown with input data (large dots).\n' % (run_start,value_col) +
                      'Time offset marked in center of input data; if there are multiple input\n' +
                      'data points in same spatial location at different times, they write on\n' + 
                      'top of each other. Just trying to do a reality check, not make a perfect plot'))
        plt.tight_layout()
        fig.savefig('%s_initial_condition.png' % value_col)
    
        # finally, write the initial condition!!!
        out_fn = os.path.join(abs_bc_dir, '%s-topini.xyz' % value_col)
        print('writing %s ...' % out_fn)
        with open(out_fn, 'wt+') as f:
    
            for i in range(len(ongrid)):
    
                f.write('%16.7f%16.7f%16.7f\n' % (xg[i], yg[i], ongrid[i]))
    
    
    
    