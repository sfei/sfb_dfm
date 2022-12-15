from __future__ import print_function
import os, sys
import six
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from stompy.spatial import interp_4d

from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
from stompy.spatial import wkb2shp, proj_utils
from stompy import utils
from stompy.model import unstructured_diffuser

import logging

from stompy.io.local import usgs_sfbay

from . import common # cache_dir

import matplotlib.pylab as plt
import utm
import xarray as xr
from scipy.interpolate import griddata

import logging
logger=logging.getLogger('initial_salinity')

# allie added temperatuere to this too....
def add_initial_salinity(run_base_dir,
                         static_dir, 
                         abs_bc_dir, 
                         old_bc_fn,
                         all_flows_unit=False):
    #static_dir_rel=os.path.relpath(static_dir,run_base_dir)
    bc_dir_rel=os.path.relpath(abs_bc_dir,run_base_dir)
    
    # Spatial salinity initial condition and friction
    lines=[]
    if not all_flows_unit: # real initial condition:
        lines+=[ "QUANTITY=initialsalinity",
                 #"FILENAME=%s/orig-saltopini.xyz"%static_dir_rel, # modified by allie
                 "FILENAME=%s/Salinity-topini.xyz"%bc_dir_rel,
                 "FILETYPE=7",
                 "METHOD=5",
                 "OPERAND=O",
                 ""]
        lines+=[ "QUANTITY=initialtemperature",                      # added by allie
                 "FILENAME=%s/Temperature-topini.xyz"%bc_dir_rel,
                 "FILETYPE=7",
                 "METHOD=5",
                 "OPERAND=O",
                 ""]

    else: #  constant 35 ppt initial condition:
        print("-=-=-=- USING 35 PPT WHILE TESTING! -=-=-=-")
        lines+=[ "QUANTITY=initialsalinity",
                 "FILENAME=constant_35ppt.xyz",
                 "FILETYPE=7",
                 "METHOD=5",
                 "OPERAND=O",
                 ""]
        orig_salt=np.loadtxt(os.path.join(static_dir,'orig-saltopini.xyz'))

        orig_salt[:,2]=35
        np.savetxt(os.path.join(run_base_dir,'constant_35ppt.xyz'),
                   orig_salt,
                   delimiter=' ')
    with open(old_bc_fn,'at') as fp:
        fp.write("\n".join(lines))

def make_initial_sal_temp(run_start, abs_init_dir, abs_bc_dir):

    ''' allie wrote this around Oct 2022 for HAB simulation, then generalized in Dec
    start time is a numpy datetime object correspondign to the start date and time of the model run
    and abs_init_dir is the absolute path to sfb_dfm/init_sal_temp, i.e. the folder where this script is found'''

    # cruise data filename -- prior to running this script, download it manually 
    # from here: https://sfbay.wr.usgs.gov/water-quality-database/
    cruise_data_fn = os.path.join(abs_init_dir,'..','inputs-static','wqdata.csv')
    
    # here's a file with the cruise station coordinates -- found this in /richmondvol1/hpcshared/inputs/stations/
    # hope they're correct!
    cruise_coords_fn = os.path.join(abs_init_dir,'..','inputs-static','station_coords.csv')
    
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




# fancier version pulls some data and extrapolates



# allie comment : damn, I didn't see the following code Rusty wrote, so I guess I might have redone this , woops!


##

def samples_to_cells(init_salt,g):
    """
    init_salt: array [N,3] samples point * {x,y,salt}  
    g: the unstructured_grid to extrapolate onto
    """
    # try again, but with the interp_4d code:
    samples=pd.DataFrame()
    samples['x']=init_salt[:,0]
    samples['y']=init_salt[:,1]
    samples['value']=init_salt[:,2]
    # doesn't really matter, though should be kept in check with alpha
    # and K_j
    samples['weight']=1e6*np.ones_like(init_salt[:,0])

    # alpha=2e-5 is too sharp
    # 5e-6 still too sharp
    # Not sure why higher values looked fine when running this directly,
    # but when it's part of the sfb_dfm.py script it needs really low
    # values of alpha.
    salt=interp_4d.weighted_grid_extrapolation(g,samples,alpha=3e-7)

    # If these fail, the extrapolation approach may be running into
    # numerical difficulties, often made worse by an alpha which is too
    # large (which might be attempting to do less smoothing).
    # decreasing alpha (which results in a smoother field) may help.
    assert np.all( np.isfinite(salt) )
    assert salt.max() < 40
    assert salt.min() > -1 # allow a bit of slop

    # use centroids as they are more predictable
    cc=g.cells_centroid()

    cc_salt=np.concatenate( ( cc, salt[:,None] ),axis=1 )
    return cc_salt


def samples_from_usgs_erddap(run_start,field='salinity'):
    """
    DEPRECATED.  Better to use the dynamically cached data from
    samples_from_usgs().

    --
    Pull some USGS transect data and return a set of 2D salinity
    samples appropriate for the given date.
    Caveats: This is currently using ERDDAP behind the scenes, which
    does not have the most recent data.  It will use a prior year
    if the date cannot be matched within 30 days.

    field: name of the field to pull from usgs_crusies.
       'salinity','temperature'
    """
    # This copy of the USGS data ends early:
    usgs_data_end=np.datetime64('2016-04-28')
    usgs_pad=np.timedelta64(30,'D')

    usgs_target=run_start

    # so we may have to grab a previous years cruise and pretend
    while usgs_target + usgs_pad > usgs_data_end:
        usgs_target -= np.timedelta64(365,'D')

    usgs_cruises=usgs_sfbay.cruise_dataset(usgs_target - usgs_pad,
                                           usgs_target + usgs_pad )

    # lame filling
    scal3d=usgs_cruises[field]
    scal2d=scal3d.mean(dim='prof_sample')
    assert scal2d.dims[0]=='date'
    scal2d_fill=utils.fill_invalid(scal2d.values,axis=0)

    scal_f=interp1d(utils.to_dnum(scal2d.date.values),
                    scal2d_fill,
                    axis=0,bounds_error=False)(utils.to_dnum(usgs_target))

    usgs_init_scal=np.c_[scal2d.x.values,scal2d.y.values,scal_f]
    return usgs_init_scal


def samples_from_usgs(run_start,field='salinity'):
    run_start=utils.to_dt64(run_start)
    pad=np.timedelta64(30,'D')
    
    df=usgs_sfbay.query_usgs_sfbay(period_start=run_start-pad,
                                   period_end=run_start+pad,
                                   cache_dir=common.cache_dir)

    # narrow to the columns we care about:
    if field=='salinity':
        dfield='Salinity'
    elif field=='temperature':
        dfield='Temperature'
    else:
        assert False,"Unknown field %s"%field
        
    df2=df.loc[:,['Date','Station Number','Depth',dfield]]
    # will have to get coordinates from elsewhere

    # depth average:
    df3=df2.groupby( ['Date','Station Number'] )[dfield].mean()

    run_start_dnum=utils.to_dnum(run_start)
    def time_interp(grp):
        dnums=[ utils.to_dnum(v[0]) for v in grp.index.values ]
        return np.interp( run_start_dnum, dnums, grp.values )

    ser4=df3.groupby('Station Number').apply(time_interp)

    lls=[ usgs_sfbay.station_number_to_lonlat(s)
          for s in ser4.index.values ]
    lls=np.array(lls)
    xys=proj_utils.mapper('WGS84','EPSG:26910')(lls)

    # glue those together to get [N,3] array, {x,y,salt}
    return np.c_[ xys, ser4.values ]

##

def samples_from_sfei_erddap(run_start,cache_dir=None):
    """
    return [N,3] array of salinity data from SFEI moorings appropriate for
    given date.  Note that this may have no data, but will be returned
    as a [0,3] array

    This version fetches and caches data from SFEI's ERDDAP server
    """
    if cache_dir is None:
        cache_dir=common.cache_dir
        
    dt_str=utils.to_datetime(run_start).strftime('%Y%m%d%H%M')
    if cache_dir is not None:
        my_cache_dir=os.path.join(cache_dir,'enviz_erddap')
        os.path.exists(my_cache_dir) or os.mkdir(my_cache_dir)

        cache_fn=os.path.join(my_cache_dir,"temp_salt-%s.csv"%dt_str)
        print("Cache fn: %s"%cache_fn)
    else:
        cache_fn=None

    if cache_fn is not None and os.path.exists(cache_fn):
        csv_data=cache_fn
    else:
        # Fetch data before/after run_start by this much
        pad=np.timedelta64(30*60,'s')
        fetch_period=[run_start-pad,run_start+pad]
        fetch_strs=[ utils.to_datetime(p).strftime('%Y-%m-%dT%H:%M:00Z')
                     for p in fetch_period ]

        # Because the table in ERDDAP is stored by sample, there is not guarantee that
        # times are increasing.  That makes access via opendap inefficient, so instead
        # specify the query to ERDDAP more directly, and grab CSV for easier human
        # readability

        # choose dataset
        base_url="http://sfbaynutrients.sfei.org/erddap/tabledap/enviz_mirror.csv"
        # choose fields to download
        params=",".join( ['stationcode','time','spcond_uS_cm','temp_C','stationname',
                          'latitude','longitude'] )
        # And the time range
        criteria="time%%3E=%s&time%%3C=%s"%tuple(fetch_strs)
        url=base_url + "?" + params + "&" + criteria

        import requests
        logging.info("Fetching SFEI data from %s"%url)
        resp=requests.get(url)

        if cache_fn is not None:
            with open(cache_fn,'wt') as fp:
                fp.write(resp.content.decode())
            csv_data=cache_fn
        else:
            csv_data=six.StringIO(resp.content.decode())
        

    # 2nd row of file has units, which we ignore.
    df=pd.read_csv(csv_data,skiprows=[1],parse_dates=['time'])

    # Could get fancier and choose the closest in time reading, or
    # interpolate.  But this is not too bad, averaging over a total of
    # 1 hour.
    dfm=df.groupby('stationcode').mean()

    # Get salinity from specific conductance
    import seawater as sw
    # specific conductance to mS/cm, and ratio to conductivityt at 35 psu, 15 degC.
    # Note that mooring data comes in already adjusted to "specific conductance
    # in uS/cm at 25 degC"
    rt=dfm['spcond_uS_cm'].values/1000. / sw.constants.c3515
    dfm['salinity']=sw.sals(rt,25.0)

    ll=np.c_[dfm.longitude.values,dfm.latitude.values]
    xy=proj_utils.mapper('WGS84','EPSG:26910')(ll)

    xys=np.c_[xy,dfm['salinity'].values]
    valid=np.isfinite(xys[:,2])
    return xys[valid,:]

def samples_from_sfei_moorings(run_start,static_dir):
    """
    return [N,3] array of salinity data from SFEI moorings appropriate for
    given date.  Note that this may have no data, but will be returned
    as a [0,3] array
    """
    # And pull some SFEI data:

    mooring_xy=[]
    mooring_salt=[]

    L2_dir='/opt/data/sfei/moored_sensors_csv/L2/'

    # tuples (<name in observation points shapefile>, <L2 data file name> )
    sfei_moorings=[
        ('ALV',"ALV_all_data_L2.csv"),
        ('SFEI_Coyote',"COY_all_data_L2.csv"),
        ('DB',"DMB_all_data_L2.csv"),
        ('SFEI_Guadalupe',"GL_all_data_L2.csv"),
        ('SFEI_Mowry',"MOW_all_data_L2.csv"),
        ('SFEI_Newark',"NW_all_data_L2.csv"),
        ('SFEI_A8Notch',"POND_all_data_L2.csv"),
        ('SMB',"SM_all_data_L2.csv")
    ]

    # lat/lon from observation-points
    # FRAGILE - FIX!
    obs_shp=wkb2shp.shp2geom(os.path.join(static_dir,"observation-points.shp"))

    for name,l2_file in sfei_moorings:
        print(name)
        fn=os.path.join(L2_dir,l2_file)
        if not os.path.exists(fn):
            logger.warning("No SFEI mooring data - will not be able to initialize LSB initial condition")
            continue
        sfei=pd.read_csv(fn,parse_dates=['Datetime','dt'],low_memory=False)
        sfei_salt=sfei['S_PSU']
        valid=~(sfei_salt.isnull())
        # limit to data within 20 days of the request
        sfei_salt_now=utils.interp_near(utils.to_dnum(run_start),
                                        utils.to_dnum(sfei.Datetime[valid]),sfei_salt[valid],
                                        max_dx=20.0)
        geom=obs_shp['geom'][ np.nonzero(obs_shp['name']==name)[0][0] ]
        xy=np.array(geom)
        if np.isfinite(sfei_salt_now):
            mooring_xy.append(xy)
            mooring_salt.append(sfei_salt_now)

    if len(mooring_xy):
        xy=np.array(mooring_xy)
        sfei_init_salt=np.c_[xy[:,0],xy[:,1],mooring_salt]
    else:
        sfei_init_salt=np.zeros( (0,3),'f8')
    return sfei_init_salt

def initial_salinity_dyn(run_base_dir,
                         mdu,
                         static_dir,
                         run_start):
    # Get some observations:
    usgs_init_salt=samples_from_usgs(run_start)

    # mooring_salt=samples_from_sfei_moorings(run_start,static_dir=static_dir)
    mooring_salt=samples_from_sfei_erddap(run_start)

    init_salt=np.concatenate( (usgs_init_salt,
                               mooring_salt) )
    ##
    g=dfm_grid.DFMGrid(os.path.join(run_base_dir,mdu['geometry','NetFile']))

    # Above here is assembling init_salt
    # Below is extrapolating -- needs g
    cc_salt = samples_to_cells(init_salt,g)


    # Because DFM is going to use some interpolation, and will not reach outside
    # the convex hull, we have to be extra cautious and throw some points out farther
    # afield.

    xys_orig=np.loadtxt(os.path.join(static_dir,'orig-saltopini.xyz'))

    combined_xys=np.concatenate( (cc_salt,xys_orig), axis=0 )

    ##
    return combined_xys


def add_initial_salinity_dyn(run_base_dir,
                             static_dir,
                             mdu,run_start):
    
    # Spatial salinity initial condition 
    lines=[]
    lines+=[ "QUANTITY=initialsalinity",
             "FILENAME=saltopini.xyz",
             "FILETYPE=7",
             "METHOD=5",
             "OPERAND=O",
             ""]
    xys=initial_salinity_dyn(run_base_dir,
                             mdu,static_dir,run_start)
    np.savetxt(os.path.join(run_base_dir,'saltopini.xyz'),
               xys,
               delimiter=' ')
    old_bc_fn=os.path.join(run_base_dir,
                           mdu['external forcing','ExtForceFile'] )
        
    with open(old_bc_fn,'at') as fp:
        fp.write("\n".join(lines))



##

# Pull and cache USGS cruise data directly

