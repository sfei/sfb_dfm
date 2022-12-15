"""
A first cut at direct precipitation and evaporation.

This depends on '/opt/data/cimis/union_city-hourly-2001-2016.nc'
which is currently downloaded through 2017-01-08 on hpc.
"""
import os
import datetime

import xarray as xr
import numpy as np
from stompy import utils

from stompy.io.local import cimis

# Last SUNTANS run had used NARR
# it's way way coarse.  Seems better to use an in-Bay climatology
# than to use NARR.
from . import common

##

def load_cimis(start_date,end_date):
    
    # undo this change rusty made because don't want to deal with putting cimis key in local environment
    # but also change the name to not limit to 2001-2016 and check that dates are included (alliek@sfei.org)
    #union_city=cimis.cimis_fetch_to_xr(stations=171,
    #                                   start_date=start_date,
    #                                   end_date=end_date,
    #                                   cache_dir=common.cache_dir)
                                       
    # union_city=xr.open_dataset('/opt/data/cimis/union_city-hourly-2001-2016.nc')
    union_city=xr.open_dataset('/opt/data/cimis/union_city-hourly.nc')
    
    # add a check that the date range is ok
    uc_start = (np.datetime64(union_city.Date.values[0][0:10]) 
                  + np.timedelta64(union_city.Date.values[0][11:13],'h')
                  + np.timedelta64(union_city.Date.values[0][13:],'m'))
    uc_end = (np.datetime64(union_city.Date.values[-1][0:10]) 
                  + np.timedelta64(union_city.Date.values[-1][11:13],'h')
                  + np.timedelta64(union_city.Date.values[-1][13:],'m'))
    if uc_start>start_date:
        raise Exception('/opt/data/cimis/union_city-hourly.nc start date is after simulation start date')
    if uc_end<end_date:
        raise Exception('/opt/data/cimis/union_city-hourly.nc end date is before simulation end date')

    # https://cals.arizona.edu/azmet/et1.htm
    # which says cool period, divide ETO by 0.7 to get pan evaporation,
    # warm period, divide by 0.6.

    temps=utils.fill_invalid(union_city.HlyAirTmp.values)
    temp_zscore=((temps-temps.mean()) / temps.std()).clip(-1,1)
    # score of 1 means warm temperature
    factors=np.interp(temp_zscore,
                      [-1,1],[0.7,0.6])
    union_city['HlyEvap']=union_city.HlyEto/factors

    union_city.time.values += np.timedelta64(8,'h')
    return union_city

# Data are in mm/hour
# relevant fields are union_city.time, union_city.HlyEvap, and union_city.HlyPrecip
# Times are adjusted from assumed PST to UTC

def add_cimis_evap_precip(run_base_dir,mdu,scale_precip=1.0,scale_evap=1.0):
    pad=np.timedelta64(5,'D')
    t_ref,t_start,t_stop=mdu.time_range()
    data = load_cimis(t_start-pad,t_stop+pad)
                                       
    old_bc_fn=os.path.join(run_base_dir,mdu['external forcing','ExtForceFile'])

    with open(old_bc_fn,'at') as fp:
        # some mentions of "rainfall", other times "rain"
        lines=["QUANTITY=rainfall",
               "FILENAME=precip_evap.tim",
               "FILETYPE=1", # uniform
               "METHOD=1", # ?? interpolate in space and time?
               "OPERAND=O",
               ""]
        fp.write("\n".join(lines))

    sel=((data.time>=t_start-pad)&(data.time<=t_stop+pad)).values
    minutes=(data.time.values[sel] - t_ref) / np.timedelta64(60,'s')

    # emma added the following two lines
    ind = np.where(np.isfinite(data.HlyPrecip.values)==False)
    data.HlyPrecip[ind] = 0.0

    # starting unit: mm/hr
    # What unit do they want? some mentions in the source of mm/hr, other times
    # mm/day.  I think they want mm/day.  Initial look at the rain values in the
    # output confirm that units are good.
    net_precip=24*(scale_precip*data.HlyPrecip.values[sel] - scale_evap*data.HlyEvap.values[sel])

    time_series = np.c_[minutes,net_precip]
    
    np.savetxt(os.path.join(run_base_dir,'precip_evap.tim'),
               time_series,
               fmt="%g")
        

        
