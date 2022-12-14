import os, sys
import numpy as np
import xarray as xr
import pandas as pd

def make_ROMS_temp_tim(abs_init_dir,abs_bc_dir,ref_date,run_start,run_stop):

    ## reference, start, and end time of the simulation
    #ref_date = np.datetime64('2003-08-01')
    #run_start = np.datetime64('2003-08-01')
    #run_stop = np.datetime64('2004-10-01')
    
    # units for *.tim file times
    units = 'm' # minutes
    
    # daily time step seems fine
    deltat = np.timedelta64(1,'D').astype(units)
    
    # create the time axis
    time = np.arange(run_start,run_stop+deltat,deltat)
    
    # convert time to day of year, since that's what the ROMS data is referenced to
    year = pd.Timestamp(run_start).year
    day = (time - np.datetime64('%d-01-01' % year))/np.timedelta64(1,'D')
    
    # read ROMS data
    data = xr.open_dataset(os.path.join(abs_init_dir,'..','inputs-static','sea_temp_ROMS.nc'))
    
    # convert time axis to time in units of tim file
    time_units = (time - ref_date)/np.timedelta64(1,units)
    
    # for each station, interpolate data onto time axis and write to tim file
    for st in range(1,7):
        data1 = data.sel(station=st)
        temp = np.interp(day, data1.time.values, data1.temperature.values)
        
        outfn = os.path.join(abs_bc_dir,'sea_temp_ROMS_000%d.tim' % st)
        with open(outfn, 'wt+') as f:
            print('writing %s' % outfn)
            for i in range(len(time_units)):
                f.write('%0.1f %0.4f\n' % (time_units[i], temp[i]))
    
    