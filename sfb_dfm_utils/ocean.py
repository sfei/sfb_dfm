import os,sys
import logging
log=logging.getLogger('sfb_dfm_utils')

import numpy as np
import xarray as xr
import pandas as pd

from stompy import (utils, filters)
from stompy.io.local import noaa_coops
import stompy.model.delft.io as dio
from . import common


# fill any missing data via linear interpolation
def fill_data(da):
    valid=np.isfinite(da.values)
    valid_times=da.time[valid]
    dt=np.median(np.diff(da.time))
    dt_gap = np.diff(da.time[valid]).max()
    if dt_gap > dt:
        log.warning("%s: gaps up to %.1f minutes"%(da.name, dt_gap/np.timedelta64(60,'s')))
    da.values = utils.fill_invalid(da.values)    

def add_ocean(run_base_dir,
              run_start,run_stop,ref_date,
              static_dir,
              grid,old_bc_fn,
              all_flows_unit=False,
              lag_seconds=0.0,
              factor=1.0):
    """
    Ocean:
    Silvia used:
        Water level data from station 46214 (apparently from Yi Chao's ROMS?)
          no spatial variation
        Maybe salinity from Yi Chao ROMS?  That's what the thesis says, but the
        actual inputs look like constant 33
    Here I'm using data from NOAA Point Reyes.
        waterlevel, water temperature from Point Reyes.
    When temperature is not available, use constant 15 degrees

    factor: a scaling factor applied to tide data to adjust amplitude around MSL.
    lag_seconds: to shift ocean boundary condition in time, a positive value 
    applying it later in time.
    """
    # get a few extra days of data to allow for transients in the low pass filter.
    pad_time=np.timedelta64(5,'D')
    
    if 1:
        if 0: # This was temporary, while NOAA had an issue with their website. 
            log.warning("TEMPORARILY USING FORT POINT TIDES")
            tide_gage="9414290" # Fort Point
        else:
            tide_gage="9415020" # Pt Reyes 

        tides_raw_fn=os.path.join(run_base_dir,'tides-%s-raw.nc'%tide_gage)
        if not os.path.exists(tides_raw_fn):
            tides_raw=noaa_coops.coops_dataset(tide_gage,
                                               run_start-pad_time,run_stop+pad_time,
                                               ["water_level","water_temperature"],
                                               days_per_request=30)

            tides_raw.to_netcdf(tides_raw_fn,engine='scipy')

    # Fort Point mean tide range is 1.248m, vs. 1.193 at Point Reyes.
    # apply rough correction to amplitude.
    # S2 phase 316.2 at Pt Reyes, 336.2 for Ft. Point.
    # 20 deg difference for a 12h tide, or 30 deg/hr, so
    # that's a lag of 40 minutes.
    # First go I got this backwards, and wound up with lags
    # at Presidio and Alameda of 4600 and 4400s.  That was
    # with lag_seconds -= 40*60.
    # Also got amplitudes 13% high at Presidio, so further correction...
    if tide_gage=="9414290":
        # 
        factor *= 1.193 / 1.248 * 1.0/1.13
        lag_seconds += 35*60.

    if 1:
        # Clean that up, fabricate salinity
        tides=xr.open_dataset(tides_raw_fn).isel(station=0)

        water_level=utils.fill_tidal_data(tides.water_level)

        # IIR butterworth.  Nicer than FIR, with minor artifacts at ends
        # 3 hours, defaults to 4th order.
        water_level[:] = filters.lowpass(water_level[:].values,
                                         utils.to_dnum(water_level.time),
                                         cutoff=3./24)

        if 1: # apply factor:
            msl= 2.152 - 1.214 # MSL(m) - NAVD88(m) for Point Reyes
            if factor!=1.0:
                log.info("Scaling tidal forcing amplitude by %.3f"%factor)
            water_level[:] = msl + factor*(water_level[:].values - msl)

        if 1: # apply lag
            if lag_seconds!=0.0:
                # sign:  if lag_seconds is positive, then I want the result
                # for time.values[0] to come from original data at time.valules[0]-lag_seconds
                if 0: # Why interpolate here? Just alter the timebase.
                    water_level[:] = np.interp( utils.to_dnum(tides.time.values),
                                                utils.to_dnum(tides.time.values)-lag_seconds/86400.,
                                                tides.water_level.values )
                else:
                    # Adjust time base directly.
                    water_level.time.values[:] = water_level.time.values + np.timedelta64(lag_seconds,'s')
            
        if 'water_temperature' not in tides:
            log.warning("Water temperature was not found in NOAA data.  Will use constant 15")
            water_temp=15+0*tides.water_level
            water_temp.name='water_temperature'
        else:
            fill_data(tides.water_temperature)
            water_temp=tides.water_temperature
                                             
        if all_flows_unit:
            print("-=-=-=- USING 35 PPT WHILE TESTING! -=-=-=-")
            salinity=35 + 0*water_level
        else:
            salinity=33 + 0*water_level
        salinity.name='salinity'
            
    if 1: # Write it all out
        # Add a stanza to FlowFMold_bnd.ext:
        src_name='Sea'
        
        src_feat=dio.read_pli(os.path.join(static_dir,'%s.pli'%src_name))[0]
        
        forcing_data=[('waterlevelbnd',water_level,'_ssh'),
                      ('salinitybnd',salinity,'_salt'),
                      ('temperaturebnd',water_temp,'_temp')]

        for quant,da,suffix in forcing_data:
            with open(old_bc_fn,'at') as fp:
                lines=["QUANTITY=%s"%quant,
                       "FILENAME=%s%s.pli"%(src_name,suffix),
                       "FILETYPE=9",
                       "METHOD=3",
                       "OPERAND=O",
                       ""]
                fp.write("\n".join(lines))

            feat_suffix=dio.add_suffix_to_feature(src_feat,suffix)
            dio.write_pli(os.path.join(run_base_dir,'%s%s.pli'%(src_name,suffix)),
                          [feat_suffix])

            # Write the data:
            columns=['elapsed_minutes',da.name]

            df=da.to_dataframe().reset_index()
            df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')

            if len(feat_suffix)==3:
                node_names=feat_suffix[2]
            else:
                node_names=[""]*len(feat_suffix[1])

            for node_idx,node_name in enumerate(node_names):
                # if no node names are known, create the default name of <feature name>_0001
                if not node_name:
                    node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

                tim_fn=os.path.join(run_base_dir,node_name+".tim")
                df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)