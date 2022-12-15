# Interface for adding COAMPS wind to a DFM run
import os
import re
import time
import requests
import datetime

import numpy as np

import xarray as xr

from stompy.spatial import field
from stompy.grid import unstructured_grid
from stompy import utils
import stompy.model.delft.io as dio

from . import local_config

## 

cache_path=os.path.join(local_config.cache_path,"coamps")

assert os.path.exists(cache_path)

# Anatomy of COAMPS filename:
#                                         nest 3?
#                                         |     HHH - hour of the output from simulation start
#                                         |     |  ??????
#                                         |     |  |     YYYYMMDDHH - start of simulation
#                                         |     |  |     |          elevation code
#                                         |     |  |     |          |    elevation
# US058GMET-GR1dyn.COAMPS-CENCOOS_CENCOOS-n3-c1_00000F0NL2017081600_0100_002000-000000wnd_utru    

missing_files= [ 'cencoos_4km/2017/2017091612/',
                 'cencoos_4km/2017/2017091700/']

def known_missing(recs):
    """
    There are some skipped dates in the COAMPS data online.  This checks
    a download URL against a list of known missing files, returning true
    if the URL is expected to be missing.
    """
    # Get a representative URL
    url=list(recs.values())[0]['url']
    for patt in missing_files:
        if re.search(patt,url):
            return True
    return False

def coamps_files(start,stop):
    """ 
    Generate urls, filenames, and dates for
    fetching or reading COAMPS data

    Tries to pull the first 12 hours of runs, but if a run is known to be missing,
    will pull later hours of an older run.
    
    returns a generator, which yields for each time step of coamps output
    {'wnd_utru':{'url=..., local=..., timestamp=...}, ...}

    """
    dataset_name="cencoos_4km"
    # round to days
    start=start.astype('M8[D]')
    stop=stop.astype('M8[D]') + np.timedelta64(1,'D')

    # The timestamps we're trying for
    target_hours=np.arange(start,stop,np.timedelta64(1,'h'))
    
    for hour in target_hours:
        day_dt=utils.to_datetime(hour)

        # Start time of the ideal run:
        run_start0=day_dt - datetime.timedelta(hours=day_dt.hour%12)

        # runs go for 48 hours, so we have a few chances to get the
        # same output timestamp
        for step_back in [0,12,24,36]:
            run_start=run_start0-datetime.timedelta(hours=step_back)

            # how many hours into this run is the target datetime?
            hour_of_run = int(round((day_dt - run_start).total_seconds() / 3600))

            run_tag="%04d%02d%02d%02d"%(run_start.year,
                                        run_start.month,
                                        run_start.day,
                                        run_start.hour)
            base_url=("http://www.usgodae.org/pub/outgoing/fnmoc/models/"
                      "coamps/calif/cencoos/cencoos_4km/%04d/%s/")%(run_start.year,run_tag)
            
            recs=dict()

            for field_name in ['wnd_utru','wnd_vtru','pres_msl']:
                if field_name in ['wnd_utru','wnd_vtru']:
                    elev_code=105 # 0001: surface?  0105: above surface  0100: pressure?
                    elev=100
                else:
                    # pressure at sea level
                    elev_code=102
                    elev=0
                url_file=("US058GMET-GR1dyn.COAMPS-CENCOOS_CENCOOS-n3-c1_"
                          "%03d"
                          "00F0NL"
                          "%s_%04d_%06d-000000%s")%(hour_of_run,run_tag,elev_code,elev,field_name)

                output_fn=os.path.join(cache_path, dataset_name, url_file)
                recs[field_name]=dict(url=base_url+url_file,
                                      local=output_fn,
                                      timestamp=hour)
            if known_missing(recs):
                continue
            yield recs
            break
        else:
            raise Exception("Couldn't find a run for date %s"%day_dt.strftime('%Y-%m-%d %H:%M'))


def fetch_coamps_wind(start,stop):
    """
    Download all COAMPS outputs between the given np.datetime64()s.
    Does not do any checking against available data, so requesting data
    before or after what is available will fail.

    Returns nothing, just downloads to cache_path
    """
    for recs in coamps_files(start,stop):
        for field_name in recs:
            rec=recs[field_name]
            output_fn=rec['local']
            if os.path.exists(output_fn):
                # print("Skip %s"%os.path.basename(output_fn))
                continue

            print("Fetch %s"%os.path.basename(output_fn))

            response=requests.get(rec['url'],stream=True)

            # Throw an error for bad status codes
            response.raise_for_status()

            with open(output_fn,'wb') as fp:
                for block in response.iter_content(1024):
                    fp.write(block)
            time.sleep(2) # be a little nice.

def coamps_press_windxy_dataset(g_target,start,stop):
    """
    Downloads COAMPS winds for the given period (see fetch_coamps_wind),
    trims to the bounds of g_target, and returns an xarray Dataset.
    """
    fetch_coamps_wind(start,stop)

    xy_min=g_target.nodes['x'].min(axis=0)
    xy_max=g_target.nodes['x'].max(axis=0)

    pad=10e3
    crop=[xy_min[0]-pad,xy_max[0]+pad,
          xy_min[1]-pad,xy_max[1]+pad]

    dss=[] 

    for recs in coamps_files(start,stop):
        timestamp=recs['wnd_utru']['timestamp']
        timestamp_dt=utils.to_datetime(timestamp)
        timestamp_str=timestamp_dt.strftime('%Y-%m-%d %H:%M')
        # use the local file dirname to get the same model subdirectory
        # i.e. cencoos_4km
        cache_fn=os.path.join(os.path.dirname(recs['pres_msl']['local']),
                              "%s.nc"%timestamp_dt.strftime('%Y%m%d%H%M'))
        if not os.path.exists(cache_fn):
            print(timestamp_str)

            # load the 3 fields:
            wnd_utru=field.GdalGrid(recs['wnd_utru']['local'])
            wnd_vtru=field.GdalGrid(recs['wnd_vtru']['local'])
            pres_msl=field.GdalGrid(recs['pres_msl']['local'])

            # Reproject to UTM: these come out as 3648m resolution, compared to 4km input.
            # Fine.  366 x 325.  Crops down to 78x95
            wnd_utru_utm=wnd_utru.warp("EPSG:26910").crop(crop)
            wnd_vtru_utm=wnd_vtru.warp("EPSG:26910").crop(crop)
            pres_msl_utm=pres_msl.warp("EPSG:26910").crop(crop)

            ds=xr.Dataset()
            ds['time']=timestamp
            x,y = wnd_utru_utm.xy()
            ds['x']=('x',),x
            ds['y']=('y',),y
            # copy, in hopes that we can free up ram more quickly
            ds['wind_u']=('y','x'), wnd_utru_utm.F.copy()
            ds['wind_v']=('y','x'), wnd_vtru_utm.F.copy()
            ds['pres']=('y','x'), pres_msl_utm.F.copy()

            ds.to_netcdf(cache_fn)
            ds.close()
        ds=xr.open_dataset(cache_fn)
        ds.load() # force load of data
        ds.close() # and close out file handles
        dss.append(ds) # so this is all in ram.

    ds=xr.concat(dss,dim='time')
    return ds

def write_coamps_press_windxy(g_target,start,stop,output_base):
    """
    Prepare and writes the atmospheric forcing.
    """
    ds=coamps_press_windxy_dataset(g_target,start,stop)
    dio.dataset_to_dfm_wind(ds,start,stop,output_base,
                            extra_header="# derived from COAMPS CENCOOS model output")

def add_coamps_to_mdu(mdu,run_base_dir,g_target,use_existing=True):
    """
    Download, trim, write, and edit FlowFM.ext to add COAMPS wind and
    atmospheric pressure forcing.
    use_existing: if all three files exist in run_base_dir, assume they're
    legit and don't recreate.
    """
    ref_date,start,stop = mdu.time_range()
    pad=np.timedelta64(2,'D')
    start=start-pad
    stop=stop+pad
    
    output_base=os.path.join(run_base_dir,'coamps4km')

    old_bc_fn = os.path.join(run_base_dir,mdu['external forcing','ExtForceFile'])

    meteo_quants=['windx','windy','atmosphericpressure']
    meteo_suffixes=['amu','amv','amp'] # must match code above

    if use_existing:
        for suffix in meteo_suffixes:
            local_path=output_base+"."+suffix
            if not os.path.exists(local_path):
                use_existing=False

    if not use_existing:
        # These need to be synchronized within this file --
        # pretty sure I meant that the filenames and choice of 
        # variables must be consistent between meteo_quants, and
        # what is populated in coamps_press_windxy_dataset
        write_coamps_press_windxy(g_target,start,stop,output_base)

    for suffix,quant in zip(meteo_suffixes,meteo_quants):
        # The filename given to DFM:
        local_path=output_base+"."+suffix
        local_file = os.path.basename(local_path)

        # and add wind to the boundary forcing
        stanza="\n".join(["QUANTITY=%s"%quant,
                          "FILENAME=%s"%local_file,
                          "FILETYPE=4",
                          "METHOD=2",
                          "OPERAND=O",
                          "\n"])
        with open(old_bc_fn,'at') as fp:
            fp.write(stanza)
                       


