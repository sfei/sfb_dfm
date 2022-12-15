import os
import logging

import numpy as np
import xarray as xr

import stompy.model.delft.io as dio
from stompy.model.delft import dfm_grid

log=logging.getLogger('sfb_dfm_utils')

DAY=np.timedelta64(86400,'s') # useful for adjusting times

def add_erddap_ludwig_wind(run_base_dir,run_start,run_stop,old_bc_fn,fallback=None):
    """
    fetch wind data, write to format supported by DFM, and append wind
    forcing stanzas to old-style DFM boundary forcing file.

    Wind data are fetched from ERDDAP, and written to DFM format

    run_base_dir: path to the run
    run_start,run_stop: target period of the run as datetime64
    old_bc_fn: path to the old-style boundary forcing file

    fallback: None, or [u,v] for constant in space/time when no ludwig available

    returns True is data was found, otherwise return False
    """
    target_filename_base=os.path.join(run_base_dir,"wind") # will have amu, amv appended

    wind_u_fn=target_filename_base+".amu"
    wind_v_fn=target_filename_base+".amv"
    if (os.path.exists(wind_u_fn) and
        os.path.exists(wind_v_fn)):
        log.info('Wind files already exist')
    else:
        data_start=run_start-1*DAY
        data_stop=run_stop+1*DAY

        url='http://sfbaynutrients.sfei.org/erddap/griddap/wind_ludwig_20170621'
        ds=xr.open_dataset(url)

        # somewhat arbitrary cutoff of at least 4 per day
        min_records=(data_stop-data_start)/DAY * 4
        avail_records=dio.dataset_to_dfm_wind(ds,data_start,data_stop,target_filename_base,
                                              min_records=min_records,
                                              extra_header="# downloaded from %s"%url)
        if avail_records<min_records:
            log.warning("Wind data not available for %s -- %s"%(run_start.astype('M8[D]'),
                                                                run_stop.astype('M8[D]')))
            return False

    # and add wind to the boundary forcing
    wind_stanza=["QUANTITY=windx",
                 "FILENAME=wind.amu",
                 "FILETYPE=4",
                 "METHOD=2",
                 "OPERAND=O",
                 "",
                 "QUANTITY=windy",
                 "FILENAME=wind.amv",
                 "FILETYPE=4",
                 "METHOD=2",
                 "OPERAND=O",
                 "\n"]
    with open(old_bc_fn,'at') as fp:
        fp.write("\n".join(wind_stanza))
    return True


def add_constant_wind(run_base_dir,mdu,wind,run_start,run_stop):
    """
    Add a constant in time and space wind field.
    """
    base_wind=xr.Dataset()
    base_wind['wind_u']=wind[0]
    base_wind['wind_v']=wind[1]

    return add_wind_dataset(mdu,base_wind)

def add_wind_dataset(mdu,base_wind):
    """
    mdu: MDUFile object (assumes that mdu.base_path has been set)
    """
    run_base_dir=mdu.base_path
    assert run_base_dir is not None

    grid_fn=mdu.filepath( ['geometry','NetFile'] )
    g=dfm_grid.DFMGrid(grid_fn)

    t_ref,run_start,run_stop=mdu.time_range()

    # manufacture a constant in time, constant in space wind field
    ds=xr.Dataset()
    if 'time' in base_wind.dims:
        ds['time']=('time',),base_wind.time
        if ds.time.values[0]>run_start:
            log.warning("In add_wind_dataset(), start time is after model start")
        if ds.time.values[-1]<run_stop:
            log.warning("In add_wind_dataset(), end time is before model end")
    else:
        DAY=np.timedelta64(1,'D')
        data_start=run_start-1*DAY
        data_stop=run_stop+1*DAY
        t=np.array([data_start, data_stop] )
        ds['time']=('time',),t

    xxyy=g.bounds()
    pad=0.1*(xxyy[1]-xxyy[0])

    if 'x' in base_wind.dims:
        ds['x']=('x',),base_wind.x
    else:
        ds['x']=('x',),np.linspace(xxyy[0]-pad,xxyy[1]+pad,2)

    if 'y' in base_wind.dims:
        ds['y']=('y',),base_wind.y
    else:
        ds['y']=('y',),np.linspace(xxyy[2]-pad,xxyy[3]+pad,3)

    # Here is the magic where xarray broadcasts the dimensions as needed to get
    # (time,y,x) wind data.
    _,_,_,new_u,new_v=xr.broadcast(ds.time,ds.y,ds.x,base_wind.wind_u,base_wind.wind_v)

    ds['wind_u']=new_u 
    ds['wind_v']=new_v

    count=dio.dataset_to_dfm_wind(ds,ds.time.values[0],ds.time.values[-1],
                                  target_filename_base=os.path.join(run_base_dir,"const_wind"),
                                  extra_header="# generated from xarray dataset input")
    assert count>0

    # and add wind to the boundary forcing
    wind_stanza=["QUANTITY=windx",
                 "FILENAME=const_wind.amu",
                 "FILETYPE=4",
                 "METHOD=2",
                 "OPERAND=O",
                 "",
                 "QUANTITY=windy",
                 "FILENAME=const_wind.amv",
                 "FILETYPE=4",
                 "METHOD=2",
                 "OPERAND=O",
                 "\n"]
    old_bc_fn=mdu.filepath(['external forcing','ExtForceFile'])

    with open(old_bc_fn,'at') as fp:
        fp.write("\n".join(wind_stanza))
    return True
