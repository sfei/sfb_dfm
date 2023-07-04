import os
import numpy as np
import xarray as xr
import pandas as pd

import stompy.model.delft.io as dio
from stompy import utils
from stompy.io.local import usgs_nwis

from . import dredge_grid

# copied Silvia's pli files to inputs-static
# even though there are 14 of these, one per node of the sea boundary,
# they all appear to have the same data, even the tidal height time series.
# Sea_temp: probably grabbed from Point Reyes?
# Sea_sal: constant 33
# Sea_0001.pli - 


def add_delta_inflow(mdu, base_dir,
                     rel_bc_dir,
                     static_dir,
                     grid,dredge_depth,
                     all_flows_unit=False,
                     temp_jersey=True,
                     temp_rio=True):
    """
    Fetch river USGS river flows, add to FlowFM_bnd.ext:
    Per Silvia's Thesis:
    Jersey: Discharge boundary affected by tides, discharge and temperature taken
    from USGS 11337190 SAN JOAQUIN R A JERSEY POINT, 0 salinity
    (Note that Dutch Slough should probably be added in here)
    Rio Vista: 11455420 SACRAMENTO A RIO VISTA, temperature from DWR station RIV.
    0 salinity.

    run_base_dir: location of the DFM inputs
    run_start,run_stop: target period for therun
    statiC_dir: path to static assets, specifically Jersey.pli and RioVista.pli
    grid: UnstructuredGrid instance, to be modified at inflow locations
    old_bc_fn: path to old-style boundary forcing file
    all_flows_unit: if True, override all flows to be 1 m3 s-1 for model diagnostics
    """
    
    # get run directory and time and forcing file info
    run_base_dir=mdu.base_path
    ref_date,run_start,run_stop = mdu.time_range()
    old_bc_fn=mdu.filepath(["external forcing","ExtForceFile"])
    
    pad=np.timedelta64(3,'D')
    
    if 1: 
        # Cache the original data from USGS, then clean it and write to DFM format
        jersey_raw_fn=os.path.join(run_base_dir,rel_bc_dir,'jersey-raw.nc')
        if not os.path.exists(jersey_raw_fn):
            if temp_jersey==True:
                jersey_raw=usgs_nwis.nwis_dataset(station="11337190",
                                                  start_date=run_start-pad,end_date=run_stop+pad,
                                                  products=[60, # "Discharge, cubic feet per second"
                                                            10], # "Temperature, water, degrees Celsius"
                                                  days_per_request=30)
                jersey_raw.to_netcdf(jersey_raw_fn,engine='scipy')
            if temp_jersey==False:
                jersey_raw=usgs_nwis.nwis_dataset(station="11337190",
                                                  start_date=run_start-pad,end_date=run_stop+pad,
                                                  products=[60], # "Discharge, cubic feet per second" 
                                                  days_per_request=30)
                jersey_raw.to_netcdf(jersey_raw_fn,engine='scipy')                

        rio_vista_raw_fn=os.path.join(run_base_dir,rel_bc_dir,'rio_vista-raw.nc')
        if not os.path.exists(rio_vista_raw_fn):
            if temp_rio==True:
                rio_vista_raw=usgs_nwis.nwis_dataset(station="11455420",
                                                     start_date=run_start-pad,end_date=run_stop+pad,
                                                     products=[60, # "Discharge, cubic feet per second"
                                                               10], # "Temperature, water, degrees Celsius"
                                                     days_per_request=30)
                rio_vista_raw.to_netcdf(rio_vista_raw_fn,engine='scipy')
            if temp_rio==False:
                rio_vista_raw=usgs_nwis.nwis_dataset(station="11455420",
                                                     start_date=run_start-pad,end_date=run_stop+pad,
                                                     products=[60], # "Discharge, cubic feet per second"
                                                     days_per_request=30)
                rio_vista_raw.to_netcdf(rio_vista_raw_fn,engine='scipy')

    if 1: # Clean and write it all out
        jersey_raw=xr.open_dataset(jersey_raw_fn)
        rio_vista_raw=xr.open_dataset(rio_vista_raw_fn)
        temp_logical = [temp_jersey, temp_rio]
        i = 0
        for src_name,source in [ ('Jersey',jersey_raw),
                                 ('RioVista',rio_vista_raw)]:
            src_feat=dio.read_pli(os.path.join(static_dir,'%s.pli'%src_name))[0]            
            dredge_grid.dredge_boundary(grid,src_feat[1],dredge_depth)

            if temp_logical[i]==True:
                # Add stanzas to FlowFMold_bnd.ext:
                for quant,suffix in [('dischargebnd','_flow'),
                                     ('salinitybnd','_salt'),
                                     ('temperaturebnd','_temp')]:
                    with open(old_bc_fn,'at') as fp:
                        lines=["QUANTITY=%s"%quant,
                               "FILENAME=%s/%s%s.pli"%(rel_bc_dir,src_name,suffix),
                               "FILETYPE=9",
                               "METHOD=3",
                               "OPERAND=O",
                               ""]
                        fp.write("\n".join(lines))
                        
                    feat_suffix=dio.add_suffix_to_feature(src_feat,suffix)
                    dio.write_pli(os.path.join(run_base_dir,rel_bc_dir,'%s%s.pli'%(src_name,suffix)),
                                  [feat_suffix])

                    # Write the data:
                    if quant=='dischargebnd':
                        da=source.stream_flow_mean_daily
                        da.values[:] *= 0.028316847 
                        da2=utils.fill_tidal_data(da)
                        if all_flows_unit:
                            da2.values[:]=1.0 
                    elif quant=='salinitybnd':
                        da2=source.stream_flow_mean_daily.copy(deep=True)
                        da2.values[:]=0.0
                    elif quant=='temperaturebnd':
                        da=source.temperature_water
                        da2=utils.fill_tidal_data(da) # maybe safer to just interpolate?
                        if all_flows_unit:
                            da2.values[:]=20.0
                            
                    df=da2.to_dataframe().reset_index()
                    df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
                    columns=['elapsed_minutes',da2.name]
                        
                    if len(feat_suffix)==3:
                        node_names=feat_suffix[2]
                    else:
                        node_names=[""]*len(feat_suffix[1])
                        
                    for node_idx,node_name in enumerate(node_names):
                        # if no node names are known, create the default name of <feature name>_0001
                        if not node_name:
                            node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

                        tim_fn=os.path.join(run_base_dir,rel_bc_dir,node_name+".tim")
                        df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)



            elif temp_logical[i]==False:
                # Add stanzas to FlowFMold_bnd.ext:
                for quant,suffix in [('dischargebnd','_flow'),
                                     ('salinitybnd','_salt'),
                                     ('temperaturebnd','_temp')]:
                    with open(old_bc_fn,'at') as fp:
                        lines=["QUANTITY=%s"%quant,
                               "FILENAME=%s/%s%s.pli"%(rel_bc_dir,src_name,suffix),
                               "FILETYPE=9",
                               "METHOD=3",
                               "OPERAND=O",
                               ""]
                        fp.write("\n".join(lines))
                        
                    feat_suffix=dio.add_suffix_to_feature(src_feat,suffix)
                    dio.write_pli(os.path.join(run_base_dir,rel_bc_dir,'%s%s.pli'%(src_name,suffix)),
                                  [feat_suffix])

                    # Write the data:
                    if quant=='dischargebnd':
                        da=source.stream_flow_mean_daily
                        da.values[:] *= 0.028316847 
                        da2=utils.fill_tidal_data(da)
                        if all_flows_unit:
                            da2.values[:]=1.0 
                    elif quant=='salinitybnd':
                        da2=source.stream_flow_mean_daily.copy(deep=True)
                        da2.values[:]=0.0
                    elif quant=='temperaturebnd':

                        # load climatological temperature
                        df_temp = pd.read_csv(os.path.join(base_dir,'delta_temperature_climatology','Delta_Climatological_Temperature.csv'))
                        day_0 = df_temp['Decimal Day (UTC)'].values
                        temp_0 = df_temp['%s Temp (oC)' % src_name].values

                        # dummy data
                        da2=source.stream_flow_mean_daily.copy(deep=True)
                        time_1 = da.time.values
                        years_1 = pd.DatetimeIndex(time_1).year.values
                        time_1_jan1 = np.array([np.datetime64('%d-01-01' % year) for year in years_1])
                        day_1 = np.floor((time_1 - time_1_jan1)/np.timedelta64(1,'D'))

                        # interpolate to get temperature
                        temp_1 = np.interp(day_1, day_0, temp_0)

                        da2.values[:] = temp_1
                            
                    df=da2.to_dataframe().reset_index()
                    df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
                    columns=['elapsed_minutes',da2.name]
                        
                    if len(feat_suffix)==3:
                        node_names=feat_suffix[2]
                    else:
                        node_names=[""]*len(feat_suffix[1])
                        
                    for node_idx,node_name in enumerate(node_names):
                        # if no node names are known, create the default name of <feature name>_0001
                        if not node_name:
                            node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

                        tim_fn=os.path.join(run_base_dir,rel_bc_dir,node_name+".tim")
                        df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)

            i+=1
