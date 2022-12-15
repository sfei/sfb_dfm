"""
Add freshwater inputs to a model based on a combination of 
gaged streams, and scaling from gaged watersheds to ungaged.

This method was used in the SUNTANS SF Bay model, but has
been largely supplanted by modeled flows from the Bay Area 
Hydrologic Model for 2000--2016.

This code requires an input shapefile giving gages, locations,
and adjustments for area.
"""


import os
import numpy as np
import xarray as xr
import pandas as pd
import logging

import stompy.model.delft.io as dio
from stompy.spatial import wkb2shp

try:
    from . import dredge_grid
    from . import local_config
    usgs_inventory_shp_fn=os.path.join(os.path.dirname(__file__),'usgs_inventory.shp')
except ImportError: # testing..
    print("TESTING")
    from sfb_dfm_utils import dredge_grid
    from sfb_dfm_utils import local_config
    usgs_inventory_shp_fn="sfb_dfm_utils/usgs_inventory.shp"
    

DAY=np.timedelta64(86400,'s') # useful for adjusting times


FT3_to_M3=0.028316847

from stompy.io.local import usgs_nwis

def write_QST_data(mdu,stn_ds,src_name,time_offset=None):
    """
    write flow, salinity and temperature time series from an xarray 
    Dataset.
    mdu: MDUFile for the run.  expects ref_date, ExtForceFile, base_path.
    src_name: sanitized name for filenames and ExtForceFile.
    time_offset: offset to apply to stn_ds.time.  Note that the sign
    here is opposite of add_sfbay_freshwater, since it is being added
    to stn_ds, instead of start_date.  
    """

    old_bc_fn=mdu.filepath( ('external forcing','ExtForceFile') )
    ref_date,run_start,run_stop=mdu.time_range()

    df=stn_ds.to_dataframe().reset_index()
    df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
    # default values:
    df['salinity']=0*df.flow_cms
    df['temperature']=20+0*df.flow_cms

    assert np.all( np.isfinite(df.flow_cms) )
    
    quant_suffix=[ ('dischargebnd','_flow') ]
    
    if mdu['physics','temperature']:
        quant_suffix.append( ('temperaturebnd','_temp') )
    if mdu['physics','salinity']:
        quant_suffix.append( ('salinitybnd','_salt') )
        
    for quantity,suffix in quant_suffix:
        lines=['QUANTITY=%s'%quantity,
               'FILENAME=%s%s.pli'%(src_name,suffix),
               'FILETYPE=9',
               'METHOD=3',
               'OPERAND=O',
               ""]
        with open(old_bc_fn,'at') as fp:
            fp.write("\n".join(lines))

        # read the pli back to know how to name the per-node timeseries
        feats=dio.read_pli(os.path.join(mdu.base_path,
                                        "%s%s.pli"%(src_name,suffix)))
        feat=feats[0] # just one polyline in the file

        if len(feat)==3:
            node_names=feat[2]
        else:
            node_names=[""]*len(feat[1])

        for node_idx,node_name in enumerate(node_names):
            # if no node names are known, create the default name of <feature name>_0001
            if not node_name:
                node_name="%s%s_%04d"%(src_name,suffix,1+node_idx)

            tim_fn=os.path.join(mdu.base_path,node_name+".tim")

            columns=['elapsed_minutes']
            if quantity=='dischargebnd':
                columns.append('flow_cms')
            elif quantity=='salinitybnd':
                columns.append('salinity')
            elif quantity=='temperaturebnd':
                columns.append('temperature')

            df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)
            # To avoid hitting the limit of open files, only write the first
            # node.  It's not actually necessary here to write more than one.
            break


def add_sfbay_freshwater(mdu,
                         flow_locations_shp,
                         grid,dredge_depth,
                         time_offset=None):
    """
    Add freshwater flows from a combination of gaged and ungaged 
    watersheds, with simple scaling between them.
    This is the approach that was used for SUNTANS runs, was replaced
    by BAHM for sfbay_dfm_v2, but is useful for periods outside 
    existing BAHM runs.

    flow_locations_shp: A shapefile with linestring giving each input
    location, 
    fields: 
      gages: A '|' separate listed of USGS gage numbers from which flow data
        will be pulled.
      area_sq_mi: watershed area for this flow.  This area is compared to the
        area in USGS inventory, in order to establish a scaling factor.
      amplify: an additional adjustment to scaling factor.
    
    grid: UnstructuredGrid to add the flows to.  Depths in this grid may be
     "dredged" down to dredge_depth in order to keep inflow locations wet.

    time_offset: pull freshwater flows from this timedelta off from the
    specified.  I.e. if your run is really 2016, but you want 2015 flows,
    specify np.timedelta64(-365,'D').

    Flows are given 0 salinity and 20degC temperature.
    """
    ref_date,run_start,run_stop=mdu.time_range()

    if time_offset is not None:
        run_start = run_start + time_offset
        run_stop = run_stop + time_offset
        ref_date = ref_date + time_offset
    else:
        time_offset=np.timedelta64(0)

    flow_features=wkb2shp.shp2geom(flow_locations_shp)

    # create a pair of bc and pli files, each including all the sources.

    # First need the observations --
    # get a list of all the gages that are referenced:
    all_gages=np.unique( np.concatenate( [gages.split('|') for gages in flow_features['gages']] ) )

    usgs_gage_cache=os.path.join( local_config.cache_path,
                                  'usgs','streamflow')
    flows_ds=usgs_nwis.nwis_dataset_collection(all_gages,
                                               start_date=run_start-5*DAY,
                                               end_date=run_stop+5*DAY,
                                               products=[60], # streamflow
                                               days_per_request='M', # monthly chunks
                                               frequency='daily', # time resolution of the data
                                               cache_dir=usgs_gage_cache)

    usgs_inventory=wkb2shp.shp2geom(usgs_inventory_shp_fn)
    station_to_area=dict( [ ("%d"%site, area)
                            for site,area
                            in zip(usgs_inventory['site_no'],
                                   usgs_inventory['drain_area']) ] )

    unique_names={}
    
    for feat_i,feat in enumerate(flow_features):
        gages=feat['gages'].split('|')
        sub_flows=flows_ds.sel(site=gages)

        featA=feat['area_sq_mi']
        gage_areas=np.array(  [ float(station_to_area[g] or 'nan')
                                for g in gages ] )

        # assume the variable name here, and that dims are [site,time],
        # and units start as cfs.

        # Weighted average of reference gages based on watershed area, and
        # data availability
        # total flow from all reference gages
        site_axis=0
        ref_cms=np.nansum( sub_flows['stream_flow_mean_daily'].values,axis=site_axis ) * FT3_to_M3
        # area represented by reference gages at each time step
        ref_area=np.sum( np.isfinite(sub_flows['stream_flow_mean_daily'].values) * gage_areas[:,None],
                         axis=site_axis )
        
        # avoid division by zero for steps missing all flows
        feat_cms=featA * ref_cms
        feat_cms[ref_area>0] /= ref_area[ ref_area>0 ]
        feat_cms[ref_area==0.0] = np.nan

        stn_ds=xr.Dataset()
        stn_ds['time']=flows_ds.time
        missing=np.isnan(feat_cms)
        if np.all(missing):
            raise Exception("Composite from gages %s has no data for period %s - %s"%(gages,
                                                                                      stn_ds.time.values[0],
                                                                                      stn_ds.time.values[-1]))
        if np.any(missing):
            logging.warning("Composite from gages %s has missing data in period %s - %s"%(gages,
                                                                                          stn_ds.time.values[0],
                                                                                          stn_ds.time.values[-1]))
            # Best guess is period average
            feat_cms[missing]=np.mean(feat_cms[~missing])
        stn_ds['flow_cms']=('time',),feat_cms

        # sanitize and trim the feature name 
        src_name=feat['name'].replace(' ','_').replace(',','_')[:13]
        if src_name in unique_names:
            serial=1
            while True:
                test_name="%s_%dser"%(src_name,serial)
                if test_name not in unique_names:
                    break
                serial+=1
            logging.warning("Source name %s duplicate - will use %s"%(src_name,test_name))
            src_name=test_name

        unique_names[src_name]=src_name

        if 1: #-- Write a PLI file
            pli_feat=(src_name,np.array(feat['geom']))

            # Write copies for flow, salinity and temperatures

            for suffix in ['_flow','_salt','_temp']:
                if suffix=='_temp' and not mdu['physics','Temperature']: # present and not 0.
                    continue
                if suffix=='_salt' and not mdu['physics','Salinity']:
                    continue

                # function to add suffix
                pli_feat_with_suffix=dio.add_suffix_to_feature(pli_feat,suffix)
                pli_fn=os.path.join(mdu.base_path,"%s%s.pli"%(src_name,suffix))
                dio.write_pli(pli_fn,[pli_feat_with_suffix])

            
            dredge_grid.dredge_boundary(grid,pli_feat[1],dredge_depth)

        if 1: #-- Write the time series and stanza in FlowFM_bnd.ext
            write_QST_data(mdu,stn_ds,src_name,time_offset=-time_offset)
