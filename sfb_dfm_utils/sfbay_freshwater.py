import os
import numpy as np
import xarray as xr
import pandas as pd

import stompy.model.delft.io as dio
from stompy import utils, filters
from stompy.io.local import usgs_nwis

from . import dredge_grid, common


DAY=np.timedelta64(86400,'s') # useful for adjusting times

def add_sfbay_freshwater(mdu,
                         rel_bc_dir,     # added directory for bc files alliek dec 2020
                         adjusted_pli_fn,
                         freshwater_dir,
                         grid,dredge_depth,
                         all_flows_unit=False,
                         time_offset=None):
    """
    Add freshwater flows from sfbay_freshwater git submodule.
    run_base_dir: location of DFM input files
    run_start,run_stop: target period for run, as np.datetime64
    ref_date: DFM reference date, as np.datetime64[D]
    adjusted_pli_fn: path to pli file to override source locations
    freshwater_dir: path to sfbay_freshwater git submodule
    grid: UnstructuredGrid instance to be modified at input locations
    old_bc_fn: path to old-style forcing input file

    time_offset: pull freshwater flows from this timedelta off from the
    specified.  I.e. if your run is really 2016, but you want 2015 flows,
    specify np.timedelta64(-365,'D').
    Slightly safer to use days than years here.
    """
    run_base_dir=mdu.base_path
    ref_date,run_start,run_stop = mdu.time_range()
    old_bc_fn=mdu.filepath(["external forcing","ExtForceFile"])
    
    if time_offset is not None:
        run_start = run_start + time_offset
        run_stop = run_stop + time_offset
        ref_date = ref_date + time_offset
        
    def write_flow_data(stn_ds,src_name,flow_scale=1.0):
        df=stn_ds.to_dataframe().reset_index()
        df['elapsed_minutes']=(df.time.values - ref_date)/np.timedelta64(60,'s')
        df['salinity']=0*df.flow_cms
        df['temperature']=20+0*df.flow_cms

        if all_flows_unit:
            df['flow_cms']=1.0+0*df.flow_cms
        else:
            df['flow_cms'] = flow_scale * df.flow_cms

        for quantity,suffix in [ ('dischargebnd','_flow'),
                                 ('salinitybnd','_salt'),
                                 ('temperaturebnd','_temp') ]:
            lines=['QUANTITY=%s'%quantity,
                   'FILENAME=%s/%s%s.pli'%(rel_bc_dir,src_name,suffix), # added rel_bc_dir alliek dec 2020
                   'FILETYPE=9',
                   'METHOD=3',
                   'OPERAND=O',
                   ""]
            with open(old_bc_fn,'at') as fp:
                fp.write("\n".join(lines))

            # read the pli back to know how to name the per-node timeseries
            feats=dio.read_pli(os.path.join(run_base_dir,rel_bc_dir, # added rel_bc_dir alliek dec 2020
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

                tim_fn=os.path.join(run_base_dir,rel_bc_dir,node_name+".tim") # added rel_bc_dir alliek dec 2020

                columns=['elapsed_minutes']
                if quantity=='dischargebnd':
                    columns.append('flow_cms')
                elif quantity=='salinitybnd':
                    columns.append('salinity')
                elif quantity=='temperaturebnd':
                    columns.append('temperature')

                df.to_csv(tim_fn, sep=' ', index=False, header=False, columns=columns)


    adjusted_features=dio.read_pli(adjusted_pli_fn)
    # Add the freshwater flows - could come from erddap, but use github submodule
    # for better control on version

    # create a pair of bc and pli files, each including all the sources.
    # exact placement will
    # be done by hand in the GUI

    full_flows_ds = xr.open_dataset(os.path.join(freshwater_dir, 'outputs', 'sfbay_freshwater.nc'))
    # period of the full dataset which will be include for this run
    sel=(full_flows_ds.time > run_start - 5*DAY) & (full_flows_ds.time < run_stop + 5*DAY)
    flows_ds = full_flows_ds.isel(time=sel)

    nudge_by_gage(flows_ds,'11169025',station='SCLARAVCc',decorr_days=20)
    nudge_by_gage(flows_ds,'11180700',station='UALAMEDA', decorr_days=20)

    if 1: # Special handling for Mowry Slough
        mowry_feat=None
        src_name="MOWRY"
        for adj_feat in adjusted_features:
            if adj_feat[0]==src_name:
                mowry_feat=adj_feat

                # Write copies for flow, salinity and temperatures
                for suffix in ['_flow','_salt','_temp']:
                    # function to add suffix
                    feat_suffix=dio.add_suffix_to_feature(mowry_feat,suffix)
                    pli_fn=os.path.join(run_base_dir,rel_bc_dir,"%s%s.pli"%(src_name,suffix)) # added rel_bc_dir alliek dec 2020
                    dio.write_pli(pli_fn,[feat_suffix])

                dredge_grid.dredge_boundary(grid,mowry_feat[1],dredge_depth)
    
    for stni in range(len(flows_ds.station)):
        stn_ds=flows_ds.isel(station=stni)

        src_name=stn_ds.station.item() # kind of a pain to get scalar values back out...

        # At least through the GUI, pli files must have more than one node.
        # Don't get too big for our britches, just stick a second node 50m east
        # if the incoming data is a point, but check for manually set locations
        # in adjusted_features
        if 1: #-- Write a PLI file
            feat=(src_name,
                  np.array( [[stn_ds.utm_x,stn_ds.utm_y],
                             [stn_ds.utm_x + 50.0,stn_ds.utm_y]] ))
            # Scan adjusted features for a match to use instead
            for adj_feat in adjusted_features:
                if adj_feat[0] == src_name:
                    feat=adj_feat
                    break
            # Write copies for flow, salinity and temperatures
            for suffix in ['_flow','_salt','_temp']:
                # function to add suffix
                feat_suffix=dio.add_suffix_to_feature(feat,suffix)
                pli_fn=os.path.join(run_base_dir,rel_bc_dir,"%s%s.pli"%(src_name,suffix)) # added rel_bc_dir alliek dec 2020
                dio.write_pli(pli_fn,[feat_suffix])

            dredge_grid.dredge_boundary(grid,feat[1],dredge_depth)

        if 1: #-- Write the time series and stanza in FlowFM_bnd.ext
            if src_name=="EBAYS" and mowry_feat is not None:
                write_flow_data(stn_ds,src_name)
                # EBAYS watershed is something like 13000 acres.
                # don't worry about scaling back EBAYS, but add in some extra
                # here for MOWRY
                write_flow_data(stn_ds,"MOWRY",flow_scale=12.8/13000)
            else:
                write_flow_data(stn_ds,src_name)
    
    full_flows_ds.close()
    
##

# Override BAHM with gage data when available
# what are possible overrides?
#  - COYOTE => 11172175 COYOTE C AB HWY 237 A MILPITAS CA
#  - SCLARAVCc => 11169025 GUADALUPE R ABV HWY 101 A SAN JOSE CA
#  - UALEMADAg => 11180700 ALAMEDA C FLOOD CHANNEL A UNION CITY CA
#  - USANLORZ => 11181040 SAN LORENZO C A SAN LORENZO CA
#  - EBAY Cc4 => 374336122095801 SAN LEANDRO C A ALVARADO ST A SAN LEANDRO CA
#  - MARINS3 => 11460000 CORTE MADERA C A ROSS CA
#  - MARINN => 11459500 NOVATO C A NOVATO CA
#  - PETALUMA => 381519122385601 PETALUMA R NR PETALUMA CA
#  - NAPA => 11458000 NAPA R NR NAPA CA


def nudge_by_gage(ds,usgs_station,station,decorr_days,period_start=None,period_end=None):
    # This slicing may be stopping one sample shy, shouldn't be a problem.
    if period_start is None:
        period_start=ds.time.values[0]
    if period_end is None:
        period_end=ds.time.values[-1]
        
    usgs_gage=usgs_nwis.nwis_dataset(usgs_station,
                                     products=[60],
                                     start_date=period_start,
                                     end_date=period_end,
                                     days_per_request='M',
                                     cache_dir=common.cache_dir)

    # Downsample to daily
    df=usgs_gage['stream_flow_mean_daily'].to_dataframe()
    df_daily=df.resample('D').mean()

    # Get the subset of BAHM data which overlaps this gage data
    time_slc=slice(np.searchsorted( ds.time, df_daily.index.values[0]),
                   1+np.searchsorted( ds.time, df_daily.index.values[-1] ) )

    bahm_subset=ds.sel(station=station).isel( time=time_slc )

    assert len(bahm_subset.time) == len(df_daily),"Maybe BAHM data doesn't cover entire period"

    errors=bahm_subset.flow_cfs - df_daily.stream_flow_mean_daily 

    # Easiest: interpolate errors over nans, apply to bahm data array.
    # the decorrelation time is tricky, though.

    # Specify a decorrelation time scale then relax from error to zero
    # over that period
    valid=np.isfinite(errors.values)
    errors_interp=np.interp( utils.to_dnum(ds.time),
                             utils.to_dnum(df_daily.index[valid]),
                             errors[valid])
    all_valid=np.zeros( len(ds.time),'f8' )
    all_valid[time_slc]=1*valid

    weights=(2*filters.lowpass_fir(all_valid,decorr_days)).clip(0,1)
    weighted_errors=weights*errors_interp

    # Does this work? 
    subset=dict(station=station)

    cfs_vals=ds.flow_cfs.loc[subset] - weighted_errors
    ds.flow_cfs.loc[subset] = cfs_vals.clip(0,np.inf) 
    ds.flow_cms.loc[subset] = 0.028316847 * ds.flow_cfs.loc[subset]

    # user feedback
    cfs_shifts=weighted_errors[time_slc]
    print("Nudge: %s => %s, shift in CFS: %.2f +- %.2f"%(usgs_station,station,
                                                         np.mean(cfs_shifts),
                                                         np.std(cfs_shifts)))

