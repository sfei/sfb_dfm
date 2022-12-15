import os
import time
import logging

import pandas as pd
import xarray as xr

import numpy as np
from . import local_config


log=logging.getLogger('hycom')

def hycom_opendap_url_for_time(t):
    """
    Different periods have different opendap urls.
    t: np.datetime64
    """
    if t>=np.datetime64("2017-10-01"):
        return 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_92.9'
    if t>=np.datetime64("2014-07-01"):
        return 'http://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_57.7'
    raise Exception("Not ready for other times")

datasets={} # url => dataset
def hycom_opendap_for_time(t):
    url=hycom_opendap_url_for_time(t)
    if url not in datasets:
        ds=xr.open_dataset(url,decode_times=False)
        del ds['tau']
        ds2=xr.decode_cf(ds)
        datasets[url]=ds2
    return datasets[url]

def fetch_range( lon_range,lat_range,time_range ):
    """
    lon_range: [lon_min,lon_max]
    lat_range: [lat_min,lat_max]
    time_range: [time_min,time_max]
    returns a list of local netcdf filenames, one per time step

    Limitations:
     * fetch daily snapshots, though source has every 3 hours
     * not ready for time ranges that span multiple HYCOM experiments.
    """
    times=pd.DatetimeIndex(start=time_range[0],end=time_range[1],freq='D')

    cache_dir=os.path.join( local_config.cache_path, 'hycom')
    os.path.exists(cache_dir) or os.mkdir(cache_dir)

    last_ds=None
    lon_slice=None
    lat_slice=None

    lon_range=np.asarray(lon_range)

    filenames=[]

    for t in times:
        time_str = t.strftime('%Y%m%d%H')
        cache_name=os.path.join( cache_dir,
                                 "%s-%.2f_%.2f_%.2f_%.2f.nc"%(time_str,
                                                              lon_range[0],lon_range[1],
                                                              lat_range[0],lat_range[1]) )
        if not os.path.exists(cache_name):
            log.info("Fetching %s"%cache_name)
            ds=hycom_opendap_for_time(t)
            if ds is not last_ds:
                # Careful to deal with +-360 possibility
                lon0=ds.lon.values[0]
                lon_slice=slice(*np.searchsorted((ds.lon.values-lon0)%360.0,(lon_range-lon0)%360.0))
                lat_slice=slice(*np.searchsorted(ds.lat.values,lat_range))
                last_ds=ds

            time_i=np.searchsorted(ds.time.values,t.to_datetime64())
            ds_region=ds.isel(lat=lat_slice,lon=lon_slice,time=time_i)
            ds_region.to_netcdf(cache_name)
            time.sleep(1.0)
        filenames.append(cache_name)
    return filenames


def annotate_grid_from_data(g,start,stop,candidate_edges=None):
    """
    Add src_idx_in,src_idx_out fields to edges for ROMS-adjacent boundary 
    edges.

    g: unstructured_grid to be annotated
    start,stop: datetime64 date range, for selecting wet cells from ROMS
    candidate_edges: if specified, only consider these edges, an array
      of edge indices
    """
    # Get the list of files
    ca_roms_files=fetch_ca_roms(start,stop)

    # Scan all of the ROMS files to find cells which are always wet
    wet=True
    for ca_roms_file in ca_roms_files:
        logging.info(ca_roms_file)
        ds=xr.open_dataset(ca_roms_file)
        one_step=np.isfinite( ds.zeta.isel(time=0).values )
        assert not np.all(~one_step), "ROMS file %s has no wet cells?!"%ca_roms_file

        wet = wet & one_step
        ds.close()

    boundary_cells=[]
    boundary_edges=[]

    # record indices into src for the src cells just inside and
    # outside each boundary edge, or -1 for non-boundary
    edge_src_idx_in=np.zeros( (g.Nedges(),2),'i4') - 1
    edge_src_idx_out=np.zeros( (g.Nedges(),2),'i4') - 1
    edge_norm_in=np.zeros( (g.Nedges(),2),'f8')*np.nan

    N=g.Nedges()
    g.edge_to_cells() # Make sure that's been populated

    centroid=g.cells_centroid()
    edge_ctr=g.edges_center()

    src=xr.open_dataset(ca_roms_files[0])

    if candidate_edges is None:
        candidate_edges=np.arange(g.Nedges())
    
    for j in candidate_edges: # range(g.Nedges()):
        if j%1000==0:
            logger.info("%d/%d"%(j,N))
        c1,c2=g.edges['cells'][j]
        if (c1>=0) and (c2>=0):
            continue # not a boundary in this grid
        if c1<0:
            cin=c2
        else:
            cin=c1
        # A bit goofy, but make this a bit more general...
        # Construct the center of the cell that would be our neighbor:
        cout_cc=2*edge_ctr[j] - centroid[cin]
        cout_cc_ll=utm2ll(cout_cc)

        lon_idx_out=utils.nearest(src.lon.values,cout_cc_ll[0]%360.0)
        lat_idx_out=utils.nearest(src.lat.values,cout_cc_ll[1])

        if not wet[lat_idx_out,lon_idx_out]:
            continue

        # good - it's a real boundary edge
        cin_cc_ll=utm2ll(centroid[cin])
        lon_idx_in=utils.nearest(src.lon.values,cout_cc_ll[0]%360.0)
        lat_idx_in=utils.nearest(src.lat.values,cout_cc_ll[1])

        boundary_edges.append(j)
        boundary_cells.append(cin)
        edge_src_idx_out[j]=(lat_idx_out,lon_idx_out) # assumes an order in src
        edge_src_idx_in[j]=(lat_idx_in,lon_idx_in) # assumes an order in src

        edge_norm_in[j,:] = utils.to_unit( centroid[cin] - cout_cc )

    g.add_edge_field('src_idx_out',edge_src_idx_out,on_exists='overwrite')
    g.add_edge_field('src_idx_in',edge_src_idx_in,on_exists='overwrite')
    g.add_edge_field('bc_norm_in',edge_norm_in,on_exists='overwrite')
    



class DFMZModel(object):
    def __init__(self,map_out,mdu):
        self.map_out=map_out
        self.mdu=mdu
        zvar=self.map_out.LayCoord_cc # could be more dynamic about this
        std_name=zvar.attrs.get('standard_name','')
        long_name=zvar.attrs.get('long_name','')
        
        if std_name=="ocean_sigma_coordinate":
            self.coord_type='sigma'
        elif std_name=="ocean_zlevel_coordinate":
            self.coord_type='zlevel'
        elif long_name.startswith('sigma layer'):
            self.coord_type='sigma'
        elif long_name.startswith('z layer'):
            self.coord_type='zlevel'
        else:
            raise Exception("Cannot decipher type of vertical coordinate system")

        self.all_bl=map_out.FlowElem_bl.values
        self.all_wd=map_out.waterdepth.isel(time=0).values
        
        if self.coord_type=='zlevel':
            print("Detected z layers!")
            self.max_depth=self.all_bl.min()
            self.n_layers=len(zvar)
            if 0: # uniform
                bounds=np.linspace(self.max_depth,0,self.n_layers+1)
            else: # exponential
                bounds=dio.exp_z_layers(mdu)
                print("Bounds are %s"%bounds)
            middles=0.5*(bounds[:-1] + bounds[1:])
            self.z_layers=middles
        elif self.coord_type=='sigma':
            self.sigma=zvar.values
            # Not ready for anything spatially variable
            assert self.sigma.ndim==1

            
    def get_dfm_z(self,c):
        if self.coord_type=='sigma':
            # This way is nice but very slow:
            # bl=map_out.FlowElem_bl.isel(nFlowElem=c).values # positive up from the z datum
            # wd=map_out.waterdepth.isel(nFlowElem=c,time=0).values
            # sigma=map_out.LayCoord_cc.values
            bl=self.all_bl[c]
            wd=self.all_wd[c]

            return self.bl + self.wd*self.sigma
        else:
            return self.z_layers
    
def set_ic_from_map_output(snap,map_file,mdu,output_fn='initial_conditions_map.nc',
                           missing=0,tol_km=10):
    """
    copy the structure of the map_file at one time step, but overwrite
    fields with ROMS snapshot data

    snap: A ROMS output file, loaded as xr.Dataset
    map_file: path to map output, must be single processor run.
    mdu: An MDUFile object, for discerning vertical coordinates
    output_fn: If specified, the updated map file is written to the given path.
    
    missing: for locations in the map file which are missing or unmatched in the ROMS
    file, set scalars to this value.  If set to None, leave those water columns
    as is in the map file.

    tol_km: if the best ROMS match is farther away than this, set it to missing

    returns a Dataset() of the updated map
    """
    dest_fields=['sa1','tem1'] # names of the fields to write to in the map file
    roms_fields=['salt','temp'] # source fields in the ROMS data

    map_in=xr.open_dataset(map_file)

    map_out=map_in.isel(time=[0])
    # there are some name clashes -- drop any coordinates attributes
    #
    for dv in map_out.data_vars:
        if 'coordinates' in map_out[dv].attrs:
            del map_out[dv].attrs['coordinates']

    # DFM reorders the cells, so read it back in.
    g_map=dfm_grid.DFMGrid(map_out)
    g_map_cc=g_map.cells_centroid()
    g_map_ll=utm2ll(g_map_cc)

    ## 
    dlon=np.median(np.diff(snap.lon.values))
    dlat=np.median(np.diff(snap.lat.values))

    roms_z=-snap.depth.values[::-1]

    snap0=snap.isel(time=0)
    # make the dimension order explicit so we can safely index
    # this via numpy below
    wet=np.isfinite(snap0.zeta.transpose('lat','lon').values)
    snap_scalars=[snap0[roms_field].transpose('lat','lon','depth').values
                  for roms_field in roms_fields]
    
    snap_lon=snap0.lon.values
    snap_lat=snap0.lat.values

    sel_time=xr.DataArray([0],dims=['time'])
    sel_cell=xr.DataArray([1000000],dims=['nFlowElem'])

    dfm_z_model=DFMZModel(map_out,mdu)

    # much faster to write straight to numpy array rather than
    # via xarray.  But for a bit more robustness, hack together
    # dynamic indexing

    dest_arrays=[map_out[fld].values
                 for fld in dest_fields]
    # dest_salt=map_out.sa1.values
    dest_idx=[]
    for dimi,dim in enumerate(map_out.sa1.dims):
        if dim=='time':
            dest_idx.append(0)
        elif dim=='laydim':
            dest_idx.append(slice(None))
        else:
            dest_idx_cell=dimi
            dest_idx.append(100000)
            
    for c in g_map.valid_cell_iter():
        if c%1000==0:
            print("%d/%d"%(c,g_map.Ncells()))
        is_missing=False
        
        loni=utils.nearest(snap_lon%360,g_map_ll[c,0]%360)
        lati=utils.nearest(snap_lat,g_map_ll[c,1])
        
        err_km=utils.haversine([snap_lon[loni],snap_lat[lati]],
                               g_map_ll[c,:])
        
        if err_km>tol_km:
            is_missing=True
        elif not bool(wet[lati,loni]):
            is_missing=True
        else:
            # grab a water column, and flip it vertically to be
            # in order of increasing, positive-up, depth, i.e.
            # bed to surface.
            # xarray = slow
            #   roms_salt=snap0.salt.isel(lat=lati,lon=loni).values[::-1]
            # numpy = fast
            roms_scalars=[ snap_scalar[lati,loni,::-1]
                           for snap_scalar in snap_scalars]
            # roms_salt=snap_salt[lati,loni,::-1]

            #valid=np.isfinite(roms_salt)
            # Assume that if the first (salt) is valid, then so are the rest (temperature)
            valid=np.isfinite(roms_scalars[0])
            if not np.any(valid):
                is_missing=True 
            else:
                dfm_z=dfm_z_model.get_dfm_z(c)

                dfm_scalars=[ np.interp(dfm_z,roms_z[valid],roms_scalar[valid])
                              for roms_scalar in roms_scalars ]
                #dfm_salt=np.interp(dfm_z,
                #                   roms_z[valid],roms_salt[valid])

        dest_idx[dest_idx_cell]=c
        if is_missing:
            if missing is None:
                continue
            else:
                for dest in dest_arrays:
                    dest[dest_idx]=missing
                # dest_salt[dest_idx]=missing
        else:
            #dest_salt[dest_idx]=dfm_salt
            for dest,dfm_scalar in zip(dest_arrays,dfm_scalars):
                dest[dest_idx]=dfm_scalar

    #map_out.sa1.values[:,:,:]=dest_salt # a little dicey
    for dest_field,dest in zip(dest_fields,dest_arrays):
        map_out[dest_field].values[:,:,:]=dest # a little dicey        
    # that was legal, and took, right? check the first one
    assert np.allclose( map_out.sa1.values, dest_arrays[0] )

    # unorm isn't that useful, and in z-layers it stalls the whole show.
    # DFM writes a bunch of warnings about vicwwu
    for bad_field in ['unorm','vicwwu']:
        if bad_field in map_out:
            del map_out[bad_field]
    
    # Does the map timestamp have to match what's in the mdu?
    # Yes.
    if output_fn is not None:
        map_out.to_netcdf(output_fn,format='NETCDF3_64BIT')
    return map_out
