"""
Utilities to download CA ROMS data to local cache, extract grids, 
and extract BC data.
"""

import time
import os
import glob
import logging

logger=logging.getLogger('ca_roms')

import numpy as np
import xarray as xr

from . import local_config

from stompy import utils, memoize

from stompy.grid import unstructured_grid
from stompy.spatial import wkb2shp
from stompy.plot import plot_utils, plot_wkb
from stompy.spatial import proj_utils, wkb2shp,field, linestring_utils
from stompy.model.delft import dfm_grid
import stompy.model.delft.io as dio

from shapely import geometry
from shapely.ops import cascaded_union


cache_path = os.path.join(local_config.cache_path,'ca_roms')

utm2ll=proj_utils.mapper('EPSG:26910','WGS84')
ll2utm=proj_utils.mapper('WGS84','EPSG:26910')


def fetch_ca_roms(start,stop):
    """
    Download the 6-hourly outputs for the given period.
    start,stop: np.datetime64 
    returns a list of paths to local files falling in that period
    
    When a file cannot be downloaded, but falls in the period, the previous
    valid file is written out with a fake name, and a variable "original_filename"
    indicates the source of the copy.  This currently does *not* handle the first
    file being invalid.

    There are a few cases of a file having no data -- this method checks for
    all of zeta being nan, in which case the file is treated as missing.

    Pads the dates out by 12-36 h, as the data files are every 6h, staggered,
    and the stop date gets truncated 1 day
    """

    start=start-np.timedelta64(12,'h')
    stop=stop+np.timedelta64(36,'h')
    
    # Ways of getting the list of files:
    # http://west.rssoffice.com:8080/thredds/catalog/roms/CA3km-nowcast/CA/catalog.xml
    # That provides xml with file names, sizes, but no metadata about period simulated.

    assert os.path.exists(cache_path)
    
    local_files=[]

    last_valid_fn=None
    
    # As a first step, construct the url by hand:
    for day in np.arange(start,stop,np.timedelta64(1,'D')):
        logger.info(day)
        for hour in [3,9,15,21]: # 6 hourly output staggered 3 hours
            ymd=utils.to_datetime(day).strftime("%Y%m%d")
            url="http://west.rssoffice.com:8080/thredds/dodsC/roms/CA3km-nowcast/CA/ca_subCA_das_%s%02d.nc"%(ymd,hour)
            base_fn=url.split('/')[-1]
            local_fn=os.path.join(cache_path,base_fn)

            if os.path.exists(local_fn):
                logger.info("%30s: exists"%base_fn)
                last_valid_fn=local_fn
            else:
                logger.info("%30s: fetching"%base_fn)
                try:
                    ds_fetch=xr.open_dataset(url)
                except OSError as exc:
                    logger.warning("%27s  FAILED"%'')
                    ds_fetch=None

                # Check for an invalid file:
                if np.all( np.isnan( ds_fetch.zeta.isel(time=0) ) ):
                    logger.warning("%27s has no wet cells?!"%"")
                    ds_fetch=None # Treat as missing

                if ds_fetch is None:
                    if last_valid_fn is None:
                        print(" -- last valid is None, can't rescue, will omit")
                        # missing files at start of period -- omit
                        # from the returned list
                        continue
                    else:
                        ds_fetch=xr.open_dataset(last_valid_fn)
                        if 'original_filename' not in ds_fetch.attrs:
                            ds_fetch.attrs['original_filename']=last_valid_fn
                        # And overwrite title, which is the only reliable timestamp
                        # in these files.
                        ds_fetch.attrs['title']='CA-%s%02d'%(ymd,hour)
                    
                if ds_fetch is not None:
                    # in the future, we could subset at this point for a possible speedup.
                    logger.info("%30s  saving"%"")
                    ds_fetch.to_netcdf(local_fn)
                    ds_fetch.close()
                    if 'original_filename' not in ds_fetch.attrs:
                        last_valid_fn=local_fn
                    
                logger.info("%30s  sleeping"%"")
                time.sleep(5) # be nice!
            local_files.append(local_fn)
    return local_files

   
# Choose a subset of that to make into the DFM domain
# picked off a map

def extract_roms_subgrid(ul_xy=(358815., 4327282.),
                         lr_xy=(633421., 4009624.)):
    """
    Loads CA ROMS 4km output
    selects a subset of that grid based on the UTM extents specified
    creates a curvilinear unstructured grid matching that subset of the ROMS
    grid
    projects to UTM, trims dry cells (based on ROMS data).
    cleans up the numbering and returns that grid, without bathymetry

    cells['lati'] and ['loni'] had been indexes to just the subset of the ROMS 
    grid.  Those are now indexes into the original roms grid.
    """

    snap_fns=glob.glob(os.path.join(cache_path,'*.nc'))
    
    snap_fns.sort()

    ds=xr.open_dataset(snap_fns[0])

    ds0=ds.isel(time=0)

    ul_ll=utm2ll(ul_xy)
    lr_ll=utm2ll(lr_xy)

    lon_range=[ul_ll[0], lr_ll[0]]
    lat_range=[lr_ll[1], ul_ll[1]]

    lat_sel=(ds0.lat.values>=lat_range[0]) & (ds0.lat.values<=lat_range[1])
    lon_sel=(ds0.lon.values>=(lon_range[0]%360)) & (ds0.lon.values<=(lon_range[1]%360))
    # starting indices 
    lat0i=np.nonzero(lat_sel)[0][0]
    lon0i=np.nonzero(lon_sel)[0][0]
    ds0_sub=ds0.isel( lat=lat_sel, lon=lon_sel )

    Lat,Lon=np.meshgrid( ds0_sub.lat.values, ds0_sub.lon.values)

    ll_mesh=np.array( [Lon,Lat] ).transpose(1,2,0)

    xy_mesh=ll2utm(ll_mesh)

    # Extract a grid:

    g=unstructured_grid.UnstructuredGrid()

    # pretty sure that the specified lon/lat are for cell centers (rho points)
    node_lat=utils.center_to_edge(ds0_sub.lat.values)
    node_lon=utils.center_to_edge(ds0_sub.lon.values)

    # Maybe 20s
    mappings = g.add_rectilinear(p0=[node_lon[0],node_lat[0]],
                                 p1=[node_lon[-1],node_lat[-1]],
                                 nx=len(node_lon),
                                 ny=len(node_lat))
                                
    # Remove dry cells and store indices to get back to lat/lon grid
    dry=np.isnan(ds0_sub.zeta.values)

    cell_values=np.zeros(g.Ncells())
    cell_lati=np.zeros(g.Ncells(),'i4')
    cell_loni=np.zeros(g.Ncells(),'i4')

    for lati,lat in enumerate(ds0_sub.lat.values):
        for loni,lon in enumerate(ds0_sub.lon.values):
            celli=mappings['cells'][loni,lati]
            # These are now referenced to the full ROMS grid
            cell_lati[celli]=lati + lat0i
            cell_loni[celli]=loni + lon0i
            if dry[lati,loni]:
                cell_values[ celli ] = 1.0

    g.add_cell_field('lati',cell_lati)
    g.add_cell_field('loni',cell_loni)

    # Reproject to UTM zone 10 meters
    g.nodes['x']=proj_utils.mapper('WGS84','EPSG:26910')(g.nodes['x'])

    # trim the dry cells:
    for c in np.nonzero(cell_values)[0]:
        g.delete_cell(c)

    # Clean up the nodes and edges related to deleted cells,
    # get everything consecutively numbered
    g.renumber_cells()
    g.make_edges_from_cells()
    g.delete_orphan_nodes()
    g.renumber()
    return g

def extract_roms_subgrid_poly(poly):
    """
    Extract a ROMS subgrid based on a polygon (in UTM coordinates)
    """
    # Start with an oversized, rectangular grid:
    g=extract_roms_subgrid( ul_xy=(poly.bounds[0],poly.bounds[3]),
                            lr_xy=(poly.bounds[2],poly.bounds[1]))
    # Then trim the fat:
    to_delete=~g.select_cells_intersecting(poly)
    for c in np.nonzero(to_delete)[0]:
        g.delete_cell(c)

    # And clean up
    g.renumber_cells()
    g.make_edges_from_cells()
    g.delete_orphan_nodes()
    g.renumber()
    return g

# 6k cells.  Not bad.
@memoize.memoize()
def coastal_dem():
    # Note - this is not the CA ROMS cache path, since it's not referencing
    # a CA ROMS dataset, but a more general bathy dataset
    fn=os.path.join(local_config.cache_path,"ngdc-etopo-merge-utm.tif")

    assert os.path.exists(fn)

    # The recipe for making that file
    # if not os.path.exists(fn):
    #     crm_ll=field.GdalGrid('ngdc-crm.tif')
    #     missing=(crm_ll.F==-32768)
    #     crm_ll.F=crm_ll.F.astype('f8')
    # 
    #     dem_etopo=field.GdalGrid('etopo1.tif')
    # 
    #     X,Y = crm_ll.XY()
    #     bad_x=X[missing]
    #     bad_y=Y[missing]
    # 
    #     fill=dem_etopo( np.c_[bad_x,bad_y] )
    #     crm_ll.F[missing] = fill
    # 
    #     dem=crm_ll.warp(t_srs='EPSG:26910',
    #                     s_srs='WGS84',
    #                     fn=fn)
    # else:
    #     dem=field.GdalGrid(fn)
    return field.GdalGrid(fn)
    

def add_coastal_bathy(g,dem=None):
    dem=dem or coastal_dem()
    
    # Add depth to the grid:
    def clip(x):
        return x.clip(-np.inf,-4)
    
    node_depth=clip( dem( g.nodes['x'] ) )
    cell_depth=clip( dem( g.cells_centroid() ) )
    edge_depth=clip( dem( g.edges_center() ) )

    # the redundant names are so it can get written to ugrid
    # where all variables have a single namespace, and when it
    # comes back in we will still have node depth called 'depth'
    g.add_node_field('depth',node_depth)
    g.add_edge_field('edge_depth',edge_depth)
    g.add_cell_field('cell_depth',cell_depth)

    return g

def infer_variable(ds,canon):
    """
    ds: xr.Dataset
    canon: 'zeta','salt','temp','u','v'
    """
    # hacky - but different models have different names
    # HYCOM at least has a good standard name for eta:
    #     standard_name:  sea_surface_elevation
    # but CA ROMS just has a long name 
    #     long_name:  Sea Surface Height
    # For now, hard code...
    if canon in ds:
        return canon 

    synonyms={'zeta':['surf_el'],
              'salt':['salinity'],
              'temp':['temperature','water_temp'],
              'u':['water_u'],
              'v':['water_v']}

    for syn in synonyms[canon]:
        if syn in ds:
            return syn
    raise None

def annotate_grid_from_data(g,coastal_files,candidate_edges=None,check_wet=False):
    """
    Add src_idx_in,src_idx_out fields to edges for edges which are 
    adjacent to active cells in the lat/lon rectilinear grid inputs 
    given in coastal_files.

    g: unstructured_grid to be annotated
    coastal_files: list of file paths, giving the snapshots of a coastal
      model (assumed to have a rectilinear lat/lon grid)
    candidate_edges: if specified, only consider these edges, an array
      of edge indices
    check_wet: if true scan the files to limit the selection of boundary edges
      to cells which are wet at all time steps
    """
    # Scan all of the ROMS files to find cells which are always wet
    wet=True

    ds=xr.open_dataset(coastal_files[0])
    eta_var=infer_variable(ds,'zeta')
    wet=np.ones( ds[eta_var].values.shape, np.bool )

    for coastal_file in coastal_files:
        logging.info(coastal_file)
        ds=xr.open_dataset(coastal_file)
        eta_values=ds[eta_var]
        if 'time' in eta_values.dims:
            eta_values=eta_values.isel(time=0)
        one_step=np.isfinite( eta_values.values )
        assert not np.all(~one_step), "Coastal model file %s has no wet cells?!"%coastal_file

        wet = wet & one_step
        ds.close()
        if not check_wet:
            logging.info("Will assume wet-cells in first time step true for eternity")
            break

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

    src=xr.open_dataset(coastal_files[0])

    if candidate_edges is None:
        # at least limit to boundary edges
        candidate_edges=np.nonzero( g.edges['cells'].min(axis=1)<0 )[0]
        # old code: arange(g.Nedges())
    
    for j in candidate_edges:
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

        lon_idx_out=utils.nearest(src.lon.values%360.0,cout_cc_ll[0]%360.0)
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
    dest_fields=[]
    roms_fields=[]
    if int(mdu['physics','salinity']):
        dest_fields.append('sa1')
        roms_fields.append('salt')
    if int(mdu['physics','temperature']):
        dest_fields.append('tem1')
        roms_fields.append('temp')

    normalize_variables(snap)

    map_in=xr.open_dataset(map_file)

    map_out=map_in.isel(time=[0])

    # there are some name clashes -- drop any coordinates attributes
    #
    for dv in map_out.data_vars:
        if 'coordinates' in map_out[dv].attrs:
            del map_out[dv].attrs['coordinates']

    if len(dest_fields)==0:
        return map_out

    # DFM reorders the cells, so read it back in.
    g_map=dfm_grid.DFMGrid(map_out)
    g_map_cc=g_map.cells_centroid()
    g_map_ll=utm2ll(g_map_cc)

    ## 
    dlon=np.median(np.diff(snap.lon.values))
    dlat=np.median(np.diff(snap.lat.values))

    roms_z=-snap.depth.values[::-1]

    if snap.time.ndim>0:
        snap0=snap.isel(time=0)
    else:
        snap0=snap

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

# Not exactly ROMS-specific, but helps with using ROMS coupling
def add_sponge_layer(mdu,run_base_dir,grid,edges,sponge_visc,background_visc,sponge_L,
                     quantity='viscosity'):
    """
    quantity: 'viscosity' or 'diffusivity'
    """
    obc_centers=grid.edges_center()[edges]

    # sponge length scale: 25000 [m] is roughly 8 cells

    sample_sets=[ np.c_[obc_centers[:,0],obc_centers[:,1],sponge_visc*np.ones(len(obc_centers))] ]

    circs=[geometry.Point(xy).buffer(sponge_L)
           for xy in obc_centers]
    obc_buff=cascaded_union(circs).boundary
    obc_buff_pnts=np.array(obc_buff)
    obc_buff_pnts_resamp=linestring_utils.downsample_linearring(obc_buff_pnts,sponge_L*0.5)

    sample_sets.append( np.c_[obc_buff_pnts_resamp[:,0],
                              obc_buff_pnts_resamp[:,1],
                              background_visc*np.ones(len(obc_buff_pnts_resamp))] )

    # And some far flung values
    x0,x1,y0,y1=grid.bounds()
    corners=np.array( [[x0-sponge_L,y0-sponge_L],
                       [x0-sponge_L,y1+sponge_L],
                       [x1+sponge_L,y1+sponge_L],
                       [x1+sponge_L,y0-sponge_L]] )

    sample_sets.append( np.c_[corners[:,0],
                              corners[:,1],
                              background_visc*np.ones(len(corners))] )

    visc_samples=np.concatenate(sample_sets,axis=0)

    np.savetxt(os.path.join(run_base_dir,'%s.xyz'%quantity),
               visc_samples)

    txt="\n".join(["QUANTITY=horizontaleddy%scoefficient"%quantity,
                   "FILENAME=%s.xyz"%quantity,
                   "FILETYPE=7",
                   "METHOD=4",
                   "OPERAND=O"
                   "\n"])
    
    old_bc_fn = os.path.join(run_base_dir,mdu['external forcing','ExtForceFile'])
    with open(old_bc_fn,'at') as fp:
        fp.write(txt)
        

def extract_data_at_boundary(coastal_files,g,boundary_edges):
    """
    coastal_files: list of file paths for lat/lon rectilinear coastal ocean model
      output (one file per time step)
    g: unstructured_grid
    boundary_edges: list/array of edge indices on g for which to extract 
      data.

    returns an xarray Dataset with cell values pulled out for each boundary edge,
      over the times represented by coastal_files.

    Also normalizes some variable names:
     water surface elevation: zeta
    """
    # Pre-extract some fields from the ROMS data
    # Might move this to ca_roms.py in the future
    # This might be getting the wrong locations.
    # This is pretty slow
    extracted=[]
    lat_da=xr.DataArray(g.edges['src_idx_out'][boundary_edges,0],dims='boundary')
    lon_da=xr.DataArray(g.edges['src_idx_out'][boundary_edges,1],dims='boundary')

    for coastal_file in coastal_files:
        ds=xr.open_dataset(coastal_file)
        ds.load()

        # Total hack - 
        # ROMS snapshots have a bogus time stamp in them.
        if 'title' in ds: # probably a ROMS file, fix the timestamp
            # timestamps appear to be wrong in the files, always
            # holding 2009-01-02.
            # use the title, which appears to be consistent with the filename
            # and the true time
            t=utils.to_dt64( datetime.datetime.strptime(ds.title,'CA-%Y%m%d%H') )
            ds.time.values[0]=t
            
        if ds.time.ndim>0: 
            # ROMS files have a 1-element time vector
            sub_ds=ds.isel(time=0)
        else:
            # HYCOM files have a scalar time value
            sub_ds=ds
        sub_ds=sub_ds.isel(lat=lat_da,lon=lon_da)
        # The above has gotten warnings along the lines of:
        # .../numeric.py:1466: VisibleDeprecationWarning: converting an array with
        #  ndim > 0 to an index will result in an error in the future
        # Yet these still check out - 
        assert np.allclose(ds.lon[lon_da].values,sub_ds.lon.values)
        assert np.allclose(ds.lat[lat_da].values,sub_ds.lat.values)

        extracted.append(sub_ds)
        ds.close()

    data_at_boundary=xr.concat(extracted,dim='time')

    normalize_variables(data_at_boundary)

    return data_at_boundary


def normalize_variables(ds):
    for v in ['zeta','salt','temp','u','v']:
        v_inf=infer_variable(ds,v)
        if v_inf is None:
            continue
        if v_inf!=v:
            # make a shallow copy with the canonical name
            ds[v]=ds[v_inf]

