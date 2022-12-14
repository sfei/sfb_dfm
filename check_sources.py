"""
Debugging the sources -- early checks do not make it clear whether the sources are
properly included or not.
"""
import os

import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import collections
import numpy as np
from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid

##
run_name="short_20120801_p16"
run_base_dir="runs/" +run_name
output_dir=os.path.join(run_base_dir,"DFM_OUTPUT_%s"%run_name)

Nproc=16
procs=range(Nproc)

all_ds=[xr.open_dataset(os.path.join(output_dir,'%s_%04d_20120801_000000_map.nc'%(run_name,proc)))
        for proc in procs]

all_g=[unstructured_grid.UnstructuredGrid.from_ugrid(ds)
       for ds in all_ds]

# Global grid so we can pull out clusters of mass
G=dfm_grid.DFMGrid('sfei_v19_net.nc')

##

# Walk through the subdomains, figure out which cell in G each
# subdomain cell belongs to.  Could be faster, but as is takes
# about 30s.
global_elems=np.zeros(g.Ncells(),'int32')-1

for proc in procs:
    print(proc)
    g=all_g[proc]
    centers=g.cells_centroid()
    for c in range(g.Ncells()):
        Gcell=G.select_cells_nearest( centers[c],inside=True )
        assert Gcell>=0
        global_elems[Gcell]=all_ds[proc].FlowElemGlobalNr.values[c]
G.add_cell_field('global_elem',global_elems)
# and establish the mapping from global element to cell index in G
# fill with invalid first
global_to_G=np.zeros(G.cells['global_elem'].max()+1,'int32') + 10*g.Ncells()
# then update valid
global_to_G[global_elems]=np.arange(G.Ncells())
##

fresh=np.zeros(G.Ncells(),'float64')

for proc in procs:
    # steps 3 and 4 should be clear of the 3h ramp up.
    salt0=all_ds[proc].sa1.isel(time=0,laydim=9).values
    saltN=all_ds[proc].sa1.isel(time=4,laydim=9).values
    dsalt=saltN-salt0
    fresh[global_to_G[all_ds[proc].FlowElemGlobalNr.values]]=dsalt

## 
plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

mask=np.abs(fresh)>0.001
coll=G.plot_cells(ax=ax,lw=0.5,values=fresh,mask=mask)
coll.set_clim([-.1,0])

e2c=G.edge_to_cells()
G.plot_edges(ax=ax,lw=1.0,color='k',mask=e2c.min(axis=1)<0)

plt.setp(coll,edgecolor='face')
plt.colorbar(coll)

import glob
for pli_fn in glob.glob(os.path.join(run_base_dir,'*.pli')):
    if '_salt' in pli_fn or '_temp' in pli_fn: continue
    feat=dio.read_pli(pli_fn)[0]
    name=feat[0].replace('_flow','')
    ax.text(feat[1][0,0],feat[1][0,1],name)
    ax.plot(feat[1][:,0],feat[1][:,1],'m-',lw=2)

##

# use [3,4] for fully-spun up boundary conditions, but earlier
# periods allow better distinction between sources.
times=[3,4] # 3600s per step

dfresh_vols=np.zeros(G.Ncells(),'float64')

for proc in procs:
    salt=all_ds[proc].sa1.isel(time=times).mean(dim='laydim').values
    fresh=(35-salt)/35.
    areas=all_g[proc].cells_area()
    depths=all_ds[proc].waterdepth.isel(time=times).values

    fresh_vol=(fresh*areas*depths)
    dfresh_vol=fresh_vol[1] - fresh_vol[0]
    dfresh_vols[global_to_G[all_ds[proc].FlowElemGlobalNr.values]]=dfresh_vol

##

# now we have the per-2d-element, per-hour, accumulation of freshwater.
plt.figure(2).clf()
fig,ax=plt.subplots(num=2)

mask=np.abs(dfresh_vols)>0.001
coll=G.plot_cells(ax=ax,lw=0.5,values=dfresh_vols,mask=mask)
coll.set_clim([0,500])

e2c=G.edge_to_cells()
G.plot_edges(ax=ax,lw=1.0,color='k',mask=e2c.min(axis=1)<0)

plt.setp(coll,edgecolor='face')
plt.colorbar(coll)

import glob
for pli_fn in glob.glob(os.path.join(run_base_dir,'*.pli')):
    if '_salt' in pli_fn or '_temp' in pli_fn: continue
    feat=dio.read_pli(pli_fn)[0]
    name=feat[0].replace('_flow','')
    ax.text(feat[1][0,0],feat[1][0,1],name)
    ax.plot(feat[1][:,0],feat[1][:,1],'m-',lw=2)


## 

# select a region based on the current axes:
zoom=ax.axis()

# iterate through grids:
cells=np.nonzero( G.cell_clip_mask(zoom) )[0]

print("freshwater volume per h: %.3f"%(dfresh_vols[cells].sum()))


##


# find contiguous areas of freshwater, count up how much freshwater mass
# is there

threshold=0.001

marked=np.zeros( G.Ncells(),'i4')
def visit(c,group):
    "DFS on cell c, adding to group when conditions are met"
    if marked[c]:
        return group

    marked[c]=1

    if np.abs(dfresh_vols[c])<threshold:
        return group
    group.append(c)
    for nbr in G.cell_to_cells(c):
        visit(nbr,group)
    return group

ax.texts=ax.texts[:88]
for c in range(G.Ncells()):
    group=visit(c,[])
    if len(group):
        gcentroid=G.cells_center()[group].mean(axis=0)
        vol=dfresh_vols[group].sum()
        ax.text(gcentroid[0],gcentroid[1],"%.2f"%vol,color='orange')

##

# G.plot_cells(mask=group,color='r',ax=ax)

# Richardson Bay source is just 614.
# EBDA is spot on!
# Mud Slough is way low
# Coyote creek is way low.
# San Jose is spot on.

# Seems like there is a problem with the unit forcing for freshwater flows
# time series for Coyote is 1.0 for both nodes.
##

# Carve out bathymetry near sources:
cells_to_dredge=[]

for pli_fn in glob.glob(os.path.join(run_base_dir,'*.pli')):
    if '_salt' in pli_fn or '_temp' in pli_fn: continue
    # A bit dirty, but use the presence of flow to signal that
    # this is used as an external flow boundary, otherwise taken
    # to be a discharge.
    feat=dio.read_pli(pli_fn)[0]
    name=feat[0].replace('_flow','')
    
    if ('_flow' in pli_fn) or ('Sea' in pli_fn):
        print("%s: external flow boundary"%name)
        feat_edges=dfm_grid.polyline_to_boundary_edges(G,feat[1])
        assert len(feat_edges)
        cells_to_dredge.append( G.edges['cells'][feat_edges].max(axis=1) )
    else:
        print("%s: discharge"%name)
        pnt=feat[1][-1,:]
        cell=G.select_cells_nearest(pnt,inside=True)
        assert cell is not None
        cells_to_dredge.append([cell])
        
cells_to_dredge=np.concatenate(cells_to_dredge)

## 

plt.figure(3).clf()
fig,ax=plt.subplots(num=3)
g.plot_edges(ax=ax,lw=0.3)

coll=g.plot_cells(mask=cells_to_dredge,ax=ax,color='orange',lw=3)
coll.set_edgecolors('orange')

##

Gdredge=dfm_grid.DFMGrid('sfei_v19_net.nc')

nodes_to_dredge=np.concatenate( [G.cell_to_nodes(c)
                                 for c in cells_to_dredge] )
nodes_to_dredge=np.unique(nodes_to_dredge)

dredge_depth=-0.5

Gdredge.nodes['depth'][nodes_to_dredge] = np.minimum(Gdredge.nodes['depth'][nodes_to_dredge],
                                                     dredge_depth)
# No cell depth in this grid, so no need to dredge that
dfm_grid.write_dfm(Gdredge,'sfei_v19_dredge_net.nc')

##

# Plot nice version of that
plt.figure(4).clf()
fig,ax=plt.subplots(num=4)

Gdredge.contourf_node_values( Gdredge.nodes['depth'],
                              np.linspace(-2,3.5,12),
                              extend='both',
                              ax=ax)
ax.axis('equal')
ax.axis((587218.08074605372, 594548.22297297686, 4144516.3193028411, 4149286.1178460922))
