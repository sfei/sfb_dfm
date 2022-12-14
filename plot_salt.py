"""
Some basic scalar plots to see how the salinity field is looking
"""
import os

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
import stompy.plot.cmap as scmap

##
run_name="wy2013"
run_base_dir="runs/" +run_name
output_dir=os.path.join(run_base_dir,"DFM_OUTPUT_%s"%run_name)

Nproc=16
procs=range(Nproc)

all_ds=[xr.open_dataset(os.path.join(output_dir,'%s_%04d_20120801_000000_map.nc'%(run_name,proc)))
        for proc in procs]

all_g=[unstructured_grid.UnstructuredGrid.from_ugrid(ds)
       for ds in all_ds]

# Global grid so we can pull out clusters of mass
G=dfm_grid.DFMGrid(os.path.join(run_base_dir,'sfei_v19_net.nc'))

##

# Walk through the subdomains, figure out which cell in G each
# subdomain cell belongs to.  Could be faster, but as is takes
# about 30s.
global_elems=np.zeros(G.Ncells(),'int32')-1

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

time=1

salt=np.zeros(G.Ncells(),'f8')

for proc in procs:
    # steps 3 and 4 should be clear of the 3h ramp up.
    one_salt=all_ds[proc].sa1.isel(time=time,laydim=9).values
    salt[global_to_G[all_ds[proc].FlowElemGlobalNr.values]]=one_salt

##

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

coll=G.plot_cells(ax=ax,lw=0.5,values=salt)
coll.set_clim([10,36])
coll.set_cmap( scmap.load_gradient('hot_desaturated.cpt') )

plt.setp(coll,edgecolor='face')
plt.colorbar(coll)

ax.xaxis.set_visible(0)
ax.yaxis.set_visible(0)
fig.tight_layout()
##

def update_salt(tidx):
    salt=np.zeros(G.Ncells(),'f8')

    for proc in procs:
        # steps 3 and 4 should be clear of the 3h ramp up.
        one_salt=all_ds[proc].sa1.isel(time=tidx,laydim=9).values
        salt[global_to_G[all_ds[proc].FlowElemGlobalNr.values]]=one_salt
    coll.set_array(salt)
    ax.texts=[]
    ax.text(0.05,0.92,str(all_ds[0].time.values[tidx]),
            transform=ax.transAxes)
    plt.draw()


##

for tidx,time in enumerate(all_ds[0].time.values):
    update_salt(tidx)
    plt.pause(0.05)
    if tidx>10:
        break

##

numlimdt=np.zeros(G.Ncells())

for proc in procs:
    # steps 3 and 4 should be clear of the 3h ramp up.
    one_lim=all_ds[proc].numlimdt.isel(time=-1).values
    numlimdt[global_to_G[all_ds[proc].FlowElemGlobalNr.values]]=one_lim
coll.set_array(np.log10(1+numlimdt))
coll.set_clim([0,6])
