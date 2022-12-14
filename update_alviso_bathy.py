import xarray as xr
import numpy as np

from stompy import utils
from stompy.grid import unstructured_grid
from stompy.model.delft import dfm_grid
from stompy.plot import plot_utils

## 
ds=xr.open_dataset('sfei_v18_net.nc')
# g=unstructured_grid.UnstructuredGrid.from_ugrid(ds)
g=dfm_grid.DFMGrid(ds)

##

plt.figure(10).clf()
fig,ax=plt.subplots(num=10)


g.plot_edges(lw=0.5,color='k')

# just need to update NetNode_z:
ncoll=g.plot_nodes(values=ds.NetNode_z)

plot_utils.cbar(ncoll)

## 
rect_1=(593511.,594720.,4136812.,4137722.)
rect_2=(589736., 591613., 4138290, 4139705)

node_xy=np.c_[ds.NetNode_x.values,ds.NetNode_y.values]

for xxyy in [rect_1,rect_2]:
    sel=utils.within_2d(node_xy,xxyy)
    ds.NetNode_z.values[sel]=-1

##

ds.to_netcdf('dredged-sfei_v18_net.nc',
             format='NETCDF3_CLASSIC')
