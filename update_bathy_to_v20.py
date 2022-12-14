# Fix some bad bathymetry values in sfei_v19_net.nc
import matplotlib.pyplot as plt
from stompy.plot import plot_utils
from shapely import geometry


from stompy.grid import unstructured_grid
from stompy.spatial import field
from stompy.model.delft import dfm_grid

##

g=dfm_grid.DFMGrid('sfei_v19_net.nc')

##

dem=field.GdalGrid('/opt/data/GIS/bathymetry/joined-40m.tif')

## 


##
if 0:
    # Manually drawn polygon around Napa, and only muck with those nodes:
    plt.figure(2).clf()
    fig,ax=plt.subplots(num=2)

    g.contourf_node_values(g.nodes['depth'],np.linspace(-10,2,20),extend='both',cmap='jet')
    ax.axis('equal')
    
    pnts=plot_utils.draw_polyline(ax=ax)
else:
    pnts=np.array([[  565555.92717692,  4214418.23511226],
                   [  566275.67385296,  4213788.45677072],
                   [  567295.31497736,  4214388.24566743],
                   [  563996.47604549,  4226863.85471885],
                   [  559707.98543407,  4227793.52750873],
                   [  559108.19653736,  4226264.06582214],
                   [  563456.66603846,  4216757.4118094 ],
                   [  565136.07494923,  4215167.97123314]])

##

poly=geometry.Polygon(pnts)
sel_nodes=g.select_nodes_intersecting(poly)
g.nodes['depth'][sel_nodes] = dem( g.nodes['x'][sel_nodes] )

##

dfm_grid.write_dfm(g,'sfei_v20_net.nc')
