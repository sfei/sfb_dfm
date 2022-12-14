"""
Write out a shapefile for the grid edges.

Run this after changing the grid.

Output written to "derived" subdirectory.
"""

from stompy.model.delft import dfm_grid

g=dfm_grid.DFMGrid('sfei_v19_net.nc')

##

g.write_edges_shp('derived/grid-edges.shp')

##

cell_depth=g.interp_node_to_cell(g.nodes['depth'])

g.write_cells_shp('derived/grid-cells.shp',
                  extra_fields=[ ('depth',cell_depth) ])

##

g.write_shore_shp('derived/grid-boundary.shp')
