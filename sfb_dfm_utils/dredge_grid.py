"""
Utility method for modifying depths on a grid to allow inflows
in poorly resolved areas which would otherwise be dry.
"""
import numpy as np
from stompy.model.delft import dfm_grid

def dredge_boundary(g,linestring,dredge_depth):
    """
    Lower bathymetry in the vicinity of external boundary, defined
    by a linestring.

    g: instance of unstructured_grid, with a node field 'depth'
    linestring: [N,2] array of coordinates
    dredge_depth: positive-up bed-level for dredged areas

    Modifies depth information in-place.
    """
    # Carve out bathymetry near sources:
    cells_to_dredge=[]

    feat_edges=dfm_grid.polyline_to_boundary_edges(g,linestring)
    assert len(feat_edges)
    cells_to_dredge=g.edges['cells'][feat_edges].max(axis=1)

    nodes_to_dredge=np.concatenate( [g.cell_to_nodes(c)
                                     for c in cells_to_dredge] )
    nodes_to_dredge=np.unique(nodes_to_dredge)

    g.nodes['depth'][nodes_to_dredge] = np.minimum(g.nodes['depth'][nodes_to_dredge],
                                                   dredge_depth)


def dredge_discharge(g,linestring,dredge_depth):
    linestring=np.asarray(linestring)
    pnt=linestring[-1,:]
    cell=g.select_cells_nearest(pnt,inside=True)
    assert cell is not None
    nodes_to_dredge=g.cell_to_nodes(cell)

    g.nodes['depth'][nodes_to_dredge] = np.minimum(g.nodes['depth'][nodes_to_dredge],
                                                   dredge_depth)
        
    
