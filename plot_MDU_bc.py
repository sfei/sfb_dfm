# -*- coding: utf-8 -*-
"""
plotting utility for DFM run. trying to take advantage of rusty's stompy utilities. 

This is a stand alone script that calls the same function as called in sfb_dfm. 
If you already have an mdu file and don't want to generate again, this is the way to go. 

@author: siennaw
"""

# MAIN INPUT: the mdu filename     
mdu_file  = r'/hpcvol1/hpcshared/Hydro_Runs/wy2018/runs/wy2018/wy2018.mdu' #<<<< input 
            # mdu file name + path
            
grid_file = r'/hpcvol1/hpcshared/hydro/WY2013/agg141-tau-lp-pass-params/flowgeom.nc' #<<<< input
            #  grid name  --> it's just for a visual, so unless grid has changed radically you can likely leave it as is)

#-------------------------------------------------
import os
abspath = os.path.abspath(__file__)
dname   = os.path.dirname(abspath)
os.chdir(dname) # Change working directory to script location
import sys
from sfb_dfm_utils import plot_mdu 

print('Now making plots...')
plot_mdu.plot_MDU(mdu_file, grid_file)
print('Done.')
