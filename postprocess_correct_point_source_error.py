''' this python script reads the path to a *.hyd file (master file for DFM DWAQ output), and makes 
a mass conservation correction to the associated binary files at the sites of POTW and tributary inflows'''


from __future__ import print_function
import argparse
import sys,os
import numpy as np
import stompy.model.delft.waq_scenario as waq

hyd_fn= os.environ.get('HYDRO_PATH')        # get path to *.hyd file from BASH enviornment variable
hydro=waq.HydroFiles(hyd_fn)                # load hydro run
hydro.adjust_boundaries_for_conservation()  # make mass correction in place to DWAQ output

