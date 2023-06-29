# -*- coding: utf-8 -*-
"""
plotting utility for DFM run. trying to take advantage of rusty's stompy utilities. 

if you point this to the MDU file it should do the rest!!! (fingers crossed)

this script works by reading in the PLI and TIM files. The .pli file has all the geometric/ geographic 
information about the boundary condition while the tim file has the time series with data. 

@author: siennaw
"""
import numpy as np 
import sys, os
import pandas as pd
import matplotlib
#matplotlib.use('agg',warn=False, force=True)    # need to do this so we can re-load pyplot with new backend
import stompy.model.delft.io as dio
from pathlib import Path
#matplotlib.pyplot.switch_backend('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages
from osgeo import ogr
import json
import matplotlib.gridspec as gridspec

# MAIN INPUT: the mdu filename     
# mdu_filename = r'/hpcvol2/Open_Bay/Hydro_model/Full_res/WY2017/wy2017_usgs_bc/wy2017_usgs_bc.mdu' #<<<< input here


def plot_MDU(mdu_filename, gridpath): 

    numnanfiles = 0
    nanprint = '\n'

    #------------- script now takes over -------------------------------------------
    mdu_filename = Path(mdu_filename)
    base_dir     = mdu_filename.parent  # The assumption is that we'll find all our bc's in the same folder as the mdu.
    folder_dir   = base_dir / 'bc_figures'
    pdf_folder_dir = base_dir / 'bc_figures_compiled'
    folder_dir.exists() or folder_dir.mkdir() 
    pdf_folder_dir.exists() or pdf_folder_dir.mkdir() 
    
    # read the grid outline
    file = ogr.Open(gridpath)
    shape = file.GetLayer(0)
    #first feature of the shapefile
    feature = shape.GetFeature(0)
    grid = json.loads(feature.ExportToJson())['geometry']['coordinates']
    
    # Open MDU, strip time information using stompy functionality
    MDU = dio.MDUFile(filename=str(mdu_filename))
    t_ref, t_start, t_stop =  MDU.time_range() 

    # get the path to the external forcing fiele
    for row in MDU.rows:
        if 'extforcefile' in row.lower():
            if row.lower().replace(' ','').split('=')[0] == 'extforcefile':
                ext_fn = row.replace(' ','').split('=')[1]
                if '#' in ext_fn:
                    ext_fn = ext_fn.split('#')[0]
        if 'tunit' in row.lower():
            if row.lower().replace(' ','').split('=')[0] == 'tunit':
                tunit = row.replace(' ','').split('=')[1]
                if '#' in tunit:
                    tunit = tunit.split('#')[0]
                tunit = tunit.lower()

    # list of types of boundary conditions in order we want them plotted
    boundary_type_list =  ['frictioncoefficient',
                           'initialsalinity',
                           'initialtemperature',
                           'secchidepth',
                           'rainfall',
                           'waterlevelbnd',
                           'discharge_salinity_temperature_sorsin',
                           'dischargebnd',
                           'salinitybnd',
                           'temperaturebnd',
                           'windx',
                           'windy',
                           'humidity_airtemperature_cloudiness']

    # read the external forcing file and build up a dictionary describing the boundary conditions
    boundary_conditions = {}
    ibc = 0
    with open(os.path.join(base_dir, ext_fn), 'r') as f:
        lines = f.readlines()
    iline = 0
    nline = len(lines)
    while iline < nline:
        if 'QUANTITY' in lines[iline]:
            btype = lines[iline].replace(' ','').strip().split('=')[1]
            iline += 1
            while not 'FILENAME' in lines[iline]:
                iline += 1
                if iline >= nline:
                    raise Exception('reached end of file without finding FILENAME')
            filename = lines[iline].replace(' ','').strip().split('=')[1]
            boundary_conditions[ibc] = {'type' : btype, 'filename' : filename}
            iline += 1
            ibc += 1
            if not btype in boundary_type_list:
                raise Exception('plot_mdu.py cannot handle boundary type %s' % btype)
        else:
            iline += 1
    nbc = len(boundary_conditions) 

    # define shared plotting functions 
    def format_xaxis (axis):
        months = mdates.MonthLocator(interval = 2)  # every other month
        fmt = mdates.DateFormatter('%b/%Y')
        axis.xaxis.set_major_locator(months)
        axis.xaxis.set_major_formatter(fmt)
        axis.set_xlim(t_ref, t_stop)
        
    def save_image(fig, name, pdf):
        fullname = folder_dir / (name + '.png')
        fig.savefig(str(fullname), dpi = 300, bbox_inches='tight')
        pdf.savefig(dpi = 300, bbox_inches='tight')
        print('Saved %s' % fullname)
        plt.close()

    # open pdf to comile plots
    pdffile = str(pdf_folder_dir / 'compiled_bc_plots.pdf')
    with PdfPages(pdffile) as pdf: 

        # loop through the different types of boundary conditions
        for boundary_type in boundary_type_list:
    
            for ibc in range(nbc):

                btype = boundary_conditions[ibc]['type']
                fname = boundary_conditions[ibc]['filename']

                if btype == boundary_type:
                
                    if btype in ['windx','windy','humidity_airtemperature_cloudiness']:
    
                        print('Scanning %s for NaN, need to add code if you want to plot it...' % fname)

                        with open(os.path.join(base_dir,fname), 'r') as f:
                            lines = f.readlines()
                            for line in lines:
                                if 'nan' in line.lower():
                                    nanprint += ('WARNING!!! NaN found in %s\n' % fname)
                                    numnanfiles+=1

                        fig, ax = plt.subplots(figsize=(8,8))
                        ax.axis((0,1,0,1))
                        ax.text(0.4,0.5, 'NOT PLOTTED')
                        ax.set_title('%s (%s)' % (fname, boundary_type))
                        save_image(fig, os.path.basename(fname).replace('.','_DOT_'), pdf)

    
                    elif btype in ['frictioncoefficient','initialtemperature','initialsalinity','secchidepth']:
    
                        df1 = pd.read_csv(os.path.join(base_dir,fname), delim_whitespace=True, header=None)
    
                        ####### check for NaN ######
                        if np.any(np.isnan(df1.values)):
                            nanprint += ('WARNING!!! NaN found in %s\n' % fname)
                            numnanfiles+=1
    
                        # make a plot
                        fig, ax = plt.subplots(tight_layout=True, figsize=(11,11))
    
                        # add the grid outline
                        for poly in grid:
                            xp = np.array(poly)[:,0]
                            yp = np.array(poly)[:,1]
                            ax.plot(xp,yp,'lightblue')
    
                        # plot the data as scatter
                        sc = ax.scatter(df1[0], df1[1], c=df1[2], cmap='jet')
                        cbar = plt.colorbar(sc, ax=ax)
                        ax.set_title('%s (%s)' % (fname, boundary_type))
                        ax.axis('off')
                        ax.axis('equal')
                        save_image(fig, os.path.basename(fname).strip('.xyz'), pdf)

                    elif btype == 'discharge_salinity_temperature_sorsin':

                        # initialize the figure and add the map
                        fig = plt.figure(tight_layout=True, figsize=(11,9))
                        gs = gridspec.GridSpec(3, 4)
                        ax1 = fig.add_subplot(gs[0,0:3])
                        ax2 = fig.add_subplot(gs[1,0:3])
                        ax3 = fig.add_subplot(gs[2,0:3])
                        map_axis = fig.add_subplot(gs[0,3])
                        for poly in grid:
                            xp = np.array(poly)[:,0]
                            yp = np.array(poly)[:,1]
                            map_axis.plot(xp,yp,'lightblue')
                        map_axis.axis('off')
                        map_axis.axis('equal')
        
                        # read the pli file, then read the tim files and plot each one
                        df = pd.read_csv(os.path.join(base_dir, fname), skiprows=2, delim_whitespace=True, header=None)
                        xs = df[0].values
                        ys = df[1].values
                        for i in range(len(df)):
                            map_axis.plot(xs[i], ys[i], 'x')
                        df1 = pd.read_csv(os.path.join(base_dir, fname.replace('.pli','.tim')), delim_whitespace=True, header=None)
                            
                        ####### check for NaN ######
                        if np.any(np.isnan(df1.values)):
                            nanprint += ('WARNING!!! NaN found in %s\n' % fname)
                            numnanfiles+=1
                            
                        # extract the time and scalar values
                        time = t_ref + df1[0].values.astype('timedelta64[%s]' % tunit)
                        flo = df1[1].values
                        sal = df1[2].values
                        tem = df1[3].values
                        ax1.plot(time, flo)
                        ax1.set_ylabel('discharge')
                        ax2.plot(time, sal)
                        ax2.set_ylabel('salinity')
                        ax3.plot(time, tem) 
                        ax3.set_ylabel('temperature')
                        ax1.set_title('%s (%s)' % (fname, boundary_type))
                        for ax in [ax1,ax2,ax3]:
                            ax.grid(b = True, alpha = 0.25)
                            format_xaxis(ax) 
                        save_image(fig, os.path.basename(fname).strip('.pli'), pdf)
    
                    elif btype == 'rainfall':
        
                        # initialize the figure, no map for this one
                        fig, ax1 = plt.subplots(tight_layout=True, figsize=(9,3))
        
                        # read the pli file, then read the tim files and plot each one
                        df1 = pd.read_csv(os.path.join(base_dir, fname), skiprows=2, delim_whitespace=True, header=None)
                        
                        ####### check for NaN ######
                        if np.any(np.isnan(df1.values)):
                            nanprint += ('WARNING!!! NaN found in %s\n' % fname)
                            numnanfiles+=1
                        
                        # extract the time and scalar values
                        time = t_ref + df1[0].values.astype('timedelta64[%s]' % tunit)
                        scalar = df1[1].values
                        ax1.plot(time, scalar) 
                        ax1.set_title('%s (%s)' % (fname, boundary_type))
                        ax1.grid(alpha = 0.25)
                        ax1.set_ylabel('precipitation - evaporation')
                        format_xaxis(ax1) 
                        save_image(fig, os.path.basename(fname).strip('.tim'), pdf)
    
                    elif btype == 'dischargebnd':
        
                        # initialize the figure and add the map
                        fig = plt.figure(tight_layout=True, figsize=(11,3))
                        gs = gridspec.GridSpec(1, 4)
                        ax1 = fig.add_subplot(gs[0:3])
                        map_axis = fig.add_subplot(gs[3])
                        for poly in grid:
                            xp = np.array(poly)[:,0]
                            yp = np.array(poly)[:,1]
                            map_axis.plot(xp,yp,'lightblue')
                        map_axis.axis('off')
                        map_axis.axis('equal')
        
                        # read the pli file, then read the flow files and sum them up
                        df = pd.read_csv(os.path.join(base_dir, fname), skiprows=2, delim_whitespace=True, header=None)
                        xs = df[0].values
                        ys = df[1].values
                        try:
                            names = df[2].values
                        except:
                            names = [os.path.basename(fname).replace('.pli','_%04d' % (i+1)) for i in range(len(df))]
                        for i in range(len(df)):
                            map_axis.plot(xs[i], ys[i], 'rx')
                            df1 = pd.read_csv(os.path.join(base_dir, os.path.dirname(fname), names[i]+'.tim'), delim_whitespace=True, header=None)
                            
                            ####### check for NaN ######
                            if np.any(np.isnan(df1.values)):
                                nanprint += ('WARNING!!! NaN found in %s\n' % fname)
                                numnanfiles+=1
                            
                            # add up the flows
                            if i==0:
                                time = t_ref + df1[0].values.astype('timedelta64[%s]' % tunit)
                                flow = df1[1].values
                                legstr = names[i] + '.tim'
                            else:
                                flow += df1[1].values
                                legstr += ' + ' + names[i] + '.tim'
                        ax1.set_title('%s (%s)' % (fname, boundary_type))
                        ax1.plot(time, flow, label=legstr)    
                        ax1.grid(alpha = 0.25)
                        ax1.legend()
                        ax1.set_ylabel('discharge')
                        format_xaxis(ax1) 
                        save_image(fig, os.path.basename(fname).strip('.pli'), pdf)
    
                    elif btype in ['salinitybnd','temperaturebnd','waterlevelbnd']:
        
                        # initialize the figure and add the map
                        fig = plt.figure(tight_layout=True, figsize=(11,3))
                        gs = gridspec.GridSpec(1, 4)
                        ax1 = fig.add_subplot(gs[0:3])
                        map_axis = fig.add_subplot(gs[3])
                        for poly in grid:
                            xp = np.array(poly)[:,0]
                            yp = np.array(poly)[:,1]
                            map_axis.plot(xp,yp,'lightblue')
                        map_axis.axis('off')
                        map_axis.axis('equal')
        
                        # read the pli file, then read the tim files and plot each one
                        df = pd.read_csv(os.path.join(base_dir, fname), skiprows=2, delim_whitespace=True, header=None)
                        xs = df[0].values
                        ys = df[1].values
                        try:
                            names = df[2].values
                        except:
                            names = [os.path.basename(fname).replace('.pli','_%04d' % (i+1)) for i in range(len(df))]
                        for i in range(len(df)):
                            map_axis.plot(xs[i], ys[i], 'x')
                            df1 = pd.read_csv(os.path.join(base_dir, os.path.dirname(fname), names[i]+'.tim'), delim_whitespace=True, header=None)
                            
                            ####### check for NaN ######
                            if np.any(np.isnan(df1.values)):
                                nanprint += ('WARNING!!! NaN found in %s\n' % fname)
                                numnanfiles+=1
                            
                            # extract the time and scalar values
                            time = t_ref + df1[0].values.astype('timedelta64[%s]' % tunit)
                            scalar = df1[1].values
                            ax1.plot(time, scalar, label=names[i] + '.tim') 
                        ax1.set_title('%s (%s)' % (fname, boundary_type))
                        ax1.grid(alpha = 0.25)
                        ax1.legend()
                        ax1.set_ylabel(btype.replace('bnd',''))
                        format_xaxis(ax1) 
                        save_image(fig, os.path.basename(fname).strip('.pli'), pdf)



                    else:

                        raise Exception('plot_mdu.py cannot handle boundary type %s...' % btype)



    if numnanfiles>0:
        print(nanprint)
    else:
        print('HOORAY! No NaNs found in input files')