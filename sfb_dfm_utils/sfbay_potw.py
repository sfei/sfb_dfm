import os
import numpy as np
import xarray as xr
from stompy import utils
import stompy.model.delft.io as dio
import six

#from . import dredge_grid

DAY=np.timedelta64(1,'D')

def add_sfbay_potw(mdu,
                   rel_src_dir, # added rel_src_dir alliek dec 2020
                   potw_dir,
                   adjusted_pli_fn,
                   grid,dredge_depth,
                   all_flows_unit=False,
                   time_offset=None,
                   write_salt=True,write_temp=True):
    """
    time_offset: shift all dates by the given offset.  To run 2016 
    with data from 2015, specify np.timedelta64(-365,'D')

    write_salt: leave as True for older DFM, and newer DFM only set to
    true when the simulation includes salinity.

    write_temp: same, but for temperature
    """
    run_base_dir=mdu.base_path
    ref_date,run_start,run_stop = mdu.time_range()
    old_bc_fn=mdu.filepath(["external forcing","ExtForceFile"])

    if time_offset is not None:
        run_start = run_start + time_offset
        run_stop = run_stop + time_offset
        ref_date = ref_date + time_offset
        
    potws=xr.open_dataset(os.path.join(potw_dir,'outputs','sfbay_delta_potw.nc'))
    adjusted_features=dio.read_pli(adjusted_pli_fn)

    # select a time subset of the flow data, starting just before the
    # simulation period, and running beyond the end:
    time_pnts = np.searchsorted(potws.time, [run_start-DAY,run_stop+DAY])
    time_pnts = time_pnts.clip(0,len(potws.time)-1)
    time_idxs=range(time_pnts[0],time_pnts[1]) # enumerate them for loops below

    with open(old_bc_fn,'at') as fp:
        for site in potws.site.values:
            # NB: site is bytes at this point
            potw=potws.sel(site=site)
            try:
                site_s=site.decode()
            except AttributeError:
                site_s=site # py2

            if site_s in ['false_sac','false_sj']:
                six.print_("(skip %s) "%site_s,end="")
                continue

            if np.isnan(potw.utm_x.values.mean()):
                # Delta POTWs are in this file, too, but not in this
                # grid.  Luckily they are easy to identify based on
                # x coordinate.
                six.print_("(skip %s -- coordinates are nan) "%site_s,end="")
                continue
            if potw.utm_x.values.mean() > 610200:
                # Delta POTWs are in this file, too, but not in this
                # grid.  Luckily they are easy to identify based on
                # x coordinate.
                six.print_("(skip %s -- too far east) "%site_s,end="")
                continue
            
            six.print_("%s "%site_s,end="")

            fp.write( ("QUANTITY=discharge_salinity_temperature_sorsin\n"
                       "FILENAME=%s/%s.pli\n"
                       "FILETYPE=9\n"
                       "METHOD=1\n"
                       "OPERAND=O\n"
                       "AREA=0 # no momentum\n"
                       "\n")%(rel_src_dir,site_s) ) # added rel_src_dir alliek dec 2020

            # Write the location - writing a single point appears to work,
            # based on how it shows up in the GUI.  Otherwise we'd have to
            # manufacture a point outside the domain.
            with open(os.path.join(run_base_dir,rel_src_dir,'%s.pli'%site_s),'wt') as pli_fp: # added rel_src_dir alliek dec 2020
                # Scan adjusted features for a match to use instead
                # This is handled slightly differently with POTWs - use the

                # put the depth at -50, should be at the bed
                feat=[site_s,
                      np.array([[potw.utm_x.values,potw.utm_y.values,-50.0]]),
                      ['']]

                for adj_feat in adjusted_features:
                    if adj_feat[0] == site_s:
                        # Merge them if the adjusted feature is more than 10 m away
                        # (to allow for some rounding in the ascii round-trip.)
                        offset=utils.dist( adj_feat[1][-1][:2] - feat[1][-1][:2] )
                        if offset > 10.0:

                            # Just add on the extra point - but may have to promote one 
                            # or the other to 3D.
                            old_geo=feat[1]
                            new_geo=adj_feat[1][-1:]
                            if old_geo.shape[1] != new_geo.shape[1]:
                                if old_geo.shape[1]<3:
                                    old_geo=np.concatenate( (old_geo,0*old_geo[:,:1]), axis=1)
                                else:
                                    # copy depth from old_geo
                                    new_geo=np.concatenate( (new_geo,
                                                             old_geo[-1,-1]+0*new_geo[:,:1]),
                                                            axis=1)

                            # if the original feature was outside the grid, then all is well,
                            # and it will show up in the GUI as a line from the original location
                            # outside the domain to the new location in the domain.
                            if grid.select_cells_nearest(old_geo[-1,:2],inside=True) is None:
                                feat[1]=np.concatenate( (old_geo,new_geo),axis=0 )
                                if len(feat)==3: # includes labels, but they don't matter here, right?
                                    feat[2].append('')
                            else:
                                # but if the original location is inside the grid, this will be interpreted
                                # as a sink-source pair, so we instead just put the single, adjusted
                                # location in.  This is done after potentially copying z-coordinate
                                # data from the original location.
                                feat[1]=new_geo
                        break

                dio.write_pli(pli_fp,[feat])

                # removed by alliek august 2024, doesn't work with new ugrid
                # format, and our updated grid already has dredged boundaries
                # to -0.5m so it is unnecessary to update at this time
                #dredge_grid.dredge_discharge(grid,feat[1],dredge_depth)
                
            with open(os.path.join(run_base_dir,rel_src_dir,'%s.tim'%site_s),'wt') as tim_fp: # added rel_src_dir alliek dec 2020
                for tidx in time_idxs:
                    tstamp_minutes = (potw.time[tidx]-ref_date) / np.timedelta64(1,'m')

                    if all_flows_unit:
                        flow_cms=1.0
                    else:
                        flow_cms=potw.flow[tidx]

                    items=[tstamp_minutes,flow_cms]
                    if write_salt:
                        items.append(0.0)
                    if write_temp:
                        items.append(20.0)

                    tim_fp.write(" ".join(["%g"%v for v in items])+"\n")

    six.print_("Done with POTWs")
