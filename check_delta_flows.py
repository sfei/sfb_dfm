# the validation plots are showing bad juju in the north.

import xarray as xr
from stompy.spatial import proj_utils
from stompy import utils

# sfb1332 is nearest the inflow boundary.

##

adcp=xr.open_dataset('/opt/data/noaa/ports/SFB1332-2013.nc')

adcp_xy=proj_utils.mapper('WGS84','EPSG:26910')([adcp.longitude,adcp.latitude])
##

output_path = "/opt/data/delft/sfb_dfm_v2/runs/wy2013a/DFM_OUTPUT_wy2013a/"	
model_files = output_path + "wy2013a_0000_20120801_000000_his.nc"
model = xr.open_dataset(model_files)


model_xy=np.c_[ model.station_x_coordinate.values,
                model.station_y_coordinate.values ]
##

# What

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

ax.plot( model.station_x_coordinate.values,
         model.station_y_coordinate.values,
         'go')

ax.plot( [adcp_xy[0]],[adcp_xy[1]],'ro')
##

dists=utils.dist( model_xy, adcp_xy )
best=np.argmin(dists)
dists[best] # yep!

##

# Rio Vista:
rio=xr.open_dataset('/opt/data/delft/sfb_dfm_v2/runs/wy2013a/rio_vista-raw.nc')

##

plt.figure(2).clf()

fig,axs=plt.subplots(3,1,num=2,sharex=True)


axs[0].plot(adcp.time,adcp.u_davg,label='U davg')
axs[1].plot(adcp.time,adcp.v_davg,label='V davg')

axs[2].plot(ds.time,
            ds.stream_flow_mean_daily,
            label='Rio new Raw')

#axs[2].plot(rio.time + np.timedelta64(8,'h'),
#            -rio.stream_flow_mean_daily,
#            label='Rio Flip/shift')
axs[0].legend()
axs[1].legend()
axs[2].legend()

axs[2].axis( (735104.96660250274,
              735107.16195219487,
              -126058.97073039225,
              127174.6378024802) )

##

# Okay -
run_start=adcp.time[0].astype('M8[D]')

run_stop=run_start + np.timedelta64(10,'D')
from stompy.io.local import usgs_nwis

import pdb
ds=usgs_nwis.nwis_dataset(station="11455420",
                          start_date=run_start,end_date=run_stop,
                          products=[60],days_per_request=30)

