"""
A first cut at direct precipitation and evaporation.
"""
import xarray as xr
from stompy import utils

# Last SUNTANS run had used NARR
# it's way way coarse.  Seems better to use an in-Bay climatology
# than to use NARR.

##

union_city=xr.open_dataset('/opt/data/cimis/union_city-hourly-2001-2016.nc')

##

temps=utils.fill_invalid(union_city.HlyAirTmp.values)
temp_zscore=((temps-temps.mean()) / temps.std()).clip(-1,1)
# score of 1 means warm temperature
factors=np.interp(temp_zscore,
                  [-1,1],[0.7,0.6])
union_city['HlyEvap']=union_city.HlyEto/factors

## 

narr_data_dir="/opt/data/suntans/spinupdated/suntans-spinup/narr-data/data/"

precip=xr.open_dataset(os.path.join(narr_data_dir,"apcp.mon.mean.nc"),decode_cf=False)
# bug in the nc file
if precip.apcp._FillValue != precip.apcp.missing_value:
    precip.apcp.attrs['_FillValue']=precip.apcp.missing_value
precip=xr.decode_cf(precip)

## 
evap  =xr.open_dataset(os.path.join(narr_data_dir,"pevap.mon.mean.nc"),decode_cf=False)
if evap.pevap._FillValue != evap.pevap.missing_value:
    evap.pevap.attrs['_FillValue']=evap.pevap.missing_value
evap=xr.decode_cf(evap)

## 
precip=precip.isel(time=precip.time>=union_city.time[0])
evap=evap.isel(time=evap.time>=union_city.time[0])
##
if 0:
    plt.figure(1).clf()
    fig,axs=plt.subplots(2,1,num=1,sharex=True,sharey=True)

    unwrap_lon=precip.lon.values % 360.0
    axs[0].pcolormesh( unwrap_lon,precip.lat, precip.apcp.isel(time=-1),vmin=0,vmax=50 )
    axs[1].pcolormesh( unwrap_lon,evap.lat, evap.pevap.isel(time=-1),vmin=-0.5,vmax=2.5 )
    
##

# find a nearby point and get a timeseries
ll=np.array( [precip.lon.values, precip.lat.values] ).transpose(1,2,0)

myll=np.array([-122.2546,37.6328]) # center of South Bay

dists=utils.dist( ll-myll )
cell=np.nonzero(dists==np.min(dists))
lat_i=cell[0][0]
lon_i=cell[1][0] # i.e. dists[lat_i,lon_i] is our man

##
# And Burlingame measured climatology, monthly accumulated, inches/month

orig_data = np.array( [1.27, 1.81, 3.60, 5.28, 6.85,  7.82,  8.42,  7.39,  5.74,  3.78,  1.98, 1.28] )

burl_ds=xr.Dataset()
burl_ds['time']=('time',),np.arange(union_city.time.values[0],
                                    union_city.time.values[-1],
                                    np.timedelta64(1,'D'))
months0=[ utils.to_datetime(t).month-1 for t in burl_ds.time.values]
# inch=>mm, month=>day
burl_ds['evap']=('time',),orig_data[months0] * 25.4 / 30.0
burl_ds['evap'].attrs['units']='mm d-1'

##

from stompy import filters

plt.figure(2).clf()
fig,(ax,ax_cumul)=plt.subplots(2,1,num=2,sharex=True)

# units: precip: kg/m2, monthly average.  Based on prior scripts, this is maybe a per day number
#   https://www.esrl.noaa.gov/psd/data/gridded/data.narr.monolevel.html#plot
#   mentions that this file is "Monthly average of Daily Accumulation"
# units: evap: kg/m2 monthly accumulated average.  Based on prior scripts, maybe a per 3h number??
# From the file metadata, "should be" mm/month.
ax.plot(utils.to_dnum(precip.time.values),
        precip.apcp.isel(x=lon_i,y=lat_i),'--',label='NARR precip')
ax.plot(utils.to_dnum(evap.time.values),
        8*evap.pevap.isel(x=lon_i,y=lat_i),label='NARR p-evap')

# data in mm, originally at hourly time scale.
ax.plot(utils.to_dnum(union_city.time),
        filters.lowpass_fir(24*union_city.HlyPrecip,10*24),
        '--',label='CIMIS precip')
ax.plot(utils.to_dnum(union_city.time),
        filters.lowpass_fir(24*union_city.HlyEto,10*24),
        label='CIMIS ETO')

ax.plot(utils.to_dnum(union_city.time),
        filters.lowpass_fir(24*union_city.HlyEvap,10*24),
        label='CIMIS Evap')

ax.plot(utils.to_dnum(burl_ds.time),burl_ds.evap,label='Burlingame')

# --

ax_cumul.plot(utils.to_dnum(precip.time.values),
              30*np.cumsum(precip.apcp.isel(x=lon_i,y=lat_i).values),
              '--',label='NARR precip')
ax_cumul.plot(utils.to_dnum(evap.time.values),
              8*30*np.cumsum(evap.pevap.isel(x=lon_i,y=lat_i).values),label='NARR p-evap')

ax_cumul.plot(utils.to_dnum(union_city.time),
              np.cumsum(utils.fill_invalid(union_city.HlyPrecip.values)),
              '--',label='CIMIS precip')
ax_cumul.plot(utils.to_dnum(union_city.time),
              np.cumsum(utils.fill_invalid(union_city.HlyEto.values)),
              label='CIMIS ETO')

ax_cumul.plot(utils.to_dnum(union_city.time),
              np.cumsum(utils.fill_invalid(union_city.HlyEvap.values)),
              label='CIMIS Evap')

ax_cumul.plot(utils.to_dnum(burl_ds.time),
              np.cumsum(burl_ds.evap.values),label='Burlingame Evap')

ax.xaxis.axis_date()

ax.legend()
ax_cumul.legend()

# cumulative precip seems to match well.
# maybe a factor of two difference in evaporation, with
# NARR twice as high

## 
def precip_datasource(self):
    factor=-(1./1000) * 1./(24*3600.) # 24h accumulation for precip
    return NarrKriged(nc_path=os.path.join(narr_data_dir,"apcp.mon.mean.nc"),
                      idxs=self.narr_idxs,
                      field='apcp',factor=factor)

def evap_datasource(self):
    # have to be a bit particular here because 104,131 is nearby, but
    # reads like an ocean point.
    # data have native units of mm/3h
    #       m/ mm      3h / sec
    factor=(1./1000) * 1./(3*3600.) # 3h accumulation for evaporation
    return NarrKriged(nc_path=os.path.join(narr_data_dir,"pevap.mon.mean.nc"),
                      idxs=self.narr_idxs,
                      field='pevap',factor=factor)

## 

# FAO - conversion from alfalfa-based ETO to open-water.
# or https://cals.arizona.edu/azmet/et1.htm
# which says cool period, divide ETO by 0.7 to get pan evaporation,
# warm period, divide by 0.6.
