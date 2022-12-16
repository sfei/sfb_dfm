import numpy as np 
import xarray as xr 
import pandas as pd 
import matplotlib.pyplot as plt 
plt.ion() 
from datetime import datetime 
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline

dat = xr.open_dataset("union_city-hourly-2001-2016.nc")
ds = dat.to_dataframe()
dates = pd.date_range('2001-02-05 01:00:00', '2019-01-29 00:00:00', freq='H')
precip = pd.Series(ds.HlyPrecip.values, index=dates)
clim_precip = precip.groupby(precip.index.dayofyear).mean()
eto = pd.Series(ds.HlyEto.values, index=dates)
clim_eto = eto.groupby(eto.index.dayofyear).mean()
solrad = pd.Series(ds.HlySolRad.values, index=dates)
clim_solrad = solrad.groupby(solrad.index.dayofyear).mean()

### fill precip data
x = np.arange(len(dates))
valid = np.where(np.isfinite(precip)==True)[0]
f = InterpolatedUnivariateSpline(x[valid], precip.values[valid], k=1)
iprecip = f(x)


### fill eto
### hourly climatology eto
dayofyear = np.asarray([d.dayofyear for d in dates])
hour = np.asarray([d.hour for d in dates])

hour_clim = np.zeros(367*24)
hour_clim_ind = np.zeros(367*24)
solrad_clim = np.zeros(367*24)
solrad_clim_ind = np.zeros(367*24)
for i in range(len(dates)):
    ind = dayofyear[i]*24 + hour[i] - 1
    if np.isfinite(eto.values[i])==True:
        hour_clim[ind] += eto.values[i]
        hour_clim_ind[ind] += 1 
    if np.isfinite(solrad.values[i])==True:
        solrad_clim[ind] += solrad.values[i]
        solrad_clim_ind[ind] += 1

hour_clim /= hour_clim_ind
solrad_clim /= solrad_clim_ind

invalid = np.where(np.isfinite(eto)==False)[0]
invalid_dates = dates[invalid]
invalid_dayofyear = np.asarray([d.dayofyear for d in invalid_dates])
invalid_hour = np.asarray([d.hour for d in invalid_dates])
invalid_ind = invalid_dayofyear*24 + invalid_hour -1
for i in range(len(invalid)):
    fill_ind = invalid_dayofyear[i]*24 + invalid_hour[i] - 1
    eto[invalid[i]] = hour_clim[fill_ind]

invalid = np.where(np.isfinite(solrad)==False)[0]
invalid_dates = dates[invalid]
invalid_dayofyear = np.asarray([d.dayofyear for d in invalid_dates])
invalid_hour = np.asarray([d.hour for d in invalid_dates])
invalid_ind = invalid_dayofyear*24 + invalid_hour -1
for i in range(len(invalid)):
    fill_ind = invalid_dayofyear[i]*24 + invalid_hour[i] - 1
    solrad[invalid[i]] = solrad_clim[fill_ind]


#### fix date formatting and input new data into existing xarray dataset
dat['HlyPrecip'] = xr.DataArray(iprecip, dims=['Date'], coords={'Date': dat.Date.values}, attrs={'units': 'mm'})
dat['HlyEto'] = xr.DataArray(eto, dims=['Date'], coords={'Date': dat.Date.values}, attrs={'units': 'mm'}) 
dat['HlySolRad'] = xr.DataArray(solrad, dims=['Date'], coords={'Date': dat.Date.values}, attrs={'units': 'mm'})

dat.to_netcdf('union_city-hourly-2001-2016-filled.nc')