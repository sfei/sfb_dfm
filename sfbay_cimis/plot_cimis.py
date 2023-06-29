import xarray as xr
import matplotlib.pylab as plt

## 

dat=xr.open_dataset('union_city-hourly.nc')

## 
fig, ax = plt.subplots(3,1,figsize=(20,12))

for i,col in enumerate(['HlyPrecip', 'HlyEto', 'HlySolRad']):

	ax[i].plot(dat['time'],dat[col].values)
	ax[i].set_ylabel(col + ' ' + dat[col].units)

ax[0].set_title('Union City CIMIS')

fig.savefig('union_city-hourly.png')