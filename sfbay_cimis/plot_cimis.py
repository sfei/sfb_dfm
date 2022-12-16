import xarray as xr
import matplotlib.pyplt as plt

## 

df=xr.open_dataset('union_city-hourly-2001-2016.nc')

## 

plt.figure(1).clf()

fig,ax=plt.subplots(num=1)

ax.plot(df.time,df['HlySolRad'])

dnums=utils.to_dnum(df.time.values)

## 

# True daily averages:
daily_rad =df['HlySolRad'].values.reshape([-1,24]).mean(axis=1)
daily_dnum=dnums.reshape([-1,24]).mean(axis=1)

ax.plot(daily_dnum,daily_rad,color='orange',lw=3)
