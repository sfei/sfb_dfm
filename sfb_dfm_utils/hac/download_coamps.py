
# This was temporarily down, but appears to be working now?
ds=xr.open_dataset('http://thredds.cencoos.org/thredds/dodsC/COAMPS_4KM_2M_AIR_TEMP.nc')

ds_x = ds.x
ds_y = ds.y

## 
x_sel=slice(*np.searchsorted(ds_x.values, [450,600] ))
y_sel=slice(*np.searchsorted(ds_y.values, [4000,4250]))

# one problem is that time starts in 2013...
ds_region=ds.isel( x=x_sel,y=y_sel,time=slice(0,30))

##




plt.clf()
ax=plt.gca()

for tidx in range(len(ds_region.time)):
    ax.cla()
    temps=ds_region.air_temperature.isel(time=tidx).values
    ax.imshow(temps,vmin=8,vmax=22,extent=[ds_region.x.min(),
                                           ds_region.x.max(),
                                           ds_region.y.min(),
                                           ds_region.y.max()],
              origin='bottom')
    plt.draw()
    plt.pause(0.1)
