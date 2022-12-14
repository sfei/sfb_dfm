import geopandas as gpd
import pandas as pd
import matplotlib.pylab as plt

# read grid boundary
gdf = gpd.read_file('../derived/grid-boundary.shp')

# read original salinity initial condition for open bay model
df1 = pd.read_csv('../inputs-static/orig-saltopini.xyz', header=None, delim_whitespace=True)
df1.columns = ['x','y','sal']

# read original salinity IC for delta model (top and bottom -- does it really use both???)
df2 = pd.read_csv('saltopini.xyz', header=None, delim_whitespace=True)
df2.columns = ['x','y','sal']
df3 = pd.read_csv('salbotini.xyz', header=None, delim_whitespace=True)
df3.columns = ['x','y','sal']

# make a nice scatter plot
fig, ax = plt.subplots(1,3,figsize=(16,4))
gdf.boundary.plot(ax=ax[0], edgecolor='k')
im = ax[0].scatter(df1.x,df1.y,c=df1.sal, vmin=0,vmax=33, cmap='jet')
ax[0].set_title('open bay IC')
gdf.boundary.plot(ax=ax[1], edgecolor='k')
im = ax[1].scatter(df2.x,df2.y,c=df2.sal, vmin=0,vmax=33, cmap='jet')
ax[1].set_title('delta IC top')
gdf.boundary.plot(ax=ax[2], edgecolor='k')
im = ax[2].scatter(df3.x,df3.y,c=df3.sal, vmin=0,vmax=33, cmap='jet')
ax[2].set_title('delta IC bottom')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
fig.suptitle('salinitiy initial conditions, original')


fig.savefig('sfb_dfm_orig-saltopini.png')

