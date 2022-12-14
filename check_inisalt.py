from stompy.model import unstructured_diffuser
from scipy.interpolate import interp1d
from stompy.model.delft import dfm_grid
import sfb_dfm_utils

## 
run_base_dir='./runs/wy2013c'
static_dir='./inputs-static'
run_start=np.datetime64('2012-08-01')

#sfb_dfm_utils.add_initial_salinity_dyn(run_base_dir,
#                                       abs_static_dir,
#                                       mdu,
#                                       run_start)

g=dfm_grid.DFMGrid(os.path.join(run_base_dir,mdu['geometry','NetFile']))

##

usgs_data_end=np.datetime64('2016-04-28')
usgs_pad=np.timedelta64(30,'D')

usgs_target=run_start

# so we may have to grab a previous years cruise and pretend
while usgs_target + usgs_pad > usgs_data_end:
    usgs_target -= np.timedelta64(365,'D')
        
usgs_cruises=usgs_sfbay.cruise_dataset(usgs_target - usgs_pad,
                                           usgs_target + usgs_pad )

# lame filling
salt3d=usgs_cruises['salinity']
salt2d=salt3d.mean(dim='prof_sample')
assert salt2d.dims[0]=='date'
salt2d_fill=utils.fill_invalid(salt2d.values,axis=0)

salt_f=interp1d(utils.to_dnum(salt2d.date.values),
                salt2d_fill,
                axis=0,bounds_error=False)(utils.to_dnum(usgs_target))
usgs_init_salt=np.c_[salt2d.x.values,salt2d.y.values,salt_f]

##

mooring_xy=[]
mooring_salt=[]

L2_dir='/opt/data/sfei/moored_sensors_csv/L2/'

# tuples (<name in observation points shapefile>, <L2 data file name> )
sfei_moorings=[
    ('ALV',"ALV_all_data_L2.csv"),
    ('SFEI_Coyote',"COY_all_data_L2.csv"),
    ('DB',"DMB_all_data_L2.csv"),
    ('SFEI_Guadalupe',"GL_all_data_L2.csv"),
    ('SFEI_Mowry',"MOW_all_data_L2.csv"),
    ('SFEI_Newark',"NW_all_data_L2.csv"),
    ('SFEI_A8Notch',"POND_all_data_L2.csv"),
    ('SMB',"SM_all_data_L2.csv")
]

# lat/lon from observation-points
# FRAGILE - FIX!
obs_shp=wkb2shp.shp2geom(os.path.join(static_dir,"observation-points.shp"))
##

for name,l2_file in sfei_moorings:
    print(name)
    # low_memory option is just to silence annoying warning
    sfei=pd.read_csv(os.path.join(L2_dir,l2_file),
                     parse_dates=['Datetime','dt'],low_memory=False)
    sfei_salt=sfei['S_PSU']
    valid=~(sfei_salt.isnull())
    sfei_salt_now=utils.interp_near(utils.to_dnum(run_start),
                                    utils.to_dnum(sfei.Datetime[valid]),sfei_salt[valid],
                                    20.0)
    geom=obs_shp['geom'][ np.nonzero(obs_shp['name']==name)[0][0] ]
    xy=np.array(geom)
    if np.isfinite(sfei_salt_now):
        mooring_xy.append(xy)
        mooring_salt.append(sfei_salt_now)

##

if len(mooring_xy):
    xy=np.array(mooring_xy)
    sfei_init_salt=np.c_[xy[:,0],xy[:,1],mooring_salt]
    init_salt=np.concatenate( (usgs_init_salt,
                               sfei_init_salt) )
else:
    init_salt=usgs_init_salt
    
## 
differ=unstructured_diffuser.Diffuser(grid=g)

for x,y,salt in init_salt:
    try:
        differ.set_dirichlet(salt,xy=[x,y],on_duplicate='skip')
    except differ.PointOutsideDomain as exc:
        print(exc)
        continue

differ.construct_linear_system()
differ.solve_linear_system(animate=False)

##
cc=g.cells_centroid()

cc_salt=np.concatenate( ( cc, differ.C_solved[:,None] ),axis=1 )

##

plt.figure(1).clf()
fig,ax=plt.subplots(num=1)

coll=g.plot_cells(ax=ax,values=cc_salt[:,2],cmap='jet')

plt.colorbar(coll)

##

from stompy.spatial import interp_4d

# try again, but with the interp_4d code:
samples=pd.DataFrame()
samples['x']=init_salt[:,0]
samples['y']=init_salt[:,1]
samples['value']=init_salt[:,2]
samples['weight']=1e6*np.ones_like(init_salt[:,0])


salt2,weights=interp_4d.weighted_grid_extrapolation(g,samples,alpha=1e-4,
                                                    return_weights=True)

##


plt.figure(2).clf()
fig,ax=plt.subplots(num=2,sharex=True,sharey=True)

coll2=g.plot_cells(ax=ax,values=salt2,cmap='jet')

scat=ax.scatter(samples.x,samples.y,40,samples.value,lw=1.,cmap='jet')
scat.set_edgecolor('k')
scat.set_clim(coll2.get_clim())

plt.colorbar(coll2,ax=ax)

##

# Something is amiss - the weights should be decaying more strongly as
# it goes out to the coastal ocean, but the weights don't actually decay
# that strongly.
alpha=1e-5
weight_col='weight'
x_col='x'
y_col='y'
value_col='value'


D=unstructured_diffuser.Diffuser(g)
D.set_decay_rate(alpha)
Dw=unstructured_diffuser.Diffuser(g)
Dw.set_decay_rate(alpha)

D.K_j[:]=10.0
Dw.K_j[:]=10.0

for i in range(len(samples)):
    rec=samples.iloc[i]
    weight=rec[weight_col]
    xy=rec[[x_col,y_col]].values
    cell=D.grid.point_to_cell(xy) or D.grid.select_cells_nearest(xy)
    D.set_flux(weight*rec[value_col],cell=cell)
    Dw.set_flux(weight,cell=cell)

D.construct_linear_system()
D.solve_linear_system(animate=False)
Dw.construct_linear_system()
Dw.solve_linear_system(animate=False)

C=D.C_solved
W=Dw.C_solved
T=C / W
salt2=T
weights=W

##

# First - are the fluxes being set in the right places?
# and maybe I need to specify no-flux BCs???
nbc_cells=[nbc[0] for nbc in D.neumann_bcs]

g.plot_cells(mask=np.array(nbc_cells),ax=axw)

# probably there is an issue with how I use alpha.
#
# alpha is decremented from the diagonal
# A[mic2,mic2] gets the sum of fluxes

# the overall equation being solved is
# A*x=b
# with x being concentrations per cell
# A*x gives the new concentrations
# K_j is for this problem constant.  Also constant edge and cell depths
# flux_per_gradient is m2/s * m * m / m * s
#   => m3
# which should get multiplied by a concentration delta, which would have
# units of mass/m3, so flux_per_gradient*(Delta-C) gives mass.
# So line 192 of unstructured_diffuser
#   A[mic1,mic1] += v1
#  v1 here is flux_per_gradient[j] / volume(v1)
#   that's saying that the change in concentrations for cell 1
#  is the gradient between cell 1 and 2, time flux_per_gradient which gives
# mass, so this is a change in concentration.
# adding alpha in there, v=alpha*dt
# A[mic,mic]= -alpha*dt
# xnew=xold - alpha*dt*xold
# so all cell values are trending toward 0 with an exponential decay rate
#  of alpha*dt

xnew=xld

# if it's a numerical problem, maybe K_j can be much higher?
# diagonal of A is already maxed out at -9.  And results are
# no better with K_j=10. (looks like changing alpha by a factor
# of 10.)

# There are some oddly large gradients in places near shorelines.
# e.g. something is funny at (561665.1446250669, 4180633.2094672034
# That's g.cell... 20067

# Dw.A.todok
## 
ic=20067
# no dirichlet, so all cells are calculated
mic=Dw.c_map[ic]

## 
Adok=Dw.A.todok()

## 
# All entries involving ic
for row in range(Adok.shape[0]):
    col=mic
    if (row,col) in Adok:
        print("A[%6d,%6d] = %g"%(row,col,Adok[row,col]))

for col in range(Adok.shape[1]):
    row=mic
    if (row,col) in Adok:
        print("A[%6d,%6d] = %g"%(row,col,Adok[row,col]))

##

nbrs=[20062, 20067,20098,39935,49897]
g.plot_cells(mask=nbrs,ax=axw)
    

# okay! 49897 is a wild card!
# yep - the grid doesn't necessarily keep all negative edges['cells'] in the second slot.
# 
