import pandas as pd
import matplotlib.pylab as plt
import numpy as np

df_jersey = pd.read_csv('nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00010=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2000-01-01&end_date=2023-07-04.txt',comment='#',skiprows=0,sep='\t').loc[1:]
df_riovista = pd.read_csv('nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00010=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2000-01-01&end_date=2023-07-04.txt',comment='#',skiprows=0,sep='\t').loc[1:]
df_avg = pd.read_csv('Delta_Climatological_Temperature.csv')

# extract Jersey temps, fill in missing with BGC project temps
time_J = df_jersey['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') 
temp_J = pd.to_numeric(df_jersey['15568_00010'].values, errors='coerce')
temp1_J = pd.to_numeric(df_jersey['164023_00010'].values, errors='coerce')
#ind = np.isnan(temp_J)
#temp_J[ind] = temp1_J[ind]

# extract Rio Vista temps, fill in missing with BGC project temps
time_R = df_riovista['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') 
temp_R = pd.to_numeric(df_riovista['16194_00010'].values, errors='coerce')
temp1_R = pd.to_numeric(df_riovista['300439_00010'].values, errors='coerce')
#ind = np.isnan(temp_R)
#temp_R[ind] = temp1_R[ind]

# cast climatalogical average onto time axis
time_CA = np.datetime64('2009-01-01') + np.arange(525600)*np.timedelta64(15,'m')
DD_CA = (time_CA - time_CA.astype('datetime64[Y]').astype('datetime64[D]'))/np.timedelta64(1,'D')
temp_R_CA = np.interp(DD_CA, df_avg['Decimal Day (UTC)'].values, 
	                        df_avg['RioVista Temp (oC)'].values)
temp_J_CA = np.interp(DD_CA, df_avg['Decimal Day (UTC)'].values, 
	                        df_avg['Jersey Temp (oC)'].values)

# plot
fig, ax = plt.subplots(2,1,figsize=(16,8.5),constrained_layout=True)
ax[0].plot(time_R, temp_R, 'b', label='Rio Vista NWIS')
ax[0].plot(time_R, temp1_R, 'k', label='Rio Vista NWIS (BGC Project)')
ax[0].plot(time_CA, temp_R_CA, 'r', label='Climatalogical Average')
ax[0].legend()
ax[0].set_ylabel('Water Temp (oC)')
ax[1].plot(time_J, temp_J, 'b', label='Jersey NWIS')
ax[1].plot(time_J, temp1_J, 'k', label='Jersey NWIS (BGC Project)')
ax[1].plot(time_CA, temp_J_CA, 'r', label='Climatalogical Average')
ax[1].legend()
ax[1].set_ylabel('Water Temp (oC)')
fig.savefig('Compare_Raw_to_Avg_Delta_Inflow_Temps.png',dpi=300)










