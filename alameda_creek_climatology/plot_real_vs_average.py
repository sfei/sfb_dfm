import pandas as pd
import matplotlib.pylab as plt
import numpy as np

df_in = pd.read_csv('nwis.waterdata.usgs.gov.txt',comment='#',skiprows=0,sep='\t').loc[1:]
df_avg = pd.read_csv('Alameda_Creek_Climatological_Temperature.csv')

# measured temps
time = df_in['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') 
temp = pd.to_numeric(df_in['14848_00010'].values, errors='coerce')

# cast climatalogical average onto time axis
time_CA = np.datetime64('2007-10-01T09:00:00.000000000') + np.arange(584470)*np.timedelta64(15,'m')
DD_CA = (time_CA - time_CA.astype('datetime64[Y]').astype('datetime64[D]'))/np.timedelta64(1,'D')
temp_CA = np.interp(DD_CA, df_avg['Decimal Day (UTC)'].values, 
	                       df_avg['Temp (oC)'].values)

# plot
fig, ax = plt.subplots(figsize=(16,6),constrained_layout=True)
ax.plot(time, temp, 'b', label='Alameda Creek Near Niles NWIS')
ax.plot(time_CA, temp_CA, 'r', label='Climatalogical Average')
ax.legend()
ax.set_ylabel('Water Temp (oC)')
fig.savefig('Compare_Raw_to_Avg_Alameda_Creek_Temps.png',dpi=300)










