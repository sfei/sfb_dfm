import pandas as pd
import matplotlib.pylab as plt
import numpy as np

df_in = pd.read_csv('nwis.waterdata.usgs.gov.txt',comment='#',skiprows=0,sep='\t').loc[1:]

year_out = 2013
time_out = np.arange(np.datetime64('%d-01-01' % year_out), np.datetime64('%d-01-01' % (year_out+1)), np.timedelta64(1,'D'))
day_out = (time_out - time_out[0])/np.timedelta64(1,'D')
nout = len(day_out)

df_out = pd.DataFrame(index=day_out)
df_out.index.name = 'Decimal Day (UTC)'

# munge the data
time = df_in['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') 
temp = pd.to_numeric(df_in['14848_00010'].values, errors='coerce')

# take only non nan data
ind = ~np.isnan(temp)
time = time[ind]
temp = temp[ind]

# find years in the dataset
years = pd.DatetimeIndex(time).year.values
time_jan1 = np.array([np.datetime64('%d-01-01' % year) for year in years])

# subtract first of the year to get day of year at each time
day = np.floor((time - time_jan1)/np.timedelta64(1,'D'))

# now interpolate
temp_out = np.nan*np.ones(nout)
for iout in range(nout):
	ind = day == day_out[iout]
	temp_out[iout] = np.mean(temp[ind])

df_out['Temp (oC)'] = temp_out

# plot
fig, ax = plt.subplots(figsize=(11,4))
ax.plot(day, temp, '.', label='data from %s to %s' % (time[0],time[-1]))
ax.plot(day_out, temp_out, label='climatological average')
ax.legend()
ax.set_ylabel('Temperature (oC)')
ax.set_xlabel('Decimal Day (UTC)')
ax.set_title('Alameda Creek Near Niles')
fig.savefig('Alameda_Creek.png')


df_out.to_csv('Alameda_Creek_Climatological_Temperature.csv')











