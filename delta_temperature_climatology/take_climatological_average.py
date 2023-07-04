import pandas as pd
import matplotlib.pylab as plt
import numpy as np

df_jersey = pd.read_csv('nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00010=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2000-01-01&end_date=2023-07-04.txt',comment='#',skiprows=0,sep='\t').loc[1:]
df_riovista = pd.read_csv('nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00010=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2000-01-01&end_date=2023-07-04.txt',comment='#',skiprows=0,sep='\t').loc[1:]

year_out = 2013
time_out = np.arange(np.datetime64('%d-01-01' % year_out), np.datetime64('%d-01-01' % (year_out+1)), np.timedelta64(1,'D'))
day_out = (time_out - time_out[0])/np.timedelta64(1,'D')
nout = len(day_out)

df_out = pd.DataFrame(index=day_out)
df_out.index.name = 'Decimal Day (UTC)'

# munge the jersey data
for river in ['Jersey', 'RioVista']:

	if river=='Jersey':
		time = df_jersey['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') 
		temp = pd.to_numeric(df_jersey['15568_00010'].values, errors='coerce')
		temp1 = pd.to_numeric(df_jersey['164023_00010'].values, errors='coerce')

	elif river=='RioVista':
		time = df_riovista['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') 
		temp = pd.to_numeric(df_riovista['16194_00010'].values, errors='coerce')
		temp1 = pd.to_numeric(df_riovista['300439_00010'].values, errors='coerce')
	
	# fill in missing temp data with alternative temperature (from BGC project)
	ind = np.isnan(temp)
	temp[ind] = temp1[ind]

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

	df_out['%s Temp (oC)' % river] = temp_out

	# plot
	fig, ax = plt.subplots(figsize=(11,4))
	ax.plot(day, temp, '.', label='data from %s to %s' % (time[0],time[-1]))
	ax.plot(day_out, temp_out, label='climatological average')
	ax.legend()
	ax.set_ylabel('Temperature (oC)')
	ax.set_xlabel('Decimal Day (UTC)')
	ax.set_title(river)
	fig.savefig('%s.png' % river)


df_out.to_csv('Delta_Climatological_Temperature.csv')











