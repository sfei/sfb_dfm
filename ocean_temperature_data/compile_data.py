import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import os, sys

maxyear=2024

# compile SF Data
for iy,year in enumerate(range(2000,maxyear)):

	fn=os.path.join('data','46026h%d.txt' % year)
	
	with open(fn,'r') as f:
		line1 = f.readline()
	
	if line1[0]=='Y':
		df = pd.read_csv(fn,delim_whitespace=True)
		df['mm'] = 0
	elif line1[0]=='#':
		names=['YYYY','MM','DD','hh','mm','WD','WSPD','GST','WVHT','DPD','APD','MWD','BAR','ATMP','WTMP','DEWP','VIS','TIDE']
		df = pd.read_csv(fn,delim_whitespace=True,comment='#',names=names)
	else:
		raise Exception('header not recognized')

	nrows = len(df)
	datetime=np.zeros(nrows,dtype='datetime64[ns]')
	for i in range(nrows):
		row=df.iloc[i]
		datetime[i] = np.datetime64('%04d-%02d-%02d %02d:%02d' % 
			(row['YYYY'],row['MM'],row['DD'],row['hh'],row['mm']))
	datetime_UTC_1 = datetime + np.timedelta64(8,'h') # time is in PST/PDT, don't worry about hour
	temp_oC_1 = df['WTMP'].values
	if iy==0:
		datetime_UTC_SF = datetime_UTC_1.copy()
		temp_oC_SF = temp_oC_1.copy()
	else:
		datetime_UTC_SF = np.concatenate((datetime_UTC_SF,datetime_UTC_1))
		temp_oC_SF = np.concatenate((temp_oC_SF,temp_oC_1)) 

# eliminate 999
ind = temp_oC_SF<900.
temp_oC_SF = temp_oC_SF[ind]
datetime_UTC_SF = datetime_UTC_SF[ind]


# compile point reyes Data
for iy,year in enumerate(range(2005,maxyear)):

	# buoy must have been down this year, no file available
	if year==2006:
		continue

	fn=os.path.join('data','pryc1h%d.txt' % year)
		
	with open(fn,'r') as f:
		line1 = f.readline()
	
	if line1[0]=='Y':
		df = pd.read_csv(fn,delim_whitespace=True)
		df['mm'] = 0
	elif line1[0]=='#':
		names=['YYYY','MM','DD','hh','mm','WD','WSPD','GST','WVHT','DPD','APD','MWD','BAR','ATMP','WTMP','DEWP','VIS','TIDE']
		df = pd.read_csv(fn,delim_whitespace=True,comment='#',names=names)
	else:
		raise Exception('header not recognized')

	nrows = len(df)
	datetime=np.zeros(nrows,dtype='datetime64[ns]')
	for i in range(nrows):
		row=df.iloc[i]
		datetime[i] = np.datetime64('%04d-%02d-%02d %02d:%02d' % 
			(row['YYYY'],row['MM'],row['DD'],row['hh'],row['mm']))
	datetime_UTC_1 = datetime + np.timedelta64(8,'h') # time is in PST/PDT, don't worry about hour
	temp_oC_1 = df['WTMP'].values
	if iy==0:
		datetime_UTC_PR = datetime_UTC_1.copy()
		temp_oC_PR = temp_oC_1.copy()
	else:
		datetime_UTC_PR = np.concatenate((datetime_UTC_PR,datetime_UTC_1))
		temp_oC_PR = np.concatenate((temp_oC_PR,temp_oC_1)) 

# eliminate 999
ind = temp_oC_PR<900.
temp_oC_PR = temp_oC_PR[ind]
datetime_UTC_PR = datetime_UTC_PR[ind]

# now create a single time axis and take daily average temperature, leaving nan where none available
datetime_UTC_daily = np.arange('2000-01-01','%d-01-01' % (maxyear+1),dtype='datetime64[D]')
ndays = len(datetime_UTC_daily)
temp_oC_PR_daily = np.nan*np.ones(ndays)
temp_oC_SF_daily = np.nan*np.ones(ndays)
for iday in range(ndays):
	ind = np.logical_and(datetime_UTC_PR>=datetime_UTC_daily[iday]-np.timedelta64(12,'h'),
		                 datetime_UTC_PR<=datetime_UTC_daily[iday]+np.timedelta64(12,'h'))
	temp_oC_PR_daily[iday] = np.mean(temp_oC_PR[ind])
	ind = np.logical_and(datetime_UTC_SF>=datetime_UTC_daily[iday]-np.timedelta64(12,'h'),
		                 datetime_UTC_SF<=datetime_UTC_daily[iday]+np.timedelta64(12,'h'))
	temp_oC_SF_daily[iday] = np.mean(temp_oC_SF[ind])

# combine temperature time series, preferring SF over Point Reyes but using Point Reyes 
# to fill in the gaps, get rid of nan data
temp_oC_combined_daily = temp_oC_SF_daily
ind = np.isnan(temp_oC_SF_daily)
temp_oC_combined_daily[ind] = temp_oC_PR_daily[ind]
ind = ~np.isnan(temp_oC_combined_daily)
temp_oC_combined_daily = temp_oC_combined_daily[ind]
datetime_UTC_combined_daily = datetime_UTC_daily[ind]

# plot
fig, ax = plt.subplots(3,1,figsize=(11,8.5),constrained_layout=True)
ax[0].plot(datetime_UTC_daily, temp_oC_SF_daily, label='San Francisco (46026)')
ax[1].plot(datetime_UTC_daily, temp_oC_PR_daily, label='Point Reyes (pryc1)')
ax[2].plot(datetime_UTC_combined_daily, temp_oC_combined_daily, label='Combined')
ax[0].set_title('Use San Francisco temps when available, fill gaps with Point Reyes')
for ax1 in ax:
	ax1.legend()
	ax1.set_xlim((datetime_UTC_daily[0],datetime_UTC_daily[-1]))
fig.savefig('ocean_temperature.png')

# output to dataframe
df_all = pd.DataFrame()
df_all['datetime_UTC'] = datetime_UTC_combined_daily
df_all['temp_oC'] = temp_oC_combined_daily

df_all.to_csv('ocean_temperature.csv')

