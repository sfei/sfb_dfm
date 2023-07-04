import sys, os
sys.path.append('pytides2-0.0.5')
from datetime import datetime
from pytides2.tide import Tide
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

''' This script creates a netcdf file with the Delta inflows, filled in using pytides where data are not available
The user must periodically update the data in the data folders and re-run this script -- the more recent sfb_dfm 
will get the delta inflow data from this repository'''

# specify a time series to interpolate on
time_out = np.arange(np.datetime64('2000-01-01'), np.datetime64('2023-01-01'), np.timedelta64(15,'m')).astype('datetime64[ns]')

# read the dayflow data
df_dayflow = pd.read_csv('data_dayflow/dayflow-results-1997-2022.csv')

# add 12 hours to dayflow time to center in the middle of the day, then add 8 hours to convert from PST to UTC
timet = pd.DatetimeIndex(df_dayflow['DATE']).values + np.timedelta64(12,'h') + np.timedelta64(8,'h') 

# create a pandas dataframe to store the resulting predicted flows
df_out = pd.DataFrame(index=time_out)
df_out.index.name = 'time (UTC)'

# loop through 2 rivers
for river in ['riovista','jersey']:

	# some params
	if river=='jersey':
		training_window = [np.datetime64('2001-03-01'), np.datetime64('2002-03-01')]
		flowt = df_dayflow['WEST'].values  * 0.0283168
		TS_param = '15560_00060'
		fn_list = ['nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2000-01-01&end_date=2004-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2005-01-01&end_date=2009-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2010-01-01&end_date=2014-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2015-01-01&end_date=2019-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11337190&legacy=1&period=&begin_date=2020-01-01&end_date=2022-12-31.txt']
	elif river=='riovista':
		training_window = [np.datetime64('2010-01-21'),np.datetime64('2013-11-26')]
		flowt = df_dayflow['RIO'].values  * 0.0283168
		TS_param = '16184_00060'
		fn_list = ['nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2000-01-01&end_date=2004-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2005-01-01&end_date=2009-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2010-01-01&end_date=2014-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2015-01-01&end_date=2019-12-31.txt',
				   'nwis.waterdata.usgs.gov_usa_nwis_uv__cb_00060=on&format=rdb&site_no=11455420&legacy=1&period=&begin_date=2020-01-01&end_date=2022-12-31.txt']


	# import decade long chunk of NWIS data 
	for ifn, fn in enumerate(fn_list):
		data1 = pd.read_csv(os.path.join('data_%s' % river, fn),comment='#',skiprows=0,sep='\t').loc[1:]
		if ifn==0:
			data = data1.copy()
		else:
			data = pd.concat([data,data1])
	time = data['datetime'].values.astype('datetime64[ns]') + np.timedelta64(8,'h') # convert PST to UTC
	flow = pd.to_numeric(data[TS_param].values, errors='coerce') * 0.0283168

	# filter out the nan data, leaving only times with good data
	ind = ~np.isnan(flow)
	time = time[ind]
	flow = flow[ind]

	## inspect the dataset to find the longest stretch without significant data gaps
	#dtmin = np.diff(time) / np.timedelta64(1,'m')
	#fig,ax = plt.subplots()
	#ax.semilogy(time[0:-1], dtmin, '.')
	#ax.grid()
	#ax.set_ylabel('time between adjacent flow measurements (min)')
	
	# for training, use the biggest continuous chunk of data available
	ind = np.logical_and(time>=training_window[0],time<=training_window[1])
	time_train = time[ind]
	flow_train = flow[ind]
	
	# interpolate dayflow flow to compute "fluctuating" componeent of the flow, for the training dataset
	tref = time_train[0]
	dtref = np.timedelta64(1,'m')
	flowt_train_i = np.interp((time_train-tref)/dtref, (timet-tref)/dtref, flowt)
	flowf_train = flow_train - flowt_train_i
	
	# interpolate dayflow flow onto the output time axis
	flowt_out = np.interp((time_out-tref)/dtref, (timet-tref)/dtref, flowt)
	
	# also get the fluctuating component of the flow for the full dataset
	flowt_i = np.interp((time-tref)/dtref, (timet-tref)/dtref, flowt)
	flowf = flow - flowt_i
	
	# convert all our time axes to datetime
	datetime_train = np.array([pd.Timestamp(t).to_pydatetime() for t in time_train])
	datetime_out = np.array([pd.Timestamp(t).to_pydatetime() for t in time_out])
	
	##Fit the tidal data to the harmonic model using Pytides
	my_tide = Tide.decompose(flowf_train, datetime_train)
	##Predict the tides using the Pytides model.
	flowf_out = my_tide.at(datetime_out)
	flow_out = flowt_out + flowf_out
	
	# now, compare the reconstructed time series to ALL the available data by interpolating onto times where data area available
	flow_predicted = np.interp((time-tref)/dtref, (time_out-tref)/dtref, flow_out)
	
	MS = 3 # marker size for plot
	
	# compute model skill, bias, rmse, r2
	num = np.sum( (flow_predicted - flow)**2 )
	den = np.sum( (np.abs(flow_predicted - flow.mean()) + np.abs(flow - flow.mean()))**2 )
	skill = 1 - num / den
	r2 = np.corrcoef(flow_predicted, flow)[0,1]**2
	bias = np.mean(flow_predicted - flow)
	rmse = np.sqrt( np.mean( (flow_predicted - flow)**2 ) )
	legstr = 'bias=%0.2f m3/s\nrmse=%0.2f m3/s\nr2=%0.3f\nskill=%0.3f' % (bias,rmse,r2,skill)
	
	fig1, ax1 = plt.subplots(1,2,gridspec_kw={'width_ratios': [4, 1]},figsize=(24,8))
	ax1[0].plot(time_out, flow_out,'.',markersize=MS, label='predicted')
	ax1[0].plot(time, flow, '.',markersize=MS, label='measured')
	ax1[0].plot(time_train, flow_train, '.', markersize=MS, label='training data')
	ax1[0].set_ylabel('flow (m3/s)')
	ax1[0].legend()
	ax1[0].set_title(river)
	ax1[1].plot(flow_predicted, flow, '.')
	ax1[1].grid()
	ax1[1].set_xlabel('predicted flow (m3/s)')
	ax1[1].set_ylabel('measured flow (m3/s)')
	ax1[1].plot([-4500,11500],[-4500,11500], label=legstr)
	ax1[1].legend()
	fig1.savefig('measured_vs_pytides_%s.png' % river)

	# add to dataframe
	df_out['%s flow (cms)' % river] = flow_out

df_out.to_csv('delta_inflows_from_pytides.csv')



	
	

	
	