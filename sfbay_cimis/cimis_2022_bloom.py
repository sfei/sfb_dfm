#!/usr/bin/env python
from __future__ import print_function
import time

import datetime
import requests

import numpy as np
import xarray as xr

from stompy import utils

## 

appKey='98b9b76d-2ce5-4210-b4ae-9274524ad28a'
url="http://et.water.ca.gov/api/data"

def cimis_json_to_xr(data):
    df=xr.Dataset()

    data2=data['Data']['Providers'][0]

    # Owner, Records, Type, Name
    df.attrs['data_owner'] = data2['Owner']
    df.attrs['station_type'] = data2['Type']
    df.attrs['service_name'] = data2['Name']

    records=data2['Records']
    # 768 records..

    if len(records)==0:
        print("No data")
        return None

    df['Date']= ( ('Date',), [ "%s %s"%(rec['Date'],rec['Hour'])
                               for rec in records] )

    def cnv(s):
        try:
            return float(s)
        except (ValueError,TypeError):
            return np.nan

    for field in ['HlyAirTmp', 'HlyEto','HlyNetRad',
                  'HlyPrecip','HlyRelHum','HlyResWind',
                  'HlySolRad','HlyWindDir','HlyWindSpd']:
        # use zip to transpose the lists
        qc,value = zip( *[ (rec[field]['Qc'],cnv(rec[field]['Value']))
                           for rec in records ] )
        df[field]=( ('Date',), np.array(value) )
        df[field+'_qc'] = ( ('Date',), np.array(qc) )

        df[field].attrs['units']=records[0][field]['Unit']

    df.attrs['station_num']=int(records[0]['Station'])

    return df

def cimis_fetch_station_metadata(station,df=None):
    df=df or xr.Dataset()
    station=171
    req=requests.get("http://et.water.ca.gov/api/station/%s"%station,
                     headers=dict(Accept='application/json'))
    station_meta=req.json()

    # add station metadata to attrs:
    stn=data['Stations'][0]

    df.attrs['elevation'] = float(stn['Elevation'])
    df.attrs['is_active'] = stn['IsActive']
    df.attrs['station_name']=stn['Name']
    lat=float(stn['HmsLatitude'].split('/')[1]) #  u"37\xba35'56N / 37.598758"
    lon=float(stn['HmsLongitude'].split('/')[1])

    df.attrs['latitude']=lat
    df.attrs['longitude']=lon
    return df


def cimis_fetch_to_xr(stations, # Union City
                      start_date,end_date,
                      fields=None,
                      station_meta=True):
    if fields is None:
        fields=['hly-air-tmp','hly-eto',
                'hly-net-rad','hly-precip',
                'hly-rel-hum','hly-res-wind',
                'hly-sol-rad','hly-wind-dir',
                'hly-wind-spd']

    if isinstance(stations,int) or isinstance(stations,str):
        stations=[stations]

    stations=[str(s) for s in stations] 

    req=requests.get(url,params=dict(appKey=appKey,
                                     targets=",".join(stations), 
                                     startDate=start_date,
                                     endDate=end_date,
                                     unitOfMeasure='M', # metric please
                                     dataItems=",".join( fields ) ))
    df=cimis_json_to_xr(req.json())

    if station_meta:
        if len(stations)!=1:
            print("Can only record station metadata for a single station")
        else:
            cimis_fetch_station_metadata(stations[0],df=df)

    return df

## 

## 

#2/5/2001 is start of record for union city
start_date = '2022-04-01'
end_date = '2022-10-28'

#period=[np.datetime64('2001-02-01'),
#        np.datetime64('2022-08-19')]


#period_dns=[utils.to_dnum(d) for d in period]
# on richmond, this returns period_dns = [730517.0, 738386.0]

span_days=30

# create lists of start and end dates to break data into chunks
ts = np.datetime64(start_date,'D')
te = np.datetime64(end_date,'D')
t1 = ts
t2 = t1 + np.timedelta64(span_days,'D')
time_windows = [(t1,t2)]
while t2<te:
    t1 = t2
    t2 = t1 + np.timedelta64(span_days,'D')
    if t2>te:
        t2=te
    time_windows.append((t1,t2))

dfs={}

### note when we run 
counter = 0
for time_window in time_windows:
    
    counter = counter + 1

    win_start = time_window[0]
    win_end = time_window[1]

    print("Fetching %s"%win_start)
    df=cimis_fetch_to_xr(171,win_start,win_end,station_meta=False)


    # on richmond fails with ConnectionError: HTTPConnectionPool(host='et.water.ca.gov', port=80): Max retries exceeded with url: /api/data?startDate=2001-02-01&targets=171&dataItems=hly-air-tmp%2Chly-eto%2Chly-net-rad%2Chly-precip%2Chly-rel-hum%2Chly-res-wind%2Chly-sol-rad%2Chly-wind-dir%2Chly-wind-spd&unitOfMeasure=M&appKey=98b9b76d-2ce5-4210-b4ae-9274524ad28a&endDate=2001-03-03 (Caused by NewConnectionError('<requests.packages.urllib3.connection.HTTPConnection object at 0x7f0ec0560eb8>: Failed to establish a new connection: [Errno -3] Temporary failure in name resolution',))

    dfs[counter]=df
    print("Delay...")
    time.sleep(1.5) # be kind to others
    print("Done")

## 

dnums=list(dfs.keys())
dnums.sort()


df_seq=[ dfs[dnum] for dnum in dnums ]

df_full = xr.concat( df_seq, dim='Date' )

# annoying convention where midnight is 2400.
dates=[ (datetime.datetime.strptime(s.replace('2400','0000'),
                                    '%Y-%m-%d %H%M')
         + datetime.timedelta(days=int(s.endswith('2400'))))
       for s in df_full.Date.values]

df_full['time']=( ('Date',), dates )

## 

#sel=utils.select_increasing(df.time.values)
#df_full=df_full.isel(Date=sel)

## 

df_full.to_netcdf('union_city-hourly-2022_bloom.nc')
