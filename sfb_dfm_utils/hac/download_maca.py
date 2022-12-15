"""
Early development processing of MACA and LOCA datasets to synthesize 
heat flux parameters, following the methods of Vroom et al.

Since we can do reasonably well interpolating observations, this is
on hold for the moment.
"""

import requests

## 
base_args=dict(temporal='all',accept='netcdf',point='false')

region_args=dict(north=38.7345,
                 south=37.1642,
                 west=-123.4607,
                 east=-121.1311)


# See
# http://maca.northwestknowledge.net/data_portal.php
# for details, including how to access other periods.
base_url=('http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/'
          'macav2metdata_%s_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc' )

# [ (label_for_var_in_url, name for var in parameters, label for local use), ...]
maca_vars=[ ('huss','specific_humidity','huss'),
            ('rhsmax','relative_humidity','rhsmax'),
            ('rhsmin','relative_humidity','rhsmin'),
            ('rsds','surface_downwelling_shortwave_flux_in_air','rsds'),
            ('tasmin','air_temperature','tasmin'),
            ('tasmax','air_temperature','tasmax') ]

for label,var_name,name in maca_vars:
    save_fn="rcp45_2011_2015-%s.nc"%name

    if os.path.exists(save_fn):
        continue
    
    params=dict(base_args)
    params.update(region_args)
    params['var']=var_name
    response=requests.get( base_url%label, params=params, stream=True)
    
    # thanks you https://stackoverflow.com/questions/14114729/save-a-large-file-using-the-python-requests-library
    # Throw an error for bad status codes
    response.raise_for_status()

    with open(save_fn, 'wb') as fp:
        for block in response.iter_content(1024):
            fp.write(block)
##

# That effectively fetches these files:
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/macav2metdata_huss_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc?&var=specific_humidity&north=38.7345&south=37.1642&west=-123.4607&east=-121.1311&temporal=all&accept=netcdf&point=false 
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/macav2metdata_rhsmax_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc?&var=relative_humidity&north=38.7345&south=37.1642&west=-123.4607&east=-121.1311&temporal=all&accept=netcdf&point=false 
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/macav2metdata_rhsmin_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc?&var=relative_humidity&north=38.7345&south=37.1642&west=-123.4607&east=-121.1311&temporal=all&accept=netcdf&point=false 
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/macav2metdata_rsds_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc?&var=surface_downwelling_shortwave_flux_in_air&north=38.7345&south=37.1642&west=-123.4607&east=-121.1311&temporal=all&accept=netcdf&point=false 
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/macav2metdata_tasmin_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc?&var=air_temperature&north=38.7345&south=37.1642&west=-123.4607&east=-121.1311&temporal=all&accept=netcdf&point=false 
# http://thredds.northwestknowledge.net:8080/thredds/ncss/grid/MACAV2/bcc-csm1-1/macav2metdata_tasmax_bcc-csm1-1_r1i1p1_rcp45_2011_2015_CONUS_daily.nc?&var=air_temperature&north=38.7345&south=37.1642&west=-123.4607&east=-121.1311&temporal=all&accept=netcdf&point=false 

##

# Carry out the processing described in Julia Vroom's temperature paper
Patm_kPa=101.3

# Goal is a 5km gridded dataset of relative humidity, air temper, cloudiness, and wind.
# Wind we already have
# Relative humidity is derived from a combination of MACA daily average specific humidity
# and air temperature.
