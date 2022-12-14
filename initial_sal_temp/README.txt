October 2022 Allie wrote these scripts to make initial salinity and temperature initial conditions 
from USGS Peterson cruise data...

How to use this stuff?

1. Hand download USGS Peterson cruise data *.csv file from here: https://sfbay.wr.usgs.gov/water-quality-database/
and put it in this folder. It will be called wqdata.csv. Make sure it is from the year matching the start date for the 
simulation

2. sfb_dfm.py should set things up automaticaly from there, should complain if you forgot to download wqdata.csv from the right time period
