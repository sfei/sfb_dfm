sfb_dfm.py automatically downloads ocean surface temperatures from the point reyes buoy, the same one it uses for tidal water level boundary condition

however, sometimes there are gaps in the temperature data so in this folder
we curate backup data from the NDBC san francisco buoy, and sfb_dfm.py uses
this to fill in any gaps

get more data here:
https://www.ndbc.noaa.gov/station_history.php?station=46026