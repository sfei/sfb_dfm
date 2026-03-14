
import geopandas as gpd

import matplotlib.pylab as plt
from shapely.geometry import Point
import pandas as pd

import utm

# Ari at SFEI says we have wave data at Hayward from 01/2022-01/2024
# (the Hayward buoy is already in our list of observation points)

# Babak from Integral says Timeline of measurement is roughly Oct 2021-April 2022:

grid = gpd.read_file('sfei_v25_straightened_net_OUTLINE.shp')

obs = gpd.read_file('observation_points_v2.shp')

integral = gpd.read_file('SFBuoys_shp.shp')

integral_babak = {
"Integral_Edens_Landing": (565150.025, 4165573.093, '2021-10-01', '2022-04-30'),
"Integral_Hunters_Point": (554445.096, 4172810.370, '2021-10-01', '2022-04-30'),
"Integral_Richmond": (550271.981, 4199565.911, '2021-10-01', '2022-04-30'),
"Integral_Petaluma": (550062.468, 4208031.705, '2021-10-01', '2022-04-30')}

# Sam from Integral :

# The project consisted of three field deployments: Summer (July - August 2018), 
# Winter (January - February 2019), and Spring (April - May 2019). Each field 
# campaign involved deploying four platforms containing EFML instruments at 
# various locations in South San Francisco Bay. The data set includes measurements 
# taken by a Vectrino Profiler (boundary layer turbulence and sediment flux profiles), 
# acoustic Doppler velocimeters (turbulent momentum and sediment fluxes), acoustic 
# Doppler current profilers (mean current profiles), pressure loggers (wave statistics), 
# CTDs (salinity, temperature, pressure), and a LISST (suspended sediment particle size 
# distributions). We encourage the free use of this data set, provided that it is 
# properly cited. The platform IDs and locations for each deployment are as follows:

# Summer

 

# P1: 
x1, y1, dum1, dum2 = utm.from_latlon(37.58745,-122.18530)
s1 = '2018-07-01'
e1 = '2018-08-31'

# P2
x2, y2, dum1, dum2 = utm.from_latlon(37.58728 , -122.17167)
s2 = '2018-07-01'
e2 = '2018-08-31'

#P3: 
x3, y3, dum1, dum2 = utm.from_latlon(37.58550 , -122.23141) 
s3 = '2018-07-01'
e3 = '2018-08-31'

#P4a: 
x4a, y4a, dum1, dum2 = utm.from_latlon(37.58681 , -122.21182)
s4a = '2018-07-01'
e4a = '2018-08-31'

integral_summer_xy = {'Integral_Summer_P1' : (x1,y1,s1,e1), 'Integral_Summer_P2' : (x2,y2,s2,e2), 'Integral_Summer_P3' : (x3,y3,s3,e3), 'Integral_Summer_P4a' : (x4a,y4a,s4a,e4a)}
 
#Winter

#P1: 
x1, y1, dum1, dum2 = utm.from_latlon(37.58744, -122.18534)
s1 = '2019-01-01'
e1 = '2019-02-28'

#P2: 
x2, y2, dum1, dum2 = utm.from_latlon(37.58728, -122.17167)
s2 = '2019-01-01'
e2 = '2019-02-28'

#P3: 
x3, y3, dum1, dum2 = utm.from_latlon(37.58550, -122.23141)
s3 = '2019-01-01'
e3 = '2019-02-28'

#P4b: 
x4b, y4b, dum1, dum2 = utm.from_latlon(37.5613, -122.18530)
s4b = '2019-01-01'
e4b = '2019-02-28'

integral_winter_xy = {'Integral_Winter_P1' : (x1,y1,s1,e1), 'Integral_Winter_P2' : (x2,y2,s2,e2), 'Integral_Winter_P3' : (x3,y3,s3,e3), 'Integral_Winter_P4b' : (x4b,y4b,s4b,e4b)}

#Spring

#P1: 
x1, y1, dum1, dum2 = utm.from_latlon(37.58742 , -122.18555)
s1 = '2019-04-01' 
e1 = '2019-05-31'

#P2: 
x2, y2, dum1, dum2 = utm.from_latlon(37.58728 , -122.17167)
s2 = '2019-04-01'
e2 = '2019-05-31'

#P3: 
x3, y3, dum1, dum2 = utm.from_latlon(37.58550 , -122.23141)
s3 = '2019-04-01'
e3 = '2019-05-31'

#P4b: 
x4b, y4b, dum1, dum2 = utm.from_latlon(37.56128 , -122.18537)
s4b = '2019-04-01'
e4b = '2019-05-31'

integral_spring_xy = {'Integral_Spring_P1' : (x1,y1,s1,e1), 'Integral_Spring_P2' : (x2,y2,s2,e2), 'Integral_Spring_P3' : (x3,y3,s3,e3), 'Integral_Spring_P4b' : (x4b,y4b,s4b,e4b)}

usgs_latlon = {'CHC13CRC': (38.010787, -122.492974, '2013-12-03', '2014-02-25'),
 'CHC13CRD': (38.008795, -122.487951, '2013-12-02', '2014-02-25'),
 'CHC13MFW': (38.01618, -122.48846, '2013-12-04', '2014-02-11'),
 'CHC13MFE': (38.01214, -122.48266, '2013-12-04', '2014-02-11'),
 'CHC13CHM': (38.0542, -122.37591, '2013-12-04', '2014-02-25'),
 'CHC13CHT': (38.03522, -122.37678, '2013-12-04', '2014-02-25'),
 'CHC13M2M': (38.04094, -122.44654, '2013-12-03', '2014-02-24'),
 'CHC13S1T': (38.01332, -122.46468, '2013-12-03', '2014-02-26'),
 'CHC13N1T': (38.07596, -122.44707, '2013-12-04', '2014-02-26'),
 'CHC13M1T': (38.04575, -122.46277, '2013-12-03', '2014-02-25'),
 'CHC13M2T': (38.04187, -122.44679, '2013-12-03', '2014-02-25'),
 'CHC13TCD': (38.00989, -122.48476, '2013-12-04', '2014-02-11'),
 'CHC13TCC': (38.01401, -122.49005, '2013-12-04', '2014-02-11'),
 'CHC14DX1': (38.00946, -122.48347, '2015-01-15', '2015-01-28'),
 'CHC14CRD': (38.00922, -122.48609, '2014-12-10', '2015-01-28'),
 'CHC14DX5': (38.00901, -122.48414, '2015-01-15', '2015-01-28'),
 'CHC14PD5': (38.009, -122.48412, '2014-12-19', '2014-12-29'),
 'CHC14PD4': (38.00919, -122.48391, '2014-12-19', '2014-12-29'),
 'CHC14DX2': (38.00939, -122.48354, '2015-01-15', '2015-01-28'),
 'CHC14DX3': (38.0093, -122.48371, '2015-01-15', '2015-01-28'),
 'CHC14PD1': (38.00946, -122.48347, '2014-12-19', '2014-12-29'),
 'CHC14PD3': (38.0093, -122.48371, '2014-12-19', '2014-12-29'),
 'CHC14MFD': (38.00988, -122.48307, '2014-12-03', '2015-02-03'),
 'CHC14SPD': (38.00959, -122.48327, '2014-12-03', '2014-12-29'),
 'CHC14CHT': (38.03533172607422,
  -122.37660217285156,
  '2014-12-02',
  '2015-02-02'),
 'CHC14M1T': (38.04570007324219,
  -122.46263885498047,
  '2014-12-02',
  '2015-02-02'),
 'CHC14M2T': (38.046021, -122.447121, '2014-12-02', '2015-04-13'),
 'CHC14N1T': (38.07622, -122.44732, '2014-12-03', '2015-02-02'),
 'CHC14S1T': (38.01326, -122.46476, '2014-12-03', '2015-02-02'),
 'CHC14TDH': (38.00922, -122.48609, '2014-12-05', '2015-01-28'),
 'CHC14TDL': (38.00981, -122.48486, '2014-12-05', '2015-01-28'),
 'CHC16DX1': (38.0095, -122.48349, '2016-05-31', '2016-06-09'),
 'CHC16DX3': (38.00936, -122.48366, '2016-05-31', '2016-06-09'),
 'CHC16DX4': (38.00919, -122.48394, '2016-05-31', '2016-06-09'),
 'CHC16DXA': (38.00941, -122.4836, '2016-05-31', '2016-06-09'),
 'CHC16MFD': (38.00982, -122.48302, '2016-05-05', '2016-06-23'),
 'CHC16SPD': (38.01204, -122.47959, '2016-05-05', '2016-06-23'),
 'CHC16DX2': (38.00948, -122.48352, '2016-05-31', '2016-06-09'),
 'CHC16CRD': (38.008842, -122.487946, '2016-05-07', '2016-06-22'),
 'CHC16CHT': (38.035191, -122.375191, '2016-05-05', '2016-06-21'),
 'CHC16MFT': (38.01204, -122.47959, '2016-05-04', '2016-06-21'),
 'CHC16M2M': (38.046379, -122.447327, '2016-05-04', '2016-06-21'),
 'CHC16M2T': (38.04607, -122.44705, '2016-05-04', '2016-06-22'),
 'CHC16MLT': (38.04618835449219,
  -122.44737243652344,
  '2016-05-04',
  '2016-06-22'),
 'CHC16N1T': (38.07616, -122.44709, '2016-05-04', '2016-06-21'),
 'CHC16S1T': (38.013351, -122.464447, '2016-05-04', '2016-06-21'),
 'CHC16LSS': (38.00988, -122.484772, '2016-06-21', '2016-06-21'),
 'CHC16TDL': (38.0098, -122.48489, '2016-05-05', '2016-06-21'),
 'SPA14N2M': (38.055241, -122.432137, '2014-02-27', '2014-06-23'),
 'SPA14N2T': (38.05524, -122.43214, '2014-02-27', '2014-03-15'),
 'SPA14S1T': (38.013321, -122.464371, '2014-02-27', '2014-06-23'),
 'SPA14N1T': (38.076069, -122.447006, '2014-02-27', '2014-06-23'),
 'SPB14N1T': (38.076271, -122.447304, '2014-06-24', '2014-09-16'),
 'SPB14N2M': (38.050091, -122.429939, '2014-06-24', '2014-09-16'),
 'SPB14N2T': (38.05002, -122.43033, '2014-06-24', '2014-09-16'),
 'SPB14S1T': (38.01334, -122.464691, '2014-06-24', '2014-09-16'),
 'SPC14N1T': (38.076012, -122.446739, '2014-09-18', '2014-12-01'),
 'SPC14N2M': (38.056019, -122.431358, '2014-09-18', '2014-12-01'),
 'SPC14N2T': (38.056019, -122.431358, '2014-09-18', '2014-12-01'),
 'SPC14S1T': (38.01329, -122.46476, '2014-09-18', '2014-12-01'),
 'SPD15N1T': (38.07618, -122.44728, '2015-02-04', '2015-04-04'),
 'SPD15N2T': (38.04998, -122.430367, '2015-02-04', '2015-04-14'),
 'SPD15S1T': (38.013271, -122.464729, '2015-02-04', '2015-04-14'),
 'ERO19PCT': (38.034271, -122.378571, '2019-06-11', '2019-08-19'),
 'ERO19PBM': (38.08363, -122.391098, '2019-06-11', '2019-08-19'),
 'ERO19PHT': (38.08323, -122.39112, '2019-06-11', '2019-08-19'),
 'ERO19PST': (38.07914, -122.393059, '2019-06-11', '2019-08-19'),
 'ERO19GBM': (38.116711, -122.040718, '2019-06-11', '2019-08-20'),
 'ERO19GHT': (38.11639, -122.0408, '2019-06-11', '2019-08-20'),
 'ERO19GST': (38.114189, -122.045982, '2019-06-11', '2019-08-20')}


usgs_xy = {}
for key in usgs_latlon.keys():

	x,y,dum1,dum2=utm.from_latlon(usgs_latlon[key][0],usgs_latlon[key][1])
	usgs_xy[key] = (x,y,usgs_latlon[key][2],usgs_latlon[key][3])

# compile all the new points in a list and use to create a dataframe
point_list = []
name_list = []
start_list = []
stop_list = []
for dict in [usgs_xy, integral_spring_xy, integral_winter_xy, integral_summer_xy, integral_babak]:
	for key in dict.keys():
		name_list.append(key)
		point_list.append(Point((dict[key][0],dict[key][1])))
		start_list.append(dict[key][2])
		stop_list.append(dict[key][3])

obs_plus = gpd.GeoDataFrame(geometry=point_list)
obs_plus['name'] = name_list
obs_plus['start_date'] = start_list
obs_plus['stop_date'] = stop_list

obs_new = pd.concat([obs,obs_plus],ignore_index=True)

ind = obs_new['name']=='Hayward'
obs_new.loc[ind,'start_date'] = '2022-01-01'
obs_new.loc[ind,'stop_date'] = '2024-01-01'


obs_new.to_file('observation-points.shp')


