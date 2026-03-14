
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
"Integral_Edens_Landing": (565150.025, 4165573.093),
"Integral_Hunters_Point": (554445.096, 4172810.370),
"Integral_Richmond": (550271.981, 4199565.911),
"Integral_Petaluma": (550062.468, 4208031.705)}

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

# P2
x2, y2, dum1, dum2 = utm.from_latlon(37.58728 , -122.17167)

#P3: 
x3, y3, dum1, dum2 = utm.from_latlon(37.58550 , -122.23141)

#P4a: 
x4a, y4a, dum1, dum2 = utm.from_latlon(37.58681 , -122.21182)

integral_summer_xy = {'Integral_Summer_P1' : (x1,y1), 'Integral_Summer_P2' : (x2,y2), 'Integral_Summer_P3' : (x3,y3), 'Integral_Summer_P4a' : (x4a,y4a)}
 
#Winter

#P1: 
x1, y1, dum1, dum2 = utm.from_latlon(37.58744, -122.18534)

#P2: 
x2, y2, dum1, dum2 = utm.from_latlon(37.58728, -122.17167)

#P3: 
x3, y3, dum1, dum2 = utm.from_latlon(37.58550, -122.23141)

#P4b: 
x4b, y4b, dum1, dum2 = utm.from_latlon(37.5613, -122.18530)

integral_winter_xy = {'Integral_Winter_P1' : (x1,y1), 'Integral_Winter_P2' : (x2,y2), 'Integral_Winter_P3' : (x3,y3), 'Integral_Winter_P4b' : (x4b,y4b)}

#Spring

#P1: 
x1, y1, dum1, dum2 = utm.from_latlon(37.58742 , -122.18555)

#P2: 
x2, y2, dum1, dum2 = utm.from_latlon(37.58728 , -122.17167)

#P3: 
x3, y3, dum1, dum2 = utm.from_latlon(37.58550 , -122.23141)

#P4b: 
x4b, y4b, dum1, dum2 = utm.from_latlon(37.56128 , -122.18537)

integral_spring_xy = {'Integral_Spring_P1' : (x1,y1), 'Integral_Spring_P2' : (x2,y2), 'Integral_Spring_P3' : (x3,y3), 'Integral_Spring_P4b' : (x4b,y4b)}

usgs_latlon = {'CHC13CRC' : (38.010787, -122.492974),
'CHC13CRD' : (38.008795, -122.487951),
'CHC13MFW' : (38.01618, -122.48846),
'CHC13MFE' : (38.01214, -122.48266),
'CHC13CHM' : (38.0542, -122.37591),
'CHC13CHT' : (38.03522, -122.37678),
'CHC13M2M' : (38.04094, -122.44654),
'CHC13S1T' : (38.01332, -122.46468),
'CHC13N1T' : (38.07596, -122.44707),
'CHC13M1T' : (38.04575, -122.46277),
'CHC13M2T' : (38.04187, -122.44679),
'CHC13TCD' : (38.00989, -122.48476),
'CHC13TCC' : (38.01401, -122.49005),
'CHC14DX1' : (38.00946, -122.48347),
'CHC14CRD' : (38.00922, -122.48609),
'CHC14DX5' : (38.00901, -122.48414),
'CHC14PD5' : (38.009, -122.48412),
'CHC14PD4' : (38.00919, -122.48391),
'CHC14DX2' : (38.00939, -122.48354),
'CHC14DX3' : (38.0093, -122.48371),
'CHC14PD1' : (38.00946, -122.48347),
'CHC14PD3' : (38.0093, -122.48371),
'CHC14MFD' : (38.00988, -122.48307),
'CHC14SPD' : (38.00959, -122.48327),
'CHC14CHT' : (38.03533172607422, -122.37660217285156),
'CHC14M1T' : (38.04570007324219, -122.46263885498047),
'CHC14M2T' : (38.046021, -122.447121),
'CHC14N1T' : (38.076220, -122.447320),
'CHC14S1T' : (38.013260, -122.464760),
'CHC14TDH' : (38.009220, -122.486090),
'CHC14TDL' : (38.009810, -122.484860),
'CHC16DX1' : (38.009500, -122.483490),
'CHC16DX3' : (38.009360, -122.483660),
'CHC16DX4' : (38.009190, -122.483940),
'CHC16DXA' : (38.009410, -122.483600),
'CHC16MFD' : (38.009820, -122.483020),
'CHC16SPD' : (38.012040, -122.479590),
'CHC16DX2' : (38.009480, -122.483520),
'CHC16CRD' : (38.008842, -122.487946),
'CHC16CHT' : (38.035191, -122.375191),
'CHC16MFT' : (38.012040, -122.479590),
'CHC16M2M' : (38.046379, -122.447327),
'CHC16M2T' : (38.046070, -122.447050),
'CHC16MLT' : (38.04618835449219, -122.44737243652344),
'CHC16N1T' : (38.076160, -122.447090),
'CHC16S1T' : (38.013351, -122.464447),
'CHC16LSS' : (38.009880, -122.484772),
'CHC16TDL' : (38.009800, -122.484890),
'SPA14N2M' : (38.055241, -122.432137),
'SPA14N2T' : (38.055240, -122.432140),
'SPA14S1T' : (38.013321, -122.464371),
'SPA14N1T' : (38.076069, -122.447006),
'SPB14N1T' : (38.076271, -122.447304),
'SPB14N2M' : (38.050091, -122.429939),
'SPB14N2T' : (38.050020, -122.430330),
'SPB14S1T' : (38.013340, -122.464691),
'SPC14N1T' : (38.076012, -122.446739),
'SPC14N2M' : (38.056019, -122.431358),
'SPC14N2T' : (38.056019, -122.431358),
'SPC14S1T' : (38.013290, -122.464760),
'SPD15N1T' : (38.076180, -122.447280),
'SPD15N2T' : (38.049980, -122.430367),
'SPD15S1T' : (38.013271, -122.464729),
'ERO19PCT' : (38.034271, -122.378571),
'ERO19PBM' : (38.083630, -122.391098),
'ERO19PHT' : (38.083230, -122.391120),
'ERO19PST' : (38.079140, -122.393059),
'ERO19GBM' : (38.116711, -122.040718),
'ERO19GHT' : (38.116390, -122.040800),
'ERO19GST' : (38.114189, -122.045982)}

usgs_xy = {}
for key in usgs_latlon.keys():

	x,y,dum1,dum2=utm.from_latlon(usgs_latlon[key][0],usgs_latlon[key][1])
	usgs_xy[key] = (x,y)

# compile all the new points in a list and use to create a dataframe
point_list = []
name_list = []
for dict in [usgs_xy, integral_spring_xy, integral_winter_xy, integral_summer_xy, integral_babak]:
	for key in dict.keys():
		name_list.append(key)
		point_list.append(Point(dict[key]))

obs_plus = gpd.GeoDataFrame(geometry=point_list)
obs_plus['name'] = name_list

obs_new = pd.concat([obs,obs_plus],ignore_index=True)


obs_new.to_file('observation-points.shp')


