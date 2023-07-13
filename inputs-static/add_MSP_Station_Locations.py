import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

# read original observation points
gdf = gpd.read_file('observation-points_original.shp')

# read points to be added
df = pd.read_csv('2022_MSP_Station_Locations.csv')

# add those points
df.rename(columns={"Station Name": "name"}, inplace=True)
for i in range(len(df)):
	df.loc[i,'name'] = df.loc[i,'name'].replace(' ','_')
	df.loc[i,'geometry'] = Point(df.loc[i,'UTM_E'], df.loc[i,'UTM_N'])
gdf1 = gpd.GeoDataFrame(df[['name','geometry']])
gdf1.crs = gdf.crs
gdf = pd.concat([gdf,gdf1])

# save
gdf.to_file('observation-points.shp')

