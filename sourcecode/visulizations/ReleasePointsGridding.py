import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color
import pandas as pd
import h3
from shapely.geometry.polygon import Polygon
from matplotlib.ticker import EngFormatter


home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()
print(mask_ds)
# get the corner points to plot on the map
x = mask_ds['glamf']
print(x[0])
y = mask_ds['gphif']
print(y[0])

# get the mask values of the corner points
c = mask_ds['tmask'][:]

# 1. seed points from the shapefiles and station locations
seed_points = pd.read_csv(home_folder + 'Nemo_H3Release_LatLon_Res5.csv')
lats = seed_points['Latitudes']
lons = seed_points['Longitudes']
res_id = seed_points['Res3_HexId']

uni_lats = lats.unique()
uni_lons = lons.unique()
print(len(seed_points['Latitudes']), len(uni_lats), len(uni_lons))

# get points in a specific region
gom_lats_index = np.where(np.logical_and(lats < 40, lats > 5))
final_res3 = res_id.loc[gom_lats_index].unique()


fig = plt.figure(dpi=150)
ax = plt.axes()
plt.tick_params(axis='both', which='major', labelsize=10)
colormap = color.ListedColormap(['gainsboro', 'white'])

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

ax.scatter(lons, lats, s=0.1, c='blue')
ax.yaxis.set_major_formatter(EngFormatter(unit=u"°N"))
ax.xaxis.set_major_formatter(EngFormatter(unit=u"°W"))


def plot_hex(hex, c):
    polygons = h3.h3_set_to_multi_polygon([hex], geo_json=False)
    p = Polygon(polygons[0][0])
    y, x = p.exterior.xy
    ax.plot(x, y, color=c, linewidth=0.5)


for h in final_res3:
    plot_hex(h, 'red')

plt.show()