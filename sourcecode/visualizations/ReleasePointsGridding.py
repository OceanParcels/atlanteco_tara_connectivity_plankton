import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color
import h3
from shapely.geometry.polygon import Polygon
from matplotlib.ticker import EngFormatter
from sourcecode.core import connectivityplots as cp
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task12/'

# 1. seed points from the shapefiles and station locations
seed_points = np.load(home_folder + 'H3_Res5_release_points.npz')
lats = seed_points['Latitude']
lons = seed_points['Longitude']
child_res = 5
parent_res = 3
res_id = pd.Series([h3.h3_to_parent(h3.geo_to_h3(lat, lon, child_res), parent_res) for lat, lon in
                    zip(seed_points['Latitude'], seed_points['Longitude'])])

# get points in a specific region
gom_lats_index = np.where(np.logical_and(lats < 31, lats > 15))[0]
final_res3 = res_id.loc[gom_lats_index].values
custom_size = 15
fig = plt.figure(dpi=300, figsize=(16, 7))
# ax = plt.axes()
ax = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND, color='gainsboro')
ax.add_feature(cfeature.BORDERS, linewidth=0.1)
ax.add_feature(cfeature.COASTLINE, linewidth=0.1)
gl = ax.gridlines(draw_labels=True)
gl.xlines = False
gl.ylines = False
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': custom_size, 'color': 'k'}
gl.ylabel_style = {'size': custom_size, 'color': 'k'}
gl.xlocator = mticker.FixedLocator([-85, -90, -95, -100])
gl.ylocator = mticker.FixedLocator([20, 25, 30])
# x, y, c = cp.load_mask_file(home_folder + 'GLOB16L98_mesh_mask_atlantic.nc')
# plt.tick_params(axis='both', which='major', labelsize=10)
# colormap = color.ListedColormap(['gainsboro', 'white'])
#
# # remove the first row and first column from the glamf/gphif to access points enclosed in the center
# ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

ax.scatter(lons[gom_lats_index], lats[gom_lats_index], s=0.2, c='blue')
ax.yaxis.set_major_formatter(EngFormatter(unit=u"°"))
ax.xaxis.set_major_formatter(EngFormatter(unit=u"°"))
ax.set_xlim(-100, -80)
ax.set_ylim(16, 31)


def plot_hex(hex, c):
    polygons = h3.h3_set_to_multi_polygon([hex], geo_json=False)
    p = Polygon(polygons[0][0])
    y, x = p.exterior.xy
    ax.plot(x, y, color=c, linewidth=0.5)


for h in final_res3:
    plot_hex(h, 'red')

# plt.show()
plt.savefig(home_folder + 'GulfofMexico_gridding_new.pdf', bbox_inches='tight',
            pad_inches=0.2)
