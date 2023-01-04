import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.colors as clr
import h3
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task12/'
stations = pd.read_csv(home_folder + 'AtlanticStations.csv', header=0)
lon = stations['Longitude']
lat = stations['Latitude']
code = stations['Station']

model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]
custom_size = 10
fig = plt.figure(figsize=(8, 6), dpi=500)
# ax = plt.axes()
ax = plt.axes(projection=ccrs.PlateCarree())
# ax.add_feature(cfeature.LAND, color='gainsboro')
# ax.add_feature(cfeature.COASTLINE, linewidth=0.1)
gl = ax.gridlines(draw_labels=True)
gl.xlines = False
gl.ylines = False
gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': custom_size, 'color': 'k'}
gl.ylabel_style = {'size': custom_size, 'color': 'k'}

colormap = clr.ListedColormap(['whitesmoke', 'lightskyblue'])
plt.tick_params(axis='both', which='major', labelsize=10)
# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap, label='Ocean model domain')

seed_points = np.load(home_folder + 'H3_Res5_release_points.npz')
release_lats = seed_points['Latitude']
release_lons = seed_points['Longitude']

ax.scatter(release_lons, release_lats, c='gold', s=0.1, alpha=0.7)
st = ax.scatter(lon, lat, c='r', s=5, label='Sample stations')

custom_lines = [Line2D([0], [0], color='lightskyblue', lw=4),
                Line2D([0], [0], color='gold', lw=4)]
plt.legend(custom_lines, ['Ocean model domain', 'Particles released'], loc='lower right')
# https://matplotlib.org/3.3.3/tutorials/intermediate/legend_guide.html#multiple-legends-on-the-same-axes

for i in range(len(lon)):
    xy = (lon[i] + 0.5, lat[i])
    ax.annotate('%s' % code[i], xy=xy, textcoords='data', bbox=dict(boxstyle='square,pad=5', fc='none', ec='none'))

# plt.show()

print('saving file')
plt.savefig(home_folder + "StationsDomain.jpeg", bbox_inches='tight',
            pad_inches=0.5)
