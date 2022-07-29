import h3
import numpy as np
import xarray as xr
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.colors as color

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task10P/'
seed_point_resolution = 5

# For test region: Gulf of Mexico
# geoJSON = {'type': 'Polygon',
#            'coordinates': [[[0, -98], [0, -70], [30, -70], [30, -98]]]}
geoJSON = {'type': 'Polygon',
           'coordinates': [[[-90, -100], [-90, 20], [90, 20],
                            [90, -100]]]}  # number of hexagons in the region with grid resolution 5: 669644

hexagons = list(h3.polyfill(geoJSON, seed_point_resolution))
print("number of hexagons in the region with grid resolution 5:", len(hexagons))

centroids = np.array([h3.h3_to_geo(hex) for hex in hexagons])
full_lons = centroids[:, 1]
full_lats = centroids[:, 0]

# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()
mask_lon = mask_ds['glamf'].values
mask_lat = mask_ds['gphif'].values
mask_land = mask_ds['tmask'].values[:, 0, :, :]

fig = plt.figure()
ax = plt.axes()
colormap = color.ListedColormap(['grey', 'gainsboro'])

ax.scatter(full_lons, full_lats, s=0.2)

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(mask_lon[0], mask_lat[0], mask_land[0, 1:, 1:], cmap=colormap, alpha=0.5)

plt.show()

# From Mikael Kandoorp's: https://github.com/OceanParcels/Global_Analysis_Mikael/blob/main/create_release_uniform_h3.py
# Interpolate the release points onto this True/False mask using nearest neighbor.
# This mask can now be used to filter out points on land (by using ~land_val_release)
land_val_release = griddata((mask_lon.ravel(), mask_lat.ravel()), mask_land.ravel(),
                            (full_lons, full_lats), method='nearest')
bool_mask = land_val_release.astype(bool)

fig = plt.figure()
ax = plt.axes()
colormap = color.ListedColormap(['gainsboro', 'grey'])

# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(mask_lon[0], mask_lat[0], mask_land[0, 1:, 1:], cmap=colormap)

release_lons = full_lons[bool_mask]
release_lats = full_lats[bool_mask]

ax.scatter(release_lons, release_lats, s=0.3)
plt.show()

# export the lat-lon for release points to npz
np.savez_compressed('/nethome/manra003/analysis/paper01/H3_Res5_release_points.npz',
                    Longitude=release_lons,
                    Latitude=release_lats)
