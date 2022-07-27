import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import matplotlib.colors as clr
import h3

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
hex_res = 4
stations = pd.read_excel(home_folder + 'AllStations_Tara.xls', header=1)
lon = stations['Longitude']
lat = stations['Latitude']
code = stations['Station']

data = np.load(home_folder + 'Full_connectivity_output/Stations_min-T_connectivity_nan_TR2deg.npz', allow_pickle=True)
atlantic_codes = data['codes']

oa_series = pd.Series(atlantic_codes).str.contains('OA')
oa = atlantic_codes[oa_series]
sur_series = pd.Series(atlantic_codes).str.contains('SUR')
sur = atlantic_codes[sur_series]

common, oa_ind, ar2_ind = np.intersect1d(code, oa, return_indices=True)
common, sur_ind, ar2_ind = np.intersect1d(code, sur, return_indices=True)

all_indices = np.append(sur_ind, oa_ind)
lats = lat[all_indices].values
lons = lon[all_indices].values
# export stations list to file
ds = pd.DataFrame({'Code': code[all_indices].values,
                   'Latitude': lats,
                   'Longitude': lons,
                   'H3Id_res4': [h3.geo_to_h3(lat, lon, hex_res) for lat, lon in zip(lats, lons)]})
ds.sort_values(by='Latitude').to_csv(home_folder + 'Tara_Stations_hexId_res4.csv', index=False)
# Get mask land mask
# we need to coordinates file to access the corner points - glamf/gphif
model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

fig = plt.figure(figsize=(8, 6))
ax = plt.axes()
colormap = clr.ListedColormap(['whitesmoke', 'lightskyblue'])
plt.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel("Longitude(°)")
ax.set_ylabel("Latitude(°)")
# remove the first row and first column from the glamf/gphif to access points enclosed in the center
ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
# plt.show()
ax.scatter(lon[oa_ind], lat[oa_ind], c='b', s=5)
ax.scatter(lon[sur_ind], lat[sur_ind], c='r', s=5)

# for i in range(len(lon)):
#     xy = (lon[i], lat[i])
#     ax.annotate('(%s)' % code[i], xy=xy, textcoords='data')
#     print('(%s)' % code[i])
print('saving file')
plt.savefig(home_folder + "TaraStations.jpeg", bbox_inches='tight',
            pad_inches=0.5, dpi=300)
# coords = {"Station": code, "Latitudes": lat, "Longitudes": lon}
# df = pd.DataFrame(coords, columns=["Station", "Latitudes", "Longitudes"])
# df.to_csv(home_folder + 'Tara_stations.csv', index=False)
