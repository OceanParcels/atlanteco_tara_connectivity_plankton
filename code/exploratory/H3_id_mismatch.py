import h3
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from shapely.geometry.polygon import Polygon
import numpy as np

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'

temp_ds = xr.open_dataset(home_folder + 'FullTara_Res5_TS_1Jul2016_dt600_z200.nc').load()
lons, lats = temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values
temp_ds.close()

indexes = [376694, 376701, 12369]

res = 3
hex_local = np.empty(3, dtype='object')

p_lons = lons[indexes]
p_lats = lats[indexes]

for i in range(len(indexes)):
    ind = indexes[i]
    hex_local[i] = h3.geo_to_h3(lats[ind], lons[ind], res)
    print('index:', ind, lats[ind], lons[ind], hex_local[i])

model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'

mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()

# get the corner points to plot on the map
x = mask_ds['glamf']
y = mask_ds['gphif']

# get the mask values of the corner points
c = mask_ds['tmask'][:]

fig = plt.figure()
ax = plt.axes()
colormap = clr.ListedColormap(['gainsboro', 'white'])

ax.pcolormesh(x[0, 9:1500, 749:], y[0, 9:1500, 749:], c[0, 0, 10:1500, 750:], cmap=colormap)

plt.scatter(p_lons, p_lats, c='black', s=5)


def plot_hex(hex, c):
    polygons = h3.h3_set_to_multi_polygon([hex], geo_json=False)
    p = Polygon(polygons[0][0])
    y, x = p.exterior.xy
    ax.plot(x, y, color=c, linewidth=1)


hex_lorenz = ['83d066fffffffff', '83d066fffffffff', '83ee13fffffffff']

for h in hex_lorenz:
    plot_hex(str(h), 'red')

for h in hex_local:
    plot_hex(str(h), 'green')

plt.show()
### RESOLUTION 3
# from lorenz:
# Location 1: 376694 -43.248012286677714 19.941899283553155 83d066fffffffff
# Location 2: 376701 -42.98395329985144 19.949365337539927 83d066fffffffff
# Location 3: 12369 -73.84536843918157 -21.25173268157056 83ee13fffffffff

# from local:
# Location 1: 376694 -43.248012286677714 19.941899283553155 83d15bfffffffff
# Location 2: 376701 -42.98395329985144 19.949365337539927 83d15bfffffffff
# Location 3: 12369 -73.84536843918157 -21.25173268157056 83ee12fffffffff

### RESOLUTION 4
# From lorenz
# index: 376694 -43.248012286677714 19.941899283553155 84d15b3ffffffff
# index: 376701 -42.98395329985144 19.949365337539927 84d066dffffffff
# index: 12369 -73.84536843918157 -21.25173268157056 84ee13dffffffff
#
# from local
# index: 376694 -43.248012286677714 19.941899283553155 84d15b3ffffffff
# index: 376701 -42.98395329985144 19.949365337539927 84d066dffffffff
# index: 12369 -73.84536843918157 -21.25173268157056 84ee13dffffffff


### RESOLUTION 5
# form lorenz:
# index: 376694 -43.248012286677714 19.941899283553155 85d15b27fffffff
# index: 376701 -42.98395329985144 19.949365337539927 85d066dbfffffff
# index: 12369 -73.84536843918157 -21.25173268157056 85ee13dbfffffff

# from Local
# index: 376694 -43.248012286677714 19.941899283553155 85d15b27fffffff
# index: 376701 -42.98395329985144 19.949365337539927 85d066dbfffffff
# index: 12369 -73.84536843918157 -21.25173268157056 85ee13dbfffffff
