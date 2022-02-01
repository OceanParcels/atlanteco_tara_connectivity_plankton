"""
Program to generate a master hex list- unique resolution 3 hex ids of the study region- 
In total= 8244 Hex ids
"""

import xarray as xr
import h3
from glob import glob
import numpy as np

home_folder = '/nethome/manra003/sim_out/'
sim_depth = 0
hex_res = 3 


def get_hex_id(lons, lats):
    return np.array([h3.geo_to_h3(y, x, hex_res) for x, y in zip(lons, lats)])


# prepare a master hex list from a random file from the final dataset
temp_ds = xr.open_dataset(np.random.choice(glob(home_folder + 'tara{0}m/FullTara_Res5_TS_1*'.format(sim_depth)))).load()
master_all_hex_t0 = get_hex_id(temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values)
temp_ds.close()
no_particles = len(master_all_hex_t0)
master_uni_hex = np.unique(master_all_hex_t0)
np.save('/nethome/manra003/data/MasterHexList_Res3.npy', master_uni_hex)