"""
Program to generate a master hex list- unique resolution 3 hex ids of the study region- 
Updated method- now using h3_to_parent method
In total= 8243 Hex ids
"""

import xarray as xr
import h3
from glob import glob
import numpy as np

home_folder = '/nethome/manra003/sim_out/'
sim_depth = 100
parent_hex_res = 3
child_hex_res = 5

def get_parent_id(hex_list, res):
    return np.array([h3.h3_to_parent(str(h), res) for h in hex_list])


def get_hex_id(lons, lats,res):
    return np.array([h3.geo_to_h3(y, x, res) for x, y in zip(lons, lats)])


# prepare a master hex list from a random file from the final dataset
temp_ds = xr.open_dataset(np.random.choice(glob(home_folder + 'tara{0}m/FullTara_Res5_TS_1*'.format(sim_depth)))).load()

res5_children_hex = get_hex_id(temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values, child_hex_res)
temp_ds.close()
parent_all_hex = get_parent_id(res5_children_hex, parent_hex_res)
parent_uni_hex = np.unique(parent_all_hex)

np.save('/nethome/manra003/data/MasterHexList_Res3.npy', parent_uni_hex)

np.save('/nethome/manra003/data/AllRes5Children.npy', res5_children_hex)