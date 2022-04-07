"""
To compute the connectivity matrix and min/max temperature/salinity per connection
from all the 120 simulations- 10 simulations for each month
using MasterHexList_Res3 for mapping
"""
import xarray as xr
from glob import glob
import numpy as np
import pandas as pd
from time import time
import sourcecode.core.matrixhelper as mxh
from scipy.sparse import csr_matrix
import os
import sys

home_folder = '/nethome/manra003/sim_out/'
data_folder = '/nethome/manra003/data/'
export_folder = '/nethome/manra003/analysis/paper01/depths/'

SIM_PER_MONTH = 10
parent_res = 3
child_res = 5


def compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth, output_path):
    t_mon1 = time()
    files = sorted(glob(home_folder + 'tara{0}m/FullTara_Res5_TS_1{1}*_dt600_z{0}.nc'.format(sim_depth, mon)))
    assert len(files) == SIM_PER_MONTH

    trans_array = np.zeros((no_grids, no_grids + 2))
    min_temp_array = np.full((no_grids, no_grids + 2), -999, dtype='float')
    max_temp_array = np.full((no_grids, no_grids + 2), 999, dtype='float')
    min_sal_array = np.full((no_grids, no_grids + 2), -999, dtype='float')
    max_sal_array = np.full((no_grids, no_grids + 2), 999, dtype='float')
    delete_count = 0

    for file in files:
        ds = xr.open_dataset(file)
        assert np.all(np.round(ds['z'][:, -1].values) == sim_depth)
        assert mxh.check_default_values(ds)

        invalid_indices = mxh.get_invalid_trajectories(ds)
        delete_count += len(invalid_indices)
#         print("invalid indices count:", len(invalid_indices))

        trans_array, min_temp_array, max_temp_array, min_sal_array, max_sal_array = mxh.get_monthly_matrix(
            ds, invalid_indices, trans_array, min_temp_array, max_temp_array,
            min_sal_array, max_sal_array, no_grids, parent_res, child_res, hex_indices, map_h3_to_mat)

    print("Total invalid trajectories removed: ", delete_count)

    min_temp_array[min_temp_array == -999] = 0
    max_temp_array[max_temp_array == 999] = 0
    min_sal_array[min_sal_array == -999] = 0
    max_sal_array[max_sal_array == 999] = 0

    assert (max_temp_array - min_temp_array).all() >= 0
    assert (max_sal_array - min_sal_array).all() >= 0

    mon_trans_matrix = csr_matrix(trans_array)
    print('Min/Max SUM of Transitions: {0} / {1}'.format(np.min(mon_trans_matrix.data), np.max(mon_trans_matrix.data)))

    mon_min_temp_matrix = csr_matrix(min_temp_array)
    print('Min/Max SUM of Minimum Temperature: {0} / {1}'.format(np.min(mon_min_temp_matrix.data),
                                                                 np.max(mon_min_temp_matrix.data)))

    mon_max_temp_matrix = csr_matrix(max_temp_array)
    print('Min/Max SUM of Maximum Temperature: {0} / {1}'.format(np.min(mon_max_temp_matrix.data),
                                                                 np.max(mon_max_temp_matrix.data)))

    mon_min_sal_matrix = csr_matrix(min_sal_array)
    print('Min/Max SUM of Minimum Salinity: {0} / {1}'.format(np.min(mon_min_sal_matrix.data),
                                                              np.max(mon_min_sal_matrix.data)))

    mon_max_sal_matrix = csr_matrix(max_sal_array)
    print('Min/Max SUM of Maximum Salinity: {0} / {1}'.format(np.min(mon_max_sal_matrix.data),
                                                              np.max(mon_max_sal_matrix.data)))

    # verify before exporting data
    # order of saving data is same for all fields
    assert np.array_equal(mon_trans_matrix.indices, mon_min_temp_matrix.indices)
    assert np.array_equal(mon_max_temp_matrix.indices, mon_max_sal_matrix.indices)
    assert np.array_equal(mon_trans_matrix.indptr, mon_max_temp_matrix.indptr)
    assert np.array_equal(mon_min_temp_matrix.indptr, mon_min_sal_matrix.indptr)

    print(
        "-------------------------------\nMonth: %s- \nTotalNumber of connections: %d" % (mon, mon_trans_matrix.sum()))
    new_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-2])[0]
    print('new particles: ', np.sum(mon_trans_matrix.data[new_index]))
    del_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-1])[0]
    print('deleted particles: ', np.sum(mon_trans_matrix.data[del_index]))

    # export all matrices to npz file
    np.savez_compressed(output_path + 'CSR_{0}.npz'.format(mon), transprob=mon_trans_matrix.data,
                        mintemp=mon_min_temp_matrix.data, maxtemp=mon_max_temp_matrix.data, minsal=mon_min_sal_matrix.data,
                        maxsal=mon_max_sal_matrix.data, indices=mon_trans_matrix.indices, indptr=mon_trans_matrix.indptr)

    t_mon2 = time()
    print("analysis time: ", t_mon2 - t_mon1)
    print("-------------------------------")


def main():
    args = sys.argv
    assert len(args) == 2
    sim_depth = np.int32(args[1])
    assert 0 <= sim_depth <= 500

    master_uni_hex = np.load(data_folder + 'MasterHexList_Res3.npy').tolist()
    assert len(master_uni_hex) == 8243

    no_particles = len(np.load(data_folder + 'AllRes5Children.npy'))
    assert no_particles == 377583

    no_grids = len(master_uni_hex)

    hex_indices = mxh.main_hex_list(master_uni_hex)
    mat_indices = np.arange(0, len(hex_indices))
    map_h3_to_mat = pd.Series(index=hex_indices, data=mat_indices)

    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    path_dir = export_folder + 't{0}m/'.format(sim_depth)
    os.makedirs(path_dir, exist_ok=True)
    [compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth, path_dir) for mon in
     months]


if __name__ == '__main__':
    main()
