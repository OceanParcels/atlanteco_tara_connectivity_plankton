"""
To compute the monthly averaged transition matrix from all the 120 simulations- 
10 simulations for each month
using MasterHexList_Res3 for mapping
"""
import xarray as xr
from glob import glob
import numpy as np
import pandas as pd
from time import time
from sklearn.preprocessing import normalize
import sourcecode.core.matrixhelper as mxh
import os
import sys

home_folder = '/nethome/manra003/sim_out/'
data_folder = '/nethome/manra003/data/'
export_folder = '/nethome/manra003/data/'

SIM_PER_MONTH = 10
parent_res = 3
child_res = 5

# set option to 1 for normalized/binary TM, 2 for Sum of transitions.
option = 1


def compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_particles, no_grids, sim_depth, output_path):
    t_mon1 = time()
    # files = sorted(glob(home_folder + 'tara_data/FullAtlantic_2D_01{0}*_1month.nc'.format(mon)))
    files = sorted(glob(home_folder + 'tara{0}m/FullTara_Res5_TS_1{1}*_dt600_z{0}.nc'.format(sim_depth, mon)))
    assert len(files) == SIM_PER_MONTH
    trans_array = np.empty(0)
    min_temp_array = np.empty(0)
    max_temp_array = np.empty(0)
    min_sal_array = np.empty(0)
    max_sal_array = np.empty(0)
    rows_array = np.empty(0)
    cols_array = np.empty(0)
    delete_count = 0

    for file in files:
        ds = xr.open_dataset(file)
        assert np.all(np.round(ds['z'][:, -1].values) == sim_depth)
        assert mxh.check_default_values(ds)

        invalid_indices = mxh.get_invalid_trajectories(ds)
        delete_count += len(invalid_indices)
        print("invalid indices count:", len(invalid_indices))

        trans_array, rows_array, cols_array, min_temp_array, max_temp_array, min_sal_array, max_sal_array = mxh.get_monthly_tm(
            ds, invalid_indices, trans_array, rows_array, cols_array, min_temp_array, max_temp_array,
            min_sal_array, max_sal_array, no_grids, parent_res, child_res, hex_indices, map_h3_to_mat)

    print("Total invalid trajectories removed: ", delete_count)

    # collate entries for same row and column pair
    mon_trans_matrix = mxh.get_coo_matrix(trans_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Transitions: {0} / {1}'.format(np.min(mon_trans_matrix.data), np.max(mon_trans_matrix.data)))

    mon_min_temp_matrix = mxh.get_coo_matrix(min_temp_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Minimum Temperature: {0} / {1}'.format(np.min(mon_min_temp_matrix.data),
                                                                 np.max(mon_min_temp_matrix.data)))

    mon_max_temp_matrix = mxh.get_coo_matrix(max_temp_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Maximum Temperature: {0} / {1}'.format(np.min(mon_max_temp_matrix.data),
                                                                 np.max(mon_max_temp_matrix.data)))

    mon_min_sal_matrix = mxh.get_coo_matrix(min_sal_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Minimum Salinity: {0} / {1}'.format(np.min(mon_min_sal_matrix.data),
                                                              np.max(mon_min_sal_matrix.data)))

    mon_max_sal_matrix = mxh.get_coo_matrix(max_sal_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Maximum Salinity: {0} / {1}'.format(np.min(mon_max_sal_matrix.data),
                                                              np.max(mon_max_sal_matrix.data)))

    # verify before exporting data
    # order of saving data is same for all fields
    assert np.array_equal(mon_trans_matrix.indices, mon_min_temp_matrix.indices)
    assert np.array_equal(mon_max_temp_matrix.indices, mon_max_sal_matrix.indices)
    assert np.array_equal(mon_trans_matrix.indptr, mon_max_temp_matrix.indptr)
    assert np.array_equal(mon_min_temp_matrix.indptr, mon_min_sal_matrix.indptr)

    print("-------------------------------\nMonth: %s- \nTotalNumber of particles: %d" % (mon, mon_trans_matrix.sum()))
    assert mon_trans_matrix.sum() == no_particles * SIM_PER_MONTH - delete_count
    print("recorded transition pairs: ", len(mon_trans_matrix.data))
    new_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-2])[0]
    print('new particles: ', np.sum(mon_trans_matrix.data[new_index]))
    del_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-1])[0]
    print('deleted particles: ', np.sum(mon_trans_matrix.data[del_index]))

    # # perform row normalization for transitions and confirm order
    # norm_matrix = normalize(mon_trans_matrix, 'l1', axis=1, copy=True)
    # assert np.array_equal(norm_matrix.indptr, mon_min_temp_matrix.indptr)
    # assert np.array_equal(norm_matrix.indices, mon_max_sal_matrix.indices)

    # create binary matrix from transitional data and confirm order
    bin_matrix = mxh.binary_matrix(mon_trans_matrix)
    assert np.array_equal(bin_matrix.indptr, mon_min_temp_matrix.indptr)
    assert np.array_equal(bin_matrix.indices, mon_max_sal_matrix.indices)

    # Store Average or total number of transitions
    # compute the average min and max T/S for each grid cell
    def get_avg_field_per_grid(data, f_type, field):
        avg_field = data / mon_trans_matrix.data
        print('Min/Max average {0} {1}: {2} / {3}'.format(f_type, field, np.min(avg_field), np.max(avg_field)))
        return avg_field

    # Set option
    if option == 1:
        avg_min_temp_per_grid = get_avg_field_per_grid(mon_min_temp_matrix.data, 'minimum', 'temperature')
        avg_max_temp_per_grid = get_avg_field_per_grid(mon_max_temp_matrix.data, 'maximum', 'temperature')
        avg_min_sal_per_grid = get_avg_field_per_grid(mon_min_sal_matrix.data, 'minimum', 'salinity')
        avg_max_sal_per_grid = get_avg_field_per_grid(mon_max_sal_matrix.data, 'maximum', 'salinity')
        # export all matrices to npz file
        np.savez_compressed(output_path + 'CSR_{0}.npz'.format(mon), transprob=bin_matrix.data,
                            mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
                            maxsal=avg_max_sal_per_grid, indices=bin_matrix.indices, indptr=bin_matrix.indptr)
    elif option == 2:
        np.savez_compressed(output_path + 'Sum_CSR_{0}.npz'.format(mon),
                            transprob=mon_trans_matrix.data,
                            mintemp=mon_min_temp_matrix.data,
                            maxtemp=mon_max_temp_matrix.data,
                            minsal=mon_min_sal_matrix.data,
                            maxsal=mon_max_sal_matrix.data,
                            indices=mon_trans_matrix.indices,
                            indptr=mon_trans_matrix.indptr)
    else:
        raise ValueError('option value is incorrect')

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
    [compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_particles, no_grids, sim_depth, path_dir) for mon in
     months]


if __name__ == '__main__':
    main()
