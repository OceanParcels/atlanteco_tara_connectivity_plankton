import xarray as xr
from glob import glob
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.utils import resample
import matrixhelper as mxh
import os
from time import time

data_folder = '/nethome/manra003/data/'
home_folder = '/nethome/manra003/sim_out/'
export_folder = '/nethome/manra003/atlanteco_tara_connectivity_plankton/data/matrices/'
NEW = 'new'
DEL = 'deleted'
SIM_PER_MONTH = 10
child_res = 5
parent_res = 3
sim_depth = 0


def get_all_matrices_for_month(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth, sample_set, path_dir):
    tf1 = time()
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

#         invalid_indices = mxh.get_invalid_trajectories(ds)
#         delete_count += len(invalid_indices)
#         print("invalid indices count:", len(invalid_indices))

        t0 = mxh.get_hexid_from_parent(ds['lon'][:, 0].values, ds['lat'][:, 0].values, child_res,
                                       parent_res)
        t1 = mxh.get_hexid_from_parent(ds['lon'][:, -1].values, ds['lat'][:, -1].values, child_res,
                                       parent_res)
        # Not checking this for now in bootstrapping
        # def get_valid_data(field_name, loc):
        #     return np.delete(ds[field_name][:, loc].values, invalid_indices)
        #
        # hex_t0 = mxh.get_hexid_from_parent(get_valid_data('lon', 0), get_valid_data('lat', 0), child_res,
        #                                         parent_res)
        # # assert np.array_equal(hex_t0, master_all_hex_t0)
        # hex_t1 = mxh.get_hexid_from_parent(get_valid_data('lon', -1), get_valid_data('lat', -1), child_res,
        #                                         parent_res)

        hex_t0 = t0[sample_set]
        hex_t1 = t1[sample_set]

        # mask hex ids that are new
        hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, NEW)
        # mask hex ids in hex_t1_new that were deleted during the simulation
        max_time = np.max(ds['time'][:, -1].values)
        # hex_t1_new = np.where(ds['time'][:, -1].values < np.max(ds['time'][:, -1].values), DEL, hex_t1_new)
        # hex_t1_new = np.where(ds['time'][sample_set, -1].values < max_time, DEL, hex_t1_new)
        hex_t1_new = np.where(np.take(ds['time'][:, -1].values, indices=sample_set) < max_time, DEL, hex_t1_new)

        rows = map_h3_to_mat[hex_t0].values
        cols = map_h3_to_mat[hex_t1_new].values
        transitions = np.ones((len(hex_t0)))

        t_matrix = mxh.get_coo_matrix(transitions, rows, cols, no_grids)
        trans_array = np.append(trans_array, t_matrix.data)
        rows_array = np.append(rows_array, t_matrix.row)
        cols_array = np.append(cols_array, t_matrix.col)

        # get min and max temperature data
        # min_temperature, max_temperature = ds['min_temp'][:, -1].values, ds['max_temp'][:, -1].values
        min_temperature, max_temperature = np.take(ds['min_temp'][:, -1].values, indices=sample_set), \
                                           np.take(ds['max_temp'][:, -1].values, indices=sample_set)

        min_temp_matrix = mxh.get_coo_matrix(min_temperature, rows, cols, no_grids)
        max_temp_matrix = mxh.get_coo_matrix(max_temperature, rows, cols, no_grids)
        min_temp_array = np.append(min_temp_array, min_temp_matrix.data)
        max_temp_array = np.append(max_temp_array, max_temp_matrix.data)

        # get min and max salinity data
        # min_salinity, max_salinity = ds['min_sal'][:, -1].values, ds['max_sal'][:, -1].values
        min_salinity, max_salinity = np.take(ds['min_sal'][:, -1].values, indices=sample_set), \
                                     np.take(ds['max_sal'][:, -1].values, indices=sample_set)

        min_sal_matrix = mxh.get_coo_matrix(min_salinity, rows, cols, no_grids)
        max_sal_matrix = mxh.get_coo_matrix(max_salinity, rows, cols, no_grids)
        min_sal_array = np.append(min_sal_array, min_sal_matrix.data)
        max_sal_array = np.append(max_sal_array, max_sal_matrix.data)

    # collate entries for same row and column pair
    mon_trans_matrix = mxh.get_coo_matrix(trans_array, rows_array, cols_array, no_grids).tocsr()
    # print('Min/Max SUM of Transitions: {0} / {1}'.format(np.min(mon_trans_matrix.data), np.max(mon_trans_matrix.data)))

    mon_min_temp_matrix = mxh.get_coo_matrix(min_temp_array, rows_array, cols_array, no_grids).tocsr()
    # print('Min/Max SUM of Minimum Temperature: {0} / {1}'.format(np.min(mon_min_temp_matrix.data),
    #                                                              np.max(mon_min_temp_matrix.data)))

    mon_max_temp_matrix = mxh.get_coo_matrix(max_temp_array, rows_array, cols_array, no_grids).tocsr()
    # print('Min/Max SUM of Maximum Temperature: {0} / {1}'.format(np.min(mon_max_temp_matrix.data),
    #                                                              np.max(mon_max_temp_matrix.data)))

    mon_min_sal_matrix = mxh.get_coo_matrix(min_sal_array, rows_array, cols_array, no_grids).tocsr()
    # print('Min/Max SUM of Minimum Salinity: {0} / {1}'.format(np.min(mon_min_sal_matrix.data),
    #                                                           np.max(mon_min_sal_matrix.data)))

    mon_max_sal_matrix = mxh.get_coo_matrix(max_sal_array, rows_array, cols_array, no_grids).tocsr()
    # print('Min/Max SUM of Maximum Salinity: {0} / {1}'.format(np.min(mon_max_sal_matrix.data),
    #                                                           np.max(mon_max_sal_matrix.data)))

    # verify before exporting data
    # order of saving data is same for all fields
    assert np.array_equal(mon_trans_matrix.indices, mon_min_temp_matrix.indices)
    assert np.array_equal(mon_max_temp_matrix.indices, mon_max_sal_matrix.indices)
    assert np.array_equal(mon_trans_matrix.indptr, mon_max_temp_matrix.indptr)
    assert np.array_equal(mon_min_temp_matrix.indptr, mon_min_sal_matrix.indptr)

    # print("-------------------------------\nMonth: %s- \nTotalNumber of particles: %d" % (mon, mon_trans_matrix.sum()))
    assert mon_trans_matrix.sum() == len(sample_set) * SIM_PER_MONTH

    # create binary matrix from transitional data and confirm order
    bin_matrix = csr_matrix((np.ones(len(mon_trans_matrix.data)), mon_trans_matrix.indices, mon_trans_matrix.indptr),
                            shape=mon_trans_matrix.shape)
    assert np.array_equal(bin_matrix.indptr, mon_min_temp_matrix.indptr)
    assert np.array_equal(bin_matrix.indices, mon_max_sal_matrix.indices)

    def get_avg_field_per_grid(data, f_type, field):
        avg_field = data / mon_trans_matrix.data
        # print('Min/Max average {0} {1}: {2} / {3}'.format(f_type, field, np.min(avg_field), np.max(avg_field)))
        return avg_field
    
    avg_min_temp_per_grid = get_avg_field_per_grid(mon_min_temp_matrix.data, 'minimum', 'temperature')
    avg_max_temp_per_grid = get_avg_field_per_grid(mon_max_temp_matrix.data, 'maximum', 'temperature')
    avg_min_sal_per_grid = get_avg_field_per_grid(mon_min_sal_matrix.data, 'minimum', 'salinity')
    avg_max_sal_per_grid = get_avg_field_per_grid(mon_max_sal_matrix.data, 'maximum', 'salinity')
    # export all matrices to npz file
    np.savez_compressed(path_dir + 'Bin_CSR_{0}_z{1}.npz'.format(mon, sim_depth),
                        transprob=bin_matrix.data,
                        mintemp=avg_min_temp_per_grid, 
                        maxtemp=avg_max_temp_per_grid, 
                        minsal=avg_min_sal_per_grid,
                        maxsal=avg_max_sal_per_grid, 
                        indices=bin_matrix.indices, 
                        indptr=bin_matrix.indptr)
    
    tf2 = time()
    print(mon, ":", tf2 - tf1)


def main():
    master_uni_hex = np.load(data_folder + 'MasterHexList_Res3.npy').tolist()
    assert len(master_uni_hex) == 8243
    no_grids = len(master_uni_hex)

    no_particles = len(np.load(data_folder+'AllRes5Children.npy'))
    assert no_particles == 377583

    hex_indices = np.append(master_uni_hex, (NEW, DEL))
    mat_indices = np.arange(0, len(hex_indices))
    map_h3_to_mat = pd.Series(index=hex_indices, data=mat_indices)

    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    sample_size = [5000, 10000, 50000, 100000, 200000, 300000]
    
    for size in sample_size:
        t1 = time()
        # test with less number of particles without replacement
        sample_set = resample(np.arange(0, no_particles), replace=False, n_samples=size)
        sample_set = np.sort(sample_set)
        path_dir = export_folder + 'Boot_Sample/Monthly/size_{0}/'.format(size)

        os.makedirs(path_dir, exist_ok=True)
        [get_all_matrices_for_month(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth,
                                    sample_set, path_dir) for mon in months]
        t2 = time()
        print('Sample size {0} completed in time: {1}'.format(size, t2 - t1))


if __name__ == '__main__':
    main()
