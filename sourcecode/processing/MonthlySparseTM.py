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
# from sklearn.preprocessing import normalize
import os
import sys

home_folder = '/nethome/manra003/sim_out/'
data_folder = '/nethome/manra003/analysis/paper01/'
export_folder = '/nethome/manra003/analysis/paper01/depths/'
# home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'
# data_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'
# export_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'

SIM_PER_MONTH = 10
NEW = 'new'
parent_res = 3
child_res = 5


def compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth, output_path):
    t_mon1 = time()
    files = sorted(glob(home_folder + 'tara{0}m/FullTara_Res5_TS_1{1}*_dt600_z{0}.zarr'.format(sim_depth, mon)))
    print(files)
    assert len(files) == SIM_PER_MONTH

    trans_array = np.empty(0, dtype=np.int32)
    rows_array = np.empty(0, dtype=np.int32)
    cols_array = np.empty(0, dtype=np.int32)

    min_Mintemp_array = np.full((no_grids, no_grids + 1), 999, dtype=np.float32)
    max_Mintemp_array = np.full((no_grids, no_grids + 1), -999, dtype=np.float32)
    min_Maxtemp_array = np.full((no_grids, no_grids + 1), 999, dtype=np.float32)
    max_Maxtemp_array = np.full((no_grids, no_grids + 1), -999, dtype=np.float32)
    Mintemp_array = np.empty(0, dtype=np.float32)
    Maxtemp_array = np.empty(0, dtype=np.float32)

    min_Minsal_array = np.full((no_grids, no_grids + 1), 999, dtype=np.float32)
    max_Minsal_array = np.full((no_grids, no_grids + 1), -999, dtype=np.float32)
    min_Maxsal_array = np.full((no_grids, no_grids + 1), 999, dtype=np.float32)
    max_Maxsal_array = np.full((no_grids, no_grids + 1), -999, dtype=np.float32)
    Minsal_array = np.empty(0, dtype=np.float32)
    Maxsal_array = np.empty(0, dtype=np.float32)

    delete_count = 0

    for file in files:
        ds = xr.open_dataset(file)
        nan_z = ds['z'][:, -1].values
        
        assert np.all(np.round(nan_z[~np.isnan(nan_z)]) == sim_depth)
        assert mxh.check_default_values(ds)

        invalid_indices = mxh.get_invalid_trajectories(ds)
        delete_count += len(invalid_indices)

        def get_valid_data(field_name, loc):
            return np.delete(ds[field_name][:, loc].values, invalid_indices)

        hex_t0 = mxh.get_hexid_from_parent(get_valid_data('lon', 0), get_valid_data('lat', 0), child_res,
                                           parent_res)
        # assert np.array_equal(hex_t0, master_all_hex_t0)
        hex_t1 = mxh.get_hexids(get_valid_data('lon', -1), get_valid_data('lat', -1), parent_res)

        # mask hex ids that are new
        hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, NEW)
        # mask hex ids in hex_t1_new that were deleted during the simulation
        # hex_t1_new = np.where(get_valid_data('time', -1) < np.max(ds['time'][:, -1].values), DEL, hex_t1_new)

        rows = map_h3_to_mat[hex_t0].values
        cols = map_h3_to_mat[hex_t1_new].values

        transitions = np.ones((len(hex_t0)))

        # reduces duplicate connections between grid pairs,
        # returns sum of number of connections between grid pairs
        t_matrix = mxh.get_coo_matrix(transitions, rows, cols, no_grids)
        trans_array = np.append(trans_array, t_matrix.data)
        rows_array = np.append(rows_array, t_matrix.row)
        cols_array = np.append(cols_array, t_matrix.col)

        # region: get min and max temperature and salinity data
        t_prop1 = time()
        min_temperature, max_temperature = get_valid_data('min_temp', -1), get_valid_data('max_temp', -1)
        min_salinity, max_salinity = get_valid_data('min_sal', -1), get_valid_data('max_sal', -1)

        property_df = pd.DataFrame({'rows': rows,
                                    'cols': cols,
                                    'min_t': min_temperature,
                                    'max_t': max_temperature,
                                    'min_s': min_salinity,
                                    'max_s': max_salinity})

        analysis_df = property_df.groupby(['rows', 'cols']).aggregate({'min_t': ['min', 'max'],
                                                                       'max_t': ['min', 'max'],
                                                                       'min_s': ['min', 'max'],
                                                                       'max_s': ['min', 'max']}).reset_index()

        r = analysis_df['rows'].values
        c = analysis_df['cols'].values
        # Min max property extraction
        # Example:minimum of minimum T and maximum of minimum T
        min_Mintemp_array[r, c] = np.minimum(min_Mintemp_array[r, c], analysis_df.min_t['min'].values)
        max_Mintemp_array[r, c] = np.maximum(max_Mintemp_array[r, c], analysis_df.min_t['max'].values)

        min_Maxtemp_array[r, c] = np.minimum(min_Maxtemp_array[r, c], analysis_df.max_t['min'].values)
        max_Maxtemp_array[r, c] = np.maximum(max_Maxtemp_array[r, c], analysis_df.max_t['max'].values)

        min_Minsal_array[r, c] = np.minimum(min_Minsal_array[r, c], analysis_df.min_s['min'].values)
        max_Minsal_array[r, c] = np.maximum(max_Minsal_array[r, c], analysis_df.min_s['max'].values)

        min_Maxsal_array[r, c] = np.minimum(min_Maxsal_array[r, c], analysis_df.max_s['min'].values)
        max_Maxsal_array[r, c] = np.maximum(max_Maxsal_array[r, c], analysis_df.max_s['max'].values)

        min_temp_matrix = mxh.get_coo_matrix(min_temperature, rows, cols, no_grids)
        max_temp_matrix = mxh.get_coo_matrix(max_temperature, rows, cols, no_grids)
        Mintemp_array = np.append(Mintemp_array, min_temp_matrix.data)
        Maxtemp_array = np.append(Maxtemp_array, max_temp_matrix.data)

        min_sal_matrix = mxh.get_coo_matrix(min_salinity, rows, cols, no_grids)
        max_sal_matrix = mxh.get_coo_matrix(max_salinity, rows, cols, no_grids)
        Minsal_array = np.append(Minsal_array, min_sal_matrix.data)
        Maxsal_array = np.append(Maxsal_array, max_sal_matrix.data)
        t_prop2 = time()
        print('time taken- ', t_prop2 - t_prop1)
        # endregion
    print("Total invalid trajectories removed: ", delete_count)

    min_Mintemp_array[min_Mintemp_array == 999] = 0
    max_Mintemp_array[max_Mintemp_array == -999] = 0
    min_Maxtemp_array[min_Maxtemp_array == 999] = 0
    max_Maxtemp_array[max_Maxtemp_array == -999] = 0
    min_Minsal_array[min_Minsal_array == 999] = 0
    max_Minsal_array[max_Minsal_array == -999] = 0
    min_Maxsal_array[min_Maxsal_array == 999] = 0
    max_Maxsal_array[max_Maxsal_array == -999] = 0

    mon_trans_matrix = mxh.get_coo_matrix(trans_array, rows_array, cols_array, no_grids).tocsr()
    print('Range of Transitions: {0} / {1}'.format(np.min(mon_trans_matrix.data), np.max(mon_trans_matrix.data)))

    # region: extract temperature to sparse matrices
    mon_min_mintemp_matrix = csr_matrix(min_Mintemp_array)
    print('Range of Minimum of Minimum Temperature: {0} / {1}'.format(np.min(mon_min_mintemp_matrix.data),
                                                                      np.max(mon_min_mintemp_matrix.data)))
    mon_max_mintemp_matrix = csr_matrix(max_Mintemp_array)
    print('Range of Maximum of Minimum Temperature: {0} / {1}'.format(np.min(mon_max_mintemp_matrix.data),
                                                                      np.max(mon_max_mintemp_matrix.data)))
    mon_min_maxtemp_matrix = csr_matrix(min_Maxtemp_array)
    print('Range of Minimum of Maximum Temperature: {0} / {1}'.format(np.min(mon_min_maxtemp_matrix.data),
                                                                      np.max(mon_min_maxtemp_matrix.data)))
    mon_max_maxtemp_matrix = csr_matrix(max_Maxtemp_array)
    print('Range of Minimum of Maximum Temperature: {0} / {1}'.format(np.min(mon_max_maxtemp_matrix.data),
                                                                      np.max(mon_max_maxtemp_matrix.data)))
    mon_min_temp_matrix = mxh.get_coo_matrix(Mintemp_array, rows_array, cols_array, no_grids).tocsr()
    print('Sum Minimum Temperature: {0} / {1}'.format(np.min(mon_min_temp_matrix.data),
                                                      np.max(mon_min_temp_matrix.data)))
    mon_max_temp_matrix = mxh.get_coo_matrix(Maxtemp_array, rows_array, cols_array, no_grids).tocsr()
    print('Sum Maximum Temperature: {0} / {1}'.format(np.min(mon_max_temp_matrix.data),
                                                      np.max(mon_max_temp_matrix.data)))
    # endregion

    # region: extract salinity to sparse matrices
    mon_min_minsal_matrix = csr_matrix(min_Minsal_array)
    print('Range of Minimum of Minimum Salinity: {0} / {1}'.format(np.min(mon_min_minsal_matrix.data),
                                                                   np.max(mon_min_minsal_matrix.data)))
    mon_max_minsal_matrix = csr_matrix(max_Minsal_array)
    print('Range of Maximum of Minimum Salinity: {0} / {1}'.format(np.min(mon_max_minsal_matrix.data),
                                                                   np.max(mon_max_minsal_matrix.data)))
    mon_min_maxsal_matrix = csr_matrix(min_Maxsal_array)
    print('Range of Minimum of Maximum Salinity: {0} / {1}'.format(np.min(mon_min_maxsal_matrix.data),
                                                                   np.max(mon_min_maxsal_matrix.data)))
    mon_max_maxsal_matrix = csr_matrix(max_Maxsal_array)
    print('Range of Minimum of Maximum Salinity: {0} / {1}'.format(np.min(mon_max_maxsal_matrix.data),
                                                                   np.max(mon_max_maxsal_matrix.data)))
    mon_min_sal_matrix = mxh.get_coo_matrix(Minsal_array, rows_array, cols_array, no_grids).tocsr()
    print('Sum Minimum Salinity: {0} / {1}'.format(np.min(mon_min_sal_matrix.data),
                                                   np.max(mon_min_sal_matrix.data)))
    mon_max_sal_matrix = mxh.get_coo_matrix(Maxsal_array, rows_array, cols_array, no_grids).tocsr()
    print('Sum Maximum Salinity: {0} / {1}'.format(np.min(mon_max_sal_matrix.data),
                                                   np.max(mon_max_sal_matrix.data)))
    # endregion
    # verify before exporting data
    # order of saving data is same for all fields
    assert np.array_equal(mon_trans_matrix.indices, mon_min_mintemp_matrix.indices)
    assert np.array_equal(mon_min_maxtemp_matrix.indices, mon_max_minsal_matrix.indices)
    assert np.array_equal(mon_trans_matrix.indptr, mon_max_maxtemp_matrix.indptr)
    assert np.array_equal(mon_min_maxtemp_matrix.indptr, mon_min_minsal_matrix.indptr)

    print(
        "-------------------------------\nMonth: %s- \nTotalNumber of connections: %d" % (mon, mon_trans_matrix.sum()))

    # perform row normalization for transitions and confirm order
    # norm_matrix = normalize(mon_trans_matrix, 'l1', axis=1, copy=True)
    # assert np.array_equal(norm_matrix.indptr, mon_min_temp_matrix.indptr)
    # assert np.array_equal(norm_matrix.indices, mon_max_sal_matrix.indices)

    # Store Average or total number of transitions
    # compute the average min and max T/S for each grid cell
    def get_avg_field_per_grid(data, f_type, field):
        avg_field = data / mon_trans_matrix.data
        print('Min/Max average {0} {1}: {2} / {3}'.format(f_type, field, np.min(avg_field), np.max(avg_field)))
        return avg_field

    new_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-1])[0]
    print('new particle locations: ', np.sum(mon_trans_matrix.data[new_index]))
    # del_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-1])[0]
    # print('deleted particles: ', np.sum(mon_trans_matrix.data[del_index]))

    avg_min_temp_per_grid = get_avg_field_per_grid(mon_min_temp_matrix.data, 'minimum', 'temperature')
    avg_max_temp_per_grid = get_avg_field_per_grid(mon_max_temp_matrix.data, 'maximum', 'temperature')
    avg_min_sal_per_grid = get_avg_field_per_grid(mon_min_sal_matrix.data, 'minimum', 'salinity')
    avg_max_sal_per_grid = get_avg_field_per_grid(mon_max_sal_matrix.data, 'maximum', 'salinity')

    # export all matrices to npz file
    np.savez_compressed(output_path + 'CSR_{0}z{1}.npz'.format(mon, sim_depth),
                        transprob=mon_trans_matrix.data,
                        min_mintemp=mon_min_mintemp_matrix.data,
                        max_mintemp=mon_max_mintemp_matrix.data,
                        min_maxtemp=mon_min_maxtemp_matrix.data,
                        max_maxtemp=mon_max_maxtemp_matrix.data,
                        min_minsal=mon_min_minsal_matrix.data,
                        max_minsal=mon_max_minsal_matrix.data,
                        min_maxsal=mon_min_maxsal_matrix.data,
                        max_maxsal=mon_max_maxsal_matrix.data,
                        avg_min_temp=avg_min_temp_per_grid,
                        avg_max_temp=avg_max_temp_per_grid,
                        avg_min_sal=avg_min_sal_per_grid,
                        avg_max_sal=avg_max_sal_per_grid,
                        indices=mon_trans_matrix.indices,
                        indptr=mon_trans_matrix.indptr)

    t_mon2 = time()
    print("analysis time: ", t_mon2 - t_mon1)
    print("-------------------------------")


def main():
    args = sys.argv
    assert len(args) == 2
    sim_depth = np.int32(args[1])
    assert 0 <= sim_depth <= 500

    if parent_res == 3:
        master_uni_hex = np.load(data_folder + 'H3_Res3_MasterHexList.npz')['Res3_HexId'].tolist()
        assert len(master_uni_hex) == 8191
    elif parent_res == 2:
        master_uni_hex = np.load(data_folder + 'H3_Res2_MasterHexList.npz')['Res2_HexId'].tolist()
        assert len(master_uni_hex) == 1260
    elif parent_res == 4:
        master_uni_hex = np.load(data_folder + 'H3_Res4_MasterHexList.npz')['Res4_HexId'].tolist()
        assert len(master_uni_hex) == 54947
    else:
        raise ValueError('check parent_res')
        
    no_particles = len(np.load(data_folder + 'H3_Res5_release_points.npz')['Latitude'])
    assert no_particles == 375570

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
