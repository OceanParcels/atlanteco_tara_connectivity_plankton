import xarray as xr
import h3
from glob import glob
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from time import time
from sklearn.preprocessing import normalize

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/'
export_folder = home_folder + 'Task5_trahms/ProcessedTM/'
NEW = 'new'
DEL = 'deleted'
SIM_PER_MONTH = 10
hex_res = 3


def get_hex_id(lons, lats):
    return np.array([h3.geo_to_h3(y, x, hex_res) for x, y in zip(lons, lats)])


def get_coo_matrix(array, rows, cols, no_grids):
    matrix = coo_matrix((array, (rows, cols)), shape=(no_grids, no_grids + 2))
    matrix.sum_duplicates()
    return matrix


def get_all_matrices_for_month(mon, master_all_hex_t0, hex_indices, map_h3_to_mat, no_particles, no_grids):
    t_mon1 = time()
    # files = sorted(glob(home_folder + 'tara_data/FullAtlantic_2D_01{0}*_1month.nc'.format(mon)))
    files = sorted(glob(home_folder + 'Task4/tara_res5_01/FullTara_Res5_TS_{0}*_dt600.nc'.format(mon)))
    trans_array = np.empty(0)
    min_temp_array = np.empty(0)
    max_temp_array = np.empty(0)
    min_sal_array = np.empty(0)
    max_sal_array = np.empty(0)
    rows_array = np.empty(0)
    cols_array = np.empty(0)

    for file in files:
        ds = xr.open_dataset(file)
        hex_t0 = get_hex_id(ds['lon'][:, 0].values, ds['lat'][:, 0].values)
        assert np.array_equal(hex_t0, master_all_hex_t0)
        hex_t1 = get_hex_id(ds['lon'][:, -1].values, ds['lat'][:, -1].values)

        # mask hex ids that are new
        hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, NEW)
        # mask hex ids in hex_t1_new that were deleted during the simulation
        hex_t1_new = np.where(ds['time'][:, -1].values < np.max(ds['time'][:, -1].values), DEL, hex_t1_new)

        rows = map_h3_to_mat[hex_t0].values
        cols = map_h3_to_mat[hex_t1_new].values
        transitions = np.ones((len(hex_t0)))

        t_matrix = get_coo_matrix(transitions, rows, cols, no_grids)
        trans_array = np.append(trans_array, t_matrix.data)
        rows_array = np.append(rows_array, t_matrix.row)
        cols_array = np.append(cols_array, t_matrix.col)

        # get min and max temperature data
        min_temperature, max_temperature = ds['min_temp'][:, -1].values, ds['max_temp'][:, -1].values
        min_temp_matrix = get_coo_matrix(min_temperature, rows, cols, no_grids)
        max_temp_matrix = get_coo_matrix(max_temperature, rows, cols, no_grids)
        min_temp_array = np.append(min_temp_array, min_temp_matrix.data)
        max_temp_array = np.append(max_temp_array, max_temp_matrix.data)

        # get min and max salinity data
        min_salinity, max_salinity = ds['min_sal'][:, -1].values, ds['max_sal'][:, -1].values
        min_sal_matrix = get_coo_matrix(min_salinity, rows, cols, no_grids)
        max_sal_matrix = get_coo_matrix(max_salinity, rows, cols, no_grids)
        min_sal_array = np.append(min_sal_array, min_sal_matrix.data)
        max_sal_array = np.append(max_sal_array, max_sal_matrix.data)

    # collate entries for same row and column pair
    mon_trans_matrix = get_coo_matrix(trans_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Transitions: {0} / {1}'.format(np.min(mon_trans_matrix.data), np.max(mon_trans_matrix.data)))

    mon_min_temp_matrix = get_coo_matrix(min_temp_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Minimum Temperature: {0} / {1}'.format(np.min(mon_min_temp_matrix.data),
                                                                 np.max(mon_min_temp_matrix.data)))

    mon_max_temp_matrix = get_coo_matrix(max_temp_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Maximum Temperature: {0} / {1}'.format(np.min(mon_max_temp_matrix.data),
                                                                 np.max(mon_max_temp_matrix.data)))

    mon_min_sal_matrix = get_coo_matrix(min_sal_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Minimum Salinity: {0} / {1}'.format(np.min(mon_min_sal_matrix.data),
                                                              np.max(mon_min_sal_matrix.data)))

    mon_max_sal_matrix = get_coo_matrix(max_sal_array, rows_array, cols_array, no_grids).tocsr()
    print('Min/Max SUM of Maximum Salinity: {0} / {1}'.format(np.min(mon_max_sal_matrix.data),
                                                              np.max(mon_max_sal_matrix.data)))

    # verify before exporting data
    # order of saving data is same for all fields
    assert np.array_equal(mon_trans_matrix.indices, mon_min_temp_matrix.indices)
    assert np.array_equal(mon_max_temp_matrix.indices, mon_max_sal_matrix.indices)
    assert np.array_equal(mon_trans_matrix.indptr, mon_max_temp_matrix.indptr)
    assert np.array_equal(mon_min_temp_matrix.indptr, mon_min_sal_matrix.indptr)

    print("-------------------------------\nMonth: %s- \nTotalNumber of particles: %d" % (mon, mon_trans_matrix.sum()))
    assert mon_trans_matrix.sum() == no_particles * SIM_PER_MONTH
    print("recorded transition pairs: ", len(mon_trans_matrix.data))
    new_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-2])[0]
    print('new particles: ', np.sum(mon_trans_matrix.data[new_index]))
    del_index = np.where(mon_trans_matrix.indices == map_h3_to_mat[-1])[0]
    print('deleted particles: ', np.sum(mon_trans_matrix.data[del_index]))

    # perform row normalization for transitions and confirm order
    norm_matrix = normalize(mon_trans_matrix, 'l1', axis=1, copy=True)
    assert np.array_equal(norm_matrix.indptr, mon_min_temp_matrix.indptr)
    assert np.array_equal(norm_matrix.indices, mon_max_sal_matrix.indices)

    # Store Average or totoal number of transitions
    # compute the average min and max T/S for each grid cell
    def get_avg_field_per_grid(data, f_type, field):
        avg_field = data / mon_trans_matrix.data
        print('Min/Max average {0} {1}: {2} / {3}'.format(f_type, field, np.min(avg_field), np.max(avg_field)))
        return avg_field

    # Set option
    option = 2
    if option == 1:
        avg_min_temp_per_grid = get_avg_field_per_grid(mon_min_temp_matrix.data, 'minimum', 'temperature')
        avg_max_temp_per_grid = get_avg_field_per_grid(mon_max_temp_matrix.data, 'maximum', 'temperature')
        avg_min_sal_per_grid = get_avg_field_per_grid(mon_min_sal_matrix.data, 'minimum', 'salinity')
        avg_max_sal_per_grid = get_avg_field_per_grid(mon_max_sal_matrix.data, 'maximum', 'salinity')
        # export all matrices to npz file
        np.savez_compressed(export_folder + 'CSR_{0}.npz'.format(mon), transprob=norm_matrix.data,
                            mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
                            maxsal=avg_max_sal_per_grid, indices=norm_matrix.indices, indptr=norm_matrix.indptr)
    elif option == 2:
        np.savez_compressed(export_folder + 'Sum_CSR_{0}.npz'.format(mon),
                            transprob=mon_trans_matrix.data,
                            mintemp=mon_min_temp_matrix.data,
                            maxtemp=mon_max_temp_matrix.data,
                            minsal=mon_min_sal_matrix.data,
                            maxsal=mon_max_sal_matrix.data,
                            indices=mon_trans_matrix.indices,
                            indptr=mon_trans_matrix.indptr)
    else:
        raise ValueError('option value is incorrect')
    # np.savez_compressed(export_folder + 'CSR_{0}.npz'.format(mon), transprob=norm_matrix.data,
    #                     indices=norm_matrix.indices, indptr=norm_matrix.indptr)
    t_mon2 = time()
    print("analysis time: ", t_mon2 - t_mon1)


def main():
    # prepare a master hex list from a random file from the final dataset
    temp_ds = xr.open_dataset(np.random.choice(glob(home_folder + 'Task4/tara_res5_01/FullTara_Res5_TS_*'))).load()
    # temp_ds = xr.open_dataset(np.random.choice(glob(home_folder + 'tara_data/FullAtlantic_2D_01*'))).load()
    master_all_hex_t0 = get_hex_id(temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values)
    temp_ds.close()
    no_particles = len(master_all_hex_t0)
    master_uni_hex = np.unique(master_all_hex_t0)

    # master_uni_hex = np.load(home_folder + 'Task5_trahms/MasterHexList.npy').tolist()
    # no_particles = 377583

    no_grids = len(master_uni_hex)

    hex_indices = np.append(master_uni_hex, (NEW, DEL))
    mat_indices = np.arange(0, len(hex_indices))
    map_h3_to_mat = pd.Series(index=hex_indices, data=mat_indices)

    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    [get_all_matrices_for_month(mon, master_all_hex_t0, hex_indices, map_h3_to_mat, no_particles, no_grids) for mon in
     months]


if __name__ == '__main__':
    main()