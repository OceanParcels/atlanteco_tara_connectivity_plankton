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
export_folder = '/nethome/manra003/atlanteco_tara_connectivity_plankton/data/matrices/'

NEW = 'new'
DEL = 'deleted'
SIM_PER_MONTH = 10
parent_res = 3
child_res = 5

# set option to 1 for normalized TM, 2 for Sum of transitions.
option = 1


def compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_particles, no_grids, sim_depth):
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

        def get_valid_data(field_name, loc):
            return np.delete(ds[field_name][:, loc].values, invalid_indices)

        hex_t0 = mxh.get_hexid_from_parent(get_valid_data('lon', 0), get_valid_data('lat', 0), child_res,
                                                    parent_res)
        # assert np.array_equal(hex_t0, master_all_hex_t0)
        hex_t1 = mxh.get_hexid_from_parent(get_valid_data('lon', -1), get_valid_data('lat', -1), child_res,
                                                    parent_res)

        # mask hex ids that are new
        hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, NEW)
        # mask hex ids in hex_t1_new that were deleted during the simulation
        hex_t1_new = np.where(get_valid_data('time', -1) < np.max(ds['time'][:, -1].values), DEL, hex_t1_new)

        rows = map_h3_to_mat[hex_t0].values
        cols = map_h3_to_mat[hex_t1_new].values
        transitions = np.ones((len(hex_t0)))

        t_matrix = mxh.get_coo_matrix(transitions, rows, cols, no_grids)
        trans_array = np.append(trans_array, t_matrix.data)
        rows_array = np.append(rows_array, t_matrix.row)
        cols_array = np.append(cols_array, t_matrix.col)

        # get min and max temperature data
        min_temperature, max_temperature = get_valid_data('min_temp', -1), get_valid_data('max_temp', -1)
        print('Min/Max of Minimum Temperature: {0} / {1}'.format(np.min(min_temperature), np.max(min_temperature)))
        print('Min/Max of Maximum Temperature: {0} / {1}'.format(np.min(max_temperature), np.max(max_temperature)))

        min_temp_matrix = mxh.get_coo_matrix(min_temperature, rows, cols, no_grids)
        max_temp_matrix = mxh.get_coo_matrix(max_temperature, rows, cols, no_grids)
        min_temp_array = np.append(min_temp_array, min_temp_matrix.data)
        max_temp_array = np.append(max_temp_array, max_temp_matrix.data)

        # get min and max salinity data
        min_salinity, max_salinity = get_valid_data('min_sal', -1), get_valid_data('max_sal', -1)
        print('Min/Max of Minimum Salinity: {0} / {1}'.format(np.min(min_salinity), np.max(min_salinity)))
        print('Min/Max of Maximum Salinity: {0} / {1}'.format(np.min(max_salinity), np.max(max_salinity)))

        min_sal_matrix = mxh.get_coo_matrix(min_salinity, rows, cols, no_grids)
        max_sal_matrix = mxh.get_coo_matrix(max_salinity, rows, cols, no_grids)
        min_sal_array = np.append(min_sal_array, min_sal_matrix.data)
        max_sal_array = np.append(max_sal_array, max_sal_matrix.data)

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
    if option == 1:
        avg_min_temp_per_grid = get_avg_field_per_grid(mon_min_temp_matrix.data, 'minimum', 'temperature')
        avg_max_temp_per_grid = get_avg_field_per_grid(mon_max_temp_matrix.data, 'maximum', 'temperature')
        avg_min_sal_per_grid = get_avg_field_per_grid(mon_min_sal_matrix.data, 'minimum', 'salinity')
        avg_max_sal_per_grid = get_avg_field_per_grid(mon_max_sal_matrix.data, 'maximum', 'salinity')
        # export all matrices to npz file
        np.savez_compressed(export_folder + 't{0}m/CSR_{1}.npz'.format(sim_depth, mon), transprob=norm_matrix.data,
                            mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
                            maxsal=avg_max_sal_per_grid, indices=norm_matrix.indices, indptr=norm_matrix.indptr)
    elif option == 2:
        np.savez_compressed(export_folder + 't{0}m/Sum_CSR_{1}.npz'.format(sim_depth, mon),
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

    master_uni_hex = np.load('/nethome/manra003/data/MasterHexList_Res3.npy').tolist()
    assert len(master_uni_hex) == 8243

    no_particles = len(np.load('/nethome/manra003/data/AllRes5Children.npy'))
    assert no_particles == 377583

    no_grids = len(master_uni_hex)

    hex_indices = np.append(master_uni_hex, (NEW, DEL))
    mat_indices = np.arange(0, len(hex_indices))
    map_h3_to_mat = pd.Series(index=hex_indices, data=mat_indices)

    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    #     months = np.array(['Jan'])
    output_path = export_folder + 't{0}m/'.format(sim_depth)
    os.makedirs(output_path, exist_ok=True)
    [compute_transition_matrix(mon, hex_indices, map_h3_to_mat, no_particles, no_grids, sim_depth) for mon in
     months]


if __name__ == '__main__':
    main()
