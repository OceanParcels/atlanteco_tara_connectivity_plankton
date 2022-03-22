import xarray as xr
from glob import glob
import numpy as np
import pandas as pd
from sklearn.utils import resample
import sourcecode.core.matrixhelper as mxh
import os
from time import time

data_folder = '/nethome/manra003/data/'
home_folder = '/nethome/manra003/sim_out/'
export_folder = '/nethome/storage/shared/oceanparcels/output_data/data_Darshika/TaraC/'
NEW = 'new'
DEL = 'deleted'
SIM_PER_MONTH = 10
child_res = 5
parent_res = 3
sim_depth = 0


def get_all_matrices_for_month(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth, no_particles, sample_set,
                               path_dir):
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

    for file in files:
        ds = xr.open_dataset(file)
        assert np.all(np.round(ds['z'][:, -1].values) == sim_depth)
        assert mxh.check_default_values(ds)

        invalid_indices = mxh.get_invalid_trajectories(ds)
        # here: also add particles not selected for sample set
        invalid_indices = np.union1d(invalid_indices, np.delete(range(no_particles), sample_set))

        trans_array, rows_array, cols_array, min_temp_array, max_temp_array, min_sal_array, max_sal_array = mxh.get_monthly_tm(
            ds, invalid_indices, trans_array, rows_array, cols_array, min_temp_array, max_temp_array,
            min_sal_array, max_sal_array, no_grids, parent_res, child_res, hex_indices, map_h3_to_mat)

    # collate entries for same row and column pair
    mon_trans_matrix = mxh.get_coo_matrix(trans_array, rows_array, cols_array, no_grids).tocsr()
    mon_min_temp_matrix = mxh.get_coo_matrix(min_temp_array, rows_array, cols_array, no_grids).tocsr()
    mon_max_temp_matrix = mxh.get_coo_matrix(max_temp_array, rows_array, cols_array, no_grids).tocsr()
    mon_min_sal_matrix = mxh.get_coo_matrix(min_sal_array, rows_array, cols_array, no_grids).tocsr()
    mon_max_sal_matrix = mxh.get_coo_matrix(max_sal_array, rows_array, cols_array, no_grids).tocsr()

    # verify before exporting data
    assert np.array_equal(mon_trans_matrix.indices, mon_min_temp_matrix.indices)
    assert np.array_equal(mon_max_temp_matrix.indices, mon_max_sal_matrix.indices)
    assert np.array_equal(mon_trans_matrix.indptr, mon_max_temp_matrix.indptr)
    assert np.array_equal(mon_min_temp_matrix.indptr, mon_min_sal_matrix.indptr)
    assert mon_trans_matrix.sum() == len(sample_set) * SIM_PER_MONTH

    # create binary matrix from transitional data and confirm order
    bin_matrix = mxh.binary_matrix(mon_trans_matrix)
    assert np.array_equal(bin_matrix.indptr, mon_min_temp_matrix.indptr)
    assert np.array_equal(bin_matrix.indices, mon_max_sal_matrix.indices)

    avg_min_temp_per_grid = mxh.avg_field_per_grid(mon_min_temp_matrix.data, mon_trans_matrix.data)
    avg_max_temp_per_grid = mxh.avg_field_per_grid(mon_max_temp_matrix.data, mon_trans_matrix.data)
    avg_min_sal_per_grid = mxh.avg_field_per_grid(mon_min_sal_matrix.data, mon_trans_matrix.data)
    avg_max_sal_per_grid = mxh.avg_field_per_grid(mon_max_sal_matrix.data, mon_trans_matrix.data)
    # export all matrices to npz file
    np.savez_compressed(path_dir + 'CSR_{0}_z{1}.npz'.format(mon, sim_depth),
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

    no_particles = len(np.load(data_folder + 'AllRes5Children.npy'))
    assert no_particles == 377583

    hex_indices = np.append(master_uni_hex, (NEW, DEL))
    mat_indices = np.arange(0, len(hex_indices))
    map_h3_to_mat = pd.Series(index=hex_indices, data=mat_indices)

    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    sample_size = [5000, 10000, 50000, 100000, 200000, 300000]

    init_state = 1
    end_state = 10

    for size in sample_size:
        t1 = time()

        for state in np.arange(init_state, end_state + 1):
            # test with less number of particles without replacement
            sample_set = resample(np.arange(0, no_particles), replace=False, n_samples=size, random_state=state)
            sample_set = np.sort(sample_set)
            path_dir = export_folder + 'Boot_Sample/Size_{0}/Monthly/State_{1}/'.format(size, state)

            os.makedirs(path_dir, exist_ok=True)
            [get_all_matrices_for_month(mon, hex_indices, map_h3_to_mat, no_grids, sim_depth, no_particles,
                                        sample_set, path_dir) for mon in months]

        t2 = time()
        print('Sample size {0} completed in time: {1}'.format(size, t2 - t1))


if __name__ == '__main__':
    main()
