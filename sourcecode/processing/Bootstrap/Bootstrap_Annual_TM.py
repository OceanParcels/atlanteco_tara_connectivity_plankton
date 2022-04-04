import numpy as np
from scipy.sparse import csr_matrix, save_npz
import os
from time import time
import sourcecode.core.matrixhelper as mxh

home_folder = '/nethome/manra003/data/'

depth = 0


def get_empty_array(mons, no_grids):
    return np.zeros((no_grids, no_grids + 2, mons))


def get_dense_matrix(matrix_data, indices, indptr, no_grids):
    return csr_matrix((matrix_data, indices, indptr), shape=(no_grids, no_grids + 2)).todense()


def compute_bin_matrix(matrix, no_grids, path):
    avg_matrix = np.mean(matrix, axis=2)
    sp_matrix = csr_matrix(avg_matrix)
    bin_matrix = mxh.binary_matrix(sp_matrix)
    print(np.min(bin_matrix), np.max(bin_matrix))
    print(np.min(bin_matrix[:, :no_grids]), np.max(bin_matrix[:, :no_grids]))
    # cannot perform the follwoing assert- as some grids are not selected at all in the bootstrap random sampling.
    # assert np.sum(bin_matrix, axis=1).all() == 1
    save_npz(path + 'Annual_Avg_FullAdjacency_csr.npz', bin_matrix)
    save_npz(path + 'Annual_Avg_DomainAdjacency_csr.npz', bin_matrix[:, :no_grids])


def compute_grid_nnz_average(matrix, nnz_count_matrix, matrix_type, path):
    sum_matrix = np.sum(matrix, axis=2)
    avg_matrix = sum_matrix / nnz_count_matrix
    print(matrix_type, np.nanmin(avg_matrix), np.nanmax(avg_matrix))
    save_npz(path + 'Annual_Avg_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))


def main():
    master_uni_hex = np.load(home_folder + 'MasterHexList_Res3.npy').tolist()
    assert len(master_uni_hex) == 8243
    no_grids = len(master_uni_hex)
    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    month_count = len(months)

    # Saving format: from MonthlySparseTM.py
    # np.savez_compressed(export_folder + 'CSR_{0}.npz'.format(mon), transprob=norm_matrix.data,
    #                         mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
    #                         maxsal=avg_max_sal_per_grid, indices=mon_trans_matrix.indices,
    #                         indptr=mon_trans_matrix.indptr)

    sample_size = [5000, 10000, 50000, 100000, 200000, 300000]
    init_state = 21
    end_state = 30
    
    for size in sample_size:
        t1 = time()
        
        for state in np.arange(init_state, end_state + 1):
    
            input_path = home_folder + 'Boot_Sample/Size_{0}/Monthly/State_{1}/'.format(size, state)
            output_path = home_folder + 'Boot_Sample/Size_{0}/Annual/State_{1}/'.format(size, state)
            os.makedirs(output_path, exist_ok=True)
            monthly_full_trans_prob = get_empty_array(month_count, no_grids)
            monthly_full_min_temp = get_empty_array(month_count, no_grids)
            monthly_full_max_temp = get_empty_array(month_count, no_grids)
            monthly_full_min_sal = get_empty_array(month_count, no_grids)
            monthly_full_max_sal = get_empty_array(month_count, no_grids)

            for mon, i in zip(months, range(len(months))):
                with np.load(input_path + 'CSR_{0}_z{1}.npz'.format(mon, depth)) as data:
                    monthly_full_trans_prob[:, :, i] = get_dense_matrix(data['transprob'], data['indices'], data['indptr'],
                                                                    no_grids)
                    monthly_full_min_temp[:, :, i] = get_dense_matrix(data['mintemp'], data['indices'], data['indptr'],
                                                                  no_grids)
                    monthly_full_max_temp[:, :, i] = get_dense_matrix(data['maxtemp'], data['indices'], data['indptr'],
                                                                  no_grids)
                    monthly_full_min_sal[:, :, i] = get_dense_matrix(data['minsal'], data['indices'], data['indptr'],
                                                                 no_grids)
                    monthly_full_max_sal[:, :, i] = get_dense_matrix(data['maxsal'], data['indices'], data['indptr'],
                                                                 no_grids)

            # Compute normalized avg or sum of transitions
#         compute_grid_annual_average(monthly_full_trans_prob, no_grids, output_path)
            # compute bin matrix -new approach
            compute_bin_matrix(monthly_full_trans_prob, no_grids, output_path)
        
            # remaining without the new column(8245) and deleted(8246)
            nnz_count_matrix = np.count_nonzero(monthly_full_trans_prob[:, :no_grids, :], axis=2)
            nnz_count_matrix[nnz_count_matrix == 0] = 1
            compute_grid_nnz_average(monthly_full_min_temp[:, :no_grids, :], nnz_count_matrix, 'MinTemperature',
                                 output_path)
            compute_grid_nnz_average(monthly_full_max_temp[:, :no_grids, :], nnz_count_matrix, 'MaxTemperature',
                                 output_path)
            compute_grid_nnz_average(monthly_full_min_sal[:, :no_grids, :], nnz_count_matrix, 'MinSalinity',
                                 output_path)
            compute_grid_nnz_average(monthly_full_max_sal[:, :no_grids, :], nnz_count_matrix, 'MaxSalinity',
                                 output_path)
        t2 = time()
        print('Sample size {0} completed in time: {1}'.format(size, t2 - t1))

if __name__ == '__main__':
    main()
