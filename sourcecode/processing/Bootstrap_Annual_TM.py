import numpy as np
from scipy.sparse import csr_matrix, save_npz
import os
from time import time
import sourcecode.core.matrixhelper as mxh
import sys


data_folder= '/nethome/manra003/analysis/paper01/'

def get_csr_matrix(matrix_data, indices, indptr, no_grids):
    return csr_matrix((matrix_data, indices, indptr), shape=(no_grids, no_grids + 1))


def compute_binary_matrix(matrix, no_grids, depth, output_path):
    matrix[matrix > 0] = 1
    save_npz(output_path + 'Annual_Binary_DomainAdjacency_z{0}_csr.npz'.format(depth), csr_matrix(matrix[:, :no_grids]))


def save_matrix(matrix, matrix_type, depth, output_path):
    print(matrix_type, np.nanmin(matrix), np.nanmax(matrix), np.nanmean(matrix))
    save_npz(output_path + 'Annual_{0}_z{1}_csr.npz'.format(matrix_type,depth), csr_matrix(matrix))


def main():
    args = sys.argv
    assert len(args) == 4
    sim_depth = np.int32(args[1])
    assert 0 <= sim_depth <= 500
    
    init_state = np.int32(args[2])
    end_state = np.int32(args[3])
    
    master_uni_hex = np.load(data_folder + 'H3_Res3_MasterHexList.npz')['Res3_HexId'].tolist()
    assert len(master_uni_hex) == 8191
    no_grids = len(master_uni_hex)
    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    month_count = len(months)

    # Saving format: from MonthlySparseTM.py
    # np.savez_compressed(export_folder + 'CSR_{0}.npz'.format(mon), transprob=norm_matrix.data,
    #                         mintemp=avg_min_temp_per_grid, maxtemp=avg_max_temp_per_grid, minsal=avg_min_sal_per_grid,
    #                         maxsal=avg_max_sal_per_grid, indices=mon_trans_matrix.indices,
    #                         indptr=mon_trans_matrix.indptr)

    sample_size = [5000, 10000, 50000, 100000, 200000, 300000]
    
    for size in sample_size:
        t1 = time()
        
        for state in np.arange(init_state, end_state+1):
            print('State ID: ',state)
            input_path = data_folder + 'Boot_Sample/Size_{0}/Monthly/State_{1}/'.format(size, state)
            output_path = data_folder + 'Boot_Sample/Size_{0}/Annual/State_{1}/'.format(size, state)
            os.makedirs(output_path, exist_ok=True)
            monthly_sum_trans = np.zeros((no_grids, no_grids + 1), dtype=np.int32)
            monthly_min_mintemp = np.full((no_grids, no_grids + 1), 999, dtype=np.float32)
            monthly_max_mintemp = np.full((no_grids, no_grids + 1), -999, dtype=np.float32)
            monthly_min_maxtemp = np.full((no_grids, no_grids + 1), 999, dtype=np.float32)
            monthly_max_maxtemp = np.full((no_grids, no_grids + 1), -999, dtype=np.float32)

            # export all matrices to npz file
            # np.savez_compressed(output_path + 'CSR_{0}z{1}.npz'.format(mon, sim_depth),
            #                     transprob=mon_trans_matrix.data,
            #                     min_mintemp=mon_min_mintemp_matrix.data,
            #                     max_mintemp=mon_max_mintemp_matrix.data,
            #                     min_maxtemp=mon_min_maxtemp_matrix.data,
            #                     max_maxtemp=mon_max_maxtemp_matrix.data,
            #                     min_minsal=mon_min_minsal_matrix.data,
            #                     max_minsal=mon_max_minsal_matrix.data,
            #                     min_maxsal=mon_min_maxsal_matrix.data,
            #                     max_maxsal=mon_max_maxsal_matrix.data,
            #                     avg_min_temp=avg_min_temp_per_grid,
            #                     avg_max_temp=avg_max_temp_per_grid,
            #                     avg_min_sal=avg_min_sal_per_grid,
            #                     avg_max_sal=avg_max_sal_per_grid,
            #                     indices=mon_trans_matrix.indices,
            #                     indptr=mon_trans_matrix.indptr)

            for mon, i in zip(months, range(len(months))):
                with np.load(input_path + 'CSR_{0}z{1}.npz'.format(mon, sim_depth)) as data:
                    mon_arr = get_csr_matrix(data['transprob'], data['indices'], data['indptr'], no_grids).todense()
                    mon_arr[mon_arr > 0] = 1
                    monthly_sum_trans = monthly_sum_trans + mon_arr

                    mon_arr = get_csr_matrix(data['min_mintemp'], data['indices'], data['indptr'], no_grids).todense()
                    mon_arr[mon_arr == 0] = 999
                    monthly_min_mintemp = np.minimum(monthly_min_mintemp, mon_arr)

                    mon_arr = get_csr_matrix(data['max_mintemp'], data['indices'], data['indptr'], no_grids).todense()
                    mon_arr[mon_arr == 0] = -999
                    monthly_max_mintemp = np.maximum(monthly_max_mintemp, mon_arr)

                    mon_arr = get_csr_matrix(data['min_maxtemp'], data['indices'], data['indptr'], no_grids).todense()
                    mon_arr[mon_arr == 0] = 999
                    monthly_min_maxtemp = np.minimum(monthly_min_maxtemp, mon_arr)

                    mon_arr = get_csr_matrix(data['max_maxtemp'], data['indices'], data['indptr'], no_grids).todense()
                    mon_arr[mon_arr == 0] = -999
                    monthly_max_maxtemp = np.maximum(monthly_max_maxtemp, mon_arr)

                    
            monthly_min_mintemp[monthly_min_mintemp == 999] = 0
            monthly_max_mintemp[monthly_max_mintemp == -999] = 0
            monthly_min_maxtemp[monthly_min_maxtemp == 999] = 0
            monthly_max_maxtemp[monthly_max_maxtemp == -999] = 0

            # verify the shape of all files is same
            assert np.array_equal(monthly_sum_trans.nonzero(), monthly_max_mintemp.nonzero())
            assert np.array_equal(monthly_min_mintemp.nonzero(), monthly_max_mintemp.nonzero())
            assert np.array_equal(monthly_min_mintemp.nonzero(), monthly_min_maxtemp.nonzero())
            assert np.array_equal(monthly_max_maxtemp.nonzero(), monthly_min_maxtemp.nonzero())

            # remaining without the new column(8244) and deleted(8245)
            save_matrix(monthly_min_mintemp[:, :no_grids], 'min_MinTemperature', sim_depth, output_path)
            save_matrix(monthly_max_mintemp[:, :no_grids], 'max_MinTemperature', sim_depth, output_path)
            save_matrix(monthly_min_maxtemp[:, :no_grids], 'min_MaxTemperature', sim_depth, output_path)
            save_matrix(monthly_max_maxtemp[:, :no_grids], 'max_MaxTemperature', sim_depth, output_path)

            compute_binary_matrix(monthly_sum_trans, no_grids, sim_depth, output_path)
            assert np.array_equal(monthly_sum_trans.nonzero(), monthly_max_mintemp.nonzero())
        
        t2 = time()
        print('Sample size {0} completed in time: {1}'.format(size, t2 - t1))

if __name__ == '__main__':
    main()
