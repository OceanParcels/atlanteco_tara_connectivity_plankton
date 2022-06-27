import numpy as np
from scipy.sparse import csr_matrix, save_npz
import sys

home_folder = '/nethome/manra003/atlanteco_tara_connectivity_plankton/data/matrices/'
# home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'


def get_csr_matrix(matrix_data, indices, indptr, no_grids):
    return csr_matrix((matrix_data, indices, indptr), shape=(no_grids, no_grids + 2))


def compute_binary_matrix(matrix, no_grids):
    matrix[matrix > 0] = 1
    save_npz(home_folder + 'Annual_Binary_FullAdjacency_csr.npz', csr_matrix(matrix))
    save_npz(home_folder + 'Annual_Binary_DomainAdjacency_csr.npz', csr_matrix(matrix[:, :no_grids]))


def compute_grid_nnz_average(matrix, nnz_count_matrix, matrix_type):
    avg_matrix = matrix / nnz_count_matrix
    print(matrix_type, np.nanmin(avg_matrix), np.nanmax(avg_matrix), np.nanmean(avg_matrix))
    save_npz(home_folder + 'Annual_Avg_{0}_csr.npz'.format(matrix_type), csr_matrix(avg_matrix))


def main():
    args = sys.argv
    assert len(args) == 2
    sim_depth = np.int32(args[1])
    assert 0 <= sim_depth <= 500

    global home_folder
    home_folder = home_folder + 't{0}m/'.format(sim_depth)

    no_grids = len(np.load('/nethome/manra003/data/MasterHexList_Res3.npy').tolist())
    # no_grids = len(np.load('/Users/dmanral/Desktop/Analysis/TARA/Task7D/MasterHexList_Res3.npy').tolist())
    months = np.array(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    monthly_sum_trans = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)
    monthly_sum_min_mintemp = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)
    monthly_sum_max_mintemp = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)
    monthly_sum_min_maxtemp = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)
    monthly_sum_max_maxtemp = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)
    monthly_sum_avg_mintemp = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)
    monthly_sum_avg_maxtemp = csr_matrix((no_grids, no_grids + 2), dtype=np.float32)

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
        with np.load(home_folder + 'CSR_{0}z{1}.npz'.format(mon, sim_depth)) as data:
            monthly_sum_trans = monthly_sum_trans + get_csr_matrix(data['transprob'], data['indices'],
                                                                   data['indptr'], no_grids)
            monthly_sum_min_mintemp = monthly_sum_min_mintemp + get_csr_matrix(data['min_mintemp'], data['indices'],
                                                                               data['indptr'], no_grids)
            monthly_sum_max_mintemp = monthly_sum_max_mintemp + get_csr_matrix(data['max_mintemp'], data['indices'],
                                                                               data['indptr'], no_grids)
            monthly_sum_min_maxtemp = monthly_sum_min_maxtemp + get_csr_matrix(data['min_maxtemp'], data['indices'],
                                                                               data['indptr'], no_grids)
            monthly_sum_max_maxtemp = monthly_sum_max_maxtemp + get_csr_matrix(data['max_maxtemp'], data['indices'],
                                                                               data['indptr'], no_grids)
            monthly_sum_avg_mintemp = monthly_sum_avg_mintemp + get_csr_matrix(data['avg_min_temp'], data['indices'],
                                                                               data['indptr'], no_grids)
            monthly_sum_avg_maxtemp = monthly_sum_avg_maxtemp + get_csr_matrix(data['avg_max_temp'], data['indices'],
                                                                               data['indptr'], no_grids)

    compute_binary_matrix(monthly_sum_trans.todense(), no_grids)
    # remaining without the new column(8245) and deleted(8246)
    # compute_grid_annual_average(monthly_full_trans_prob[:, :, :no_grids], nnz_count_matrix, 'DomainAdjacency')
    nnz_count_matrix = monthly_sum_trans.todense()[:, :no_grids]
    nnz_count_matrix[nnz_count_matrix == 0] = 1
    compute_grid_nnz_average(monthly_sum_min_mintemp.todense()[:, :no_grids], nnz_count_matrix, 'min_MinTemperature')
    compute_grid_nnz_average(monthly_sum_max_mintemp.todense()[:, :no_grids], nnz_count_matrix, 'max_MinTemperature')
    compute_grid_nnz_average(monthly_sum_min_maxtemp.todense()[:, :no_grids], nnz_count_matrix, 'min_MaxTemperature')
    compute_grid_nnz_average(monthly_sum_max_maxtemp.todense()[:, :no_grids], nnz_count_matrix, 'max_MaxTemperature')
    compute_grid_nnz_average(monthly_sum_avg_mintemp.todense()[:, :no_grids], nnz_count_matrix, 'avg_MinTemperature')
    compute_grid_nnz_average(monthly_sum_avg_maxtemp.todense()[:, :no_grids], nnz_count_matrix, 'avg_MaxTemperature')


if __name__ == '__main__':
    main()
