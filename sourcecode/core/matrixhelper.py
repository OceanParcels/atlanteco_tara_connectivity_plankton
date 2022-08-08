import numpy as np
import h3
from scipy.sparse import coo_matrix, csr_matrix

NEW = 'new'


def main_hex_list(m_hex):
    """
    hex list needed to process the mapping of the transitions.
    :param hex:
    :return:
    """
    return np.append(m_hex, (NEW))


def get_hexid_from_parent(lons, lats, child_res, parent_res):
    child_hex = np.array([h3.geo_to_h3(y, x, child_res) for x, y in zip(lons, lats)])
    return np.array([h3.h3_to_parent(str(h), parent_res) for h in child_hex])


def get_coo_matrix(array, rows, cols, no_grids):
    matrix = coo_matrix((array, (rows, cols)), shape=(no_grids, no_grids + 1))
    matrix.sum_duplicates()
    return matrix


def check_default_values(ds):
    """
    Quick validation to confirm the default values of 999 and -999 were updated in the simulation
    :param ds:
    :return: true if no invalid values found
    """

    if -999 in ds['max_temp'][:, -1].values or 999 in ds['max_temp'][:, -1].values or \
            -999 in ds['min_temp'][:, -1].values or 999 in ds['min_temp'][:, -1].values or \
            -999 in ds['max_sal'][:, -1].values or 999 in ds['max_sal'][:, -1].values or \
            -999 in ds['min_sal'][:, -1].values or 999 in ds['min_sal'][:, -1].values:
        return False
    return True


def get_invalid_trajectories(ds):
    """
    Particles that never moved or those that got stuck on land during simulation-
    :param ds:
    :return: indices of invalid trajectories
    """

    # static points analysis- points that did not change their position at all
    static_pts_index = \
        np.where(np.logical_and((ds['lat'][:, 0] == ds['lat'][:, -1]), (ds['lon'][:, 0] == ds['lon'][:, -1])))[0]
    print("Static Points count: ", len(static_pts_index))

    # all the deleted points get assigned values of nan for all the fields, hence need to remove it
    # using latest version of parcels deployment on lorenz-master branch- ~6-7 August 2022
    deleted_pts_index = np.where(np.logical_or(np.isnan(ds['lat'][:, -1]), np.isnan(ds['lon'][:, -1])))[0]
    print("Deleted Points count: ", len(deleted_pts_index))
    
    # points with 0 fields of temp or salinity
    maxsal_zero_index = np.where(ds['max_sal'][:, -1] == 0)[0]
    print("Zero MAX Salinity count: ", len(maxsal_zero_index))
    minsal_zero_index = np.where(ds['min_sal'][:, -1] == 0)[0]
    print("Zero MIN Salinity count: ", len(minsal_zero_index))
    # only possible scenario is when a particle will get a lower salinity-when stuck

    maxtemp_zero_index = np.where(ds['max_temp'][:, -1] == 0)[0]
    print("Zero MAX Temperature count: ", len(maxtemp_zero_index))
    mintemp_zero_index = np.where(ds['min_temp'][:, -1] == 0)[0]
    print("Zero MIN Temperature count: ", len(mintemp_zero_index))
    # here both scenarios are possible, particle with >0 min_temperature and <0 max_temp
    # it is possible to find unique scenarios when maz_temp was below 0 and hex
    # when it gets stuck the max_temp gets updated to zero

    assert np.array_equal(static_pts_index, maxsal_zero_index)

    minsalonly = np.setdiff1d(minsal_zero_index, static_pts_index)
    mintemponly = np.setdiff1d(mintemp_zero_index, static_pts_index)
    maxtemponly = np.setdiff1d(maxtemp_zero_index, static_pts_index)

    join_arrays = np.concatenate((minsalonly, mintemponly, maxtemponly, static_pts_index, deleted_pts_index))
    return np.unique(join_arrays)


# def get_monthly_tm(ds, invalid_indices, trans_array, rows_array, cols_array, min_temp_array, max_temp_array,
#                    min_sal_array, max_sal_array, no_grids, parent_res, child_res, hex_indices, map_h3_to_mat
#                    ):
#     def get_valid_data(field_name, loc):
#         return np.delete(ds[field_name][:, loc].values, invalid_indices)

#     hex_t0 = get_hexid_from_parent(get_valid_data('lon', 0), get_valid_data('lat', 0), child_res,
#                                    parent_res)
#     # assert np.array_equal(hex_t0, master_all_hex_t0)
#     hex_t1 = get_hexid_from_parent(get_valid_data('lon', -1), get_valid_data('lat', -1), child_res,
#                                    parent_res)

#     # mask hex ids that are new
#     hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, NEW)
#     # mask hex ids in hex_t1_new that were deleted during the simulation
#     hex_t1_new = np.where(get_valid_data('time', -1) < np.max(ds['time'][:, -1].values), DEL, hex_t1_new)

#     rows = map_h3_to_mat[hex_t0].values
#     cols = map_h3_to_mat[hex_t1_new].values
#     transitions = np.ones((len(hex_t0)))

#     t_matrix = get_coo_matrix(transitions, rows, cols, no_grids)
#     trans_array = np.append(trans_array, t_matrix.data)
#     rows_array = np.append(rows_array, t_matrix.row)
#     cols_array = np.append(cols_array, t_matrix.col)

#     # get min and max temperature data
#     min_temperature, max_temperature = get_valid_data('min_temp', -1), get_valid_data('max_temp', -1)

#     min_temp_matrix = get_coo_matrix(min_temperature, rows, cols, no_grids)
#     max_temp_matrix = get_coo_matrix(max_temperature, rows, cols, no_grids)
#     min_temp_array = np.append(min_temp_array, min_temp_matrix.data)
#     max_temp_array = np.append(max_temp_array, max_temp_matrix.data)

#     # get min and max salinity data
#     min_salinity, max_salinity = get_valid_data('min_sal', -1), get_valid_data('max_sal', -1)

#     min_sal_matrix = get_coo_matrix(min_salinity, rows, cols, no_grids)
#     max_sal_matrix = get_coo_matrix(max_salinity, rows, cols, no_grids)
#     min_sal_array = np.append(min_sal_array, min_sal_matrix.data)
#     max_sal_array = np.append(max_sal_array, max_sal_matrix.data)

#     return trans_array, rows_array, cols_array, min_temp_array, max_temp_array, min_sal_array, max_sal_array


# def get_monthly_matrix(ds, invalid_indices, trans_array, min_temp_array, max_temp_array,
#                        min_sal_array, max_sal_array, no_grids, parent_res, child_res, hex_indices, map_h3_to_mat
#                        ):
#     def get_valid_data(field_name, loc):
#         return np.delete(ds[field_name][:, loc].values, invalid_indices)

#     hex_t0 = get_hexid_from_parent(get_valid_data('lon', 0), get_valid_data('lat', 0), child_res,
#                                    parent_res)
#     # assert np.array_equal(hex_t0, master_all_hex_t0)
#     hex_t1 = get_hexid_from_parent(get_valid_data('lon', -1), get_valid_data('lat', -1), child_res,
#                                    parent_res)

#     # mask hex ids that are new
#     hex_t1_new = np.where(np.isin(hex_t1, hex_indices), hex_t1, NEW)
#     # mask hex ids in hex_t1_new that were deleted during the simulation
#     hex_t1_new = np.where(get_valid_data('time', -1) < np.max(ds['time'][:, -1].values), DEL, hex_t1_new)

#     rows = map_h3_to_mat[hex_t0].values
#     cols = map_h3_to_mat[hex_t1_new].values

#     trans_array[rows, cols] = 1

#     # get min and max temperature data
#     min_t = np.full(min_temp_array.shape, -999, dtype='float')
#     max_t = np.full(max_temp_array.shape, 999, dtype='float')

#     min_temperature, max_temperature = get_valid_data('min_temp', -1), get_valid_data('max_temp', -1)
#     min_t[rows, cols] = min_temperature
#     min_temp_array = np.where(min_temp_array < min_t, min_t, min_temp_array)
#     max_t[rows, cols] = max_temperature
#     max_temp_array = np.where(max_temp_array > max_t, max_t, max_temp_array)

#     # get min and max salinity data
#     min_s = np.full(min_sal_array.shape, -999, dtype='float')
#     max_s = np.full(max_sal_array.shape, 999, dtype='float')

#     min_salinity, max_salinity = get_valid_data('min_temp', -1), get_valid_data('max_temp', -1)
#     min_s[rows, cols] = min_salinity
#     min_sal_array = np.where(min_sal_array < min_s, min_s, min_sal_array)
#     max_s[rows, cols] = max_salinity
#     max_sal_array = np.where(max_sal_array > max_s, max_s, max_sal_array)

#     return trans_array, min_temp_array, max_temp_array, min_sal_array, max_sal_array


def binary_matrix(matrix):
    """
    get binary matrix : replacing all exisitng connections with 1
    :param matrix:
    :return:
    """
    return csr_matrix((np.ones(len(matrix.data)), matrix.indices, matrix.indptr),
                      shape=matrix.shape)


# def avg_field_per_grid(f_data, t_data, f_type, field):
#     avg_field = f_data / t_data
#     print('Min/Max average {0} {1}: {2} / {3}'.format(f_type, field, np.min(avg_field), np.max(avg_field)))
#     return avg_field
