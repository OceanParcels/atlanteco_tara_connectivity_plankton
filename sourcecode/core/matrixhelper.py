import numpy as np
import h3
from scipy.sparse import coo_matrix


def get_hexid_from_parent(lons, lats, child_res, parent_res):
    child_hex = np.array([h3.geo_to_h3(y, x, child_res) for x, y in zip(lons, lats)])
    return np.array([h3.h3_to_parent(str(h), parent_res) for h in child_hex])


def get_coo_matrix(array, rows, cols, no_grids):
    matrix = coo_matrix((array, (rows, cols)), shape=(no_grids, no_grids + 2))
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

    join_arrays = np.concatenate((minsalonly, mintemponly, maxtemponly, static_pts_index))
    return np.unique(join_arrays)
