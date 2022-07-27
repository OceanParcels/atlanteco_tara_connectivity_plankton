import pandas as pd
import numpy as np
from sourcecode.core import adjacencygraph as ag
from sourcecode.core import connectivityhelper as ch
import os

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task9B/'
dataset = 'NoConstraints'  # 'sample_constraints'  # 'NoConstraints'  '2011_lombard_forams'
out_folder = home_folder + 'Connectivities/{0}/'.format(dataset)

hex_res = 3
depth = 500


def get_stationcode_hexes_mapping():
    stations = pd.read_csv(home_folder + 'Tara_Stations_hexId_sorted.csv', header=0, index_col=0)
    return stations.index.values, stations['H3Id_res3'].values


def full_connectivity(species, min_accept_temp, max_accept_temp, temp_constraint_range, len_stations,
                      final_stations_master_indices, final_stations_code, width_type):
    domain_adjacency_file = 't{0}m/Annual_Binary_DomainAdjacency_z{0}_csr.npz'.format(depth)
    # to add constraints to the connectivity

    if width_type == 'average':
        min_temp_file = home_folder + 't{0}m/Annual_avg_MinTemperature_z{0}_csr.npz'.format(depth)
        max_temp_file = home_folder + 't{0}m/Annual_avg_MaxTemperature_z{0}_csr.npz'.format(depth)
    # elif width_type == 'narrow':
    #     min_temp_file = home_folder + 't{0}m/Annual_max_MinTemperature_csr.npz'.format(depth)
    #     max_temp_file = home_folder + 't{0}m/Annual_min_MaxTemperature_csr.npz'.format(depth)
    elif width_type == 'broad':
        min_temp_file = home_folder + 't{0}m/Annual_min_MinTemperature_z{0}_csr.npz'.format(depth)
        max_temp_file = home_folder + 't{0}m/Annual_max_MaxTemperature_z{0}_csr.npz'.format(depth)
    elif width_type == 'passive':
        min_temp_file = None
        max_temp_file = None
    else:
        raise ValueError("width type is not recognized")
    print(width_type)
    # create graph
    if ~np.isnan(temp_constraint_range):
        print('Temp range: ', temp_constraint_range)
        atlantic_graph = ag.create_temp_range_graph(home_folder + domain_adjacency_file,
                                                    min_temp_file,
                                                    max_temp_file,
                                                    temp_constraint_range)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        print('Min/Max Temp: ', min_accept_temp, max_accept_temp, species)
        atlantic_graph = ag.create_temp_min_max_graph(home_folder + domain_adjacency_file,
                                                      min_temp_file,
                                                      max_temp_file,
                                                      min_accept_temp, max_accept_temp)

    else:
        atlantic_graph = ag.create_simple_graph(home_folder + domain_adjacency_file)
    # get min_T paths for all pairs- forward and backward- 2 d matrix.
    min_T_matrix = np.empty((len_stations, len_stations))
    min_T_matrix[:] = np.NAN

    station_count = len(final_stations_master_indices)
    nnz_nan_count = 0
    zero_count = 0

    for i in range(station_count):
        s_idx = final_stations_master_indices[i]
        for j in range(i, station_count):
            d_idx = final_stations_master_indices[j]
            if s_idx != d_idx:
                f_path = ag.get_shortest_path(atlantic_graph, s_idx, d_idx)
                # f_path = ag.get_most_probable_path(atlantic_graph, s_idx, d_idx)
                # f_path = ag.minimum_time_path_malley2021(atlantic_graph, s_idx, d_idx)
                if f_path:
                    min_T_matrix[i][j] = len(f_path) - 1  # (time in months)
                    # min_T_matrix[i][j] = ag.get_time_from_most_probable_path(atlantic_graph, f_path)[-1]
                    nnz_nan_count += 1

                b_path = ag.get_shortest_path(atlantic_graph, d_idx, s_idx)
                # b_path = ag.get_most_probable_path(atlantic_graph, d_idx, s_idx)
                # b_path = ag.minimum_time_path_malley2021(atlantic_graph, d_idx, s_idx)
                if b_path:
                    min_T_matrix[j][i] = len(b_path) - 1  # (time in months)
                    # min_T_matrix[j][i] = ag.get_time_from_most_probable_path(atlantic_graph, b_path)[-1]
                    nnz_nan_count += 1
            else:
                if ag.check_if_edge_exists(atlantic_graph, s_idx, d_idx):
                    if i != j:
                        min_T_matrix[j][i] = 0
                        zero_count += 1
                    min_T_matrix[i][j] = 0
                    zero_count += 1

    print('maximum time (months): ', np.nanmax(min_T_matrix))
    print('Non zero-non NAN counter: ', nnz_nan_count)
    print('zero counter: ', zero_count)
    print("NAN count: ", np.count_nonzero(np.isnan(min_T_matrix)))
    print("non_zero count: ", np.count_nonzero(min_T_matrix))
    print("zero count: ", len(np.where(min_T_matrix == 0)[0]), "\n************************")
    assert nnz_nan_count + np.count_nonzero(np.isnan(min_T_matrix)) + len(np.where(min_T_matrix == 0)[0]) == \
           np.square(len_stations)
    path = out_folder + 't{0}m/{1}/'.format(depth, width_type)
    os.makedirs(path, exist_ok=True)
    np.savez_compressed(path + 'Stations_minT_connectivity_{0}z_{1}_{2}.npz'.format(depth, species, width_type),
                        codes=final_stations_code, matrix=min_T_matrix)


def main():
    # Get all sorted stations- station codes between -100 and 20 Longitude
    stations_code, stations_hex = get_stationcode_hexes_mapping()

    master_hex_ids = ch.get_all_grids_hex_ids(home_folder + 'MasterHexList_Res3.npy')

    # map station to master hex (some stations lie in the same hex- same connectivity)
    domain_mask = np.in1d(stations_hex, master_hex_ids)
    final_stations_code = stations_code[domain_mask]
    final_stations_hex = stations_hex[domain_mask]

    final_stations_master_indices = np.array([master_hex_ids.index(h) for h in final_stations_hex])
    # not use this as it sorts the results and return unique: np.in1d(master_hex_ids,stations_hex).nonzero()[0]

    # for no constraints- try any width setting- only makes use of Binary connectivity matrix.
    full_connectivity('NoConstraints', np.NAN, np.NAN, np.NAN, len(final_stations_hex), final_stations_master_indices,
                      final_stations_code, 'passive')

    # species_info = pd.read_csv(home_folder + dataset + '.csv',
    #                            delimiter=';|,', keep_default_na=True, header=0, engine='python')
    #
    # print("AVERAGE WIDTH\n-----------------")
    # [full_connectivity(entry['Species'], entry['MinTemp'], entry['MaxTemp'], entry['TRange'], len(final_stations_hex),
    #                    final_stations_master_indices, final_stations_code, 'average')
    #  for index, entry in species_info.iterrows()]
    #
    # # print("-----------------\nNARROW WIDTH\n-----------------")
    # # Not working - min_maxtemp-max_mintemp <0 hence, cannot work.
    # # [full_connectivity(entry['Species'], entry['MinTemp'], entry['MaxTemp'], np.NAN, len(final_stations_hex),
    # #                    final_stations_master_indices, final_stations_code, 'narrow')
    # #  for index, entry in species_info.iterrows()]
    #
    # print("-----------------\nBROAD WIDTH\n-----------------")
    # [full_connectivity(entry['Species'], entry['MinTemp'], entry['MaxTemp'], entry['TRange'], len(final_stations_hex),
    #                    final_stations_master_indices, final_stations_code, 'broad')
    #  for index, entry in species_info.iterrows()]


if __name__ == '__main__':
    main()
