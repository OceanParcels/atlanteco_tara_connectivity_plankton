import pandas as pd
import numpy as np
from code import adjacencygraph as ag
from code import connectivityhelper as ch

# home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task6_Sens/Resolution2/'
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
data_folder = home_folder + 'Full_connectivity_output/2005Zaric/'

hex_res = 3


# def get_stationcode_hexes_mapping():
#     stations = pd.read_excel(home_folder + 'AllS  tations_Tara.xls', header=1, index_col=0)
#     filtered_stations = stations[(stations['Longitude'] > -100) & (stations['Longitude'] < 20)]
#     # sort stations in order of their latitudes
#     sorted_stations = filtered_stations.sort_values(by=['Latitude'])
#     hex_ids = np.array([h3.geo_to_h3(lat, lon, hex_res) for lat, lon in
#                         zip(sorted_stations['Latitude'], sorted_stations['Longitude'])])
#     return sorted_stations.index.values, hex_ids


def get_stationcode_hexes_mapping():
    stations = pd.read_csv(home_folder + 'Tara_Stations_hexId_sorted.csv', header=0, index_col=0)
    # sort stations in order of their latitudes
    # sorted_stations = stations.sort_values(by=['Latitude'])
    # hex_ids = np.array([h3.geo_to_h3(lat, lon, hex_res) for lat, lon in
    #                         zip(sorted_stations['Latitude'], sorted_stations['Longitude'])])
    # return sorted_stations.index.values, hex_ids
    # sorted_stations = stations.sort_values(by=['Latitude'], ascending=False)

    return stations.index.values, stations['H3Id_res3'].values


# def get_all_grids_hex_ids():
# return pd.read_csv(home_folder + 'MasterHexList_Res2.csv', header=None).values
#     # hexId and matrix index mapping from a random simulation output file from the output dataset
# with xr.open_dataset(np.random.choice(glob(home_folder + '/tara_res5_01/FullTara_Res5_TS_*'))).load() as temp_ds:
# with xr.open_dataset(
#         '/Users/dmanral/Desktop/Analysis/TARA/Task4/tara_res5_01/FullTara_Res5_TS_Aug2015_dt600.nc').load() as temp_ds:
#     master_all_hex_t0 = np.array(
#         [h3.geo_to_h3(y, x, hex_res) for x, y in zip(temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values)])
#     master_uni_hex, counts = np.unique(master_all_hex_t0, return_counts=True)
#     return master_uni_hex.tolist()


def fullconnectivity(species, mintemp, maxtemp, len_stations, final_stations_master_indices, final_stations_code):
    # to add constraints to the connectivity
    # domain_adjacency_file = 'ProcessedTM/Annual_Avg_Prob_Cutoff_pt001_csr.npz'
    # domain_adjacency_file = 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz'
    domain_adjacency_file = 'ProcessedTM/Annual_Avg_DomainAdjacency_csr.npz'

    temp_constraint_range = np.NaN
    # min_accept_temp = np.NAN
    # max_accept_temp = np.NAN
    min_accept_temp = mintemp
    max_accept_temp = maxtemp
    prob_cutoff = np.NAN
    trahms2021 = False
    minimum_time_omalley2021 = False

    if minimum_time_omalley2021:
        t_ratio_file = home_folder + 'ProcessedTM/MinTRatio_FullAdjacency_csr.npz'
    else:
        t_ratio_file = None

    # create graph
    if ~np.isnan(temp_constraint_range):
        print('Temp range: ', temp_constraint_range)
        atlantic_graph = ag.create_temp_range_graph(home_folder + domain_adjacency_file,
                                                    home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                    home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                    temp_constraint_range,
                                                    t_ratio_file)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        print('Min/Max Temp: ', min_accept_temp, max_accept_temp)
        atlantic_graph = ag.create_temp_min_max_graph(home_folder + domain_adjacency_file,
                                                      home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                      home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                      min_accept_temp, max_accept_temp,
                                                      t_ratio_file)
    elif ~np.isnan(prob_cutoff):
        print('Probability Cutoff Value: ', prob_cutoff)
        atlantic_graph = ag.create_prob_filtered_graph(home_folder + domain_adjacency_file,
                                                       home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                                       home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                                       prob_cutoff)
    elif trahms2021:
        atlantic_graph = ag.create_full_graph(home_folder + domain_adjacency_file,
                                              home_folder + 'ProcessedTM/Annual_SUM_MinTemperature_csr.npz',
                                              home_folder + 'ProcessedTM/Annual_SUM_MaxTemperature_csr.npz')

    else:
        atlantic_graph = ag.create_full_graph(home_folder + domain_adjacency_file,
                                              home_folder + 'ProcessedTM/Annual_Avg_MinTemperature_csr.npz',
                                              home_folder + 'ProcessedTM/Annual_Avg_MaxTemperature_csr.npz',
                                              t_ratio_file)
        # atlantic_graph = ag.create_simple_graph(home_folder + domain_adjacency_file,
        #                                         t_ratio_file)

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

    print('maximum time: ', np.nanmax(min_T_matrix))
    print('Non zero-non NAN counter: ', nnz_nan_count)
    print('zero counter: ', zero_count)
    print("NAN count: ", np.count_nonzero(np.isnan(min_T_matrix)))
    print("non_zero count: ", np.count_nonzero(min_T_matrix))
    print("zero count: ", len(np.where(min_T_matrix == 0)[0]))
    assert nnz_nan_count + np.count_nonzero(np.isnan(min_T_matrix)) + len(np.where(min_T_matrix == 0)[0]) == \
           np.square(len_stations)
    np.savez_compressed(
        data_folder + 'Stations_MinT_connectivity_{0}.npz'.format(
            species),
        codes=final_stations_code, matrix=min_T_matrix)


def main():
    species_info = pd.read_csv(data_folder + '2005_zaric_forams.csv',
                               delimiter=';|,', header=0)
    # Get all sorted stations- station codes between -100 and 20 Longitude
    stations_code, stations_hex = get_stationcode_hexes_mapping()

    master_hex_ids = ch.get_all_grids_hex_ids(home_folder + 'MasterHexList.npy')
    # master_hex_ids = pd.read_csv(home_folder + 'MasterHexList_Res2.csv', header=None).values.tolist()
    # master_hex_ids = get_all_grids_hex_ids()
    # map station to master hex (some stations lie in the same hex- same connectivity)
    domain_mask = np.in1d(stations_hex, master_hex_ids)
    final_stations_code = stations_code[domain_mask]
    final_stations_hex = stations_hex[domain_mask]

    final_stations_master_indices = np.array([master_hex_ids.index(h) for h in final_stations_hex])
    # not use this as it sorts the results and return unique: np.in1d(master_hex_ids,stations_hex).nonzero()[0]

    [fullconnectivity(entry['Species'], entry['MinTemp'], entry['MaxTemp'], len(final_stations_hex),
                      final_stations_master_indices, final_stations_code) for index, entry in species_info.iterrows()]


if __name__ == '__main__':
    main()
