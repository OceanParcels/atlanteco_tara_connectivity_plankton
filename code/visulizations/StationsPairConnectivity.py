import numpy as np

from code.core import adjacencygraph as ag
from code.core import connectivityhelper as ch
from code.core import connectivityplots as cp

hex_res = 3
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'


# home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/'
# data_folder = '/Users/dmanral/Desktop/Analysis/UvA/'


# home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task5_trahms/'


# def get_all_grids_hex_ids():
# hexId and matrix index mapping from a random simulation output file from the output dataset
# with xr.open_dataset(np.random.choice(glob(home_folder + '/tara_res5_01/FullTara_Res5_TS_*'))).load() as temp_ds:
#     # with xr.open_dataset(np.random.choice(glob(home_folder + 'tara_data/FullAtlantic_2D_01*'))).load() as temp_ds:
#     master_all_hex_t0 = np.array(
#         [h3.geo_to_h3(y, x, hex_res) for x, y in zip(temp_ds['lon'][:, 0].values, temp_ds['lat'][:, 0].values)])
#     master_uni_hex, counts = np.unique(master_all_hex_t0, return_counts=True)
#     return master_uni_hex.tolist()
# master_uni_hex = np.load(home_folder + 'MasterHexList.npy')
# return master_uni_hex.tolist()


def most_likely_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                     d_index, destination_code):
    path_to_compute = 'Most_Likely'
    forward_path = ag.get_most_probable_path(atlantic_graph, s_index, d_index)
    if forward_path:
        print('f_path=', len(forward_path) - 1)
        f_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, forward_path)
        cp.plot_time_vs_transitions(master_grids_list, forward_path, f_time_laps, source_code, destination_code)
    else:
        print("Path not found from {0} to {1}".format(source_code, destination_code))
    backward_path = ag.get_most_probable_path(atlantic_graph, d_index, s_index)
    if backward_path:
        print('b_path=', len(backward_path) - 1)
        b_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, backward_path)
        cp.plot_time_vs_transitions(master_grids_list, backward_path, b_time_laps, destination_code, source_code)
    else:
        print("Path not found from {0} to {1}".format(destination_code, source_code))
    if forward_path and backward_path:
        cp.plot_paths(mask_lons, mask_lats, mask_value, master_grids_list, forward_path, backward_path, f_time_laps,
                      b_time_laps, source_code, destination_code, path_to_compute)


def single_min_T_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                      d_index, destination_code):
    path_to_compute = 'Minimum Time'
    forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
    # print('f-path probable time:', ag.get_time_from_most_probable_path(atlantic_graph, forward_path)[-1])
    f_time_laps = np.arange(0, len(forward_path), 1) / 12
    print(f_time_laps[-1])
    backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index)
    # print('b-path probable time:', ag.get_time_from_most_probable_path(atlantic_graph, backward_path)[-1])
    b_time_laps = np.arange(0, len(backward_path), 1) / 12
    print(b_time_laps[-1])
    cp.plot_paths(mask_lons, mask_lats, mask_value, master_grids_list, forward_path, backward_path, f_time_laps,
                  b_time_laps, source_code, destination_code, path_to_compute)


def subset_min_T_paths(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                       d_index, destination_code):
    path_to_compute = 'Shortest_Paths'
    forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
    forward_paths = ag.get_shortest_paths_subset(atlantic_graph, s_index, d_index, len(forward_path))
    # forward_paths = ag.get_all_shortest_paths(atlantic_graph, s_index, d_index, len(forward_path))
    # for p in forward_paths:
    #     ag.get_time_from_most_probable_path(atlantic_graph, p)
    backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index)
    backward_paths = ag.get_shortest_paths_subset(atlantic_graph, d_index, s_index, len(backward_path))
    # backward_paths = ag.get_all_shortest_paths(atlantic_graph, d_index, s_index, len(backward_path))
    cp.plot_shortest_paths_subset(mask_lons, mask_lats, mask_value, master_grids_list, forward_paths, backward_paths,
                                  source_code, destination_code)


def main():
    # parameters from arguments?? pairs(66SUR,OA002)(84SUR,210SUR)
    # pairs(66SUR,OA002)(84SUR,210SUR)(78SUR,150SUR)
    source_code = '66SUR'
    destination_code = 'OA002'
    # domain_adjacency_file = 't0m/Annual_Avg_Prob_Cutoff_pt001_csr.npz'
    domain_adjacency_file = 't0m/Annual_Avg_DomainAdjacency_csr.npz'
    # domain_adjacency_file = 'Boot_377583/Annual/State9/Annual_Avg_DomainAdjacency_csr.npz'

    temp_constraint_range = np.NaN
    min_accept_temp = np.NaN
    max_accept_temp = np.NaN
    prob_cutoff = np.NAN
    trahms_2021 = False
    minimum_time_omalley2021 = False

    if minimum_time_omalley2021:
        t_ratio_file = home_folder + 't0m/MinTRatio_FullAdjacency_csr.npz'
    else:
        t_ratio_file = None

    master_grids_list = ch.get_all_grids_hex_ids(home_folder + 'MasterHexList_Res3.npy')

    # def get_stationcode_hexes_mapping():
    #     stations = pd.read_csv(data_folder + 'Stations.csv', header=0, index_col=1)
    #     # sort stations in order of their latitudes
    #     sorted_stations = stations.sort_values(by=['Latitude'], ascending=False)
    #     return sorted_stations.loc[source_code].H3Id, sorted_stations.loc[destination_code].H3Id
    #
    # s_hex, d_hex = get_stationcode_hexes_mapping()

    s_hex, d_hex = ch.get_station_hexes_from_code(home_folder + 'AllStations_Tara.xls', hex_res, source_code,
                                                  destination_code)
    # source_code = 'Gulf Of Mexico'
    # destination_code = 'Dutch Coast'
    # s_hex = h3.geo_to_h3(25.254817, -90.689833, hex_res) #middle of gom
    # d_hex = h3.geo_to_h3(52.035329, 4.168174, hex_res)
    try:
        s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
    except KeyError:
        print("Source/destination code not present in the domain. Recheck values")
        raise

    # graph where the min-max 'temperature range' b/w grids is restricted
    if ~np.isnan(temp_constraint_range):
        atlantic_graph = ag.create_temp_range_graph(home_folder + domain_adjacency_file,
                                                    home_folder + 't0m/Annual_Avg_MinTemperature_csr.npz',
                                                    home_folder + 't0m/Annual_Avg_MaxTemperature_csr.npz',
                                                    temp_constraint_range, t_ratio_file)
    # graph where connections with min and max temperature values are within the min-max values provided by the user
    # (eg. a thermal range for a species)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        atlantic_graph = ag.create_temp_min_max_graph(home_folder + domain_adjacency_file,
                                                      home_folder + 't0m/Annual_Avg_MinTemperature_csr.npz',
                                                      home_folder + 't0m/Annual_Avg_MaxTemperature_csr.npz',
                                                      min_accept_temp, max_accept_temp, t_ratio_file)
    elif ~np.isnan(prob_cutoff):
        atlantic_graph = ag.create_prob_filtered_graph(home_folder + domain_adjacency_file,
                                                       home_folder + 't0m/Annual_Avg_MinTemperature_csr.npz',
                                                       home_folder + 't0m/Annual_Avg_MaxTemperature_csr.npz',
                                                       prob_cutoff)
    elif trahms_2021:
        atlantic_graph = ag.create_full_graph(home_folder + domain_adjacency_file,
                                              home_folder + 't0m/Annual_SUM_MinTemperature_csr.npz',
                                              home_folder + 't0m/Annual_SUM_MaxTemperature_csr.npz',
                                              t_ratio_file)
    # elif minimumTime_malley2021:
    #     atlantic_graph = ag.create_minimum_time_graph(home_folder + domain_adjacency_file,
    #                                                   home_folder + 't0m/MinTRatio_FullAdjacency_csr.npz')
    # simple graph, without any temperature boundaries
    else:
        # atlantic_graph = ag.create_full_graph(home_folder + domain_adjacency_file,
        #                                       home_folder + 't0m/Annual_Avg_MinTemperature_csr.npz',
        #                                       home_folder + 't0m/Annual_Avg_MaxTemperature_csr.npz',
        #                                       t_ratio_file)
        atlantic_graph = ag.create_simple_graph(home_folder + domain_adjacency_file, t_ratio_file)
    print('Graph ready')

    mask_lons, mask_lats, mask_value = cp.load_mask_file(home_folder + 'GLOB16L98_mesh_mask_atlantic.nc')

    if minimum_time_omalley2021:
        forward_path = ag.minimum_time_path_malley2021(atlantic_graph, s_index, d_index)
        f_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, forward_path)
        cp.plot_time_vs_transitions(master_grids_list, forward_path, f_time_laps, source_code, destination_code)
        backward_path = ag.minimum_time_path_malley2021(atlantic_graph, d_index, s_index)
        b_time_laps = ag.get_time_from_most_probable_path(atlantic_graph, backward_path)
        cp.plot_time_vs_transitions(master_grids_list, backward_path, b_time_laps, destination_code, source_code)
        cp.plot_paths(mask_lons, mask_lats, mask_value, master_grids_list, forward_path, backward_path, f_time_laps,
                      b_time_laps, source_code, destination_code, 'Minimum Travel Time')
    else:
        # most_likely_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
        #                  d_index, destination_code)
        single_min_T_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                          d_index, destination_code)
        subset_min_T_paths(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                           d_index, destination_code)


if __name__ == '__main__':
    main()

