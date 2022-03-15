"""
Tool to obtain connectivty path and timescales between a pair of Tara stations
User must confirm the transition matrix file in use and that source and destination stations that exists in the list
Else you can give your stations locations too (needs some code edits to make this transition smooth)
"""

import numpy as np

from sourcecode.core import adjacencygraph as ag
from sourcecode.core import connectivityhelper as ch
from sourcecode.visulizations import connectivityplots as cp

hex_res = 3
data_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'
depth = 100


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
    forward_paths = ag.get_shortest_paths_subset(atlantic_graph, s_index, d_index)
    backward_paths = ag.get_shortest_paths_subset(atlantic_graph, d_index, s_index)
    cp.plot_shortest_paths_subset(mask_lons, mask_lats, mask_value, master_grids_list, forward_paths, backward_paths,
                                  source_code, destination_code)


def main():
    source_code = '66SUR'
    destination_code = 'OA002'
    domain_adjacency_file = 't{0}m/Annual_Avg_DomainAdjacency_csr.npz'.format(depth)

    temp_constraint_range = np.NaN
    min_accept_temp = 11.85
    max_accept_temp = 28.85
    minimum_time_omalley2021 = False

    if minimum_time_omalley2021:
        t_ratio_file = data_folder + 't0m/MinTRatio_FullAdjacency_csr.npz'
    else:
        t_ratio_file = None

    master_grids_list = ch.get_all_grids_hex_ids(data_folder + 'MasterHexList_Res3.npy')

    s_hex, d_hex = ch.get_station_hexes_from_code(data_folder + 'AllStations_Tara.xls', hex_res, source_code,
                                                  destination_code)
    try:
        s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
    except KeyError:
        print("Source/destination sourcecode not present in the domain. Recheck values")
        raise

    # graph where the min-max 'temperature range' b/w grids is restricted
    if ~np.isnan(temp_constraint_range):
        atlantic_graph = ag.create_temp_range_graph(data_folder + domain_adjacency_file,
                                                    data_folder + 't{0}m/Annual_Avg_MinTemperature_csr.npz'.format(
                                                        depth),
                                                    data_folder + 't{0}m/Annual_Avg_MaxTemperature_csr.npz'.format(
                                                        depth),
                                                    temp_constraint_range, t_ratio_file)
    # graph where connections with min and max temperature values are within the min-max values provided by the user
    # (eg. a thermal range for a species)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        atlantic_graph = ag.create_temp_min_max_graph(data_folder + domain_adjacency_file,
                                                      data_folder + 't{0}m/Annual_Avg_MinTemperature_csr.npz'.format(
                                                          depth),
                                                      data_folder + 't{0}m/Annual_Avg_MaxTemperature_csr.npz'.format(
                                                          depth),
                                                      min_accept_temp, max_accept_temp, t_ratio_file)
    # simple graph, without any temperature boundaries
    else:
        atlantic_graph = ag.create_simple_graph(data_folder + domain_adjacency_file, t_ratio_file)
    print('Graph ready')

    mask_lons, mask_lats, mask_value = cp.load_mask_file(data_folder + 'GLOB16L98_mesh_mask_atlantic.nc')

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
        single_min_T_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                          d_index, destination_code)
        subset_min_T_paths(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                           d_index, destination_code)


if __name__ == '__main__':
    main()
