"""
Tool to obtain connectivty path and timescales between a pair of Tara stations
User must confirm the transition matrix file in use and that source and destination stations that exists in the list
Else you can give your stations locations too (needs some code edits to make this transition smooth)
"""

import numpy as np
import pandas as pd
from sourcecode.core import adjacencygraph as ag, connectivityplots as cp
from sourcecode.core import connectivityhelper as ch

hex_res = 3
data_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task11/'
depth = 50
Tara = False


def single_min_T_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                      d_index, destination_code):
    path_to_compute = 'Minimum Time'
    forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index, True)
    # print('f-path probable time:', ag.get_time_from_most_probable_path(atlantic_graph, forward_path)[-1])
    f_time_laps = np.arange(0, len(forward_path), 1) / 12
    print(f_time_laps[-1])
    backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index, True)
    # print('b-path probable time:', ag.get_time_from_most_probable_path(atlantic_graph, backward_path)[-1])
    b_time_laps = np.arange(0, len(backward_path), 1) / 12
    print(b_time_laps[-1])
    # cp.plot_paths(mask_lons, mask_lats, mask_value, master_grids_list, forward_path, backward_path, f_time_laps,
    #               b_time_laps, source_code, destination_code, path_to_compute)


def subset_min_T_paths(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                       d_index, destination_code, path_count, depth):
    forward_paths = ag.get_shortest_paths_subset(atlantic_graph, s_index, d_index, path_count)
    backward_paths = ag.get_shortest_paths_subset(atlantic_graph, d_index, s_index, path_count)
    cp.plot_shortest_paths_subset(mask_lons, mask_lats, mask_value, master_grids_list, forward_paths, backward_paths,
                                  source_code, destination_code, depth)


def main():
    source_code = 1
    destination_code = 9
    domain_adjacency_file = 't{0}m/Annual_Binary_DomainAdjacency_z{0}_csr.npz'.format(depth)

    temp_constraint_range = np.NaN
    min_accept_temp = np.nan
    max_accept_temp = np.nan
    width_type = 'passive'
    path_count = 1000
    master_grids_list = np.load(data_folder + 'H3_Res3_MasterHexList.npz')['Res3_HexId'].tolist()

    if Tara:
        s_hex, d_hex = ch.get_station_hexes_from_code(data_folder + 'AllStations_Tara.xls', hex_res, source_code,
                                                      destination_code)
    else:
        stations = pd.read_csv(data_folder + 'AtlanticStations.csv', header=0)
        src = stations.loc[stations['Station'] == source_code]
        des = stations.loc[stations['Station'] == destination_code]

        s_hex = ch.get_hexids(src['Latitude'], src['Longitude'], hex_res)
        d_hex = ch.get_hexids(des['Latitude'], des['Longitude'], hex_res)
    try:
        s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
    except KeyError:
        print("Source/destination sourcecode not present in the domain. Recheck values")
        raise

    if width_type == 'average':
        min_temp_file = data_folder + 't{0}m/Annual_avg_MinTemperature_z{0}_csr.npz'.format(depth)
        max_temp_file = data_folder + 't{0}m/Annual_avg_MaxTemperature_z{0}_csr.npz'.format(depth)
    elif width_type == 'broad':
        min_temp_file = data_folder + 't{0}m/Annual_min_MinTemperature_z{0}_csr.npz'.format(depth)
        max_temp_file = data_folder + 't{0}m/Annual_max_MaxTemperature_z{0}_csr.npz'.format(depth)
    elif width_type == 'passive':
        min_temp_file = None
        max_temp_file = None
    else:
        raise ValueError("width type is not recognized")
    print(width_type)
    # graph where the min-max 'temperature range' b/w grids is restricted
    if ~np.isnan(temp_constraint_range):
        atlantic_graph = ag.create_temp_range_graph(data_folder + domain_adjacency_file,
                                                    min_temp_file,
                                                    max_temp_file,
                                                    temp_constraint_range, None)
    # graph where connections with min and max temperature values are within the min-max values provided by the user
    # (eg. a thermal range for a species)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        atlantic_graph = ag.create_temp_min_max_graph(data_folder + domain_adjacency_file,
                                                      min_temp_file,
                                                      max_temp_file,
                                                      min_accept_temp, max_accept_temp, None)
    # simple graph, without any temperature boundaries
    else:
        atlantic_graph = ag.create_simple_graph(data_folder + domain_adjacency_file, None)
    print('Graph ready')

    mask_lons, mask_lats, mask_value = cp.load_mask_file(data_folder + 'GLOB16L98_mesh_mask_atlantic.nc')

    # single_min_T_path(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
    #                   d_index, destination_code)
    subset_min_T_paths(atlantic_graph, mask_lons, mask_lats, mask_value, master_grids_list, s_index, source_code,
                       d_index, destination_code, path_count, depth)


if __name__ == '__main__':
    main()
