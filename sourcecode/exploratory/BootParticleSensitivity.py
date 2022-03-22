"""
Program to test how sensitive the connectivity time is to the number of released particles:
Primary test based on surface connectivity without constraints
Input- adjacency files from bootstrap runs
Bootstrap runs without replacement (less than maximum number of particles).
"""
import pandas as pd
import numpy as np
from sourcecode.core import adjacencygraph as ag
from sourcecode.core import connectivityhelper as ch
import matplotlib.pyplot as plt

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'
hex_res = 3
sample_size = [5000, 10000, 50000, 100000, 200000, 300000, 377583]

master_grids_list = ch.get_all_grids_hex_ids(home_folder + 'MasterHexList_Res3.npy')
src_stations = ['80SUR', 'OA024', '71SUR', '84SUR']
des_stations = ['68SUR', 'OA245', 'OA022', '210SUR']
stations_pair = pd.read_excel(home_folder + 'AllStations_Tara.xls', header=1, index_col=0)

s_index = np.zeros((len(src_stations)), dtype='int')
d_index = np.zeros((len(src_stations)), dtype='int')

for s, d, i in zip(src_stations, des_stations, range(len(src_stations))):
    s_hex, d_hex = ch.get_station_hexid(stations_pair, s, hex_res), ch.get_station_hexid(stations_pair, d, hex_res)
    s_index[i], d_index[i] = master_grids_list.index(s_hex), master_grids_list.index(d_hex)

size_list = list()
f_time_list = list()
b_time_list = list()
src_list = list()
des_list = list()
# load graph for a given boot sample size
for size in sample_size:
    atlantic_graph = ag.create_simple_graph(home_folder + 'boot/Annual_Avg_DomainAdjacency_csr_{0}.npz'.format(size),
                                            None)

    # get connectivity between set of stations
    for s_id, d_id, i in zip(s_index, d_index, range(len(src_stations))):
        src_list.append(src_stations[i])
        des_list.append(des_stations[i])
        size_list.append(size)
        forward_path = ag.get_shortest_path(atlantic_graph, s_id, d_id)
        if forward_path:
            f_time = len(forward_path) - 1
            f_time_list.append(f_time)
        else:
            f_time_list.append(-1)
        backward_path = ag.get_shortest_path(atlantic_graph, d_id, s_id)
        if backward_path:
            b_time = len(backward_path) - 1
            b_time_list.append(b_time)
        else:
            b_time_list.append(-1)

# create dataframe to store all the information
data_frame = pd.DataFrame({'Source': src_list,
                           'Destination': des_list,
                           'Sample_size': size_list,
                           'Fw_time': f_time_list,
                           'Bw_time': b_time_list})

f"processing completed. time to prepare plots"
# for each pair of station
# plot connectivity time sensitivity to number of particles
fig, ax = plt.subplots(2, len(src_stations), sharex=True)
fig.suptitle("Minimum connectivity time sensitivity to number of particles")
fig.supxlabel('Bootstrap sample size')
fig.supylabel('Connectivity time in months')


def plot_pair(index, color):
    ds = data_frame[(data_frame['Source'] == src_stations[index]) & (data_frame['Destination'] == des_stations[index])]
    ax[0, index].set(title='{0} to {1}'.format(src_stations[index], des_stations[index]))
    ax[0, index].plot(ds['Sample_size'], ds['Fw_time'], linestyle='--', marker='o', color=color)

    ax[1, index].set(title='{0} to {1}'.format(des_stations[index], src_stations[index]))
    ax[1, index].plot(ds['Sample_size'], ds['Bw_time'], linestyle='--', marker='o', color=color)
    ax[1, index].set_xticks(ds['Sample_size'][1:])
    ax[1, index].set_xticklabels(ds['Sample_size'][1:], rotation=90, ha='center')


[plot_pair(i, c) for i, c in zip(range(len(src_stations)), ['r', 'b', 'orange', 'green'])]
plt.show()

# sampled points may not be necessarily equally distributed- select points per grid cell? too much?

# these samples are obtained from a signle random state, so maybe choose 300,000 particles set and create multiple samples
