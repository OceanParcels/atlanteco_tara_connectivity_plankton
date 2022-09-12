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


home_folder = '/nethome/manra003/analysis/paper01/'
hex_res = 3

master_grids_list = np.load(home_folder + 'H3_Res3_MasterHexList.npz')['Res3_HexId'].tolist()

stations = pd.read_csv(home_folder + 'AtlanticStations.csv', header=0)
lon = stations['Longitude']
lat = stations['Latitude']

s = 1
d = 9

source = stations.loc[stations['Station']==s]
destination = stations.loc[stations['Station']==d]

s_hex, d_hex = ch.get_hexids(source['Latitude'].values[0], source['Longitude'].values[0], hex_res), ch.get_hexids(destination['Latitude'].values[0], destination['Longitude'].values[0], hex_res)

s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
print(s_index,d_index)

class Ensemble:

    def __init__(self, size, state, fT, fP, bT, bP):
        self.sample_size = size
        self.ensemble_state = state
        self.f_min_time = fT
        self.f_min_path = fP
        self.b_min_time = bT
        self.b_min_path = bP

    def to_dict(self):
        return {
            'Sample_size': self.sample_size,
            'State': self.ensemble_state,
            'F-minT': self.f_min_time,
            'F-minP': self.f_min_path,
            'B-minT': self.b_min_time,
            'B-minP': self.b_min_path,

        }

states_count = 100
sample_size = [5000, 10000, 50000, 100000, 200000, 300000]


for size in sample_size:
    ensemble_list = list()
    for state in range(1, states_count + 1):

        atlantic_graph = ag.create_simple_graph(
            home_folder + 'Boot_Sample/Size_{0}/Annual/State_{1}/Annual_Binary_DomainAdjacency_z0_csr.npz'.format(size,
                                                                                                                state),
            None)
            
        forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
        if forward_path:
            f_time = len(forward_path) - 1
        else:
            f_time = -1
        backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index)
        if backward_path:
            b_time = len(backward_path) - 1
        else:
            b_time = -1

        ensemble_list.append(Ensemble(size, state, f_time, forward_path, b_time, backward_path))
    # export to dataframe and save file
    pd.DataFrame.from_records([e.to_dict() for e in ensemble_list]).to_csv(
        home_folder + 'Boot_Sample/outputs/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(s, d, size, states_count))