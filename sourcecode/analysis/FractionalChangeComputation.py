import numpy as np
import pandas as pd

dataset = '2011_Lombard_Species'
depth = 0
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/{0}/depth{1}m/'.format(dataset, depth)

# np.savez_compressed(home_folder + 'Stations_min-T_connectivity.npz', codes=final_stations_code, matrix=min_T_matrix)
data = np.load(
    '/Users/dmanral/Desktop/Analysis/TARA/Task4/Full_connectivity_output/' + 'Stations_min-T_connectivity_nan.npz',
    allow_pickle=True)
codes = data['codes']
original_matrix = data['matrix']
print('maximum time: ', np.nanmax(original_matrix))



species_info = pd.read_csv('/Users/dmanral/Desktop/Analysis/TARA/Task7D/{0}/'.format(dataset) + '2011_lombard_forams.csv',
                           delimiter=';|,', header=0)
max_time = np.zeros(len(species_info))
mean_fraction = np.empty(len(species_info))
max_fraction = np.empty(len(species_info))

for index, entry in species_info.iterrows():
    data_new = np.load(home_folder + 'Stations_MinT_connectivity_{0}.npz'.format(entry['Species']),
                       allow_pickle=True)

    new_codes = data_new['codes']
    new_matrix = data_new['matrix']
    max_time[index] = np.nanmax(new_matrix)

    # assure order is same
    assert np.array_equal(codes, new_codes)

    # compute fraction
    fraction = (new_matrix - original_matrix) / original_matrix
    mean_fraction[index] = np.round(np.nanmean(fraction) * 100, 2)
    max_fraction[index] = np.round(np.nanmax(fraction) * 100, 2)
    print(max_time[index], mean_fraction[index], max_fraction[index])

species_info['MaxTimeMonths'] = max_time
species_info['MeanFractionalChange'] = mean_fraction
species_info['MaxFractionalChange'] = max_fraction

species_info.to_csv(home_folder + 'FractionalChange.csv')
