import numpy as np
import pandas as pd

dataset = 'sample_constraints'  # 'sample_constraints'# '2011_lombard_forams'
depth = 50
work_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task11/'
home_folder = work_folder + 'Connectivities/{0}/t{1}m/'.format(dataset, depth)
width_type = 'broad'

# np.savez_compressed(home_folder + 'Stations_min-T_connectivity.npz', codes=final_stations_code, matrix=min_T_matrix)
data = np.load(work_folder + 'Connectivities/Stations_minT_connectivity_0z_NoConstraints_passive.npz',
               allow_pickle=True)
codes = data['codes']
original_matrix = data['matrix']
print('maximum time: ', np.nanmax(original_matrix))

species_info = pd.read_csv(work_folder + dataset + '.csv',
                           delimiter=';|,', header=0)
max_time = np.zeros(len(species_info))
mean_fraction = np.empty(len(species_info))
max_fraction = np.empty(len(species_info))
min_fraction = np.empty(len(species_info))

for index, entry in species_info.iterrows():
    data_new = np.load(
        home_folder + '{2}/Stations_minT_connectivity_{0}z_{1}_{2}.npz'.format(depth, entry['Species'], width_type),
        allow_pickle=True)

    new_codes = data_new['codes']
    new_matrix = data_new['matrix']
    max_time[index] = np.nanmax(new_matrix)

    # assure order is same
    assert np.array_equal(codes, new_codes)

    # compute fraction
    fraction = (new_matrix - original_matrix) / original_matrix
    mean_fraction[index] = np.round(np.nanmean(fraction) * 100, 2)
    min_fraction[index] = np.round(np.nanmin(fraction) * 100, 2)
    max_fraction[index] = np.round(np.nanmax(fraction) * 100, 2)
    print(max_time[index], mean_fraction[index],min_fraction[index], max_fraction[index])

species_info['MaxTimeMonths'] = max_time
species_info['MeanFractionalChange'] = mean_fraction
species_info['MinFractionalChange'] = min_fraction
species_info['MaxFractionalChange'] = max_fraction

species_info.to_csv(home_folder + 'FractionalChange_{0}_{1}z.csv'.format(width_type, depth))
