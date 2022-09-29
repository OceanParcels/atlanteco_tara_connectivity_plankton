import pandas as pd
import numpy as np
import h3

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task12/'

stations = pd.read_csv(home_folder + 'AllStations_Tara.csv', header=0)

atlantic_lon = stations['Longitude'].values
atlantic_lat = stations['Latitude'].values

stations_hex_list = [h3.geo_to_h3(y, x, 3) for x, y in zip(atlantic_lon.ravel(), atlantic_lat.ravel())]
stations_hex_res2 = [h3.geo_to_h3(y, x, 2) for x, y in zip(atlantic_lon.ravel(), atlantic_lat.ravel())]
stations_hex_res4 = [h3.geo_to_h3(y, x, 4) for x, y in zip(atlantic_lon.ravel(), atlantic_lat.ravel())]
print(len(stations_hex_list))

full_hex = np.load(home_folder + 'H3_Res3_MasterHexList.npz')['Res3_HexId'].tolist()
missing_hex_mask = np.in1d(stations_hex_list, full_hex)

output_df = pd.DataFrame(data={'Code': stations['Station'][missing_hex_mask].values,
                               'Latitude': stations['Latitude'][missing_hex_mask].values,
                               'Longitude': stations['Longitude'][missing_hex_mask].values,
                               'Res3HexId': pd.Series(stations_hex_list)[missing_hex_mask].values,
                               'Res2HexId': pd.Series(stations_hex_res2)[missing_hex_mask].values,
                               'Res4HexId': pd.Series(stations_hex_res4)[missing_hex_mask].values
                               })
output_df.sort_values('Latitude', ascending=False).to_csv(home_folder + 'TaraStationsHexIdMapping.csv', index=0)
