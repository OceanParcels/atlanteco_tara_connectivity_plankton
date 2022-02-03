import h3
import pandas as pd
import numpy as np


def get_station_hexid(stations, code, hex_res):
    try:
        lat, lon = stations.loc[code]['Latitude'], stations.loc[code]['Longitude']
        return h3.geo_to_h3(lat, lon, hex_res)
    except KeyError:
        print("Incorrect source/destination code provided. Recheck values")
        raise


def get_station_hexes_from_code(file, hex_res, s_code, d_code):
    stations = pd.read_excel(file, header=1, index_col=0)
    return get_station_hexid(stations, s_code, hex_res), get_station_hexid(stations, d_code, hex_res)


def get_all_grids_hex_ids(file):
    master_uni_hex = np.load(file)  # for npy
    return master_uni_hex.tolist()
