"""
Program to test a single simulation output
Set month, year and depth to select the file (of a given name format)
1. check particles that remained stuck from time 0 of the simulations
2. check particles that get stuck during a simulation
3. remove these particles and plot Min/Max of Temperature and Salinity
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'


def main():
    mon = 'Jan'
    year = 2012
    depth = 0

    file = home_folder + 'FullTara_Res5_TS_1{0}{1}_dt600_z{2}.nc'.format(mon, year, depth)

    ds = xr.open_dataset(file)

    t_start, t_final = (ds['time'][0, 0].values, ds['time'][0, 1].values)
    print(t_start, t_final)

    lons = ds['lon'][:, 0].values
    lats = ds['lat'][:, 0].values
    flons = ds['lon'][:, -1].values
    flats = ds['lat'][:, -1].values

    # static points analysis- points that did not change their position at all
    static_pts_index = \
        np.where(np.logical_and((ds['lat'][:, 0] == ds['lat'][:, -1]), (ds['lon'][:, 0] == ds['lon'][:, -1])))[0]
    print("Static Points count: ", len(static_pts_index))

    # points deleted
    # points with 0 fields of temp or salinity
    maxsal_zero_index = np.where(ds['max_sal'][:, -1] == 0)[0]
    print("Zero MAX Salinity count: ", len(maxsal_zero_index))
    minsal_zero_index = np.where(ds['min_sal'][:, -1] == 0)[0]
    print("Zero MIN Salinity count: ", len(minsal_zero_index))
    # only possible scenario is when a particle will get a lower salinity-when stuck

    maxtemp_zero_index = np.where(ds['max_temp'][:, -1] == 0)[0]
    print("Zero MAX Temperature count: ", len(maxtemp_zero_index))
    mintemp_zero_index = np.where(ds['min_temp'][:, -1] == 0)[0]
    print("Zero MIN Temperature count: ", len(mintemp_zero_index))
    # here both scenarios are possible, particle with >0 min_temperature and <0 max_temp
    # it is possible to find unique scenarios when maz_temp was below 0 and hex
    # when it gets stuck the max_temp gets updated to zero

    assert np.array_equal(static_pts_index, maxsal_zero_index)

    min_sal_only = np.setdiff1d(minsal_zero_index, static_pts_index)
    min_temp_only = np.setdiff1d(mintemp_zero_index, static_pts_index)
    max_temp_only = np.setdiff1d(maxtemp_zero_index, static_pts_index)

    join_arrays = np.concatenate((min_sal_only, min_temp_only, max_temp_only, static_pts_index))
    delete_indexes = np.unique(join_arrays)

    model_mask_file = home_folder + 'GLOB16L98_mesh_mask_atlantic.nc'
    mask_ds = xr.open_dataset(model_mask_file, decode_times=False).load()
    # get the corner points to plot on the map
    x = mask_ds['glamf']
    y = mask_ds['gphif']
    # get the mask values of the corner points
    c = mask_ds['tmask'][:]

    def plot_points(field_indices, field_name):
        fig = plt.figure()
        ax = plt.axes()
        colormap = clr.ListedColormap(['gainsboro', 'white'])
        plt.suptitle(
            "Highlighting land points(black) and trajectories that get stuck at land with {0} as 0.0".format(
                field_name))
        ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

        plt.scatter(lons[static_pts_index], lats[static_pts_index], s=1, c='black', label='Invalid depth')
        plt.scatter(lons[field_indices], lats[field_indices], s=10, c='orange', alpha=0.6, label='Initial location')
        plt.scatter(flons[field_indices], flats[field_indices], s=25, c='seagreen', alpha=0.4, label='Final location')
        plt.legend()
        plt.show()

    # plot all the particles locations -to be deleted
    plot_points(min_sal_only, 'minimum salinity')
    plot_points(min_temp_only, 'minimum temperature')
    plot_points(max_temp_only, 'maximum temperature')

    def plot_field(field, data, colorm):
        fig = plt.figure()
        ax = plt.axes()
        colormap = clr.ListedColormap(['gainsboro', 'white'])

        ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)

        plt.scatter(np.delete(lons, delete_indexes), np.delete(lats, delete_indexes), c=data, cmap=colorm, s=0.5)
        # plt.scatter(lons, lats, c=data, cmap=colorm, norm=clr.LogNorm(), s=0.5)

        plt.title('{} for {} {} \n min={:.4f}, max={:.4f}'.format(field, mon, year, np.min(data), np.max(data)),
                  fontsize=7)
        cbar = plt.colorbar()
        cbar.set_label(field)
        plt.show()

    plot_field('Minimum Temperature', np.delete(ds['min_temp'][:, 1].values, delete_indexes), plt.cm.Blues.reversed())
    plot_field('Maximum Temperature', np.delete(ds['max_temp'][:, 1].values, delete_indexes), plt.cm.Blues.reversed())

    plot_field('Minimum Salinity', np.delete(ds['min_sal'][:, 1].values, delete_indexes), plt.cm.Greens.reversed())
    plot_field('Maximum Salinity', np.delete(ds['max_sal'][:, 1].values, delete_indexes), plt.cm.Greens.reversed())


if __name__ == '__main__':
    main()
