import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sourcecode.core.connectivityplots as cp

# home_folder = '/nethome/manra003/analysis/paper01/'
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task12/'
states_count = 100
sample_size = np.array((5000, 10000, 50000, 100000, 200000, 300000))
custom_font_size = 15
custom_dpi = 300
source = 1
destination = 9
state = 100
max_particles = 375570

SRC_DES_FOUND = False
s_lon, s_lat = None, None
d_lon, d_lat = None, None

full_forward_time = 2.5
full_backward_time = 1.83
# for each pair of station
# plot connectivity time sensitivity to number of particles
fig, ax = plt.subplots(ncols=2, nrows=2, sharex=True,
                       dpi=custom_dpi, figsize=(16, 7))
fig.suptitle("Minimum connectivity time from 100 ensembles")
fig.supxlabel('Bootstrap sample size', fontsize=custom_font_size)

count_paths = np.zeros((len(sample_size)))

for index, size in enumerate(sample_size):
    df = pd.read_csv(
        home_folder + 'BootEx/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(source,
                                                                                    destination,
                                                                                    size,
                                                                                    states_count))
    # forward paths
    f_times = df['F-minT'].where(df['F-minT'] != -1) / 12
    ax[0, 0].scatter(np.full_like(f_times, index), f_times, c='r', alpha=0.5, marker='o')
    f_count = df['F-minT'].where(df['F-minT'] != -1).count()
    ax[1, 0].bar(index, f_count, width=0.2, color='r')

    # backward paths
    b_times = df['B-minT'].where(df['B-minT'] != -1) / 12
    ax[0, 1].scatter(np.full_like(b_times, index), b_times, c='b', alpha=0.5, marker='o')
    b_count = df['B-minT'].where(df['B-minT'] != -1).count()
    ax[1, 1].bar(index, b_count, width=0.2, color='b')

ax[0, 0].set_title('From Station {0} to {1}'.format(source, destination), fontsize=custom_font_size)
ax[0, 1].set_title('From Station {0} to {1}'.format(destination, source), fontsize=custom_font_size)

ax[1, 0].xaxis.set_ticks(range(len(sample_size)))
ax[1, 0].xaxis.set_ticklabels(sample_size, fontsize=custom_font_size)
ax[1, 0].set_ylim(0, 100)
y_labels = np.array((0, 20, 40, 60, 80, 100))
# ax[1, 0].yaxis.set_ticks(range(len(y_labels)))
ax[1, 0].yaxis.set_ticklabels(y_labels, fontsize=custom_font_size)
ax[1, 1].xaxis.set_ticks(range(len(sample_size)))
ax[1, 1].xaxis.set_ticklabels(sample_size, fontsize=custom_font_size)
ax[1, 1].set_ylim(0, 100)
# ax[1, 1].yaxis.set_ticks(range(len(y_labels)))
ax[1, 1].yaxis.set_ticklabels(y_labels, fontsize=custom_font_size)
ax[0, 0].axhline(full_forward_time, c='r', linestyle='--', alpha=0.6)
ax[0, 1].axhline(full_backward_time, c='b', linestyle='--', alpha=0.6)

t_labels = np.arange(1, 5)
ax[0, 0].set_yticks(t_labels)
ax[0, 0].yaxis.set_ticklabels(t_labels, fontsize=custom_font_size)

ax[0, 0].set_ylabel('Minimum connectivity time \n(years)', fontsize=custom_font_size)
ax[0, 1].set_yticks(t_labels)
ax[0, 1].yaxis.set_ticklabels(t_labels, fontsize=custom_font_size)
ax[1, 0].set_ylabel('Paths returned (%)', fontsize=custom_font_size)
plt.subplots_adjust(hspace=0.11, wspace=0.11)
plt.savefig(home_folder + 'BootEx/Particles_sens_Station1-930nov2022.pdf', bbox_inches='tight',
            pad_inches=0.5)
# plt.show()

# region: plot-connetivity paths

master_hex_ids = np.load(home_folder + 'H3_Res3_MasterHexList.npz')['Res3_HexId']


def plot_path(r_index, c_index, color, paths):
    ax[r_index, c_index].add_feature(cfeature.LAND, color='gainsboro')
    ax[r_index, c_index].set_extent([-100, 20, -78, 80])
    for p in paths:
        int_path = np.fromstring(p.replace('[', '').replace(']', ''), dtype=int, sep=',')
        if int_path.size != 0:
            lats, lons = cp.masterhex_to_latlon(master_hex_ids, int_path)
            ax[r_index, c_index].plot(lons, lats, color=color, linestyle='solid', alpha=0.3)
            ax[r_index, c_index].scatter(lons, lats, color=color, s=0.4)
            global SRC_DES_FOUND, s_lat, s_lon, d_lat, d_lon
            if not SRC_DES_FOUND:
                s_lon, s_lat = lons[0], lats[0]
                d_lon, d_lat = lons[-1], lats[-1]
                SRC_DES_FOUND = True


fig, ax = plt.subplots(figsize=(16, 7), ncols=len(sample_size), nrows=2,
                       sharex=True, subplot_kw={'projection': ccrs.PlateCarree()}, dpi=custom_dpi)

for index, size in enumerate(sample_size):
    df = pd.read_csv(
        home_folder + 'BootEx/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(source,
                                                                                    destination,
                                                                                    size,
                                                                                    states_count))
    plot_path(0, index, 'r', df['F-minP'].values)
    plot_path(1, index, 'b', df['B-minP'].values)
    ax[0, index].set_title('n = ' + str(size), fontsize=custom_font_size)
    ax[1, index].set_title('n = ' + str(size), fontsize=custom_font_size)
    ax[0, index].text(s_lon, s_lat, source,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'), fontsize=custom_font_size)
    ax[0, index].text(d_lon, d_lat, destination,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'), fontsize=custom_font_size)
    ax[1, index].text(s_lon, s_lat, source,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'), fontsize=custom_font_size)
    ax[1, index].text(d_lon, d_lat, destination,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'),
                      fontsize=custom_font_size)

ax[0, 0].set_ylabel('Minimum connectivity paths from station {0} to {1}'.format(source, destination))
ax[1, 0].set_ylabel('Minimum connectivity paths from station {0} to {1}'.format(destination, source))
plt.subplots_adjust(top=1, bottom=0.01, hspace=0.25, wspace=0.15)
plt.savefig(home_folder + 'BootEx/Min_paths_Station1-9_30nov2022.pdf', bbox_inches='tight',
            pad_inches=0.5)
# plt.show()
# endregion
