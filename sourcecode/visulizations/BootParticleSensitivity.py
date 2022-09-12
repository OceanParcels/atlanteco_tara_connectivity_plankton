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
fig, ax = plt.subplots(ncols=1, nrows=2, sharex=True,
                       dpi=100, figsize=(10, 6))  # , subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle("Minimum connectivity time sensitivity to number of particles")
fig.supxlabel('Bootstrap sample size')
# fig.supylabel('Minimum connectivity time (years)')

count_paths = np.zeros((len(sample_size)))

# pd.read_csv(home_folder + 'Boot_Sample/outputs/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(s, d, size, states_count))
for index, size in enumerate(sample_size):
    df = pd.read_csv(
        home_folder + 'BootEx/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(source,
                                                                                    destination,
                                                                                    size,
                                                                                    states_count))
    time = df['F-minT'].where(df['F-minT'] != -1) / 12

    ax[0].scatter(np.full_like(time, index), time, c='r', alpha=0.5, marker='o')
    count = df['F-minT'].where(df['F-minT'] != -1).count()
    ax[1].bar(index, count, width=0.2, color='r')

ax[0].set(title='Station {0} to {1}'.format(source, destination))

ax[0].xaxis.set_ticks(range(len(sample_size)))
ax[0].xaxis.set_ticklabels(sample_size)
ax[1].xaxis.set_ticks(range(len(sample_size)))
ax[1].xaxis.set_ticklabels(sample_size)
# ax[1].scatter(df['Sample_size'], df['B-minT'].where(df['B-minT'] != -1) / 12, s=10, c='b', alpha=0.5, marker='o')
# ax[1].scatter(df['Sample_size'], df['B-minT'].where((df['B-minT'] == -1)) + 1, s=10, edgecolors='b',
#               alpha=0.5, marker='o', facecolors='none')
# ax[1].set(title='Station {0} to {1}'.format(destination, source))


# ax[0].axhline(full_forward_time, c='r', linestyle='--', alpha=0.6)
# ax[1].axhline(full_backward_time, c='b', linestyle='--', alpha=0.6)
ax[0].set_yticks(np.arange(2, 5))
# ax[1].set_yticks(np.arange(0, 5))
# ax[1].set_xticks(sample_size[1:])
# ax[1].set_xticklabels(sample_size[1:], rotation=90, ha='center')
# plt.savefig(home_folder + 'BootEx/Particles_sens_Station1-9.pdf', bbox_inches='tight',
#             pad_inches=0.5)
plt.show()
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
                       sharex=True, subplot_kw={'projection': ccrs.PlateCarree()}, dpi=100)

for index, size in enumerate(sample_size):
    df = pd.read_csv(
        home_folder + 'BootEx/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(source,
                                                                                    destination,
                                                                                    size,
                                                                                    states_count))
    plot_path(0, index, 'r', df['F-minP'].values)
    plot_path(1, index, 'b', df['B-minP'].values)
    ax[0, index].set_title('n = ' + str(size))
    ax[1, index].set_title('n = ' + str(size))
    ax[0, index].text(s_lon, s_lat, source,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'), fontsize=10)
    ax[0, index].text(d_lon, d_lat, destination,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'), fontsize=10)
    ax[1, index].text(s_lon, s_lat, source,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'), fontsize=10)
    ax[1, index].text(d_lon, d_lat, destination,
                      bbox=dict(facecolor='white', alpha=0.8, pad=0.2, edgecolor='none'),
                      fontsize=10)

# plt.subplots_adjust(wspace=0, hspace=0)
ax[0, 0].set_ylabel('Minimum connectivity paths from station {0} to {1}'.format(source, destination))
ax[1, 0].set_ylabel('Minimum connectivity paths from station {0} to {1}'.format(destination, source))
plt.subplots_adjust(top=1, bottom=0.01, hspace=0.25, wspace=0.15)
# plt.savefig(home_folder + 'BootEx/Min_paths_Station1-9.pdf', bbox_inches='tight',
#             pad_inches=0.5)
plt.show()
# endregion
