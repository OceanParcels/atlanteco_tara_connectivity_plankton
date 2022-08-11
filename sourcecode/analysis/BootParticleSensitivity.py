import pandas as pd
import matplotlib.pyplot as plt

# home_folder = '/nethome/manra003/analysis/paper01/'
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task11/BootEx/'
states_count = 100
sample_size = [5000, 10000, 50000, 100000, 200000, 300000]

source = 1
destination = 9
state = 100
max_particles = 375570

fp = pd.read_csv(home_folder + 'FullPaths.csv')
sd_pair = fp.loc[(fp['Source'] == source) & (fp['Destination'] == destination)]
# for each pair of station
# plot connectivity time sensitivity to number of particles
fig, ax = plt.subplots(2, 1, sharex=True)
# fig = plt.figure()
fig.suptitle("Minimum connectivity time sensitivity to number of particles")
fig.supxlabel('Bootstrap sample size')
fig.supylabel('Connectivity time in years')
# ax[0].set_xlim(0, max_particles)
# ax[1].set_xlim(0, max_particles)

# def plot_pair(index, color):
#     ds = data_frame[(data_frame['Source'] == src_stations[index]) & (data_frame['Destination'] == des_stations[index])]
#     ax[0, index].set(title='{0} to {1}'.format(src_stations[index], des_stations[index]))
#     ax[0, index].plot(ds['Sample_size'], ds['Fw_time'], linestyle='--', marker='o', color=color)
#
#     ax[1, index].set(title='{0} to {1}'.format(des_stations[index], src_stations[index]))
#     ax[1, index].plot(ds['Sample_size'], ds['Bw_time'], linestyle='--', marker='o', color=color)
#     ax[1, index].set_xticks(ds['Sample_size'][1:])
#     ax[1, index].set_xticklabels(ds['Sample_size'][1:], rotation=90, ha='center')
#
# #
# [plot_pair(i, c) for i, c in zip(range(len(src_stations)), ['r', 'b', 'orange', 'green'])]
# plt.show()


# pd.read_csv(home_folder + 'Boot_Sample/outputs/EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(s, d, size, states_count))
for size in sample_size:
    df = pd.read_csv(
        home_folder + 'EnsemblePaths_S{0}_D{1}_size{2}_en_{3}_z0.csv'.format(source,
                                                                             destination,
                                                                             size,
                                                                             states_count))
    ax[0].scatter(df['Sample_size'], df['F-minT'] / 12, s=5, c='r', alpha=0.5, marker='o')
    ax[1].scatter(df['Sample_size'], df['B-minT'] / 12, s=5, c='b', alpha=0.5, marker='o')

ax[0].scatter(max_particles, sd_pair['F-minT'], s=7, c='r', marker='^')
ax[1].scatter(max_particles, sd_pair['B-minT'], s=7, c='b', marker='^')
plt.show()
