import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

depth1 = 0
depth2 = 100
species = "G_bulloides"

base_path = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/2011_Lombard_Species/'
home_folder1 = base_path + 'depth{0}m/'.format(depth1)
home_folder2 = base_path + 'depth{0}m/'.format(depth2)

output_folder = base_path + 'output/'

# save format
# np.savez_compressed(home_folder + 'Stations_min-T_connectivity.npz', codes=final_stations_code, matrix=min_T_matrix)
data = np.load(home_folder1 + 'Stations_MinT_connectivity_{0}.npz'.format(species), allow_pickle=True)
codes = data['codes']
original_matrix = data['matrix']
print('maximum time: ', np.nanmax(original_matrix))

data_new = np.load(home_folder2 + 'Stations_MinT_connectivity_{0}.npz'.format(species),
                   allow_pickle=True)
# data_new = np.load(home_folder + 'Stations_min-T_connectivity_nan_TR2deg.npz', allow_pickle=True)

new_codes = data_new['codes']
new_matrix = data_new['matrix']
print('maximum time: ', np.nanmax(new_matrix))

# assure order is same
assert np.array_equal(codes, new_codes)

# compute fraction
fraction = (new_matrix - original_matrix) / original_matrix * 100
# fraction = new_matrix - original_matrix

avg = np.nanmean(fraction)
min_fraction, max_fraction = np.nanmin(fraction), np.nanmax(fraction)
print(min_fraction, max_fraction, avg)

# Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
# fig = plt.figure()
fig = plt.figure(figsize=(16, 14), dpi=100)
plt.margins(0, 0)
ax = plt.gca()
ax.set_title("Fractional change = (z{0} - z{1}) / z{1} * 100 | Average: {2}%".format(depth2, depth1, np.round(avg, 2)),
             pad=70)
plt.suptitle(
    "Fractional change in minimum connectivity time between all stations ({0}m vs {1}m depth)".format(depth1, depth2))
ax.set_xlabel("Destination")
ax.set_ylabel("Source", labelpad=20)
ax.set_xticks(np.arange(len(codes)))
ax.set_yticks(np.arange(len(codes)))
ax.set_xticklabels(np.arange(1, len(codes) + 1))
ax.set_yticklabels(np.arange(1, len(codes) + 1))
# ax.set_xticklabels(codes)
# ax.set_yticklabels(codes)
plt.tick_params(axis='both', which='major', labelsize=12)

# Y axis labels on top
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
plt.xticks(rotation=0)

# Turn spines off and create white grid.
ax.spines[:].set_visible(False)
ax.set_xticks(np.arange(len(codes) + 1) - .5, minor=True)
ax.set_yticks(np.arange(len(codes) + 1) - .5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
ax.tick_params(which="minor", bottom=False, left=False)

divnorm = colors.TwoSlopeNorm(vmin=min_fraction, vcenter=0., vmax=max_fraction)
plt.imshow(fraction, cmap=plt.cm.coolwarm, norm=divnorm)

cbar = plt.colorbar(orientation='vertical')

cbar.set_label('Fractional change (%)')
# plt.show()
plt.savefig(output_folder + "FractionalChange_z{0}m_vs_z{1}m_{2}.pdf".format(depth1, depth2,species), bbox_inches='tight',
            pad_inches=0.2)
