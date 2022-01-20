import numpy as np
import matplotlib.pyplot as plt

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task4/Full_connectivity_output/'

# np.savez_compressed(home_folder + 'Stations_min-T_connectivity.npz', codes=final_stations_code, matrix=min_T_matrix)
data = np.load(home_folder + 'Stations_min-T_connectivity_nan.npz', allow_pickle=True)
codes = data['codes']
original_matrix = data['matrix']
print('maximum time: ', np.nanmax(original_matrix))
genus = "Globoturborotalita"
species = 'rubescens'

data_new = np.load(home_folder + '2021_Fabio_Species/Stations_MinT_connectivity_{0}_{1}.npz'.format(genus, species),
                   allow_pickle=True)
# data_new = np.load(home_folder + 'Stations_min-T_connectivity_nan_TR2deg.npz', allow_pickle=True)

new_codes = data_new['codes']
new_matrix = data_new['matrix']
print('maximum time: ', np.nanmax(new_matrix))

# assure order is same
assert np.array_equal(codes, new_codes)

# compute fraction
fraction = (new_matrix - original_matrix) / original_matrix
avg = np.nanmean(fraction)
print(avg)

# Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
# fig = plt.figure()
fig = plt.figure(figsize=(16, 14), dpi=100)
plt.margins(0, 0)
ax = plt.gca()
ax.set_title("Average fractional change: {0}%".format(np.round(avg * 100, 2)), pad=30)
plt.suptitle("Fractional change in Connectivity Time between all stations (Passive vs {0}. {1})".format(genus, species))
ax.set_xlabel("Destination")
ax.set_ylabel("Source", labelpad=20)
ax.set_xticks(np.arange(len(codes)))
ax.set_yticks(np.arange(len(codes)))
ax.set_xticklabels(codes)
ax.set_yticklabels(codes)
plt.tick_params(axis='both', which='major', labelsize=7)

# Y axis labels on top
ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
plt.xticks(rotation=90)
# plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

# Turn spines off and create white grid.
ax.spines[:].set_visible(False)
ax.set_xticks(np.arange(len(codes) + 1) - .5, minor=True)
ax.set_yticks(np.arange(len(codes) + 1) - .5, minor=True)
ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
ax.tick_params(which="minor", bottom=False, left=False)

plt.imshow(fraction * 100, cmap=plt.cm.cool)
cbar = plt.colorbar(orientation='vertical')
# plt.clim(0, 3.5)
cbar.set_label('Fractional change (%)')
plt.show()
# plt.savefig(home_folder + "2021_Fabio_Species/Plots/Fraction_{0}_{1}.pdf".format(genus, species), bbox_inches='tight', pad_inches=0.2)
