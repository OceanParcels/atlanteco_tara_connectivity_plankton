import numpy as np
import matplotlib.pyplot as plt

depth = 0
species = "G_Sacculifer"
dataset = '2011_lombard_forams'
width_type = 'broad'
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task8E/Connectivities/{0}/t{1}m/'.format(dataset, depth)

data = np.load(home_folder + '{2}/Stations_minT_connectivity_{0}z_{1}_{2}.npz'.format(depth, species, width_type),
               allow_pickle=True)
codes = data['codes']
con_matrix = data['matrix']
print('maximum time: ', np.nanmax(con_matrix))
final_matrix = con_matrix / 12
print(round(np.nanmin(final_matrix), 2), round(np.nanmax(final_matrix), 2))

# Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
fig = plt.figure(figsize=(16, 14), dpi=200)
plt.margins(0, 0)
ax = plt.gca()
ax.set_title(
    "Minimum Connectivity Time ({4}) between all stations at depth {0}m- {1}\n minimum={2}, maximum={3} years".
        format(depth, species, round(np.nanmin(final_matrix), 2), round(np.nanmax(final_matrix), 2), width_type),
    pad=30)  # min/max Temperature: 15$^\circ$C to 32.54$^\circ$C Temperature adaptation rate: 2$^\circ$C per month

ax.set_xlabel("Destination stations", labelpad=30)
ax.set_ylabel("Source stations", labelpad=30)
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

plt.imshow(final_matrix, cmap='inferno_r')
cbar = plt.colorbar(orientation='vertical')

cbar.set_label('Minimum connectivity time (Years)', size=20)
cbar.ax.tick_params(labelsize=20)
# plt.show()
plt.savefig(home_folder + "{2}/Stations_minT_connectivity_z{0}_{1}_{2}.pdf".format(depth, species, width_type),
            bbox_inches='tight',
            pad_inches=0.2)
