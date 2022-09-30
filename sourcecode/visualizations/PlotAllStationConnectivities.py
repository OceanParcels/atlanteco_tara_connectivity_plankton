import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

depth = 0

dataset = 'NoConstraints'  # 'sample_constraints' 2011_lombard_forams NoConstraints
width_type = 'passive'
work_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task12/'
home_folder = work_folder + 'Connectivities/{0}/t{1}m/'.format(dataset, depth)


def plot_connectivity(species):
    data = np.load(home_folder + '{2}/Stations_minT_connectivity_{0}z_{1}_{2}.npz'.format(depth, species, width_type),
                   allow_pickle=True)
    codes = data['codes']
    con_matrix = data['matrix']
    print('maximum time: ', np.nanmax(con_matrix))
    final_matrix = con_matrix / 12
    print(round(np.nanmin(final_matrix), 2), round(np.nanmax(final_matrix), 2))

    # Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    fig = plt.figure(figsize=(16, 14), dpi=300)
    plt.margins(0, 0)
    ax = plt.gca()
    ax.set_title(
        "Minimum Connectivity Time ({4}) between all stations at depth {0}m- {1}\n minimum={2}, maximum={3} years".
            format(depth, species, round(np.nanmin(final_matrix), 2), round(np.nanmax(final_matrix), 2), width_type),
        pad=70, fontsize=20)
    ax.set_xlabel("Destination stations", labelpad=30, fontsize=20)
    ax.set_ylabel("Source stations", labelpad=30, fontsize=20)
    ax.set_xticks(np.arange(len(codes)))
    ax.set_yticks(np.arange(len(codes)))
    ax.set_xticklabels(codes)
    ax.set_yticklabels(codes)
    plt.tick_params(axis='both', which='major', labelsize=20)

    # Y axis labels on top
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    # plt.xticks(rotation=90)
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
    plt.savefig(home_folder + "Stations_minT_connectivity_z{0}_{1}_{2}.png".format(depth, species, width_type),
                bbox_inches='tight',
                pad_inches=0.2)


if dataset == 'NoConstraints':
    plot_connectivity('NoConstraints')
    exit(0)

species_info = pd.read_csv(work_folder + dataset + '.csv',
                           delimiter=';|,', keep_default_na=True, header=0, engine='python')
for species in species_info.Species:
    plot_connectivity(species)
