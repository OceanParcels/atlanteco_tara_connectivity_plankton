import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd

depth1 = 0
depth2 = 0
species = 'AP_Type1'
width_type1 = 'average'
width_type2 = 'broad'
home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task9B/Connectivities/'
dataset = 'sample_constraints'  # sample_constraints 2011_lombard_forams


# depth2 = 0
# home_folder1 = base_path + 'depth{0}m/'.format(depth1)
# home_folder2 = base_path + 'depth{0}m/'.format(depth2)


def compute_fractional_change(base_matrix, matrix2):
    # compute fraction
    fraction = (matrix2 - base_matrix) / base_matrix * 100
    # fraction = new_matrix - original_matrix
    return fraction


def plot_change(fraction, species, codes, depth1, depth2, width_type1, width_type2):
    # Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    # fig = plt.figure()

    avg = np.nanmean(fraction)
    min_fraction, max_fraction = np.nanmin(fraction), np.nanmax(fraction)
    print(min_fraction, max_fraction, avg)
    fig = plt.figure(figsize=(16, 14), dpi=100)
    plt.margins(0, 0)
    ax = plt.gca()
    ax.set_title("minimum: {0}%, maximum: {1}%, average: {2}%".format(round(min_fraction, 2), round(max_fraction, 2),
                                                                      np.round(avg, 2)),
                 pad=70, size=20)
    plt.suptitle("Fractional change ({0}_{1}z-{2}_{3}z)/{2}_{3}z".format(width_type2, depth2, width_type1, depth1))

    ax.set_xlabel("Destination stations")
    ax.set_ylabel("Source stations", labelpad=20)
    ax.set_xticks(np.arange(len(codes)))
    ax.set_yticks(np.arange(len(codes)))
    ax.set_xticklabels(np.arange(1, len(codes) + 1))
    ax.set_yticklabels(np.arange(1, len(codes) + 1))
    ax.set_xticklabels(codes)
    ax.set_yticklabels(codes)
    plt.tick_params(axis='both', which='major', labelsize=7)

    # Y axis labels on top
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    plt.xticks(rotation=90)

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(len(codes) + 1) - .5, minor=True)
    ax.set_yticks(np.arange(len(codes) + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

    # divnorm = colors.TwoSlopeNorm(vmin=min_fraction, vcenter=0, vmax=max_fraction)
    # plt.imshow(fraction, cmap=plt.cm.coolwarm, norm=divnorm)
    # else:
    plt.imshow(fraction, cmap=plt.cm.Wistia)

    cbar = plt.colorbar(orientation='vertical')
    cbar.set_label('Fractional change (%)', size=20)
    cbar.ax.tick_params(labelsize=20)
    # plt.show()
    plt.savefig(
        home_folder + "{4}/t{0}m/FractionalChange_z{0}{1}_z{2}{3}_{5}.png".format(depth2, width_type2, depth1,
                                                                                  width_type1, dataset, species),
        bbox_inches='tight',
        pad_inches=0.2)


def main():
    data1 = np.load(
        home_folder + '{0}/t{1}m/{2}/Stations_minT_connectivity_{1}z_{3}_{2}.npz'.format(dataset, depth1, width_type1,
                                                                                         species), allow_pickle=True)
    data2 = np.load(
        home_folder + '{0}/t{1}m/{2}/Stations_minT_connectivity_{1}z_{3}_{2}.npz'.format(dataset, depth2, width_type2,
                                                                                         species), allow_pickle=True)
    codes = data1['codes']
    base_matrix = data1['matrix']
    data2_matrix = data2['matrix']
    print('Data1, Data2 maximum time: ', np.nanmax(base_matrix), np.nanmax(data2_matrix))

    fr = compute_fractional_change(base_matrix, data2_matrix)
    plot_change(fr, species, codes, depth1, depth2, width_type1, width_type2)


if __name__ == '__main__':
    main()
