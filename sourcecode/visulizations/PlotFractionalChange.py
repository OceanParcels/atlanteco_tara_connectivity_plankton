import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import re

depth = 0
dataset = 'sample_constraints'  # 'sample_constraints'  # '2011_lombard_forams'
width_type = 'broad'
work_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task12/'
base_path = work_folder + 'Connectivities/{0}/t{1}m/'.format(dataset, depth)
label_size = 25


def compute_fractional_change(original_matrix, species, codes):
    data_new = np.load(
        base_path + 'Stations_minT_connectivity_{0}z_{1}_{2}.npz'.format(depth, species, width_type),
        allow_pickle=True)

    new_codes = data_new['codes']
    new_matrix = data_new['matrix']
    print('maximum time: ', np.nanmax(new_matrix))

    # assure order is same
    assert np.array_equal(codes, new_codes)

    # compute fraction
    fraction = (new_matrix - original_matrix) / original_matrix * 100
    # fraction = new_matrix - original_matrix
    return fraction


def plot_change(fraction, species, codes, min_limit, limit):
    # Source: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py
    # fig = plt.figure()

    avg = np.nanmean(fraction)
    min_fraction, max_fraction = np.nanmin(fraction), np.nanmax(fraction)
    print(min_fraction, max_fraction, avg)
    fig = plt.figure(figsize=(16, 14), dpi=300)
    plt.margins(0, 0)
    ax = plt.gca()
    ax.set_title("minimum: {0}%, maximum: {1}%".format(round(min_fraction, 2), round(max_fraction, 2)),
                 pad=70, fontsize=label_size)
    plt.suptitle("Average Fractional change for {2} at depth {0}m: {1}%".format(depth, np.round(avg, 2), species),
                 fontsize=label_size)

    ax.set_xlabel("Destination stations", labelpad=30, fontsize=label_size)
    ax.set_ylabel("Source stations", labelpad=30, fontsize=label_size)
    ax.set_xticks(np.arange(len(codes)))
    ax.set_yticks(np.arange(len(codes)))
    ax.set_xticklabels(codes)
    ax.set_yticklabels(codes)
    plt.tick_params(axis='both', which='major', labelsize=label_size)

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

    if min_limit:
        divnorm = colors.TwoSlopeNorm(vmin=min_limit, vcenter=0, vmax=limit)
        plt.imshow(fraction, cmap=plt.cm.coolwarm, norm=divnorm)
    else:
        plt.imshow(fraction, cmap=plt.cm.plasma_r)
    # plt.clim(0, limit)
    cbar = plt.colorbar(orientation='vertical')

    cbar.set_label('Fractional change (%)', size=label_size)
    cbar.ax.tick_params(labelsize=20)
    # plt.show()
    if min_limit:
        plt.savefig(
            base_path + "FractionalChange_z{0}m_{1}_{2}_Fr{3}-{4}_g.png".format(depth, species, width_type, min_limit,
                                                                                limit),
            bbox_inches='tight',
            pad_inches=0.2)
    else:
        plt.savefig(
            base_path + "FractionalChange_z{0}m_{1}_{2}_g.png".format(depth, species, width_type),
            bbox_inches='tight',
            pad_inches=0.2)


def main():
    # save format
    # np.savez_compressed(home_folder + 'Stations_min-T_connectivity.npz', codes=final_stations_code, matrix=min_T_matrix)
    data = np.load(
        work_folder + 'Connectivities/Stations_minT_connectivity_0z_NoConstraints_passive.npz',
        allow_pickle=True)
    codes = data['codes']
    original_matrix = data['matrix']
    print('maximum time: ', np.nanmax(original_matrix))
    TR_regexp = re.compile(r'TN_*')
    AP_regexp = re.compile(r'AP_*')
    species_info = pd.read_csv(work_folder + dataset + '.csv',
                               delimiter=';|,', keep_default_na=True, header=0)
    for index, entry in species_info.iterrows():
        fr = compute_fractional_change(original_matrix, entry['Species'], codes)
        if TR_regexp.search(entry['Species']):
            min_limit, limit = -45, 125
        elif AP_regexp.search(entry['Species']):
            min_limit, limit = -40, 390
        else:
            min_limit, limit = None, None
        plot_change(fr, entry['Species'], codes, min_limit, limit)


if __name__ == '__main__':
    main()
