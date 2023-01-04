import h3
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature

custom_text_size = 22


def load_mask_file(file):
    mask_ds = xr.open_dataset(file, decode_times=False).load()
    x = mask_ds['glamf']
    y = mask_ds['gphif']
    c = mask_ds['tmask'][:]
    return x, y, c


def masterhex_to_latlon(master_hex_ids, path):
    centers = [h3.h3_to_geo(str(master_hex_ids[ind])) for ind in path]
    lats = [x1[0] for x1 in centers]
    lons = [x1[1] for x1 in centers]
    return lats, lons


def plot_paths(x, y, c, master_hex_ids, forward_path, backward_path, f_time_laps, b_time_laps, s_code, d_code,
               path_to_compute):
    fig = plt.figure()
    ax = plt.axes()
    colormap = clr.ListedColormap(['grey', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    # ax.set_title('{0} Path and Time estimate between station {1} and {2}'.format(path_to_compute, s_code, d_code))
    ax.set_title('{0} Path and Time estimate from {1} to {2}'.format(path_to_compute, s_code, d_code))

    # https://stackoverflow.com/questions/48520393/filling-shapefile-polygons-with-a-color-in-matplotlib
    def plot_path(time_laps, path, cmap, s, d, linecolor):
        # time_laps = np.append(0, time_laps)
        norm = plt.Normalize(vmin=0, vmax=time_laps[-1])
        patches = []

        lats, lons = masterhex_to_latlon(master_hex_ids, path)
        ax.plot(lons, lats, color=linecolor, linestyle='solid', alpha=0.3)

        for i in range(len(path)):
            color = cmap(norm(time_laps[i]))
            polygons = h3.h3_set_to_multi_polygon([master_hex_ids[path[i]]], geo_json=True)
            patches.append(Polygon(polygons[0][0], True, color=color))

        pc = PatchCollection(patches, match_original=True, edgecolor=None, linewidths=None, zorder=10)
        ax.add_collection(pc)
        sm1 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm1.set_array(time_laps)
        cb1 = fig.colorbar(sm1, ax=ax, orientation='vertical')
        cb1.set_label('{0} to {1} path time laps in years'.format(s, d))

    plot_path(b_time_laps, backward_path, plt.cm.Wistia.reversed(), d_code, s_code, 'sandybrown')
    plot_path(f_time_laps, forward_path, plt.cm.winter, s_code, d_code, 'deepskyblue')
    plt.show()


def plot_shortest_paths_subset(x, y, c, master_hex_ids, f_paths, b_paths, s_code, d_code, depth):
    fig, ax = plt.subplots(1, subplot_kw={'projection': ccrs.PlateCarree()}, dpi=300, figsize=(14, 12))
    ax.add_feature(cfeature.LAND, color='gainsboro')
    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 20, 'color': 'k'}
    gl.ylabel_style = {'size': 20, 'color': 'k'}
    ax.set_extent([-95, 20, -78, 85])
    # colormap = clr.ListedColormap(['gainsboro', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    # ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    # ax.set_title('Sample of Minimum Time Paths between station {0} and {1}'.format(s_code, d_code))
    ax.set_title('Random sample of Minimum Time Paths from {0} to {1}'.format(s_code, d_code))

    def plot_path(path, color, s, d):
        lats, lons = masterhex_to_latlon(master_hex_ids, path)
        ax.plot(lons, lats, color=color, linestyle='solid', alpha=0.3)
        ax.scatter(lons, lats, color=color, s=0.4)
        ax.text(lons[0], lats[0], s, bbox=dict(facecolor='white', alpha=0.7, pad=0.2, edgecolor='none'),
                fontsize=custom_text_size)
        ax.text(lons[-1], lats[-1], d, bbox=dict(facecolor='white', alpha=0.7, pad=0.2, edgecolor='none'),
                fontsize=custom_text_size)

    # [plot_path(paths[i], colors(i)) for i in range(len(paths))]
    [plot_path(f_paths[i], 'red', s_code, d_code) for i in range(len(f_paths))]
    [plot_path(b_paths[i], 'blue', d_code, s_code) for i in range(len(b_paths))]

    textstr = '\n'.join((
        r'[Time, paths count] at depth %d m' % (depth),
        r'%s to %s: [%.2f years, %d]' % (s_code, d_code, round((len(f_paths[0]) - 1) / 12, 2), len(f_paths)),
        r'%s to %s: [%.2f years, %d]' % (d_code, s_code, round((len(b_paths[0]) - 1) / 12, 2), len(b_paths))))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='square', facecolor='lemonchiffon', alpha=0.8)

    # place a text box in upper left in axes coords
    # 0.58 for 0m, 0.56 for 50m, 0.55 for 100m,
    ax.text(0.08, 0.085, textstr, transform=ax.transAxes, fontsize=custom_text_size,
            verticalalignment='center', bbox=props)

    # plt.show()
    plt.savefig(
        '/Users/dmanral/Desktop/Analysis/TARA/Task12/PairPlots/{0}_{1}_z{2}_passive.png'.format(s_code, d_code,
                                                                                                        depth),
        bbox_inches='tight',
        pad_inches=0.2)


def plot_time_vs_transitions(master_hex_ids, path, time_laps, s, d):
    lats = masterhex_to_latlon(master_hex_ids, path)[0]
    ax = plt.axes()
    plt.plot(lats, time_laps, marker='*', linestyle='--', linewidth='0.5')
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Time (in years)")
    ax.set_title("Time laps for path from {0} to {1}".format(s, d))
    # plt.xticks(np.arange(min(np.round(lats, 0)), max(np.round(lats, 0)), 10))

    # fig = plt.figure()
    # ax1 = fig.add_subplot(111, label='1')
    # # ax2 = fig.add_subplot(111, label='2')
    # plt.plot(steps, time_laps, marker='*', linestyle='--')
    # ax1.set_xlabel("Transition Steps")
    # ax1.set_ylabel("Time (years)")
    #
    # def forward(x):
    #     return np.interp(x, steps, lats)
    #
    # def inverse(x):
    #     return np.interp(x, lats, steps)
    #
    # secax = ax1.secondary_xaxis('top', functions=(forward, inverse))
    # secax.set_xlabel('Latitude')
    plt.show()


def plots_transitions(x, y, c, master_hex_ids, grid_index, direction, s_ids, probs, min_temps, max_temps):
    fig = plt.figure()
    ax = plt.axes()
    colormap = clr.ListedColormap(['gainsboro', 'white'])
    # remove the first row and first column from the glamf/gphif to access points enclosed in the center
    ax.pcolormesh(x[0], y[0], c[0, 0, 1:, 1:], cmap=colormap)
    print("base map ready")
    ax.set_title('{0} {1} connections to node {2}'.format(len(s_ids), direction, grid_index))

    grid_lat, grid_lon = h3.h3_to_geo(master_hex_ids[grid_index])
    plt.plot(grid_lon, grid_lat, c="red")
    connected_lats, connected_lons = masterhex_to_latlon(master_hex_ids, s_ids)
    plt.scatter(connected_lons, connected_lats, c='blue')
    [plt.plot([grid_lon, connected_lons[i]], [grid_lat, connected_lats[i]], linestyle='--', linewidth='0.5') for i in
     range(len(connected_lats))]

    for i in range(len(connected_lats)):
        xy = (connected_lons[i], connected_lats[i])
        # ax.annotate('(%.4f, %.2f, %.2f)' % (probs[i], min_temps[i], max_temps[i]), xy=xy, textcoords='data')
        ax.annotate('(%.4f)' % probs[i], xy=xy, textcoords='data')
    plt.show()
