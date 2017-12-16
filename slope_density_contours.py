from __future__ import division

import numpy as np

import sys
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

sys.path.append(slopedir)
import plot_maps as pm

def make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, diam_bin_min, diam_bin_max):

    # get diam bin indices
    diam_bin = str(diam_bin_min) + 'to' + str(diam_bin_max)
    diam_bin_idx = pm.get_diam_idx(color_arr, diam_bin)

    # make figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{log(Density)}$', fontsize=16)

    x = slope_arr_color[diam_bin_idx]
    y = density_arr_color[diam_bin_idx]

    ax.scatter(x, np.log10(y), s=5, c=color_arr[diam_bin_idx], alpha=0.4, edgecolors='none')
    ax.set_ylim(-8, 0.5)
    ax.set_xlim(0, 35)

    # draw contours
    # make sure the arrays dont have NaNs
    slope_fin_idx = np.where(np.isfinite(x))[0]
    density_fin_idx = np.where(np.isfinite(y))[0]
    fin_idx = np.intersect1d(slope_fin_idx, density_fin_idx)

    xp = x[fin_idx]
    yp = y[fin_idx]

    counts, xbins, ybins = np.histogram2d(xp, np.log10(yp), bins=25, normed=False)
    print np.min(counts), np.max(counts)
    levels_to_plot, cb_lw = get_levels_to_plot(diam_bin)
    c = ax.contour(counts.transpose(), levels=levels_to_plot, \
        extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
        cmap=cm.copper, linestyles='solid', linewidths=2, interpolation='None', zorder=10)

    # plot colorbar inside figure
    cbaxes = inset_axes(ax, width='3%', height='50%', loc=7, bbox_to_anchor=[-0.05, -0.2, 1, 1], bbox_transform=ax.transAxes)
    cb = plt.colorbar(c, cax=cbaxes, ticks=[min(levels_to_plot), max(levels_to_plot)], orientation='vertical')
    cb.ax.get_children()[0].set_linewidths(cb_lw)

    # add text on figure to indicate diameter bin
    diambinbox = TextArea(str(diam_bin_min) + ' to ' + str(diam_bin_max) + ' km', textprops=dict(color='k', size=14))
    anc_diambinbox = AnchoredOffsetbox(loc=2, child=diambinbox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.6, 0.1),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_diambinbox)

    # add ticks and grid
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    # save the figure
    fig.savefig(slope_extdir + 'slope_v_density_withcontour_' + diam_bin + 'km.png', dpi=300, bbox_inches='tight')
    #fig.savefig(slope_extdir + 'slope_v_density_withcontour_' + diam_bin + 'km.eps', dpi=300, bbox_inches='tight')

    #plt.show()

    return None

def get_levels_to_plot(diam_bin):

    if diam_bin == '1to2':
        levels_to_plot = [10, 50, 200, 500, 675]
        cb_lw = 34.0
    elif diam_bin == '2to3':
        levels_to_plot = [3, 10, 50, 100, 200]
        cb_lw = 34.0
    elif diam_bin == '3to4':
        levels_to_plot = [3, 10, 25, 50, 120]
        cb_lw = 34.0
    elif diam_bin == '4to5':
        levels_to_plot = [3, 12, 25, 43, 58]
        cb_lw = 34.0
    elif diam_bin == '5to6':
        levels_to_plot = [2, 10, 20, 30, 70]
        cb_lw = 34.0
    elif diam_bin == '6to7':
        levels_to_plot = [1.9, 10, 20, 30, 70]
        cb_lw = 34.0
    elif diam_bin == '7to8':
        levels_to_plot = [2, 8, 20, 35, 50]
        cb_lw = 34.0
    elif diam_bin == '8to9':
        levels_to_plot = [2, 8, 20, 35, 50]
        cb_lw = 34.0
    elif diam_bin == '9to10':
        levels_to_plot = [2, 6, 10, 20, 30]
        cb_lw = 34.0
    elif diam_bin == '10to15':
        levels_to_plot = [5, 10, 20, 50, 100, 175]
        cb_lw = 27.0
    elif diam_bin == '15to20':
        levels_to_plot = [3, 10, 20, 50, 100, 150]
        cb_lw = 27.0
    elif diam_bin == '20to25':
        levels_to_plot = [5, 20, 50, 100, 300, 550]
        cb_lw = 27.0
    elif diam_bin == '25to30':
        levels_to_plot = [5, 20, 50, 100, 250, 600]
        cb_lw = 27.0
    elif diam_bin == '30to35':
        levels_to_plot = [5, 20, 50, 100, 250, 600]
        cb_lw = 27.0

    return levels_to_plot, cb_lw

if __name__ == '__main__':

    # read arrays
    density_arr_color = np.load(slope_extdir + 'density_arr_color.npy')
    slope_arr_color = np.load(slope_extdir + 'slope_arr_color.npy')
    color_arr = np.load(slope_extdir + 'color_arr.npy')

    # do the actual plotting
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 1, 2)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 2, 3)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 3, 4)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 4, 5)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 5, 6)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 6, 7)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 7, 8)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 8, 9)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 9, 10)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 10, 15)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 15, 20)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 20, 25)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 25, 30)
    make_plot_diam_bin_with_contour(density_arr_color, slope_arr_color, color_arr, 30, 35)

    sys.exit(0)