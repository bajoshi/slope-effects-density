from __future__ import division

import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel

import sys
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText
import matplotlib as mpl

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'
taffy_dir = home + '/Desktop/ipac/taffy/'

sys.path.append(slopedir)
sys.path.append(taffy_dir + 'codes/')
import plot_maps as pm
import vel_channel_map as vcm

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

    ax.scatter(x, np.log10(y), s=5, c=color_arr[diam_bin_idx], alpha=0.3, edgecolors='none')
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
    # smooth counts to get smoother contours
    kernel = Gaussian2DKernel(stddev=1.4)
    counts = convolve(counts, kernel, boundary='extend')

    print np.min(counts), np.max(counts)
    levels_to_plot, cb_lw, vmin = get_levels_to_plot(diam_bin)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=max(levels_to_plot))

    c = ax.contour(counts.transpose(), levels=levels_to_plot, \
        extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
        cmap=cm.viridis, linestyles='solid', linewidths=2, \
        interpolation='None', zorder=10, norm=norm)

    # plot colorbar inside figure
    cbaxes = inset_axes(ax, width='3%', height='52%', loc=7, bbox_to_anchor=[-0.05, -0.2, 1, 1], bbox_transform=ax.transAxes)
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
        levels_to_plot = [10, 100, 200, 335, 450]
        cb_lw = 35.0
        vmin = -150
    elif diam_bin == '2to3':
        levels_to_plot = [3, 20, 50, 80, 120]
        cb_lw = 35.0
        vmin = -40
    elif diam_bin == '3to4':
        levels_to_plot = [3, 15, 30, 50, 68]
        cb_lw = 35.0
        vmin = -20
    elif diam_bin == '4to5':
        levels_to_plot = [3, 10, 15, 24, 34]
        cb_lw = 35.0
        vmin = -8
    elif diam_bin == '5to6':
        levels_to_plot = [2, 10, 18, 25, 37]
        cb_lw = 35.0
        vmin = -10
    elif diam_bin == '6to7':
        levels_to_plot = [2, 10, 18, 25, 37]
        cb_lw = 35.0
        vmin = -10
    elif diam_bin == '7to8':
        levels_to_plot = [2, 8, 13, 18, 24]
        cb_lw = 35.0
        vmin = -5
    elif diam_bin == '8to9':
        levels_to_plot = [2, 8, 13, 18, 23]
        cb_lw = 35.0
        vmin = -5
    elif diam_bin == '9to10':
        levels_to_plot = [1, 4, 7, 9, 12]
        cb_lw = 35.0
        vmin = -2
    elif diam_bin == '10to15':
        levels_to_plot = [5, 20, 40, 65, 83, 95]
        cb_lw = 28.0
        vmin = -20
    elif diam_bin == '15to20':
        levels_to_plot = [3, 12, 30, 50, 68, 80]
        cb_lw = 28.0
        vmin = -20
    elif diam_bin == '20to25':
        levels_to_plot = [5, 50, 140, 200, 250, 280]
        cb_lw = 28.0
        vmin = -50
    elif diam_bin == '25to30':
        levels_to_plot = [5, 50, 140, 200, 250, 280]
        cb_lw = 28.0
        vmin = -50
    elif diam_bin == '30to35':
        levels_to_plot = [5, 50, 140, 200, 250, 280]
        cb_lw = 28.0
        vmin = -60

    return levels_to_plot, cb_lw, vmin

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