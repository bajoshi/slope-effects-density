from __future__ import division

import numpy as np
from astropy.convolution import convolve, Gaussian2DKernel

import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText
import matplotlib.gridspec as gridspec

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

    print "Min and max point number density values in bins", str("{:.3}".format(np.min(counts))), str("{:.3}".format(np.max(counts)))
    levels_to_plot, cb_lw, vmin = get_levels_to_plot(diam_bin, plottype='cumulative')
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
    plt.cla()
    plt.clf()
    plt.close()

    return None

def get_levels_to_plot(diam_bin, plottype):

    if plottype == 'cumulative':
        cumulative = True
        nooverlap = False
        nvalue = False
    elif plottype == 'nooverlap':
        nooverlap = True
        cumulative = False
        nvalue = False
    elif plottype == 'Nvalue':
        nvalue = True
        cumulative = False
        nooverlap = False

    if diam_bin == '1to2':
        levels_to_plot = [8, 80, 200, 450, 800]
        if nvalue:
            levels_to_plot = [20, 150, 600, 1200, 1600]
        cb_lw = 35.0
        vmin = -150
    elif diam_bin == '2to3':
        levels_to_plot = [3, 20, 50, 150, 300]
        if nvalue:
            levels_to_plot = [5, 80, 300, 600, 1100]
        cb_lw = 35.0
        vmin = -40
    elif diam_bin == '3to4':
        if cumulative:
            levels_to_plot = [3, 15, 30, 50, 68]
        elif nooverlap:
            levels_to_plot = [3, 15, 50, 85, 130]
        elif nvalue:
            levels_to_plot = [5, 90, 350, 600, 1000]
        cb_lw = 35.0
        vmin = -20
    elif diam_bin == '4to5':
        if cumulative:
            levels_to_plot = [3, 10, 15, 24, 34]
        elif nooverlap:
            levels_to_plot = [3, 10, 40, 80, 125]
        elif nvalue:
            levels_to_plot = [5, 90, 350, 600, 950]
        cb_lw = 35.0
        vmin = -8
    elif diam_bin == '5to6':
        if cumulative:
            levels_to_plot = [2, 10, 18, 25, 37]
        elif nooverlap:
            levels_to_plot = [2, 10, 22, 45, 62]
        elif nvalue:
            levels_to_plot = [5, 90, 350, 600, 950]
        cb_lw = 35.0
        vmin = -10
    elif diam_bin == '6to7':
        levels_to_plot = [2, 10, 30, 60, 100]
        if nvalue:
            levels_to_plot = [5, 80, 350, 600, 900]
        cb_lw = 35.0
        vmin = -10
    elif diam_bin == '7to8':
        if cumulative:
            levels_to_plot = [2, 8, 13, 18, 24]
        elif nooverlap:
            levels_to_plot = [2, 8, 13, 25, 30]
        elif nvalue:
            levels_to_plot = [5, 80, 350, 600, 900]
        cb_lw = 35.0
        vmin = -5
    elif diam_bin == '8to9':
        if cumulative:
            levels_to_plot = [2, 8, 13, 18, 23]
        elif nooverlap:
            levels_to_plot = [2, 8, 13, 30, 42]
        elif nvalue:
            levels_to_plot = [5, 80, 350, 600, 875]
        cb_lw = 35.0
        vmin = -5
    elif diam_bin == '9to10':
        levels_to_plot = [1, 6, 20, 40, 70]
        if nvalue:
            levels_to_plot = [5, 80, 350, 600, 850]
        cb_lw = 35.0
        vmin = -2
    elif diam_bin == '10to15':
        levels_to_plot = [4, 20, 40, 65, 90, 125]
        if nvalue:
            levels_to_plot = [5, 30, 150, 420, 700, 850]
        cb_lw = 28.0
        vmin = -20
    elif diam_bin == '15to20':
        levels_to_plot = [3, 12, 30, 50, 68, 80]
        if nvalue:
            levels_to_plot = [5, 30, 150, 420, 680, 760]
        cb_lw = 28.0
        vmin = -20
    elif diam_bin == '20to25':
        if cumulative:
            levels_to_plot = [5, 50, 140, 200, 250, 280]
        elif nooverlap:
            levels_to_plot = [3, 15, 50, 80, 120, 145]
        elif nvalue:
            levels_to_plot = [5, 30, 150, 400, 620, 680]
        cb_lw = 28.0
        vmin = -50
    elif diam_bin == '25to30':
        if cumulative:
            levels_to_plot = [5, 50, 140, 200, 250, 280]
        elif nooverlap:
            levels_to_plot = [3, 15, 90, 140, 200, 260]
        elif nvalue:
            levels_to_plot = [5, 30, 150, 300, 450, 550]
        cb_lw = 28.0
        vmin = -50
    elif diam_bin == '30to35':
        levels_to_plot = [5, 50, 140, 200, 250, 280]
        if nvalue:
            levels_to_plot = [3, 20, 70, 140, 220, 290]
        cb_lw = 28.0
        vmin = -60

    return levels_to_plot, cb_lw, vmin

def make_cumulative_plots():

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

    return None

def make_no_overlap_Nvalue_plots(density_arr, slope_arr, diam_bin_min, diam_bin_max, color, plottype):

    # make figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=16)
    ax.set_ylabel(r'$\mathrm{log(Density)}$', fontsize=16)

    x = slope_arr
    y = density_arr

    ax.scatter(x, np.log10(y), s=5, c=color, alpha=0.3, edgecolors='none')
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

    print "Min and max point number density values in bins", str("{:.3}".format(np.min(counts))), str("{:.3}".format(np.max(counts)))
    diam_bin = str(diam_bin_min) + 'to' + str(diam_bin_max)
    levels_to_plot, cb_lw, vmin = get_levels_to_plot(diam_bin, plottype=plottype)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=max(levels_to_plot))

    c = ax.contour(counts.transpose(), levels=levels_to_plot, \
        extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
        cmap=cm.viridis, linestyles='solid', linewidths=2, \
        zorder=10, norm=norm)

    # plot colorbar inside figure
    cbaxes = inset_axes(ax, width='3%', height='52%', loc=7, bbox_to_anchor=[-0.05, -0.2, 1, 1], bbox_transform=ax.transAxes)
    cb = plt.colorbar(c, cax=cbaxes, ticks=[min(levels_to_plot), max(levels_to_plot)], orientation='vertical')
    cb.ax.get_children()[0].set_linewidths(cb_lw)

    # add ticks and grid
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    # save the figure
    if plottype == 'Nvalue':
        # add text on figure to indicate diameter bin
        diambinbox = TextArea(r"$\mathrm{N \geq\ }$" + str(diam_bin_min) + r'$\mathrm{\, km}$', textprops=dict(color='k', size=14))
        anc_diambinbox = AnchoredOffsetbox(loc=2, child=diambinbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.6, 0.1),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_diambinbox)
        fig.savefig(slope_extdir + 'slope_v_density_withcontour_Nvalue' + str(diam_bin_min) + 'km.png', dpi=300, bbox_inches='tight')

    elif plottype == 'nooverlap':
        # add text on figure to indicate diameter bin
        diambinbox = TextArea(str(diam_bin_min) + r'$\mathrm{\ to\ }$' + str(diam_bin_max) + r'$\mathrm{\,km}$', \
            textprops=dict(color='k', size=14))
        anc_diambinbox = AnchoredOffsetbox(loc=2, child=diambinbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.6, 0.1),\
                                             bbox_transform=ax.transAxes, borderpad=0.0)
        ax.add_artist(anc_diambinbox)

        fig.savefig(slope_extdir + 'slope_v_density_withcontour_nooverlap' + diam_bin + 'km.png', dpi=300, bbox_inches='tight')

    plt.cla()
    plt.clf()
    plt.close()

    return None

def call_no_overlap_plots():

    # read in all arrays
    density_diambin_1_2 = np.load(slope_extdir + 'density_diambin_1_2.npy')
    density_diambin_2_3 = np.load(slope_extdir + 'density_diambin_2_3.npy')
    density_diambin_3_4 = np.load(slope_extdir + 'density_diambin_3_4.npy')
    density_diambin_4_5 = np.load(slope_extdir + 'density_diambin_4_5.npy')
    density_diambin_5_6 = np.load(slope_extdir + 'density_diambin_5_6.npy')
    density_diambin_6_7 = np.load(slope_extdir + 'density_diambin_6_7.npy')
    density_diambin_7_8 = np.load(slope_extdir + 'density_diambin_7_8.npy')
    density_diambin_8_9 = np.load(slope_extdir + 'density_diambin_8_9.npy')
    density_diambin_9_10 = np.load(slope_extdir + 'density_diambin_9_10.npy')
    density_diambin_10_15 = np.load(slope_extdir + 'density_diambin_10_15.npy')
    density_diambin_15_20 = np.load(slope_extdir + 'density_diambin_15_20.npy')
    density_diambin_20_25 = np.load(slope_extdir + 'density_diambin_20_25.npy')
    density_diambin_25_30 = np.load(slope_extdir + 'density_diambin_25_30.npy')
    density_diambin_30_35 = np.load(slope_extdir + 'density_diambin_30_35.npy')

    slope_diambin_1_2 = np.load(slope_extdir + 'slope_diambin_1_2.npy')
    slope_diambin_2_3 = np.load(slope_extdir + 'slope_diambin_2_3.npy')
    slope_diambin_3_4 = np.load(slope_extdir + 'slope_diambin_3_4.npy')
    slope_diambin_4_5 = np.load(slope_extdir + 'slope_diambin_4_5.npy')
    slope_diambin_5_6 = np.load(slope_extdir + 'slope_diambin_5_6.npy')
    slope_diambin_6_7 = np.load(slope_extdir + 'slope_diambin_6_7.npy')
    slope_diambin_7_8 = np.load(slope_extdir + 'slope_diambin_7_8.npy')
    slope_diambin_8_9 = np.load(slope_extdir + 'slope_diambin_8_9.npy')
    slope_diambin_9_10 = np.load(slope_extdir + 'slope_diambin_9_10.npy')
    slope_diambin_10_15 = np.load(slope_extdir + 'slope_diambin_10_15.npy')
    slope_diambin_15_20 = np.load(slope_extdir + 'slope_diambin_15_20.npy')
    slope_diambin_20_25 = np.load(slope_extdir + 'slope_diambin_20_25.npy')
    slope_diambin_25_30 = np.load(slope_extdir + 'slope_diambin_25_30.npy')
    slope_diambin_30_35 = np.load(slope_extdir + 'slope_diambin_30_35.npy')

    # make plots with only contributions from craters that 
    # actually fall within the specified diameter bin
    make_no_overlap_Nvalue_plots(density_diambin_1_2, slope_diambin_1_2, 1, 2, 'midnightblue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_2_3, slope_diambin_2_3, 2, 3, 'blue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_3_4, slope_diambin_3_4, 3, 4, 'royalblue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_4_5, slope_diambin_4_5, 4, 5, 'dodgerblue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_5_6, slope_diambin_5_6, 5, 6, 'deepskyblue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_6_7, slope_diambin_6_7, 6, 7, 'steelblue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_7_8, slope_diambin_7_8, 7, 8, 'slateblue', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_8_9, slope_diambin_8_9, 8, 9, 'rebeccapurple', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_9_10, slope_diambin_9_10, 9, 10, 'darkcyan', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_10_15, slope_diambin_10_15, 10, 15, 'green', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_15_20, slope_diambin_15_20, 15, 20, 'olive', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_20_25, slope_diambin_20_25, 20, 25, 'goldenrod', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_25_30, slope_diambin_25_30, 25, 30, 'darkorchid', 'nooverlap')
    make_no_overlap_Nvalue_plots(density_diambin_30_35, slope_diambin_30_35, 30, 35, 'maroon', 'nooverlap')

    return None

def call_Nvalue_plots():

    # read in all arrays
    density_diambin_1 = np.load(slope_extdir + 'density_diambin_1.npy')
    density_diambin_2 = np.load(slope_extdir + 'density_diambin_2.npy')
    density_diambin_3 = np.load(slope_extdir + 'density_diambin_3.npy')
    density_diambin_4 = np.load(slope_extdir + 'density_diambin_4.npy')
    density_diambin_5 = np.load(slope_extdir + 'density_diambin_5.npy')
    density_diambin_6 = np.load(slope_extdir + 'density_diambin_6.npy')
    density_diambin_7 = np.load(slope_extdir + 'density_diambin_7.npy')
    density_diambin_8 = np.load(slope_extdir + 'density_diambin_8.npy')
    density_diambin_9 = np.load(slope_extdir + 'density_diambin_9.npy')
    density_diambin_10 = np.load(slope_extdir + 'density_diambin_10.npy')
    density_diambin_15 = np.load(slope_extdir + 'density_diambin_15.npy')
    density_diambin_20 = np.load(slope_extdir + 'density_diambin_20.npy')
    density_diambin_25 = np.load(slope_extdir + 'density_diambin_25.npy')
    density_diambin_30 = np.load(slope_extdir + 'density_diambin_30.npy')

    slope_diambin_1 = np.load(slope_extdir + 'slope_diambin_1.npy')
    slope_diambin_2 = np.load(slope_extdir + 'slope_diambin_2.npy')
    slope_diambin_3 = np.load(slope_extdir + 'slope_diambin_3.npy')
    slope_diambin_4 = np.load(slope_extdir + 'slope_diambin_4.npy')
    slope_diambin_5 = np.load(slope_extdir + 'slope_diambin_5.npy')
    slope_diambin_6 = np.load(slope_extdir + 'slope_diambin_6.npy')
    slope_diambin_7 = np.load(slope_extdir + 'slope_diambin_7.npy')
    slope_diambin_8 = np.load(slope_extdir + 'slope_diambin_8.npy')
    slope_diambin_9 = np.load(slope_extdir + 'slope_diambin_9.npy')
    slope_diambin_10 = np.load(slope_extdir + 'slope_diambin_10.npy')
    slope_diambin_15 = np.load(slope_extdir + 'slope_diambin_15.npy')
    slope_diambin_20 = np.load(slope_extdir + 'slope_diambin_20.npy')
    slope_diambin_25 = np.load(slope_extdir + 'slope_diambin_25.npy')
    slope_diambin_30 = np.load(slope_extdir + 'slope_diambin_30.npy')

    # make plots with only contributions from craters that 
    # actually fall within the specified diameter bin
    make_no_overlap_Nvalue_plots(density_diambin_1, slope_diambin_1, 1, 2, 'midnightblue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_2, slope_diambin_2, 2, 3, 'blue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_3, slope_diambin_3, 3, 4, 'royalblue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_4, slope_diambin_4, 4, 5, 'dodgerblue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_5, slope_diambin_5, 5, 6, 'deepskyblue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_6, slope_diambin_6, 6, 7, 'steelblue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_7, slope_diambin_7, 7, 8, 'slateblue', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_8, slope_diambin_8, 8, 9, 'rebeccapurple', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_9, slope_diambin_9, 9, 10, 'darkcyan', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_10, slope_diambin_10, 10, 15, 'green', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_15, slope_diambin_15, 15, 20, 'olive', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_20, slope_diambin_20, 20, 25, 'goldenrod', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_25, slope_diambin_25, 25, 30, 'darkorchid', 'Nvalue')
    make_no_overlap_Nvalue_plots(density_diambin_30, slope_diambin_30, 30, 35, 'maroon', 'Nvalue')

    return None

def make_nooverlap_nvalue_grid_plots(slope_list, density_list, plottype):

    gs = gridspec.GridSpec(4,4)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure(figsize=(9,5.5))
    ax_1to2   = fig.add_subplot(gs[0, 0])
    ax_2to3   = fig.add_subplot(gs[0, 1])
    ax_3to4   = fig.add_subplot(gs[0, 2])
    ax_4to5   = fig.add_subplot(gs[0, 3])
    ax_5to6   = fig.add_subplot(gs[1, 0])
    ax_6to7   = fig.add_subplot(gs[1, 1])
    ax_7to8   = fig.add_subplot(gs[1, 2])
    ax_8to9   = fig.add_subplot(gs[1, 3])
    ax_9to10  = fig.add_subplot(gs[2, 0])
    ax_10to15 = fig.add_subplot(gs[2, 1])
    ax_15to20 = fig.add_subplot(gs[2, 2])
    ax_20to25 = fig.add_subplot(gs[2, 3])
    ax_25to30 = fig.add_subplot(gs[3, 1])
    ax_30to35 = fig.add_subplot(gs[3, 2])

    all_axes = [ax_1to2, ax_2to3, ax_3to4, ax_4to5, ax_5to6, ax_6to7, ax_7to8, ax_8to9, ax_9to10, \
    ax_10to15, ax_15to20, ax_20to25, ax_25to30, ax_30to35]
    diam_bins = ['1to2', '2to3', '3to4', '4to5', '5to6', '6to7', '7to8', '8to9', '9to10', \
    '10to15', '15to20', '20to25', '25to30', '30to35']
    color_list = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#543005','#fdbf6f',\
    '#ff7f00','#cab2d6','#bf812d','#6a3d9a','#b15928','#01665e']

    ax_25to30.set_xlabel(r'$\mathrm{Slope}$', fontsize=14)
    ax_25to30.xaxis.set_label_coords(1.05, -0.25)

    ax_9to10.set_ylabel(r'$\mathrm{log(Density)}$', fontsize=14)
    ax_9to10.yaxis.set_label_coords(-0.3, 1.05)

    for i in range(len(all_axes)):

        all_axes[i].scatter(slope_list[i], np.log10(density_list[i]), s=5, c=color_list[i], alpha=0.25, edgecolors='none')

        all_axes[i].set_ylim(-8, 0.5)
        all_axes[i].set_xlim(0, 35)

        x = slope_list[i]
        y = density_list[i]

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

        print "Min and max point number density values in bins", str("{:.3}".format(np.min(counts))), str("{:.3}".format(np.max(counts)))
        diam_bin = diam_bins[i]
        diam_bin_min = diam_bin.split('to')[0]
        diam_bin_max = diam_bin.split('to')[1]
        levels_to_plot, cb_lw, vmin = get_levels_to_plot(diam_bin, plottype=plottype)
        norm = mpl.colors.Normalize(vmin=vmin, vmax=max(levels_to_plot))

        c = all_axes[i].contour(counts.transpose(), levels=levels_to_plot, \
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
            cmap=cm.viridis, linestyles='solid', linewidths=1, \
            zorder=10, norm=norm)

        # add minor ticks and grid
        all_axes[i].minorticks_on()
        all_axes[i].tick_params('both', width=1, length=3, which='minor')
        all_axes[i].tick_params('both', width=1, length=4.7, which='major')

        if i <= 11:
            all_axes[i].set_xticklabels([])

        if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13:
            all_axes[i].set_yticklabels([])

        if plottype == 'Nvalue':
            all_axes[i].text(0.6, 0.15, r"$\mathrm{N(}$" + str(diam_bin_min) + r'$\mathrm{)}$', \
                verticalalignment='top', horizontalalignment='left', \
                transform=all_axes[i].transAxes, color='k', size=5)

        elif plottype == 'nooverlap':
            all_axes[i].text(0.6, 0.15, str(diam_bin_min) + r'$\mathrm{\ to\ }$' + str(diam_bin_max) + r'$\mathrm{\,km}$', \
                verticalalignment='top', horizontalalignment='left', \
                transform=all_axes[i].transAxes, color='k', size=5)

    ax_9to10.set_xticklabels(['0', '10', '20', '30'])
    ax_20to25.set_xticklabels(['', '10', '20', '30'])
    ax_25to30.set_xticklabels(['0', '10', '20', '30'])
    ax_30to35.set_xticklabels(['0', '10', '20', '30'])

    # save the figure
    if plottype == 'Nvalue':
        fig.savefig(slope_extdir + 'slope_v_density_withcontour_grid_Nvalue.png', dpi=300, bbox_inches='tight')
    elif plottype == 'nooverlap':
        fig.savefig(slope_extdir + 'slope_v_density_withcontour_grid_nooverlap.png', dpi=300, bbox_inches='tight')

    return None

if __name__ == '__main__':    

    # These three function calls below will make
    # individual plots for each diam range/bin.
    # If you want all the diam ranges/bins on a 
    # single figure with a grid then use the code
    # below, after these 3 lines.
    #make_cumulative_plots()
    #call_no_overlap_plots()
    #call_Nvalue_plots()

    # ---------------------------------- GRID PLOTS ---------------------------------- #
    # ---------------------------------- No Overlap ---------------------------------- #
    # read in all arrays
    density_diambin_1_2 = np.load(slope_extdir + 'density_diambin_1_2.npy')
    density_diambin_2_3 = np.load(slope_extdir + 'density_diambin_2_3.npy')
    density_diambin_3_4 = np.load(slope_extdir + 'density_diambin_3_4.npy')
    density_diambin_4_5 = np.load(slope_extdir + 'density_diambin_4_5.npy')
    density_diambin_5_6 = np.load(slope_extdir + 'density_diambin_5_6.npy')
    density_diambin_6_7 = np.load(slope_extdir + 'density_diambin_6_7.npy')
    density_diambin_7_8 = np.load(slope_extdir + 'density_diambin_7_8.npy')
    density_diambin_8_9 = np.load(slope_extdir + 'density_diambin_8_9.npy')
    density_diambin_9_10 = np.load(slope_extdir + 'density_diambin_9_10.npy')
    density_diambin_10_15 = np.load(slope_extdir + 'density_diambin_10_15.npy')
    density_diambin_15_20 = np.load(slope_extdir + 'density_diambin_15_20.npy')
    density_diambin_20_25 = np.load(slope_extdir + 'density_diambin_20_25.npy')
    density_diambin_25_30 = np.load(slope_extdir + 'density_diambin_25_30.npy')
    density_diambin_30_35 = np.load(slope_extdir + 'density_diambin_30_35.npy')

    slope_diambin_1_2 = np.load(slope_extdir + 'slope_diambin_1_2.npy')
    slope_diambin_2_3 = np.load(slope_extdir + 'slope_diambin_2_3.npy')
    slope_diambin_3_4 = np.load(slope_extdir + 'slope_diambin_3_4.npy')
    slope_diambin_4_5 = np.load(slope_extdir + 'slope_diambin_4_5.npy')
    slope_diambin_5_6 = np.load(slope_extdir + 'slope_diambin_5_6.npy')
    slope_diambin_6_7 = np.load(slope_extdir + 'slope_diambin_6_7.npy')
    slope_diambin_7_8 = np.load(slope_extdir + 'slope_diambin_7_8.npy')
    slope_diambin_8_9 = np.load(slope_extdir + 'slope_diambin_8_9.npy')
    slope_diambin_9_10 = np.load(slope_extdir + 'slope_diambin_9_10.npy')
    slope_diambin_10_15 = np.load(slope_extdir + 'slope_diambin_10_15.npy')
    slope_diambin_15_20 = np.load(slope_extdir + 'slope_diambin_15_20.npy')
    slope_diambin_20_25 = np.load(slope_extdir + 'slope_diambin_20_25.npy')
    slope_diambin_25_30 = np.load(slope_extdir + 'slope_diambin_25_30.npy')
    slope_diambin_30_35 = np.load(slope_extdir + 'slope_diambin_30_35.npy')

    density_list_nooverlap = [density_diambin_1_2, density_diambin_2_3, density_diambin_3_4, density_diambin_4_5, \
    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9, density_diambin_9_10, \
    density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, density_diambin_25_30, density_diambin_30_35]

    slope_list_nooverlap = [slope_diambin_1_2, slope_diambin_2_3, slope_diambin_3_4, slope_diambin_4_5, \
    slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, slope_diambin_8_9, slope_diambin_9_10, \
    slope_diambin_10_15, slope_diambin_15_20, slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35]

    # ---------------------------------- N value ---------------------------------- #
    # read in all arrays
    density_diambin_1 = np.load(slope_extdir + 'density_diambin_1.npy')
    density_diambin_2 = np.load(slope_extdir + 'density_diambin_2.npy')
    density_diambin_3 = np.load(slope_extdir + 'density_diambin_3.npy')
    density_diambin_4 = np.load(slope_extdir + 'density_diambin_4.npy')
    density_diambin_5 = np.load(slope_extdir + 'density_diambin_5.npy')
    density_diambin_6 = np.load(slope_extdir + 'density_diambin_6.npy')
    density_diambin_7 = np.load(slope_extdir + 'density_diambin_7.npy')
    density_diambin_8 = np.load(slope_extdir + 'density_diambin_8.npy')
    density_diambin_9 = np.load(slope_extdir + 'density_diambin_9.npy')
    density_diambin_10 = np.load(slope_extdir + 'density_diambin_10.npy')
    density_diambin_15 = np.load(slope_extdir + 'density_diambin_15.npy')
    density_diambin_20 = np.load(slope_extdir + 'density_diambin_20.npy')
    density_diambin_25 = np.load(slope_extdir + 'density_diambin_25.npy')
    density_diambin_30 = np.load(slope_extdir + 'density_diambin_30.npy')

    slope_diambin_1 = np.load(slope_extdir + 'slope_diambin_1.npy')
    slope_diambin_2 = np.load(slope_extdir + 'slope_diambin_2.npy')
    slope_diambin_3 = np.load(slope_extdir + 'slope_diambin_3.npy')
    slope_diambin_4 = np.load(slope_extdir + 'slope_diambin_4.npy')
    slope_diambin_5 = np.load(slope_extdir + 'slope_diambin_5.npy')
    slope_diambin_6 = np.load(slope_extdir + 'slope_diambin_6.npy')
    slope_diambin_7 = np.load(slope_extdir + 'slope_diambin_7.npy')
    slope_diambin_8 = np.load(slope_extdir + 'slope_diambin_8.npy')
    slope_diambin_9 = np.load(slope_extdir + 'slope_diambin_9.npy')
    slope_diambin_10 = np.load(slope_extdir + 'slope_diambin_10.npy')
    slope_diambin_15 = np.load(slope_extdir + 'slope_diambin_15.npy')
    slope_diambin_20 = np.load(slope_extdir + 'slope_diambin_20.npy')
    slope_diambin_25 = np.load(slope_extdir + 'slope_diambin_25.npy')
    slope_diambin_30 = np.load(slope_extdir + 'slope_diambin_30.npy')

    density_list_nvalue = [density_diambin_1, density_diambin_2, density_diambin_3, density_diambin_4, \
    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8, density_diambin_9, \
    density_diambin_10, density_diambin_15, density_diambin_20, density_diambin_25, density_diambin_30]

    slope_list_nvalue = [slope_diambin_1, slope_diambin_2, slope_diambin_3, slope_diambin_4, \
    slope_diambin_5, slope_diambin_6, slope_diambin_7, slope_diambin_8, slope_diambin_9, \
    slope_diambin_10, slope_diambin_15, slope_diambin_20, slope_diambin_25, slope_diambin_30]

    make_nooverlap_nvalue_grid_plots(slope_list_nvalue, density_list_nvalue, 'Nvalue')
    make_nooverlap_nvalue_grid_plots(slope_list_nooverlap, density_list_nooverlap, 'nooverlap')

    sys.exit(0)
