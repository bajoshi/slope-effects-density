from __future__ import division

import numpy as np
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from scipy.misc import factorial
from scipy.special import gamma
import cPickle

import sys
import os
import time
import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'
massive_galaxies_dir = home + "/Desktop/FIGS/massive-galaxies/"

sys.path.append(slopedir)
sys.path.append(massive_galaxies_dir + 'codes/')
import average_density_plots as avg
#import cython_util_funcs
import mag_hist as mh
import slope_utils as su

def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)

def gamma_dist(x, k, theta):
    return (x**(k-1) * np.exp(-x/theta)) / (theta**k * gamma(k))

def hist_and_fit(fig, ax, slopearr, color, callcount, plottype=None):
    """
    Simply comment out the ax.hist() line at hte end
    if you only want to plot the fits and not the 
    histograms as well.
    """

    counts, bins = np.histogram(slopearr, 35, range=[0,35], density=True)
    fitting_bins = [(bins[i] + bins[i+1])/2 for i in bins[:-1]]
    fitting_bins = np.asarray(fitting_bins)

    x_plot_arr = np.linspace(0,35,1000)

    # fitting a gamma distribution
    popt, pcov = curve_fit(gamma_dist, fitting_bins, counts, p0=[3.0, 3.0])

    label_list = ['1-1.25 km', '1.25-1.5 km', '1.5-1.75 km', '1.75-2 km',\
                  '2-2.25 km', '2.25-2.5 km', '2.5-2.75 km', '2.75-3 km',\
                  '3-3.25 km', '3.25-3.5 km', '3.5-3.75 km', '3.75-4 km',\
                  '4-4.25 km', '4.25-4.5 km', '4.5-4.75 km', '4.75-5 km']

    #ax.hist(slopearr, 35, range=[0,35], normed=True, color=color, histtype='step')
    if plottype == 'finegrid':
        ax.plot(x_plot_arr, gamma_dist(x_plot_arr, *popt), ls='-', color=color, lw=2.0, label=label_list[callcount])
    elif plottype == None:
        # get label
        if callcount == 0:
            label = '1 to 2 km'
        elif callcount == 1:
            label = '2 to 3 km'
        elif callcount == 2:
            label = '3 to 4 km'
        elif callcount == 3:
            label = '4 to 5 km'
        elif callcount == 4:
            label = '5 to 6 km'
        elif callcount == 5:
            label = '6 to 7 km'
        elif callcount == 6:
            label = '7 to 8 km'
        elif callcount == 7:
            label = '8 to 9 km'
        elif callcount == 8:
            label = '9 to 10 km'
        elif callcount == 9:
            label = '10 to 15 km'
        elif callcount == 10:
            label = '15 to 20 km'
        elif callcount == 11:
            label = '20 to 25 km'
        elif callcount == 12:
            label = '25 to 30 km'
        elif callcount == 13:
            label = '30 to 35 km'
        ax.plot(x_plot_arr, gamma_dist(x_plot_arr, *popt), ls='-', color=color, lw=2.0, label=label)

    ax.legend(loc=0, prop={'size':8})

    callcount += 1

    return fig, ax, callcount

def call_hist_and_fits():

    # read in all arrays 
    density_diambin_1_2, density_diambin_2_3, density_diambin_3_4, density_diambin_4_5, \
    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9, \
    density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, \
    density_diambin_25_30, density_diambin_30_35, slope_diambin_1_2, slope_diambin_2_3, \
    slope_diambin_3_4, slope_diambin_4_5, slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, \
    slope_diambin_8_9, slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20, \
    slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35 = avg.read_no_overlap_arrays()

    density_diambin_1, density_diambin_2, density_diambin_3, density_diambin_4, \
    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8, \
    density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, \
    density_diambin_25, density_diambin_30, slope_diambin_1, slope_diambin_2, \
    slope_diambin_3, slope_diambin_4, slope_diambin_5, slope_diambin_6, slope_diambin_7, \
    slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15, slope_diambin_20, \
    slope_diambin_25, slope_diambin_30 = avg.read_Nvalue_arrays()

    # plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel('Slope', fontsize=12)
    ax.set_ylabel('Normalized bin counts', fontsize=12)

    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#543005','#fdbf6f',\
    '#ff7f00','#cab2d6','#bf812d','#6a3d9a','#b15928','#01665e']

    # checked the histogram plotting 
    # it seems like fitting a function and just plotting that might be more presentable.
    # Also check is histtype='step' might help better than the default histtype.
    # The histogram is normalized to make sure that hte area under the curve equals 1
    callcount = 0
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_1_2, colors[0], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_2_3, colors[1], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_3_4, colors[2], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_4_5, colors[3], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_5_6, colors[4], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_6_7, colors[5], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_7_8, colors[6], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_8_9, colors[7], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_9_10, colors[8], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_10_15, colors[9], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_15_20, colors[10], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_20_25, colors[11], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_25_30, colors[12], callcount)
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_30_35, colors[13], callcount)

    # minor ticks
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    fig.savefig(slope_extdir + 'slope_histogram_fits_and_hist.png', dpi=300, bbox_inches='tight')

    return None

def call_hist_and_fits_smallgrid():

    # read in all arrays 
    density_diambin_1_1p25, density_diambin_1p25_1p5, density_diambin_1p5_1p75, density_diambin_1p75_2, \
    density_diambin_2_2p25, density_diambin_2p25_2p5, density_diambin_2p5_2p75, density_diambin_2p75_3, \
    density_diambin_3_3p25, density_diambin_3p25_3p5, density_diambin_3p5_3p75, density_diambin_3p75_4, \
    density_diambin_4_4p25, density_diambin_4p25_4p5, density_diambin_4p5_4p75, density_diambin_4p75_5, \
    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9, \
    density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, \
    density_diambin_25_30, density_diambin_30_35, \
    slope_diambin_1_1p25, slope_diambin_1p25_1p5, slope_diambin_1p5_1p75, slope_diambin_1p75_2, \
    slope_diambin_2_2p25, slope_diambin_2p25_2p5, slope_diambin_2p5_2p75, slope_diambin_2p75_3, \
    slope_diambin_3_3p25, slope_diambin_3p25_3p5, slope_diambin_3p5_3p75, slope_diambin_3p75_4, \
    slope_diambin_4_4p25, slope_diambin_4p25_4p5, slope_diambin_4p5_4p75, slope_diambin_4p75_5, \
    slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, \
    slope_diambin_8_9, slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20, \
    slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35 = avg.read_no_overlap_arrays()

    density_diambin_1, density_diambin_2, density_diambin_3, density_diambin_4, \
    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8, \
    density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, \
    density_diambin_25, density_diambin_30, slope_diambin_1, slope_diambin_2, \
    slope_diambin_3, slope_diambin_4, slope_diambin_5, slope_diambin_6, slope_diambin_7, \
    slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15, slope_diambin_20, \
    slope_diambin_25, slope_diambin_30 = avg.read_Nvalue_arrays()

    # plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel('Slope', fontsize=12)
    ax.set_ylabel('Normalized bin counts', fontsize=12)

    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00',\
    '#cab2d6','#6a3d9a','#ffff99','#b15928','#878787','#c51b7d','#35978f','#b2182b']

    # checked the histogram plotting
    # it seems like fitting a function and just plotting that might be more presentable.
    # Also check is histtype='step' might help better than the default histtype.
    # The histogram is normalized to make sure that hte area under the curve equals 1
    callcount = 0
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_1_1p25, colors[0], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_1p25_1p5, colors[1], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_1p5_1p75, colors[2], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_1p75_2, colors[3], callcount, 'finegrid')

    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_2_2p25, colors[4], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_2p25_2p5, colors[5], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_2p5_2p75, colors[6], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_2p75_3, colors[7], callcount, 'finegrid')

    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_3_3p25, colors[8], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_3p25_3p5, colors[9], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_3p5_3p75, colors[10], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_3p75_4, colors[11], callcount, 'finegrid')

    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_4_4p25, colors[12], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_4p25_4p5, colors[13], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_4p5_4p75, colors[14], callcount, 'finegrid')
    fig, ax, callcount = hist_and_fit(fig, ax, slope_diambin_4p75_5, colors[15], callcount, 'finegrid')

    # minor ticks
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    fig.savefig(slope_extdir + 'slope_histogram_fits_smallrange_finegrid.png', dpi=300, bbox_inches='tight')

    return None

def get_crater_diams(crater_diam_m_arr, crater_ids, crater_ids_arr, total_craters):

    crater_diams = np.zeros(total_craters)

    for i in range(total_craters):
        current_id = crater_ids[i]
        current_id_idx = np.where(crater_ids_arr == current_id)[0][0]

        crater_diams[i] = crater_diam_m_arr[current_id_idx]

    return crater_diams

if __name__ == '__main__':

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    #call_hist_and_fits()
    call_hist_and_fits_smallgrid()
    sys.exit(0)

    ### Make diameter histograms for all slopes < 5 and all slopes >= 5 ###
    # Read in arrays 
    crater_ids_arr = np.load(slope_extdir + 'crater_ids_arr.npy')
    crater_x_arr = np.load(slope_extdir + 'crater_x_arr.npy')
    crater_y_arr = np.load(slope_extdir + 'crater_y_arr.npy')
    crater_diam_m_arr = np.load(slope_extdir + 'crater_diam_m_arr.npy')

    # get unique ids
    crater_ids = np.unique(crater_ids_arr)
    total_craters = len(crater_ids)

    # create and initialize crater diams corresponding to unique craters
    #crater_diams = cython_util_funcs.get_crater_diams(crater_diam_m_arr, crater_ids, crater_ids_arr, total_craters)
    crater_diams = get_crater_diams(crater_diam_m_arr, crater_ids, crater_ids_arr, total_craters)

    # Read in slopes
    slopemap_path = slope_extdir + 'hf_full_slopemap_clipped.txt'
    su.raster_to_numpy(slopemap_path)
    slope_arr = np.load(slopemap_path.replace('.txt', '.npy'))
    slope_arr = slope_arr.ravel()

    # read in crater ids associated with each pixel
    with open(slope_extdir + 'pix_crater_id_fastcomp.pkl', 'rb') as crater_id_file:
        crater_id_in_pix_arr = cPickle.load(crater_id_file)

    # now loop over all pixels and 
    # save the diams to two separate arrays 
    # depending on the slope of the pixel

    # first create the diam arrays
    total_pixels = len(crater_id_in_pix_arr)
    diam_lowslopes = []
    diam_highslopes = []

    slope_thresh = 10.0  # define the slope threshold dividing the low and high slopes

    for i in range(total_pixels):

        if (i % 100000) == 0.0:
            print '\r',
            print "At pixel number:",'{0:.2e}'.format(i),\
            "; time taken up to now:",'{0:.2f}'.format((time.time() - start)/60),"minutes.",
            sys.stdout.flush()

        current_crater_ids = crater_id_in_pix_arr[i]
        current_slope = slope_arr[i]

        total_current_craters = len(current_crater_ids)

        if total_current_craters == 0:
            continue

        elif total_current_craters == 1:

            current_diam_idx = np.where(crater_ids == current_crater_ids[0])[0]
            current_diam = float(crater_diams[current_diam_idx] / 1e3)  # converting to km

            if current_slope < slope_thresh:
                diam_lowslopes.append(current_diam)
            elif current_slope >= slope_thresh:
                diam_highslopes.append(current_diam)

        elif total_current_craters > 1:

            for j in range(total_current_craters):

                current_diam_idx = np.where(crater_ids == current_crater_ids[0])[0]
                current_diam = float(crater_diams[current_diam_idx] / 1e3)  # converting to km

                if current_slope < slope_thresh:
                    diam_lowslopes.append(current_diam)
                elif current_slope >= slope_thresh:
                    diam_highslopes.append(current_diam)

    # convert diameter lists in low and high slope to numpy arrays
    diam_lowslopes = np.asarray(diam_lowslopes).ravel()
    diam_highslopes = np.asarray(diam_highslopes).ravel()

    ### ----------------------------------------- Plotting ----------------------------------------- ###
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel('Crater Diameter [km]', fontsize=14)

    # define nicer colors
    myblue = mh.rgb_to_hex(0, 100, 180)
    myred = mh.rgb_to_hex(214, 39, 40)  # tableau 20 red

    ax.hist(diam_highslopes, 70, range=[0,35], color=myred, alpha=0.7, zorder=10)
    ax.hist(diam_lowslopes, 70, range=[0,35], color=myblue, alpha=0.7, zorder=11)

    ax.set_xticklabels(ax.get_xticks().tolist(), size=12)
    ax.set_yticklabels(ax.get_yticks().tolist(), size=12)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True, alpha=0.4)

    ax.text(0.55, 0.86, r'$\mathrm{Diameters\ at\ Slope<}$' + str(int(slope_thresh)) + r'$^\circ$',\
    verticalalignment='top', horizontalalignment='left', \
    transform=ax.transAxes, color=myblue, size=14)
    ax.text(0.55, 0.8, r'$\mathrm{Diameters\ at\ Slope\geq}$' + str(int(slope_thresh)) + r'$^\circ$',\
    verticalalignment='top', horizontalalignment='left', \
    transform=ax.transAxes, color=myred, size=14)

    fig.savefig(slope_extdir + 'diams_low_high_' + str(int(slope_thresh)) + '_slopes_hist.png', dpi=300, bbox_inches='tight')

    # total run time
    print "Total time taken --", (time.time() - start)/60, "minutes."
    sys.exit(0)