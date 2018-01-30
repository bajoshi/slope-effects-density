from __future__ import division

import numpy as np
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from scipy.misc import factorial
from scipy.special import gamma

import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

sys.path.append(slopedir)
import average_density_plots as avg

def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)

def gamma_dist(x, k, theta):
    return (x**(k-1) * np.exp(-x/theta)) / (theta**k * gamma(k))

def hist_and_fit(fig, ax, slopearr, color):

    counts, bins = np.histogram(slopearr, 35, range=[0,35], density=True)
    fitting_bins = [(bins[i] + bins[i+1])/2 for i in bins[:-1]]
    fitting_bins = np.asarray(fitting_bins)

    x_plot_arr = np.linspace(0,35,1000)

    # fitting a gamma distribution
    popt, pcov = curve_fit(gamma_dist, fitting_bins, counts, p0=[3.0, 3.0])

    # get label
    if color == 'midnightblue':
        label = '1 to 2 km'
    elif color == 'blue':
        label = '2 to 3 km'
    elif color == 'royalblue':
        label = '3 to 4 km'
    elif color == 'dodgerblue':
        label = '4 to 5 km'
    elif color == 'deepskyblue':
        label = '5 to 6 km'
    elif color == 'steelblue':
        label = '6 to 7 km'
    elif color == 'slateblue':
        label = '7 to 8 km'
    elif color == 'rebeccapurple':
        label = '8 to 9 km'
    elif color == 'darkcyan':
        label = '9 to 10 km'
    elif color == 'green':
        label = '10 to 15 km'
    elif color == 'olive':
        label = '15 to 20 km'
    elif color == 'goldenrod':
        label = '20 to 25 km'
    elif color == 'darkorchid':
        label = '25 to 30 km'
    elif color == 'maroon':
        label = '30 to 35 km'

    ax.hist(slopearr, 35, range=[0,35], normed=True, color=color, histtype='step')
    ax.plot(x_plot_arr, gamma_dist(x_plot_arr, *popt), ls='-', color=color, lw=2.0, label=label)

    ax.legend(loc=0, prop={'size':8})

    return fig, ax

if __name__ == '__main__':
    
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

    # checked the histogram plotting 
    # it seems like fitting a function and just plotting that might be more presentable.
    # Also check is histtype='step' might help better than the default histtype.
    # The histogram is normalized to make sure that hte area under the curve equals 1
    fig, ax = hist_and_fit(fig, ax, slope_diambin_1_2, 'midnightblue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_2_3, 'blue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_3_4, 'royalblue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_4_5, 'dodgerblue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_5_6, 'deepskyblue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_6_7, 'steelblue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_7_8, 'slateblue')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_8_9, 'rebeccapurple')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_9_10, 'darkcyan')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_10_15, 'green')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_15_20, 'olive')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_20_25, 'goldenrod')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_25_30, 'darkorchid')
    fig, ax = hist_and_fit(fig, ax, slope_diambin_30_35, 'maroon')

    # minor ticks
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    fig.savefig(slope_extdir + 'slope_histogram_fits_and_hist.png', dpi=150, bbox_inches='tight')

    sys.exit(0)