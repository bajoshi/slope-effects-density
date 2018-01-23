from __future__ import division

import numpy as np
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit

import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

def read_no_overlap_arrays():

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

    return density_diambin_1_2, density_diambin_2_3, density_diambin_3_4, density_diambin_4_5, \
    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9, \
    density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, \
    density_diambin_25_30, density_diambin_30_35, slope_diambin_1_2, slope_diambin_2_3, \
    slope_diambin_3_4, slope_diambin_4_5, slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, \
    slope_diambin_8_9, slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20, \
    slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35

def read_Nvalue_arrays():

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

    return density_diambin_1, density_diambin_2, density_diambin_3, density_diambin_4, \
    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8, \
    density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, \
    density_diambin_25, density_diambin_30, slope_diambin_1, slope_diambin_2, \
    slope_diambin_3, slope_diambin_4, slope_diambin_5, slope_diambin_6, slope_diambin_7, \
    slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15, slope_diambin_20, \
    slope_diambin_25, slope_diambin_30

def get_avg_finite_elements(density_arr, slope_arr):

    density_finidx = np.where(np.isfinite(density_arr))[0]
    slope_finidx = np.where(np.isfinite(slope_arr))[0]
    fin_idx = np.intersect1d(slope_finidx, density_finidx)

    density_avg = np.mean(density_arr[fin_idx])
    slope_avg = np.mean(slope_arr[fin_idx])

    return density_avg, slope_avg

def exp_func(x, amp, alpha, x0):
    return amp * np.power(x,-1*alpha) * np.exp(-1*x/x0)

if __name__ == '__main__':

    # read in all arrays 
    density_diambin_1_2, density_diambin_2_3, density_diambin_3_4, density_diambin_4_5, \
    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9, \
    density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, \
    density_diambin_25_30, density_diambin_30_35, slope_diambin_1_2, slope_diambin_2_3, \
    slope_diambin_3_4, slope_diambin_4_5, slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, \
    slope_diambin_8_9, slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20, \
    slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35 = read_no_overlap_arrays()

    density_diambin_1, density_diambin_2, density_diambin_3, density_diambin_4, \
    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8, \
    density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, \
    density_diambin_25, density_diambin_30, slope_diambin_1, slope_diambin_2, \
    slope_diambin_3, slope_diambin_4, slope_diambin_5, slope_diambin_6, slope_diambin_7, \
    slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15, slope_diambin_20, \
    slope_diambin_25, slope_diambin_30 = read_Nvalue_arrays()

    # get averages for these arrays
    # although there aren't any nans in the density arrays
    # there are some in the slope arrays.
    # no overlap
    density_diambin_1_2_avg, slope_diambin_1_2_avg = get_avg_finite_elements(density_diambin_1_2, slope_diambin_1_2)
    density_diambin_2_3_avg, slope_diambin_2_3_avg = get_avg_finite_elements(density_diambin_2_3, slope_diambin_2_3)
    density_diambin_3_4_avg, slope_diambin_3_4_avg = get_avg_finite_elements(density_diambin_3_4, slope_diambin_3_4)
    density_diambin_4_5_avg, slope_diambin_4_5_avg = get_avg_finite_elements(density_diambin_4_5, slope_diambin_4_5)
    density_diambin_5_6_avg, slope_diambin_5_6_avg = get_avg_finite_elements(density_diambin_5_6, slope_diambin_5_6)
    density_diambin_6_7_avg, slope_diambin_6_7_avg = get_avg_finite_elements(density_diambin_6_7, slope_diambin_6_7)
    density_diambin_7_8_avg, slope_diambin_7_8_avg = get_avg_finite_elements(density_diambin_7_8, slope_diambin_7_8)
    density_diambin_8_9_avg, slope_diambin_8_9_avg = get_avg_finite_elements(density_diambin_8_9, slope_diambin_8_9)
    density_diambin_9_10_avg, slope_diambin_9_10_avg = get_avg_finite_elements(density_diambin_9_10, slope_diambin_9_10)
    density_diambin_10_15_avg, slope_diambin_10_15_avg = get_avg_finite_elements(density_diambin_10_15, slope_diambin_10_15)
    density_diambin_15_20_avg, slope_diambin_15_20_avg = get_avg_finite_elements(density_diambin_15_20, slope_diambin_15_20)
    density_diambin_20_25_avg, slope_diambin_20_25_avg = get_avg_finite_elements(density_diambin_20_25, slope_diambin_20_25)
    density_diambin_25_30_avg, slope_diambin_25_30_avg = get_avg_finite_elements(density_diambin_25_30, slope_diambin_25_30)
    density_diambin_30_35_avg, slope_diambin_30_35_avg = get_avg_finite_elements(density_diambin_30_35, slope_diambin_30_35)

    # nvalue
    density_diambin_1_avg, slope_diambin_1_avg = get_avg_finite_elements(density_diambin_1, slope_diambin_1)
    density_diambin_2_avg, slope_diambin_2_avg = get_avg_finite_elements(density_diambin_2, slope_diambin_2)
    density_diambin_3_avg, slope_diambin_3_avg = get_avg_finite_elements(density_diambin_3, slope_diambin_3)
    density_diambin_4_avg, slope_diambin_4_avg = get_avg_finite_elements(density_diambin_4, slope_diambin_4)
    density_diambin_5_avg, slope_diambin_5_avg = get_avg_finite_elements(density_diambin_5, slope_diambin_5)
    density_diambin_6_avg, slope_diambin_6_avg = get_avg_finite_elements(density_diambin_6, slope_diambin_6)
    density_diambin_7_avg, slope_diambin_7_avg = get_avg_finite_elements(density_diambin_7, slope_diambin_7)
    density_diambin_8_avg, slope_diambin_8_avg = get_avg_finite_elements(density_diambin_8, slope_diambin_8)
    density_diambin_9_avg, slope_diambin_9_avg = get_avg_finite_elements(density_diambin_9, slope_diambin_9)
    density_diambin_10_avg, slope_diambin_10_avg = get_avg_finite_elements(density_diambin_10, slope_diambin_10)
    density_diambin_15_avg, slope_diambin_15_avg = get_avg_finite_elements(density_diambin_15, slope_diambin_15)
    density_diambin_20_avg, slope_diambin_20_avg = get_avg_finite_elements(density_diambin_20, slope_diambin_20)
    density_diambin_25_avg, slope_diambin_25_avg = get_avg_finite_elements(density_diambin_25, slope_diambin_25)
    density_diambin_30_avg, slope_diambin_30_avg = get_avg_finite_elements(density_diambin_30, slope_diambin_30)

    all_density_averages_nooverlap = np.array([density_diambin_1_2_avg, density_diambin_2_3_avg, density_diambin_3_4_avg, \
    density_diambin_4_5_avg, density_diambin_5_6_avg, density_diambin_6_7_avg, density_diambin_7_8_avg, \
    density_diambin_8_9_avg, density_diambin_9_10_avg, density_diambin_10_15_avg, \
    density_diambin_15_20_avg, density_diambin_20_25_avg, density_diambin_25_30_avg, density_diambin_30_35_avg])

    all_slope_averages_nooverlap = np.array([slope_diambin_1_2_avg, slope_diambin_2_3_avg, slope_diambin_3_4_avg, \
    slope_diambin_4_5_avg, slope_diambin_5_6_avg, slope_diambin_6_7_avg, slope_diambin_7_8_avg, \
    slope_diambin_8_9_avg, slope_diambin_9_10_avg, slope_diambin_10_15_avg, \
    slope_diambin_15_20_avg, slope_diambin_20_25_avg, slope_diambin_25_30_avg, slope_diambin_30_35_avg])

    all_density_averages_nvalue = np.array([density_diambin_1_avg, density_diambin_2_avg, density_diambin_3_avg, \
    density_diambin_4_avg, density_diambin_5_avg, density_diambin_6_avg, density_diambin_7_avg, \
    density_diambin_8_avg, density_diambin_9_avg, density_diambin_10_avg, density_diambin_15_avg, \
    density_diambin_20_avg, density_diambin_25_avg, density_diambin_30_avg])

    all_slope_averages_nvalue = np.array([slope_diambin_1_avg, slope_diambin_2_avg, slope_diambin_3_avg, \
    slope_diambin_4_avg, slope_diambin_5_avg, slope_diambin_6_avg, slope_diambin_7_avg, \
    slope_diambin_8_avg, slope_diambin_9_avg, slope_diambin_10_avg, slope_diambin_15_avg, \
    slope_diambin_20_avg, slope_diambin_25_avg, slope_diambin_30_avg])

    # plot
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$')
    ax.set_ylabel(r'$\mathrm{log(Density)}$')

    # add minor ticks and grid
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True, alpha=0.5)

    ax.scatter(all_slope_averages_nooverlap, all_density_averages_nooverlap, s=10, color='k', edgecolors='none')
    init = models.PowerLaw1D(amplitude=1, x_0=1, alpha=1)
    fit = fitting.LevMarLSQFitter()
    f = fit(init, all_slope_averages_nooverlap[:11], all_density_averages_nooverlap[:11])
    print f
    x_plot_arr = np.linspace(3,18,1000)
    ax.plot(x_plot_arr, f(x_plot_arr), ls='-', color='skyblue', lw=2)

    ax.set_xlim(3,18)
    ax.set_ylim(-0.01,0.2)

    fig.savefig(slope_extdir + 'nooverlap_averages_plot.png', dpi=150, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    # _----------------------------- Nvalue -------------------------- #
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$')
    ax.set_ylabel(r'$\mathrm{log(Density)}$')

    # add minor ticks and grid
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True, alpha=0.5)

    ax.scatter(all_slope_averages_nvalue, all_density_averages_nvalue, s=10, color='k', edgecolors='none')
    init = models.PowerLaw1D(amplitude=1, x_0=1, alpha=1)
    fit = fitting.LevMarLSQFitter()
    f = fit(init, all_slope_averages_nvalue[:9], all_density_averages_nvalue[:9])
    print f
    x_plot_arr = np.linspace(3,18,1000)
    ax.plot(x_plot_arr, f(x_plot_arr), ls='-', color='skyblue', lw=2)

    # Show residuals

    ax.set_xlim(9,14)
    ax.set_ylim(-0.01,0.1)

    fig.savefig(slope_extdir + 'nvalue_averages_plot.png', dpi=150, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()    

    # ----------------------------- Nvalue normalized -------------------------- #
    all_diam_values = np.array([1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 12.5, 17.5, 22.5, 27.5, 32.5])
    norm_values = np.power(10, 4.9) * np.power(all_diam_values, -2)
    print norm_values

    all_density_averages_nvalue /= norm_values

    print all_density_averages_nvalue

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$')
    ax.set_ylabel(r'$\mathrm{log(Density)}$')

    # add minor ticks and grid
    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True, alpha=0.5)

    ax.scatter(all_slope_averages_nvalue, all_density_averages_nvalue, s=10, color='k', edgecolors='none')
    init = models.PowerLaw1D(amplitude=1, x_0=1, alpha=1)
    fit = fitting.LevMarLSQFitter()
    f = fit(init, all_slope_averages_nvalue[:9], all_density_averages_nvalue[:9])
    print f
    x_plot_arr = np.linspace(3,18,1000)
    ax.plot(x_plot_arr, f(x_plot_arr), ls='-', color='skyblue', lw=2)

    # Show residuals

    ax.set_xlim(9,14)
    ax.set_ylim(0,2e-5)

    fig.savefig(slope_extdir + 'nvalue_norm_averages_plot.png', dpi=150, bbox_inches='tight')
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()  

    sys.exit(0)