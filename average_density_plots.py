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

    # Error on the mean
    #density_avg_error = np.std(density_arr[fin_idx]) / np.sqrt(len(density_arr[fin_idx]))
    #slope_avg_error = np.std(slope_arr[fin_idx]) / np.sqrt(len(slope_arr[fin_idx]))

    density_avg_error = np.std(density_arr[fin_idx])
    slope_avg_error = np.std(slope_arr[fin_idx])

    return density_avg, slope_avg, density_avg_error, slope_avg_error

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
    density_diambin_1_2_avg, slope_diambin_1_2_avg, density_diambin_1_2_avgerror, slope_diambin_1_2_avgerror \
    = get_avg_finite_elements(density_diambin_1_2, slope_diambin_1_2)
    density_diambin_2_3_avg, slope_diambin_2_3_avg, density_diambin_2_3_avgerror, slope_diambin_2_3_avgerror \
    = get_avg_finite_elements(density_diambin_2_3, slope_diambin_2_3)
    density_diambin_3_4_avg, slope_diambin_3_4_avg, density_diambin_3_4_avgerror, slope_diambin_3_4_avgerror \
    = get_avg_finite_elements(density_diambin_3_4, slope_diambin_3_4)
    density_diambin_4_5_avg, slope_diambin_4_5_avg, density_diambin_4_5_avgerror, slope_diambin_4_5_avgerror \
    = get_avg_finite_elements(density_diambin_4_5, slope_diambin_4_5)
    density_diambin_5_6_avg, slope_diambin_5_6_avg, density_diambin_5_6_avgerror, slope_diambin_5_6_avgerror \
    = get_avg_finite_elements(density_diambin_5_6, slope_diambin_5_6)
    density_diambin_6_7_avg, slope_diambin_6_7_avg, density_diambin_6_7_avgerror, slope_diambin_6_7_avgerror \
    = get_avg_finite_elements(density_diambin_6_7, slope_diambin_6_7)
    density_diambin_7_8_avg, slope_diambin_7_8_avg, density_diambin_7_8_avgerror, slope_diambin_7_8_avgerror \
    = get_avg_finite_elements(density_diambin_7_8, slope_diambin_7_8)
    density_diambin_8_9_avg, slope_diambin_8_9_avg, density_diambin_8_9_avgerror, slope_diambin_8_9_avgerror \
    = get_avg_finite_elements(density_diambin_8_9, slope_diambin_8_9)
    density_diambin_9_10_avg, slope_diambin_9_10_avg, density_diambin_9_10_avgerror, slope_diambin_9_10_avgerror \
    = get_avg_finite_elements(density_diambin_9_10, slope_diambin_9_10)
    density_diambin_10_15_avg, slope_diambin_10_15_avg, density_diambin_10_15_avgerror, slope_diambin_10_15_avgerror \
    = get_avg_finite_elements(density_diambin_10_15, slope_diambin_10_15)
    density_diambin_15_20_avg, slope_diambin_15_20_avg, density_diambin_15_20_avgerror, slope_diambin_15_20_avgerror \
    = get_avg_finite_elements(density_diambin_15_20, slope_diambin_15_20)
    density_diambin_20_25_avg, slope_diambin_20_25_avg, density_diambin_20_25_avgerror, slope_diambin_20_25_avgerror \
    = get_avg_finite_elements(density_diambin_20_25, slope_diambin_20_25)
    density_diambin_25_30_avg, slope_diambin_25_30_avg, density_diambin_25_30_avgerror, slope_diambin_25_30_avgerror \
    = get_avg_finite_elements(density_diambin_25_30, slope_diambin_25_30)
    density_diambin_30_35_avg, slope_diambin_30_35_avg, density_diambin_30_35_avgerror, slope_diambin_30_35_avgerror \
    = get_avg_finite_elements(density_diambin_30_35, slope_diambin_30_35)

    # nvalue
    density_diambin_1_avg, slope_diambin_1_avg, density_diambin_1_avgerror, slope_diambin_1_avgerror \
    = get_avg_finite_elements(density_diambin_1, slope_diambin_1)
    density_diambin_2_avg, slope_diambin_2_avg, density_diambin_2_avgerror, slope_diambin_2_avgerror \
    = get_avg_finite_elements(density_diambin_2, slope_diambin_2)
    density_diambin_3_avg, slope_diambin_3_avg, density_diambin_3_avgerror, slope_diambin_3_avgerror \
    = get_avg_finite_elements(density_diambin_3, slope_diambin_3)
    density_diambin_4_avg, slope_diambin_4_avg, density_diambin_4_avgerror, slope_diambin_4_avgerror \
    = get_avg_finite_elements(density_diambin_4, slope_diambin_4)
    density_diambin_5_avg, slope_diambin_5_avg, density_diambin_5_avgerror, slope_diambin_5_avgerror \
    = get_avg_finite_elements(density_diambin_5, slope_diambin_5)
    density_diambin_6_avg, slope_diambin_6_avg, density_diambin_6_avgerror, slope_diambin_6_avgerror \
    = get_avg_finite_elements(density_diambin_6, slope_diambin_6)
    density_diambin_7_avg, slope_diambin_7_avg, density_diambin_7_avgerror, slope_diambin_7_avgerror \
    = get_avg_finite_elements(density_diambin_7, slope_diambin_7)
    density_diambin_8_avg, slope_diambin_8_avg, density_diambin_8_avgerror, slope_diambin_8_avgerror \
    = get_avg_finite_elements(density_diambin_8, slope_diambin_8)
    density_diambin_9_avg, slope_diambin_9_avg, density_diambin_9_avgerror, slope_diambin_9_avgerror \
    = get_avg_finite_elements(density_diambin_9, slope_diambin_9)
    density_diambin_10_avg, slope_diambin_10_avg, density_diambin_10_avgerror, slope_diambin_10_avgerror \
    = get_avg_finite_elements(density_diambin_10, slope_diambin_10)
    density_diambin_15_avg, slope_diambin_15_avg, density_diambin_15_avgerror, slope_diambin_15_avgerror \
    = get_avg_finite_elements(density_diambin_15, slope_diambin_15)
    density_diambin_20_avg, slope_diambin_20_avg, density_diambin_20_avgerror, slope_diambin_20_avgerror \
    = get_avg_finite_elements(density_diambin_20, slope_diambin_20)
    density_diambin_25_avg, slope_diambin_25_avg, density_diambin_25_avgerror, slope_diambin_25_avgerror \
    = get_avg_finite_elements(density_diambin_25, slope_diambin_25)
    density_diambin_30_avg, slope_diambin_30_avg, density_diambin_30_avgerror, slope_diambin_30_avgerror \
    = get_avg_finite_elements(density_diambin_30, slope_diambin_30)

    ### Lump all avg value arrays together ###
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

    ### Lump all error arrays together ###
    all_density_avgerrors_nooverlap = np.array([density_diambin_1_2_avgerror, density_diambin_2_3_avgerror, density_diambin_3_4_avgerror, \
    density_diambin_4_5_avgerror, density_diambin_5_6_avgerror, density_diambin_6_7_avgerror, density_diambin_7_8_avgerror, \
    density_diambin_8_9_avgerror, density_diambin_9_10_avgerror, density_diambin_10_15_avgerror, \
    density_diambin_15_20_avgerror, density_diambin_20_25_avgerror, density_diambin_25_30_avgerror, density_diambin_30_35_avgerror])

    all_slope_avgerrors_nooverlap = np.array([slope_diambin_1_2_avgerror, slope_diambin_2_3_avgerror, slope_diambin_3_4_avgerror, \
    slope_diambin_4_5_avgerror, slope_diambin_5_6_avgerror, slope_diambin_6_7_avgerror, slope_diambin_7_8_avgerror, \
    slope_diambin_8_9_avgerror, slope_diambin_9_10_avgerror, slope_diambin_10_15_avgerror, \
    slope_diambin_15_20_avgerror, slope_diambin_20_25_avgerror, slope_diambin_25_30_avgerror, slope_diambin_30_35_avgerror])

    all_density_avgerrors_nvalue = np.array([density_diambin_1_avgerror, density_diambin_2_avgerror, density_diambin_3_avgerror, \
    density_diambin_4_avgerror, density_diambin_5_avgerror, density_diambin_6_avgerror, density_diambin_7_avgerror, \
    density_diambin_8_avgerror, density_diambin_9_avgerror, density_diambin_10_avgerror, density_diambin_15_avgerror, \
    density_diambin_20_avgerror, density_diambin_25_avgerror, density_diambin_30_avgerror])

    all_slope_avgerrors_nvalue = np.array([slope_diambin_1_avgerror, slope_diambin_2_avgerror, slope_diambin_3_avgerror, \
    slope_diambin_4_avgerror, slope_diambin_5_avgerror, slope_diambin_6_avgerror, slope_diambin_7_avgerror, \
    slope_diambin_8_avgerror, slope_diambin_9_avgerror, slope_diambin_10_avgerror, slope_diambin_15_avgerror, \
    slope_diambin_20_avgerror, slope_diambin_25_avgerror, slope_diambin_30_avgerror])

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

    colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#543005','#fdbf6f',\
    '#ff7f00','#cab2d6','#bf812d','#6a3d9a','#b15928','#01665e']

    ax.scatter(all_slope_averages_nooverlap[0], all_density_averages_nooverlap[0], \
        s=50, marker='o', color=colors[0], label='1-2' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[1], all_density_averages_nooverlap[1], \
        s=50, marker='o', color=colors[1], label='2-3' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[2], all_density_averages_nooverlap[2], \
        s=50, marker='o', color=colors[2], label='3-4' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[3], all_density_averages_nooverlap[3], \
        s=50, marker='o', color=colors[3], label='4-5' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[4], all_density_averages_nooverlap[4], \
        s=50, marker='o', color=colors[4], label='5-6' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[5], all_density_averages_nooverlap[5], \
        s=50, marker='o', color=colors[5], label='6-7' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[6], all_density_averages_nooverlap[6], \
        s=50, marker='o', color=colors[6], label='7-8' + ' km', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nooverlap[7], all_density_averages_nooverlap[7], \
        s=50, marker='o', color=colors[7], label='8-9' + ' km', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nooverlap[8], all_density_averages_nooverlap[8], \
        s=50, marker='o', color=colors[8], label='9-10' + ' km', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nooverlap[9], all_density_averages_nooverlap[9], \
        s=50, marker='o', color=colors[9], label='10-15' + ' km', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nooverlap[10], all_density_averages_nooverlap[10], 
        s=50, marker='o', color=colors[10], label='15-20' + ' km', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nooverlap[11], all_density_averages_nooverlap[11], 
        s=20, marker='x', color=colors[11], label='20-25' + ' km', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nooverlap[12], all_density_averages_nooverlap[12], 
        s=20, marker='x', color=colors[12], label='25-30' + ' km', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nooverlap[13], all_density_averages_nooverlap[13], 
        s=20, marker='x', color=colors[13], label='30-35' + ' km', \
        edgecolors='none', zorder=4)

    # fitting
    init = models.PowerLaw1D(amplitude=1, x_0=1, alpha=1)
    fit = fitting.LevMarLSQFitter()
    f = fit(init, all_slope_averages_nooverlap[:11], all_density_averages_nooverlap[:11])
    print f
    x_plot_arr = np.linspace(3,18,1000)
    ax.plot(x_plot_arr, f(x_plot_arr), ls='-', color='skyblue', lw=2)

    # Find chi2 and put the value on the plot
    chi2 = np.sum(((all_density_averages_nooverlap[:11] - f(all_slope_averages_nooverlap[:11])) / all_density_avgerrors_nooverlap[:11])**2)
    
    # text on plot
    # equation 
    ax.text(0.4, 0.86, r'$\mathrm{f(x) = A\left(\frac{x}{x_0}\right)^{-\alpha}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)
    # best fit parameters
    ax.text(0.315, 0.78, r'$\mathrm{Amplitude =\ }$' + str("{:.3}".format(f.parameters[0])), \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)
    ax.text(0.417, 0.72, r'$\mathrm{x_0 =\ }$' +  str("{:.3}".format(f.parameters[1])), \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)
    ax.text(0.428, 0.66, r'$\mathrm{\alpha =\ }$' +  str("{:.3}".format(f.parameters[2])), \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)

    # chi2
    ax.text(0.416, 0.6, r'$\mathrm{\chi^2 =\ }$' + str("{:.3}".format(chi2)), verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)

    ax.set_xlim(3,18)
    ax.set_ylim(-0.01,0.2)

    ax.legend(loc=0)

    fig.savefig(slope_extdir + 'nooverlap_averages_plot.png', dpi=300, bbox_inches='tight')
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

    ax.scatter(all_slope_averages_nvalue[0], all_density_averages_nvalue[0], \
        s=50, marker='o', color=colors[0], label='N(1)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[1], all_density_averages_nvalue[1], \
        s=50, marker='o', color=colors[1], label='N(2)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[2], all_density_averages_nvalue[2], \
        s=50, marker='o', color=colors[2], label='N(3)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[3], all_density_averages_nvalue[3], \
        s=50, marker='o', color=colors[3], label='N(4)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[4], all_density_averages_nvalue[4], \
        s=50, marker='o', color=colors[4], label='N(5)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[5], all_density_averages_nvalue[5], \
        s=50, marker='o', color=colors[5], label='N(6)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[6], all_density_averages_nvalue[6], \
        s=50, marker='o', color=colors[6], label='N(7)', \
        edgecolors='none', zorder=5)
    ax.scatter(all_slope_averages_nvalue[7], all_density_averages_nvalue[7], \
        s=50, marker='o', color=colors[7], label='N(8)', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nvalue[8], all_density_averages_nvalue[8], \
        s=50, marker='o', color=colors[8], label='N(9)', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nvalue[9], all_density_averages_nvalue[9], \
        s=20, marker='x', color=colors[9], label='N(10)', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nvalue[10], all_density_averages_nvalue[10], 
        s=20, marker='x', color=colors[10], label='N(15)', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nvalue[11], all_density_averages_nvalue[11], 
        s=20, marker='x', color=colors[11], label='N(20)', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nvalue[12], all_density_averages_nvalue[12], 
        s=20, marker='x', color=colors[12], label='N(25)', \
        edgecolors='none', zorder=4)
    ax.scatter(all_slope_averages_nvalue[13], all_density_averages_nvalue[13], 
        s=20, marker='x', color=colors[13], label='N(30)', \
        edgecolors='none', zorder=4)

    # fitting
    init = models.PowerLaw1D(amplitude=1, x_0=1, alpha=1)
    fit = fitting.LevMarLSQFitter()
    f = fit(init, all_slope_averages_nvalue[:9], all_density_averages_nvalue[:9])
    print f
    x_plot_arr = np.linspace(3,18,1000)
    ax.plot(x_plot_arr, f(x_plot_arr), ls='-', color='skyblue', lw=2)

    # Find chi2 and put the value on the plot
    chi2 = np.sum(((all_density_averages_nvalue[:9] - f(all_slope_averages_nvalue[:9])) / all_density_avgerrors_nvalue[:9])**2)
    
    # text on plot
    # equation 
    ax.text(0.4, 0.86, r'$\mathrm{f(x) = A\left(\frac{x}{x_0}\right)^{-\alpha}}$', \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)
    # best fit parameters
    ax.text(0.315, 0.78, r'$\mathrm{Amplitude =\ }$' + str("{:.3}".format(f.parameters[0])), \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)
    ax.text(0.417, 0.72, r'$\mathrm{x_0 =\ }$' +  str("{:.3}".format(f.parameters[1])), \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)
    ax.text(0.428, 0.66, r'$\mathrm{\alpha =\ }$' +  str("{:.3}".format(f.parameters[2])), \
        verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)

    # chi2
    ax.text(0.416, 0.6, r'$\mathrm{\chi^2 =\ }$' + str("{:.3}".format(chi2)), verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)

    # Show residuals

    ax.set_xlim(9,14)
    ax.set_ylim(-0.01,0.1)

    ax.legend(loc=0)

    fig.savefig(slope_extdir + 'nvalue_averages_plot.png', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()    

    sys.exit(0)

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

    fig.savefig(slope_extdir + 'nvalue_norm_averages_plot.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.clf()
    plt.cla()
    plt.close()  

    sys.exit(0)