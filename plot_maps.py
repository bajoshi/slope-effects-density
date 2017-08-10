from __future__ import division  

import numpy as np
from astropy.modeling import models, fitting

import sys
import os

import matplotlib.pyplot as plt

# Checks if os is windows or unix
if os.name == 'posix':
    home = os.getenv('HOME')  # does not have a trailing slash
    desktop = home + '/Desktop/'
    slopedir = desktop + '/slope-effects-density/'
    slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'
elif os.name == 'nt':
    #desktop = 'C:\Users\Heather\Desktop\\'
    #slopedir = desktop + '\\slope-effects-density\\'
    slopedir = 'E:\slope-effects-density\\'

sys.path.append(slopedir)
import slope_effects_density as se

def plot_hist(arr, xlow, xhigh):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # freedman - diaconis rule
    #iqr = np.std(pix_frac[high_idx], dtype=np.float64)
    #binsize = 2*iqr*np.power(len(pix_frac[high_idx]),-1/3)
    #totalbins = np.floor((max(pix_frac[high_idx]) - min(pix_frac[high_idx]))/binsize)
    #print binsize, totalbins

    if any(np.isnan(arr)):
        arr_nonan = []
        for i in range(len(arr)):
            if arr[i] is not np.nan:
                if arr[i] > 10:
                    continue
                else:
                    arr_nonan.append(arr[i])
        arr_nonan = np.asarray(arr_nonan)

        ax.hist(arr_nonan, 1000, alpha=0.5)
    else:
        ax.hist(arr, 1000, alpha=0.5)

    ax.set_xlim(xlow, xhigh)
    ax.set_yscale('log')
    plt.show()

    return None

if __name__ == '__main__':
    
    # load all arrays
    # read in products
    pix_frac = np.load(slope_extdir + 'pix_area_fraction.npy')
    crater_frac = np.load(slope_extdir + 'crater_area_frac_in_pix.npy')

    density = se.get_density(crater_frac, pix_frac, len(pix_frac))

    # read in pixel coordinates
    slope_arr = np.load(slope_extdir + '3km_slope_points.npy')

    pix_x_cen_arr = slope_arr['pix_x_cen']
    pix_y_cen_arr = slope_arr['pix_y_cen']
    slope = slope_arr['slope_val']
    # These are all 1D arrays of 4709560 elements

    # reshape to 2D arrays
    rows, columns = se.get_rows_columns(pix_x_cen_arr, pix_y_cen_arr)

    pix_frac_2d = pix_frac.reshape(rows, columns)
    crater_frac_2d = crater_frac.reshape(rows, columns)
    density_2d = density.reshape(rows, columns)
    slope_2d = slope.reshape(rows, columns)

    # choose valid points to plot
    # first, replace all -9999.0 values by NaN
    nodata_idx = np.where(density == -9999.0)[0]
    density[nodata_idx] = np.nan

    # second, for now, also ignore high density values
    density_upper_lim = 5.0
    high_idx = np.where(density > density_upper_lim)[0]
    density[high_idx] = np.nan

    # third, only select non-zero values of density
    #val_idx = np.where(density != 0.0)[0] # np.isfinite(density)

    # find frac pix and crater values for pixels that have high density values
    #idx = np.where(density > density_upper_lim)[0]
    ##print idx
    #print "Upper density limit is", density_upper_lim
    #print len(np.where(density > density_upper_lim)[0]), \
    #"pixels with density value greater than density upper limit."
    #print "Frac crater counts at high density value are \n", crater_frac[idx]
    #print "Frac pixel  areas  at high density value are \n", pix_frac[idx]
    #sys.exit(0)

    eps = 0.01
    idx = np.where((density >= 1.0 - eps) & ((density <= 1.0 + eps)))[0]
    print idx
    print crater_frac[idx]
    print crater_diam[idx]
    print pix_frac[idx]
    sys.exit(0)

    # ---------------------------------------- #
    #plot_hist(density)

    # --------------------------------------- # 
    # Fit gaussian
    gauss_init = models.Gaussian1D(amplitude=1.0, mean=5.0, stddev=10.0)
    fit_gauss = fitting.LevMarLSQFitter()
    #x_fit_arr = slope[val_idx]
    #y_fit_arr = density[val_idx]
    #gauss = fit_gauss(gauss_init, x_fit_arr, y_fit_arr)

    # plots
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.plot(pix_frac[val_idx], np.log10(density[val_idx]), 'o', color='k', markersize=2, markeredgecolor='None')
    #ax.set_ylim(-5,2.5)

    ax.plot(slope, density, 'o', color='k', markersize=2, markeredgecolor='None')
    #ax.plot(x_fit_arr, gauss(x_fit_arr), '-', color='r', lw=2)
    #ax.set_ylim(0,1)

    plt.show()

    sys.exit(0)