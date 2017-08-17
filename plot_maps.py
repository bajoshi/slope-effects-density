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
import slope_utils as su

def plot_hist(arr, xlow, xhigh):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if any(np.isnan(arr)):
        arr_nonan = []
        idx = np.isnan(arr)
        for i in range(len(arr)):
            if not idx[i]:
                arr_nonan.append(arr[i])
        arr_nonan = np.asarray(arr_nonan)

        arr = arr_nonan

    print "Using freedman - diaconis rule"
    iqr = np.std(arr, dtype=np.float64)
    binsize = 2*iqr*np.power(len(arr),-1/3)
    totalbins = np.floor((np.max(arr) - np.min(arr))/binsize)
    print "IQR, binsize and totalbins:", iqr, binsize, totalbins

    ax.hist(arr, totalbins, alpha=0.5, edgecolor='none')

    ax.set_xlim(xlow, xhigh)
    ax.set_yscale('log')
    plt.show()

    return None

def plot_image(arr, vmin, vmax):

    plt.imshow(arr, vmin=vmin, vmax=vmax, origin='lower', interpolation='nearest')
    plt.colorbar()
    plt.show()

    return None

if __name__ == '__main__':
    
    # first convert the clipped rasters (they are simple 
    # txt files) to numpy arrays. then load them in.
    pix_frac_path = slope_extdir + 'pix_area_fraction_clipped.txt'
    crater_frac_path = slope_extdir + 'crater_area_frac_in_pix_clipped.txt'
    #su.raster_to_numpy(pix_frac_path)
    #su.raster_to_numpy(crater_frac_path)

    # load all arrays
    # read in products
    pix_frac = np.load(pix_frac_path.replace('.txt', '.npy'))
    crater_frac = np.load(crater_frac_path.replace('.txt', '.npy'))

    pix_frac = pix_frac.ravel()
    crater_frac = crater_frac.ravel()

    density = se.get_density(crater_frac, pix_frac, len(pix_frac))

    # read in pixel coordinates
    slopemap_path = slope_extdir + 'hf_full_slopemap_clipped.txt'
    #su.raster_to_numpy(slopemap_path)
    slope_arr = np.load(slopemap_path.replace('.txt', '.npy'))
    slope_arr = slope_arr.ravel()


    #pix_x_cen_arr = slope_arr['pix_x_cen']
    #pix_y_cen_arr = slope_arr['pix_y_cen']
    #slope = slope_arr['slope_val']
    # These are all 1D arrays of 4709560 elements

    # reshape to 2D arrays
    """
    rows, columns = se.get_rows_columns(pix_x_cen_arr, pix_y_cen_arr)

    pix_frac_2d = pix_frac.reshape(rows, columns)
    crater_frac_2d = crater_frac.reshape(rows, columns)
    density_2d = density.reshape(rows, columns)
    slope_2d = slope.reshape(rows, columns)
    """

    # choose valid points to plot
    # first, replace all -9999.0 values by NaN
    nodata_idx = np.where(density == -9999.0)
    density[nodata_idx] = np.nan

    nodata_idx = np.where(slope_arr == -9999.0)
    slope_arr[nodata_idx] = np.nan

    # second, for now, also ignore high density values
    #density_upper_lim = 5.0
    #high_idx = np.where(density > density_upper_lim)[0]
    #density[high_idx] = np.nan

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

    #eps = 0.01
    #idx = np.where((density >= 1.0 - eps) & ((density <= 1.0 + eps)))[0]
    #print idx
    #print crater_frac[idx]
    #print crater_diam[idx]
    #print pix_frac[idx]
    #sys.exit(0)

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

    ax.plot(slope_arr, density, 'o', color='k', markersize=2, markeredgecolor='None')
    #ax.plot(x_fit_arr, gauss(x_fit_arr), '-', color='r', lw=2)
    #ax.set_ylim(0,1)

    plt.show()

    sys.exit(0)