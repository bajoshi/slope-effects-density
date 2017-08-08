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

if __name__ == '__main__':
    
    # load all arrays
    # read in products
    pix_frac = np.load(slope_extdir + 'pix_area_fraction.npy')
    crater_frac = np.load(slope_extdir + 'crater_area_frac_in_pix.npy')

    density = np.load(slope_extdir + 'density_array.npy')

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
    density_upper_lim = 1.0
    high_idx = np.where(density > density_upper_lim)[0]
    density[high_idx] = np.nan

    # third,     
    val_idx = np.where(density != 0.0)[0] # np.isfinite(density)

    # ---------------------------------------- #
    # find frac pix values for pixels that have high density values
    

    # --------------------------------------- # 
    # Fit gaussian
    gauss_init = models.Gaussian1D(amplitude=1.0, mean=5.0, stddev=10.0)
    fit_gauss = fitting.LevMarLSQFitter()
    x_fit_arr = slope[val_idx]
    y_fit_arr = density[val_idx]
    #gauss = fit_gauss(gauss_init, x_fit_arr, y_fit_arr)

    # plots
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #ax.plot(pix_frac[val_idx], np.log10(density[val_idx]), 'o', color='k', markersize=2, markeredgecolor='None')
    #ax.set_ylim(-5,2.5)

    ax.plot(slope[val_idx], density[val_idx], 'o', color='k', markersize=2, markeredgecolor='None')
    #ax.plot(x_fit_arr, gauss(x_fit_arr), '-', color='r', lw=2)
    ax.set_ylim(0,1)

    plt.show()

    sys.exit(0)