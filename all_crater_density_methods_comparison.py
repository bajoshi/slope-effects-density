from __future__ import division

import numpy as np
from scipy.stats import gaussian_kde

import sys
import os
import time
import datetime

import matplotlib.pyplot as plt
import matplotlib.cm as cm

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

sys.path.append(slopedir)
import slope_effects_density as se
import slope_utils as su

if __name__ == '__main__':
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in all density and slope rasters
    # and put them in numpy arrays.
    #### --------------------------------- Point Density --------------------------------- ####
    point_density_path = slope_extdir + 'pointdensityraster.txt'
    slope_pointdensity_path = slope_extdir + 'slopes4pointdensity.txt'

    #su.raster_to_numpy(point_density_path)
    #su.raster_to_numpy(slope_pointdensity_path)

    point_density = np.load(point_density_path.replace('.txt','.npy'))
    slope_pointdensity = np.load(slope_pointdensity_path.replace('.txt','.npy'))
    
    # turn them to 1D arrays
    point_density = point_density.ravel()
    slope_pointdensity = slope_pointdensity.ravel()

    # make sure to set the -9999.0 (i.e. NODATA values) to NaN
    pointdensity_invalid_idx = np.where(point_density == -9999.0)[0]
    slope_pointdensity_invalid_idx = np.where(slope_pointdensity == -9999.0)[0]
    invalid_idx = np.intersect1d(pointdensity_invalid_idx, slope_pointdensity_invalid_idx)

    point_density[invalid_idx] = np.nan
    slope_pointdensity[invalid_idx] = np.nan

    #### --------------------------------- Boolean Density --------------------------------- ####
    boolean_density_path = slope_extdir + 'booleandensityraster.txt'
    slope_booleandensity_path = slope_extdir + 'slopes4boolean.txt'

    #su.raster_to_numpy(boolean_density_path)
    #su.raster_to_numpy(slope_booleandensity_path)

    boolean_density = np.load(boolean_density_path.replace('.txt','.npy'))
    slope_booleandensity = np.load(slope_booleandensity_path.replace('.txt','.npy'))

    # turn them to 1D arrays
    boolean_density = boolean_density.ravel()
    slope_booleandensity = slope_booleandensity.ravel()

    # make sure to set the -9999.0 (i.e. NODATA values) to NaN
    booleandensity_invalid_idx = np.where(boolean_density == -9999.0)[0]
    slope_booleandensity_invalid_idx = np.where(slope_booleandensity == -9999.0)[0]
    invalid_idx = np.intersect1d(booleandensity_invalid_idx, slope_booleandensity_invalid_idx)

    boolean_density[invalid_idx] = np.nan
    slope_booleandensity[invalid_idx] = np.nan

    #### --------------------------------- Fuzzy Density --------------------------------- ####
    fuzzy_crater_frac = np.load(slope_extdir + 'crater_area_frac_in_pix_fastcomp_newbool.npy')
    fuzzy_pix_frac = np.load(slope_extdir + 'pix_area_fraction_clipped.npy')
    slope_fuzzydensity = np.load(slope_extdir + 'hf_full_slopemap_clipped.npy')

    slope_fuzzydensity = slope_fuzzydensity.ravel()
    fuzzy_pix_frac = fuzzy_pix_frac.ravel()
    fuzzy_crater_frac = fuzzy_crater_frac.ravel()

    fuzzy_density = se.get_density(fuzzy_crater_frac, fuzzy_pix_frac, len(fuzzy_pix_frac))

    # ------------
    # why does the boolean method have a zero in it?
    # there should only be ones and twos. Any pixels 
    # with zero (i.e. none) craters in them should 
    # have been assigned the NODATA value and should
    # not be in the raster or the corresponding numpy array.

    ######### --------------------------------- Plotting --------------------------------- #########
    # Now plot all three on the same figure

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    #ax3 = fig.add_subplot(gs[0,2])

    ax1.set_ylabel(r'$\mathrm{Density\ [Craters\, km^{-2}]}$', fontsize=15)
    ax2.set_xlabel(r'$\mathrm{Slope\ [Degrees]}$', fontsize=15)

    # add minor ticks and grid
    ax1.minorticks_on()
    ax2.minorticks_on()

    # Set limits 
    ax1.set_ylim(-0.001, 0.023)
    #ax2.set_ylim(-0.0625, 2.0625)
    #ax3.set_ylim(-0.05, 1.35)

    ax1.set_xlim(-2, 37)
    ax2.set_xlim(-2, 37)
    #ax3.set_xlim(-2, 37)

    # add text that says which method is plotted
    ax1.text(0.7, 0.95, 'Point' + '\n' + 'density', verticalalignment='top', horizontalalignment='left', \
        transform=ax1.transAxes, color='k', size=10)
    ax2.text(0.7, 0.95, 'Pixel-by-pixel' + '\n' + 'density', verticalalignment='top', horizontalalignment='left', \
        transform=ax2.transAxes, color='k', size=10)

    #ax3.text(0.535, 0.95, '(c) Fractional', verticalalignment='top', horizontalalignment='left', \
    #    transform=ax3.transAxes, color='k', size=10)
    #ax3.text(0.63, 0.91, 'Density', verticalalignment='top', horizontalalignment='left', \
    #    transform=ax3.transAxes, color='k', size=10)

    # --------- Contours --------- #
    print "Computing contours"
    print "Time taken until now --", (time.time() - start), "seconds."
    # First get rid of invalid entries
    slope_val_idx = np.where(slope_fuzzydensity != -9999.0)  # no NaNs in the slope arr # checked
    fuzzy_density_val_idx = np.where(~np.isnan(fuzzy_density))
    val_idx = reduce(np.intersect1d, (slope_val_idx, fuzzy_density_val_idx))
    slope_fuzzydensity = slope_fuzzydensity[val_idx]
    fuzzy_density = fuzzy_density[val_idx]
    # plot contours for point density
    """
    counts, xbins, ybins = np.histogram2d(slope_fuzzydensity, fuzzy_density, bins=100, normed=False)
    levels_to_plot = [20, 50, 100, 200, 500, 800, 1e3, 1e4]
    c = ax2.contour(counts.transpose(), levels=levels_to_plot, \
        extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
        cmap=cm.Blues_r, linestyles='solid', zorder=10)
    """

    xmin = min(slope_fuzzydensity)
    xmax = max(slope_fuzzydensity)
    ymin = min(fuzzy_density)
    ymax = max(fuzzy_density)
    
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([slope_fuzzydensity, fuzzy_density])
    kernel = gaussian_kde(values)
    
    density = np.reshape(kernel(positions).T, X.shape)

    print "Plotting KDE for Pixel-by-pixel density."
    print "Time taken until now --", (time.time() - start), "seconds."
    ax2.imshow(np.rot90(density), cmap=plt.cm.gist_earth_r)

    ax1.scatter(slope_pointdensity, point_density, marker='o', s=3, color='k', edgecolors='None', alpha=0.75)
    #ax2.scatter(slope_booleandensity, boolean_density, marker='o', s=5, color='k', edgecolors='None', alpha=0.75)
    ax2.scatter(slope_fuzzydensity, fuzzy_density, marker='o', s=1, color='k', edgecolors='None', alpha=0.75)

    ax1.set_aspect(1.0)
    ax2.set_aspect(1.0)

    # save figure
    fig.savefig(slope_extdir + 'all_crater_density_methods_comparison_newbool.png', dpi=300, bbox_inches='tight')

    # total run time
    print "Total time taken --", (time.time() - start)/60, "minutes."
    sys.exit(0)