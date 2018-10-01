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

    # Get rid of NaN entries for KDE
    slope_val_idx = np.where(~np.isnan(slope_pointdensity))
    point_density_val_idx = np.where(~np.isnan(point_density))
    val_idx = reduce(np.intersect1d, (slope_val_idx, point_density_val_idx))
    slope_pointdensity = slope_pointdensity[val_idx]
    point_density = point_density[val_idx]

    # KDE estimate
    xmin = min(slope_pointdensity)
    xmax = max(slope_pointdensity)
    ymin = min(point_density)
    ymax = max(point_density)

    print xmin, xmax, ymin, ymax
    
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([slope_pointdensity, point_density])
    kernel = gaussian_kde(values)
    
    density = np.reshape(kernel(positions).T, X.shape)

    print "Plotting KDE for point density."
    print "Time taken until now --", (time.time() - start), "seconds."

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{Density\ [Craters\, km^{-2}]}$', fontsize=15)
    ax.set_xlabel(r'$\mathrm{Slope\ [Degrees]}$', fontsize=15)

    ax.imshow(np.rot90(density), cmap=plt.cm.gist_stern_r, extent=[0, 35, 0, 0.02], alpha=0.6)
    ax.scatter(slope_pointdensity, point_density, marker='o', s=1, color='k', edgecolors='None', alpha=0.75)

    ax.axis('normal')
    ax.set_ylim(-0.001, 0.023)

    ax.minorticks_on()
    ax.text(0.65, 0.95, 'Point Density', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)

    fig.savefig(slope_extdir + 'point_density_newbool_kde.png', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ------------------------------------------------------------------ #
    #### --------------------------------- New Boolean Density --------------------------------- ####
    fuzzy_crater_frac = np.load(slope_extdir + 'crater_area_frac_in_pix_fastcomp_newbool.npy')
    fuzzy_pix_frac = np.load(slope_extdir + 'pix_area_fraction_clipped.npy')
    slope_fuzzydensity = np.load(slope_extdir + 'hf_full_slopemap_clipped.npy')

    slope_fuzzydensity = slope_fuzzydensity.ravel()
    fuzzy_pix_frac = fuzzy_pix_frac.ravel()
    fuzzy_crater_frac = fuzzy_crater_frac.ravel()

    fuzzy_density = se.get_density(fuzzy_crater_frac, fuzzy_pix_frac, len(fuzzy_pix_frac))

    # --------- Contours --------- #
    print "Computing contours"
    print "Time taken until now --", (time.time() - start), "seconds."
    # First get rid of invalid entries
    slope_val_idx = np.where(slope_fuzzydensity != -9999.0)  # no NaNs in the slope arr # checked
    fuzzy_density_val_idx = np.where(~np.isnan(fuzzy_density))
    val_idx = reduce(np.intersect1d, (slope_val_idx, fuzzy_density_val_idx))
    slope_fuzzydensity = slope_fuzzydensity[val_idx]
    fuzzy_density = fuzzy_density[val_idx]

    include_idx = np.where(fuzzy_density >= 1.0)
    fuzzy_density = fuzzy_density[include_idx]
    slope_fuzzydensity = slope_fuzzydensity[include_idx]

    # KDE estimate
    xmin = min(slope_fuzzydensity)
    xmax = max(slope_fuzzydensity)
    ymin = min(fuzzy_density)
    ymax = max(fuzzy_density)
    print xmin, xmax, ymin, ymax

    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([slope_fuzzydensity, fuzzy_density])
    kernel = gaussian_kde(values, bw_method=2)
    print "Bandwidth factor:", kernel.factor
    
    density = np.reshape(kernel(positions).T, X.shape)

    print "Plotting KDE for point density."
    print "Time taken until now --", (time.time() - start), "seconds."

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{Density\ [Craters\, km^{-2}]}$', fontsize=15)
    ax.set_xlabel(r'$\mathrm{Slope\ [Degrees]}$', fontsize=15)

    ax.imshow(np.rot90(density), cmap=plt.cm.gist_stern_r, extent=[0, 35, 1.0, 6.0], alpha=0.6)
    ax.scatter(slope_fuzzydensity, fuzzy_density, marker='o', s=1, color='k', edgecolors='None', alpha=0.75)

    ax.axis('normal')
    ax.set_ylim(0.5, 6.5)

    ax.minorticks_on()
    ax.text(0.65, 0.95, 'Modified Boolean Density', verticalalignment='top', horizontalalignment='left', \
        transform=ax.transAxes, color='k', size=10)

    fig.savefig(slope_extdir + 'newboolean_density_kde.png', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    sys.exit(0)