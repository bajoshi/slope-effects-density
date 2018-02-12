from __future__ import division

import numpy as np

import sys
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

sys.path.append(slopedir)
import slope_effects_density as se
import slope_utils as su

if __name__ == '__main__':
    
    # read in all density and slope rasters
    # and put them in numpy arrays.
    #### -------- Point Density -------- ####
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

    #### -------- Boolean Density -------- ####
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

    #### -------- Fuzzy Density -------- ####
    fuzzy_crater_frac = np.load(slope_extdir + 'crater_area_frac_in_pix_clipped.npy')
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
    gs = gridspec.GridSpec(1, 3)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.15, hspace=0.15)

    fig = plt.figure(figsize=(9,5))
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[0,2])

    ax1.scatter(slope_pointdensity, point_density, marker='o', s=3, color='k', edgecolors='None', alpha=0.75)
    ax2.scatter(slope_booleandensity, boolean_density, marker='o', s=5, color='k', edgecolors='None', alpha=0.75)
    ax3.scatter(slope_fuzzydensity, fuzzy_density, marker='o', s=1, color='k', edgecolors='None', alpha=0.75)

    ax1.set_ylabel('Density', fontsize=15)
    ax2.set_xlabel('Slope', fontsize=15)

    # add minor ticks and grid
    ax1.minorticks_on()
    ax1.tick_params('both', width=1, length=3, which='minor')
    ax1.tick_params('both', width=1, length=4.7, which='major')
    ax2.minorticks_on()
    ax2.tick_params('both', width=1, length=3, which='minor')
    ax2.tick_params('both', width=1, length=4.7, which='major')
    ax3.minorticks_on()
    ax3.tick_params('both', width=1, length=3, which='minor')
    ax3.tick_params('both', width=1, length=4.7, which='major')

    # Set limits 
    ax1.set_ylim(-0.001, 0.023)
    ax2.set_ylim(-0.0625, 2.0625)
    ax3.set_ylim(-0.05, 1.35)

    ax1.set_xlim(-2, 37)
    ax2.set_xlim(-2, 37)
    ax3.set_xlim(-2, 37)

    # add text that says which method is plotted
    ax1.text(0.603, 0.95, '(a) Point', verticalalignment='top', horizontalalignment='left', \
        transform=ax1.transAxes, color='k', size=10)
    ax1.text(0.7, 0.91, 'Density', verticalalignment='top', horizontalalignment='left', \
        transform=ax1.transAxes, color='k', size=10)

    ax2.text(0.6, 0.95, '(b) Boolean', verticalalignment='top', horizontalalignment='left', \
        transform=ax2.transAxes, color='k', size=10)
    ax2.text(0.7, 0.91, 'Density', verticalalignment='top', horizontalalignment='left', \
        transform=ax2.transAxes, color='k', size=10)

    ax3.text(0.535, 0.95, '(c) Fractional', verticalalignment='top', horizontalalignment='left', \
        transform=ax3.transAxes, color='k', size=10)
    ax3.text(0.63, 0.91, 'Density', verticalalignment='top', horizontalalignment='left', \
        transform=ax3.transAxes, color='k', size=10)

    # save figure
    fig.savefig(slope_extdir + 'all_crater_density_methods_comparison.png', dpi=300, bbox_inches='tight')

    sys.exit(0)