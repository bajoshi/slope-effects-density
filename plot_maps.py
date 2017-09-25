from __future__ import division  

import numpy as np
from astropy.modeling import models, fitting
import cPickle
import Polygon as pg

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

def get_diam(crater_vert_cat, crater_id):

    crater_vert = crater_vert_cat

    # Using the explicit x and y crater vertices create a polygon
    current_crater_vert_idx = np.where(crater_vert['ORIG_FID'] == crater_id)
    current_x_vert = crater_vert['x_coord_m'][current_crater_vert_idx]
    current_y_vert = crater_vert['y_coord_m'][current_crater_vert_idx]

    crater_poly = pg.Polygon(zip(current_x_vert, current_y_vert))

    # Check the diameter of the crater
    # The area will be in sq. meters. 
    diam = np.sqrt(crater_poly.area() * 4 / np.pi)
    diam /= 1000  # convert to km

    return diam

def get_idx_in_clipped_arr(idx_in_unclipped_arr):

    long_idx = idx_in_unclipped_arr

    return short_idx 

def plot_by_diam(density, slope):

    # first read in crater ids associated with each pixel
    with open(slope_extdir + 'pix_crater_id_fastcomp.pkl', 'rb') as crater_id_file:
        crater_id_in_pix_arr = cPickle.load(crater_id_file)
    # this crater id array is of a different length than the 
    # other "clipped" arrays because it cannot be clipped since
    # it is a list of lists. i.e. it has the original length 
    # of ~4.7 million pixels.

    # now read in crater diam from id and save them 
    # read crater vertices file
    crater_vert = np.genfromtxt(slope_extdir + 'CRATER_FullHF_Vertices_coords.txt', \
        dtype=None, names=True, delimiter=',')
    crater_ids = np.unique(crater_vert['ORIG_FID'])

    # Now loop over all pixels
    color_arr = []  
    # create color array for storing what color point should 
    # be depending on the crater diam(s) on that pixel 
    # also create density and slope arrays again for plotting
    # pixels sorted by crater diameters.
    density_arr_color = []
    slope_arr_color = []

    for i in range(len(crater_id_in_pix_arr)):

        current_crater_ids = crater_id_in_pix_arr[i]
        if len(current_crater_ids) == 0:
            continue

        elif len(current_crater_ids) == 1:
            current_id = current_crater_ids[0]
            current_diam = get_diam(crater_vert, crater_id)

            if (current_diam > 4) and (current_diam < 30):
                continue

            else:
                # I need the index in the clipped array 
                # because the crater id array is not the 
                # same size as the density and slope arrays
                # which are clipped by ArcGIS.
                idx_in_clipped_arr = get_idx_in_clipped_arr(i)
                # store censity and slope values
                density_arr_color.append(density[idx_in_clipped_arr])
                slope_arr_color.append(slope[idx_in_clipped_arr])

                if current_diam <= 4:
                    color_arr.append('b')

                if current_diam >= 30:
                    color_arr.append('r')

        elif len(current_crater_ids) > 1:
            for j in range(len(current_crater_ids)):


    return None

if __name__ == '__main__':

    use_point_density = False
    if use_point_density:
        density_path = slope_extdir + 'pointdensity_bull.txt'
        #su.raster_to_numpy(density_path)
        density = np.load(density_path.replace('.txt','.npy'))
        density = density.ravel()

    else:
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

    # save density raster
    #np.save(slope_extdir + 'density_arr_bool_clipped.npy', density)
    #su.numpy_to_asciiraster(slope_extdir + 'density_arr_bool_clipped.npy', (2109,1949), -3822716.7379517, -1728009.4370101)

    # choose valid points to plot
    # first, replace all -9999.0 values by NaN
    nodata_idx = np.where(density == -9999.0)
    density[nodata_idx] = np.nan

    density_zero_idx = np.where(density == 0.0)  # NaNing this out because these are pixels where there are no craters
    density[density_zero_idx] = np.nan

    nodata_idx = np.where(slope_arr == -9999.0)
    slope_arr[nodata_idx] = np.nan

    plot_by_diam(density, slope_arr)
    sys.exit(0)

    # plots
    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    ax.plot(slope_arr, density, 'o', color='k', markersize=1.5, markeredgecolor='None')

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(slope_extdir + 'density_slope_fuzzy', dpi=300, bbox_inches='tight')
    #plt.show()

    sys.exit(0)