from __future__ import division  

import numpy as np
from astropy.modeling import models, fitting
import cPickle
import Polygon as pg
from scipy.misc import factorial
from scipy.optimize import curve_fit

import sys
import os
import time
import datetime

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

def fit_gauss(fit_x_arr, fit_y_arr):

    # make sure that the fitting is done to only the finite
    # elements in the arrays i.e. No NaN or inf in the arrays
    fin_idx_x = np.where(np.isfinite(fit_x_arr))[0]
    fin_idx_y = np.where(np.isfinite(fit_y_arr))[0]
    fin_idx = np.intersect1d(fin_idx_x, fin_idx_y)

    fit_x_arr = fit_x_arr[fin_idx]
    fit_y_arr = fit_y_arr[fin_idx] 

    # Checked: Sorting the arrays before passing them 
    # to the fitting function makes no difference

    # fitting using astropy
    gauss_init = models.Gaussian1D(amplitude=0.5, mean=4.0, stddev=5.0)
    fit_gauss = fitting.LevMarLSQFitter()
    g = fit_gauss(gauss_init, fit_x_arr, fit_y_arr)

    print g.parameters
    print "amp", g.parameters[0]
    print "mean", g.parameters[1]
    print "std", g.parameters[2]

    return g

def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)


def get_diam_ref_arrays(density, slope):

    # first read in crater ids associated with each pixel
    with open(slope_extdir + 'pix_crater_id_fastcomp.pkl', 'rb') as crater_id_file:
        crater_id_in_pix_arr = cPickle.load(crater_id_file)
    # this crater id array is of a different length than the 
    # other "clipped" arrays because it cannot be clipped since
    # it is a list of lists. i.e. it has the original length 
    # of ~4.7 million pixels.

    # now read in crater diam from id and save them 
    # read crater vertices file
    crater_vert_cat = np.genfromtxt(slope_extdir + 'CRATER_FullHF_Vertices_coords.txt', \
        dtype=None, names=True, delimiter=',')
    crater_ids = np.unique(crater_vert_cat['ORIG_FID'])

    # Now loop over all pixels
    color_arr = []  
    # create color array for storing what color point should 
    # be depending on the crater diam(s) on that pixel 
    # also create density and slope arrays again for plotting
    # pixels sorted by crater diameters.
    density_arr_color = []
    slope_arr_color = []

    for i in range(100000):#len(crater_id_in_pix_arr)):

        if (i % 100000) == 0.0:
            print '\r',
            print "At pixel number:",'{0:.2e}'.format(i),\
            "; time taken up to now:",'{0:.2f}'.format((time.time() - start)/60),"minutes.",
            sys.stdout.flush()

        current_crater_ids = crater_id_in_pix_arr[i]
        if len(current_crater_ids) == 0:
            continue

        elif len(current_crater_ids) == 1:
            current_id = current_crater_ids[0]
            current_diam = get_diam(crater_vert_cat, current_id)

            if (current_diam > 4) and (current_diam < 30):
                continue

            else:
                # store density and slope values
                density_arr_color.append(density[i])
                slope_arr_color.append(slope[i])

                if current_diam <= 4:
                    color_arr.append('b')

                if current_diam >= 30:
                    color_arr.append('r')

        elif len(current_crater_ids) > 1:
            for j in range(len(current_crater_ids)):
                current_id = current_crater_ids[j]
                current_diam = get_diam(crater_vert_cat, current_id)

                if (current_diam > 4) and (current_diam < 30):
                    continue

                else:
                    # store density and slope values
                    density_arr_color.append(density[i])
                    slope_arr_color.append(slope[i])

                    if current_diam <= 4:
                        color_arr.append('b')

                    if current_diam >= 30:
                        color_arr.append('r')

    # convert to numpy arrays so you can do array ops
    density_arr_color = np.asarray(density_arr_color)
    slope_arr_color = np.asarray(slope_arr_color)
    color_arr = np.asarray(color_arr).astype(str)

    return density_arr_color, slope_arr_color, color_arr

def plot_by_diam(density, slope):

    # get arrays where the crater diam has been identified by color
    density_arr_color, slope_arr_color, color_arr = get_diam_ref_arrays(density, slope)

    # do the actual plotting
    # perhaps you could make the blue points bigger than the
    # red points simply because there are fewer blue points.
    # i.e. weighting by the size of the crater.
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    b_idx = np.where(color_arr == 'b')[0]
    r_idx = np.where(color_arr == 'r')[0]

    # fit a curve 
    # first define fitting arrays
    fit_x_arr = slope_arr_color[b_idx]
    fit_y_arr = density_arr_color[b_idx]

    gb = fit_gauss(fit_x_arr, fit_y_arr)

    # fit a poisson distribution
    # you could make a histogram and fit to a poisson distribution
    # which might be easier but first fitting a poisson distribution 
    # directly to the data
    fin_idx_x = np.where(np.isfinite(fit_x_arr))[0]
    fin_idx_y = np.where(np.isfinite(fit_y_arr))[0]
    fin_idx = np.intersect1d(fin_idx_x, fin_idx_y)

    fit_x_arr = fit_x_arr[fin_idx]
    fit_y_arr = fit_y_arr[fin_idx] 

    popt, pcov = curve_fit(poisson, fit_x_arr, fit_y_arr, p0=[1.5])
    print popt
    print pcov

    # fit to other diam bins
    fit_x_arr = slope_arr_color[r_idx]
    fit_y_arr = density_arr_color[r_idx]

    gr = fit_gauss(fit_x_arr, fit_y_arr)

    # plot the actual points
    ax.scatter(slope_arr_color[b_idx], density_arr_color[b_idx], s=8, c=color_arr[b_idx], alpha=0.4, edgecolors='none')
    ax.scatter(slope_arr_color[r_idx], density_arr_color[r_idx], s=1, c=color_arr[r_idx], alpha=0.4, edgecolors='none')

    # plot the best fit curves
    x_plot_arr = np.linspace(0,30,1000)
    ax.plot(x_plot_arr, gb(x_plot_arr), ls='-', color='skyblue', lw=2)
    ax.plot(x_plot_arr, gr(x_plot_arr), ls='-', color='pink', lw=2)
    ax.plot(x_plot_arr, poisson(x_plot_arr, *popt), ls='-', color='forestgreen', lw=2)

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    #fig.savefig(slope_extdir + 'slope_v_density_4_and_30_km.png', dpi=300, bbox_inches='tight')
    #fig.savefig(slope_extdir + 'slope_v_density_4_and_30_km.eps', dpi=300, bbox_inches='tight')
    plt.show()

    plt.clf()
    plt.cla()
    plt.close()

    return None

def plot_3d_hist(density, slope):

    # get arrays where the crater diam has been identified by color
    #density_arr_color, slope_arr_color, color_arr = get_diam_ref_arrays(density, slope)

    # get the code to remove all NaN from both arrays
    fin_idx_density = np.where(np.isfinite(density))[0]
    fin_idx_slope = np.where(np.isfinite(slope))[0]
    fin_idx = np.intersect1d(fin_idx_density, fin_idx_slope)

    density = density[fin_idx]
    slope = slope[fin_idx]

    # make 3d histogram
    # code taken from matplotlib examples
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    import matplotlib.colors as colors

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel(r"$\mathrm{Slope}$", fontsize=14, labelpad=10)
    ax.set_ylabel(r"$\mathrm{Density}$", fontsize=14, labelpad=10)
    ax.set_zlabel(r"$\mathrm{log(N)}$", fontsize=14, labelpad=10)

    x, y = slope, density
    hist, xedges, yedges = np.histogram2d(x, y, bins=(40,40)) #, range=[[0, 30], [0, 1.2]])

    # Construct arrays for the anchor positions of the bars.
    # Note: np.meshgrid gives arrays in (ny, nx) so we use 'F' to flatten xpos,
    # ypos in column-major order. For numpy >= 1.7, we could instead call meshgrid
    # with indexing='ij'.
    xpos, ypos = np.meshgrid(xedges[:-1] + xedges[1:], yedges[:-1] + yedges[1:])
    xpos = xpos.flatten('F') / 2
    ypos = ypos.flatten('F') / 2
    zpos = np.zeros_like(xpos)

    # Construct arrays with the dimensions for the bars.
    dx = xedges[1] - xedges[0]
    dy = yedges[1] - xedges[0]
    dz = hist.flatten()
    dz = np.log10(dz)

    # color using colormap 
    my_cmap = plt.get_cmap('cubehelix')
    cNorm = colors.Normalize(vmin=0, vmax=np.max(dz))
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=my_cmap)
    colorval = scalarMap.to_rgba(dz)

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colorval, zsort='average')#, edgecolor='k', linewidth=1)
    ax.view_init(elev=10, azim=22)

    #plt.colorbar(b)

    plt.show()

    return None

def make_plot(density, slope_arr):

    fig = plt.figure(figsize=(16,16))
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    ax.plot(slope_arr, density, 'o', color='k', markersize=1.5, markeredgecolor='None')

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    #fig.savefig(slope_extdir + 'density_slope', dpi=300, bbox_inches='tight')
    plt.show()

    return None

if __name__ == '__main__':

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

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

    #plot_by_diam(density, slope_arr)
    plot_3d_hist(density, slope_arr)

    # total run time
    print '\n'
    print "Total time taken --", (time.time() - start)/60, "minutes."
    sys.exit(0)
