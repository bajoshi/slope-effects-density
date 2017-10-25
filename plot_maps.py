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
import matplotlib.gridspec as gridspec

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

def get_diam_mycalc(crater_vert_cat, crater_id):

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

def get_diam(crater_vert_cat, crater_id):
    """
    This function returns the diameter computed by ArcGIS.
    Arc's diameter is favored over my calculation because 
    it has taken the correct projection into account.
    """

    crater_vert = crater_vert_cat

    current_crater_vert_idx = np.where(crater_vert['ORIG_FID'] == crater_id)
    diam = crater_vert['Diam_km'][current_crater_vert_idx][0]
    if diam < 1: print "Diameter less than 1km:", diam, "km. This should not be here."

    return diam

def fit_gauss(fit_x_arr, fit_y_arr):

    fit_y_arr = np.log10(fit_y_arr)

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

def plot_crater_diam_hist():

    # now read in crater diam from id and save them 
    # read crater vertices file
    crater_vert_cat = np.genfromtxt(slope_extdir + 'CRATER_FullHF_Vertices_coords.txt', \
        dtype=None, names=True, delimiter=',')
    crater_ids = np.unique(crater_vert_cat['ORIG_FID'])

    # loop over all craters and store 
    # their diameters in an array.
    crater_diam_array = np.empty(len(crater_ids))
    for i in range(len(crater_ids)):

        current_id = crater_ids[i]
        current_diam = get_diam(crater_vert_cat, current_id)
        crater_diam_array[i] = current_diam

    print len(np.where(crater_diam_array > 50)[0]), "craters larger than 50km"

    # plot hist of diameters
    # don't consider any craters larger than 50km
    #nolargediam_idx = np.where(crater_diam_array <= 50)[0]
    #crater_diam_array = crater_diam_array[nolargediam_idx]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ncount, edges, patches = ax.hist(crater_diam_array, 100, color='lightgray', align='mid')
    ax.axhline(y=10, ls='--', lw=2, color='k')

    # color all selected crater diam bins a different color
    max_crater_diam_bin = edges[min(np.where(ncount < 10)[0]) - 1]

    edges_plot = np.where(ncount >= 10)[0] #np.where((edges >= 0) & (edges <= max_crater_diam_bin))[0]
    patches_plot = [patches[edge_ind] for edge_ind in edges_plot]
    col = np.full(len(patches_plot), 'lightblue', dtype='|S9')
    # make sure the length of the string given in the array initialization is the same as the color name
    for c, p in zip(col, patches_plot):
        plt.setp(p, 'facecolor', c)

    ax.set_xlabel(r'$\mathrm{Crater\ Diameter\ [KM]}$', fontsize=15)
    ax.set_ylabel(r'$\mathrm{N}$', fontsize=15)
    ax.grid(True)

    ax.set_yscale('log')

    fig.savefig(slope_extdir + 'crater_diam_hist.png', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    return None

def get_ids():

    # first read in crater ids associated with each pixel
    with open(slope_extdir + 'pix_crater_id_fastcomp.pkl', 'rb') as crater_id_file:
        crater_id_in_pix_arr = cPickle.load(crater_id_file)

    # now read in crater diam from id and save them 
    # read crater vertices file
    crater_vert_cat = np.genfromtxt(slope_extdir + 'CRATER_FullHF_Vertices_coords.txt', \
        dtype=None, names=True, delimiter=',')
    crater_ids = np.unique(crater_vert_cat['ORIG_FID'])

    return crater_ids, crater_id_in_pix_arr, crater_vert_cat

def get_diam_ref_arrays(density, slope, crater_vert_cat, crater_id_in_pix_arr, start):

    # loop over all pixels
    color_arr = []  
    # create color array for storing what color point should 
    # be depending on the crater diam(s) on that pixel 
    # also create density and slope arrays again for plotting
    # pixels sorted by crater diameters.
    density_arr_color = []
    slope_arr_color = []
    pix_1d_idx_arr = []
    crater_id_diam = []
    # future:
    # this is a list of lists where the 1st list dimension has length equal to the number of diam bins
    # i.e. if there are say 2 diam bins then this is a 2 element list but each element is itself a list
    # containing pixel indices that fall into that diam bin

    for i in range(len(crater_id_in_pix_arr)):

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

            if (current_diam > 35):
                continue

            else:
                # store density and slope values
                density_arr_color.append(density[i])
                slope_arr_color.append(slope[i])

                if current_diam <= 5:
                    color_arr.append('blue')
                    pix_1d_idx_arr.append(i)
                    crater_id_diam.append(current_id)

                elif (current_diam > 5) and (current_diam <= 10):
                    color_arr.append('cyan')

                elif (current_diam > 10) and (current_diam <= 15):
                    color_arr.append('green')

                elif (current_diam > 15) and (current_diam <= 20):
                    color_arr.append('olive')

                elif (current_diam > 20) and (current_diam <= 25):
                    color_arr.append('goldenrod')

                elif (current_diam > 25) and (current_diam <= 30):
                    color_arr.append('darkorchid')

                elif (current_diam > 30) and (current_diam <= 35):
                    color_arr.append('maroon')
                    pix_1d_idx_arr.append(i)
                    crater_id_diam.append(current_id)

        elif len(current_crater_ids) > 1:
            for j in range(len(current_crater_ids)):
                current_id = current_crater_ids[j]
                current_diam = get_diam(crater_vert_cat, current_id)

                if (current_diam > 35):
                    continue

                else:
                    # store density and slope values
                    density_arr_color.append(density[i])
                    slope_arr_color.append(slope[i])

                    if current_diam <= 5:
                        color_arr.append('blue')

                    elif (current_diam > 5) and (current_diam <= 10):
                        color_arr.append('cyan')

                    elif (current_diam > 10) and (current_diam <= 15):
                        color_arr.append('green')

                    elif (current_diam > 15) and (current_diam <= 20):
                        color_arr.append('olive')

                    elif (current_diam > 20) and (current_diam <= 25):
                        color_arr.append('goldenrod')

                    elif (current_diam > 25) and (current_diam <= 30):
                        color_arr.append('darkorchid')

                    elif (current_diam > 30) and (current_diam <= 35):
                        color_arr.append('maroon')

    # convert to numpy arrays so you can do array ops
    density_arr_color = np.asarray(density_arr_color)
    slope_arr_color = np.asarray(slope_arr_color)
    color_arr = np.asarray(color_arr).astype(str)
    pix_1d_idx_arr = np.asarray(pix_1d_idx_arr)
    crater_id_diam = np.asarray(crater_id_diam)

    return density_arr_color, slope_arr_color, color_arr, pix_1d_idx_arr, crater_id_diam

def plot_by_diam(density, slope_arr, start):

    crater_ids, crater_id_in_pix_arr, crater_vert_cat = get_ids()

    # get arrays where the crater diam has been identified by color
    density_arr_color, slope_arr_color, color_arr, pix_1d_idx_arr, crater_id_diam = \
    get_diam_ref_arrays(density, slope_arr, crater_vert_cat, crater_id_in_pix_arr, start)

    # do the actual plotting
    # perhaps you could make the blue points bigger than the
    # red points simply because there are fewer blue points.
    # i.e. weighting by the size of the crater.
    b_idx = np.where(color_arr == 'blue')[0]  # '1to5'
    c_idx = np.where(color_arr == 'cyan')[0]  # '5to10'
    g_idx = np.where(color_arr == 'green')[0]  # '10to15'
    ol_idx = np.where(color_arr == 'olive')[0]  # '15to20'
    gr_idx = np.where(color_arr == 'goldenrod')[0]  # '20to25'
    do_idx = np.where(color_arr == 'darkorchid')[0]  # '25to30'
    m_idx = np.where(color_arr == 'maroon')[0]  # '30to35'

    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, b_idx, '1to5',    0.051, 1.27)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, c_idx, '5to10',   0.013, 0.051)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, g_idx, '10to15',  5.66e-3, 0.013)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, ol_idx, '15to20', 3.18e-3, 5.66e-3)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, gr_idx, '20to25', 2.04e-3, 3.18e-3)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, do_idx, '25to30', 1.42e-3, 2.04e-3)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, m_idx, '30to35',  1.04e-3, 1.42e-3)

    # Put all together
    # plot these 7 bins on a 7 panel grid in a single plot
    # and also all diam bins in the same plot
    # ----------------------------- GRID PLOT ----------------------------- #
    gs = gridspec.GridSpec(3,3)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure()
    ax_b = fig.add_subplot(gs[0, 0])
    ax_c = fig.add_subplot(gs[0, 1])
    ax_g = fig.add_subplot(gs[0, 2])
    ax_ol = fig.add_subplot(gs[1, 0])
    ax_gr = fig.add_subplot(gs[1, 1])
    ax_do = fig.add_subplot(gs[1, 2])
    ax_m = fig.add_subplot(gs[2, 1])

    all_axes = [ax_b, ax_c, ax_g, ax_ol, ax_gr, ax_do, ax_m]
    all_diam_idx = [b_idx, c_idx, g_idx, ol_idx, gr_idx, do_idx, m_idx]
    min_val_list = [0.051, 0.013, 5.66e-3, 3.18e-3, 2.04e-3, 1.42e-3, 1.04e-3]
    max_val_list = [1.27, 0.051, 0.013, 5.66e-3, 3.18e-3, 2.04e-3, 1.42e-3]

    ax_m.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax_ol.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    for i in range(len(all_axes)):

        all_axes[i].scatter(slope_arr_color[all_diam_idx[i]], density_arr_color[all_diam_idx[i]], s=5, c=color_arr[all_diam_idx[i]], alpha=0.4, edgecolors='none')
        all_axes[i].set_yscale('log')
        all_axes[i].set_ylim(1e-8, 2.0)
        all_axes[i].set_xlim(0, 30)

        all_axes[i].axhline(y=min_val_list[i], ls='--', color='steelblue')  # min value of density from biggest crater in bin
        all_axes[i].axhline(y=max_val_list[i], ls='--', color='lightblue')  # max value of density from smallest crater in bin

        all_axes[i].minorticks_on()
        all_axes[i].tick_params('both', width=1, length=3, which='minor')
        all_axes[i].tick_params('both', width=1, length=4.7, which='major')
        all_axes[i].grid(True)

        if i < 3:
            all_axes[i].set_xticklabels([])

        if i == 1 or i == 2 or i == 4 or i == 5:
            all_axes[i].set_yticklabels([])


    ax_ol.set_xticklabels(['0', '10', '20', ''])
    ax_do.set_xticklabels(['', '10', '20', '30'])

    fig.savefig(slope_extdir + 'slope_v_density_all_grid.png', dpi=300, bbox_inches='tight')
    fig.savefig(slope_extdir + 'slope_v_density_all_grid.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    # ----------------------------- ALL TOGETHER PLOT ----------------------------- #
    # Uses some of the lists from the grid plot so don't comment that one out
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    for i in range(len(all_axes)):

        ax.scatter(slope_arr_color[all_diam_idx[i]], density_arr_color[all_diam_idx[i]], s=5, c=color_arr[all_diam_idx[i]], alpha=0.4, edgecolors='none')

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(slope_extdir + 'slope_v_density_all_together.png', dpi=300, bbox_inches='tight')
    fig.savefig(slope_extdir + 'slope_v_density_all_together.eps', dpi=300, bbox_inches='tight')

    plt.clf()
    plt.cla()
    plt.close()

    return None

def make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, diam_bin_idx, diam_bin, min_val, max_val):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    ax.scatter(slope_arr_color[diam_bin_idx], density_arr_color[diam_bin_idx], s=5, c=color_arr[diam_bin_idx], alpha=0.4, edgecolors='none')
    ax.set_yscale('log')
    ax.set_ylim(1e-8, 2.0)
    ax.set_xlim(0, 30)

    ax.axhline(y=min_val, ls='--', color='steelblue')  # min value of density from biggest crater in bin
    ax.axhline(y=max_val, ls='--', color='lightblue')  # max value of density from smallest crater in bin
    # values are obtained for a pixel that is completely inside 
    # a crater and not overlapped by any other craters.

    # fit a curve 
    # first define fitting arrays
    """
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

    fit_y_arr = np.log10(fit_y_arr)

    popt, pcov = curve_fit(poisson, fit_x_arr, fit_y_arr, p0=[1.5])
    print popt
    print pcov
    """

    # plot the best fit curves
    #x_plot_arr = np.linspace(0,30,1000)
    #ax.plot(x_plot_arr, gb(x_plot_arr), ls='-', color='skyblue', lw=2)
    #ax.plot(x_plot_arr, gr(x_plot_arr), ls='-', color='pink', lw=2)
    #ax.plot(x_plot_arr, poisson(x_plot_arr, *popt), ls='-', color='forestgreen', lw=2)

    #print np.nanmin(density_arr_color[b_idx])
    #print np.nanmax(density_arr_color[b_idx])

    #check_idx = np.where((density_arr_color[r_idx] >= 1e-3) & (density_arr_color[r_idx] <= np.power(10, -2.9)))[0]
    #print check_idx
    #print len(check_idx)
    #print pix_1d_idx_arr[r_idx][check_idx]

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(slope_extdir + 'slope_v_density_' + diam_bin + 'km.png', dpi=300, bbox_inches='tight')
    fig.savefig(slope_extdir + 'slope_v_density_' + diam_bin + 'km.eps', dpi=300, bbox_inches='tight')

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
    hist, xedges, yedges = np.histogram2d(x, y, bins=(50,50)) #, range=[[0, 30], [0, 1.2]])

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

def convert_npy_array_tofits(npy_array, final_shape, savedir, final_name):
    """
    This code will convert the supplied numpy array
    to a 2D fits image file. The length of the array 
    must be consistent with the required final shape.
    The given array should be in the binary numpy format.
    """

    from astropy.io import fits

    npy_array = npy_array.reshape(final_shape)
    hdu = fits.PrimaryHDU(data=npy_array)
    hdu.writeto(savedir + final_name + '.fits', overwrite=True)

    return None

if __name__ == '__main__':

    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # uncomment the following two lines if you 
    # want to make the crater diam histogram
    #plot_crater_diam_hist()
    #sys.exit(0)

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
        su.raster_to_numpy(pix_frac_path)
        su.raster_to_numpy(crater_frac_path)

        # load all arrays
        # read in products
        pix_frac = np.load(pix_frac_path.replace('.txt', '.npy'))
        crater_frac = np.load(crater_frac_path.replace('.txt', '.npy'))

        pix_frac = pix_frac.ravel()
        crater_frac = crater_frac.ravel()

        density = se.get_density(crater_frac, pix_frac, len(pix_frac))

    # read in pixel coordinates
    slopemap_path = slope_extdir + 'hf_full_slopemap_clipped.txt'
    su.raster_to_numpy(slopemap_path)
    slope_arr = np.load(slopemap_path.replace('.txt', '.npy'))
    slope_arr = slope_arr.ravel()

    # use the following lines to convert arrays to fits files
    """
    rows = 2109
    columns = 1949
    convert_npy_array_tofits(density, (rows, columns), slope_extdir, 'density_array_raw')
    convert_npy_array_tofits(slope_arr, (rows, columns), slope_extdir, 'slope_array_raw')
    sys.exit(0)
    """

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

    plot_by_diam(density, slope_arr, start)
    #plot_3d_hist(density, slope_arr)

    # total run time
    print '\n'
    print "Total time taken --", (time.time() - start)/60, "minutes."
    sys.exit(0)
