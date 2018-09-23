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

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

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

def get_diam(crater_ids_arr, crater_diam_m_arr, crater_id):
    """
    This function returns the diameter computed by ArcGIS.
    Arc's diameter is favored over my calculation because 
    it has taken the correct projection into account.
    """

    current_crater_idx = np.where(crater_ids_arr == crater_id)[0]
    diam = crater_diam_m_arr[current_crater_idx[0]] / 1e3  # convert to km
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

    # read in crater diam from id and save them 
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

    # read required arrays
    crater_ids_arr = np.load(slope_extdir + 'crater_ids_arr.npy')
    crater_diam_m_arr = np.load(slope_extdir + 'crater_diam_m_arr.npy')

    crater_unique_ids_arr = np.unique(crater_ids_arr)

    # first read in crater ids associated with each pixel
    with open(slope_extdir + 'pix_crater_id_fastcomp.pkl', 'rb') as crater_id_file:
        crater_id_in_pix_arr = cPickle.load(crater_id_file)

    return crater_ids_arr, crater_diam_m_arr, crater_id_in_pix_arr

def assign_color(current_diam):

    if current_diam <= 2:
        return 'midnightblue'

    elif (current_diam > 2) and (current_diam <= 3):
        return 'blue'

    elif (current_diam > 3) and (current_diam <= 4):
        return 'royalblue'

    elif (current_diam > 4) and (current_diam <= 5):
        return 'dodgerblue'

    elif (current_diam > 5) and (current_diam <= 6):
        return 'deepskyblue'

    elif (current_diam > 6) and (current_diam <= 7):
        return 'steelblue'

    elif (current_diam > 7) and (current_diam <= 8):
        return 'slateblue'

    elif (current_diam > 8) and (current_diam <= 9):
        return 'rebeccapurple'

    elif (current_diam > 9) and (current_diam <= 10):
        return 'darkcyan'

    elif (current_diam > 10) and (current_diam <= 15):
        return 'green'

    elif (current_diam > 15) and (current_diam <= 20):
        return 'olive'

    elif (current_diam > 20) and (current_diam <= 25):
        return 'goldenrod'

    elif (current_diam > 25) and (current_diam <= 30):
        return 'darkorchid'

    elif (current_diam > 30) and (current_diam <= 35):
        return 'maroon'

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

                color_to_append = assign_color(current_diam)
                color_arr.append(color_to_append)

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

                    color_to_append = assign_color(current_diam)
                    color_arr.append(color_to_append)

    # convert to numpy arrays so you can do array ops
    density_arr_color = np.asarray(density_arr_color)
    slope_arr_color = np.asarray(slope_arr_color)
    color_arr = np.asarray(color_arr).astype(str)

    # save the arrays
    np.save(slope_extdir + 'density_arr_color.npy', density_arr_color)
    np.save(slope_extdir + 'slope_arr_color.npy', slope_arr_color)
    np.save(slope_extdir + 'color_arr.npy', color_arr)

    print "Arrays saved."

    return density_arr_color, slope_arr_color, color_arr

def get_diam_idx(color_arr, diam_bin):

    if diam_bin == '1to2':
        diam_bin_idx = np.where(color_arr == 'midnightblue')[0]  # '1to2'
    elif diam_bin == '2to3':
        diam_bin_idx = np.where(color_arr == 'blue')[0]  # '2to3'
    elif diam_bin == '3to4':
        diam_bin_idx = np.where(color_arr == 'royalblue')[0]  # '3to4'
    elif diam_bin == '4to5':
        diam_bin_idx = np.where(color_arr == 'dodgerblue')[0]  # '4to5'
    elif diam_bin == '5to6':
        diam_bin_idx = np.where(color_arr == 'deepskyblue')[0]  # '5to6'
    elif diam_bin == '6to7':
        diam_bin_idx = np.where(color_arr == 'steelblue')[0]  # '6to7'
    elif diam_bin == '7to8':
        diam_bin_idx = np.where(color_arr == 'slateblue')[0]  # '7to8'
    elif diam_bin == '8to9':
        diam_bin_idx = np.where(color_arr == 'rebeccapurple')[0]  # '8to9'
    elif diam_bin == '9to10':
        diam_bin_idx = np.where(color_arr == 'darkcyan')[0]  # '9to10'
    elif diam_bin == '10to15':
        diam_bin_idx = np.where(color_arr == 'green')[0]  # '10to15'
    elif diam_bin == '15to20':
        diam_bin_idx = np.where(color_arr == 'olive')[0]  # '15to20'
    elif diam_bin == '20to25':
        diam_bin_idx = np.where(color_arr == 'goldenrod')[0]  # '20to25'
    elif diam_bin == '25to30':
        diam_bin_idx = np.where(color_arr == 'darkorchid')[0]  # '25to30'
    elif diam_bin == '30to35':
        diam_bin_idx = np.where(color_arr == 'maroon')[0]  # '30to35'

    return diam_bin_idx

def make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, diam_bin_min, diam_bin_max):

    # perhaps you could make the blue points bigger than the
    # red points simply because there are fewer blue points.
    # i.e. weighting by the size of the crater.???

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{Slope}$', fontsize=18)
    ax.set_ylabel(r'$\mathrm{Density}$', fontsize=18)

    diam_bin = str(diam_bin_min) + 'to' + str(diam_bin_max)
    diam_bin_idx = get_diam_idx(color_arr, diam_bin)

    ax.scatter(slope_arr_color[diam_bin_idx], density_arr_color[diam_bin_idx], s=5, c=color_arr[diam_bin_idx], alpha=0.4, edgecolors='none')
    ax.set_yscale('log')
    ax.set_ylim(1e-8, 2.0)
    ax.set_xlim(0, 35)

    min_val = 4e6 / (np.pi * (diam_bin_max*1e3)**2)
    max_val = 4e6 / (np.pi * (diam_bin_min*1e3)**2)
    # note that the min and max values come from the biggest and smallest crater respectively

    ax.axhline(y=min_val, ls='--', color='k', lw=1)  # min value of density from biggest crater in bin
    ax.axhline(y=max_val, ls='--', color='k', lw=1)  # max value of density from smallest crater in bin
    # values are obtained for a pixel that is completely inside 
    # a crater and not overlapped by any other craters.

    # plot contours for point density
    # make sure the arrays dont have NaNs
    slope_fin_idx = np.where(np.isfinite(slope_arr_color[diam_bin_idx]))[0]
    density_fin_idx = np.where(np.isfinite(density_arr_color[diam_bin_idx]))[0]
    fin_idx = np.intersect1d(slope_fin_idx, density_fin_idx)

    slope_arr_color_plot = slope_arr_color[diam_bin_idx][fin_idx]
    density_arr_color_plot = density_arr_color[diam_bin_idx][fin_idx]

    counts, xbins, ybins = np.histogram2d(slope_arr_color_plot, density_arr_color_plot, \
        bins=25, normed=False)
    levels_to_plot = [10, 50, 200, 500, 1e3, 2e3, 5e3]
    c = ax.contour(counts.transpose(), levels=levels_to_plot, \
        extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
        colors='lime', linestyles='solid', interpolation='None', zorder=10)
    ax.clabel(c, inline=True, colors='darkgoldenrod', inline_spacing=8, \
        fontsize=6, fontweight='black', fmt='%d', lw=2, ls='-')

    """
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

    fit_y_arr = np.log10(fit_y_arr)

    popt, pcov = curve_fit(poisson, fit_x_arr, fit_y_arr, p0=[1.5])
    print popt
    print pcov

    # plot the best fit curves
    x_plot_arr = np.linspace(0,30,1000)
    ax.plot(x_plot_arr, gb(x_plot_arr), ls='-', color='skyblue', lw=2)
    ax.plot(x_plot_arr, gr(x_plot_arr), ls='-', color='pink', lw=2)
    ax.plot(x_plot_arr, poisson(x_plot_arr, *popt), ls='-', color='forestgreen', lw=2)

    print np.nanmin(density_arr_color[b_idx])
    print np.nanmax(density_arr_color[b_idx])

    check_idx = np.where((density_arr_color[r_idx] >= 1e-3) & (density_arr_color[r_idx] <= np.power(10, -2.9)))[0]
    print check_idx
    print len(check_idx)
    print pix_1d_idx_arr[r_idx][check_idx]
    """

    # add text on figure to indicate diameter bin
    diambinbox = TextArea(str(diam_bin_min) + ' to ' + str(diam_bin_max) + ' km', textprops=dict(color='k', size=10))
    anc_diambinbox = AnchoredOffsetbox(loc=2, child=diambinbox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.75, 0.1),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_diambinbox)

    # add ticks and grid
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

def read_indiv_diambin_crater_frac():

    # read in all invidual diambin crater contributions
    crater_frac_diambin_1_1p25 = np.load(slope_extdir + 'crater_frac_diambin_1_1p25_newbool.npy')
    crater_frac_diambin_1p25_1p5 = np.load(slope_extdir + 'crater_frac_diambin_1p25_1p5_newbool.npy')
    crater_frac_diambin_1p5_1p75 = np.load(slope_extdir + 'crater_frac_diambin_1p5_1p75_newbool.npy')
    crater_frac_diambin_1p75_2 = np.load(slope_extdir + 'crater_frac_diambin_1p75_2_newbool.npy')

    crater_frac_diambin_2_2p25 = np.load(slope_extdir + 'crater_frac_diambin_2_2p25_newbool.npy')
    crater_frac_diambin_2p25_2p5 = np.load(slope_extdir + 'crater_frac_diambin_2p25_2p5_newbool.npy')
    crater_frac_diambin_2p5_2p75 = np.load(slope_extdir + 'crater_frac_diambin_2p5_2p75_newbool.npy')
    crater_frac_diambin_2p75_3 = np.load(slope_extdir + 'crater_frac_diambin_2p75_3_newbool.npy')

    crater_frac_diambin_3_3p25 = np.load(slope_extdir + 'crater_frac_diambin_3_3p25_newbool.npy')
    crater_frac_diambin_3p25_3p5 = np.load(slope_extdir + 'crater_frac_diambin_3p25_3p5_newbool.npy')
    crater_frac_diambin_3p5_3p75 = np.load(slope_extdir + 'crater_frac_diambin_3p5_3p75_newbool.npy')
    crater_frac_diambin_3p75_4 = np.load(slope_extdir + 'crater_frac_diambin_3p75_4_newbool.npy')

    crater_frac_diambin_4_4p25 = np.load(slope_extdir + 'crater_frac_diambin_4_4p25_newbool.npy')
    crater_frac_diambin_4p25_4p5 = np.load(slope_extdir + 'crater_frac_diambin_4p25_4p5_newbool.npy')
    crater_frac_diambin_4p5_4p75 = np.load(slope_extdir + 'crater_frac_diambin_4p5_4p75_newbool.npy')
    crater_frac_diambin_4p75_5 = np.load(slope_extdir + 'crater_frac_diambin_4p75_5_newbool.npy')

    # ---

    crater_frac_diambin_5_6   = np.load(slope_extdir + 'crater_frac_diambin_5_6_newbool.npy')
    crater_frac_diambin_6_7   = np.load(slope_extdir + 'crater_frac_diambin_6_7_newbool.npy')
    crater_frac_diambin_7_8   = np.load(slope_extdir + 'crater_frac_diambin_7_8_newbool.npy')
    crater_frac_diambin_8_9   = np.load(slope_extdir + 'crater_frac_diambin_8_9_newbool.npy')
    crater_frac_diambin_9_10  = np.load(slope_extdir + 'crater_frac_diambin_9_10_newbool.npy')
    crater_frac_diambin_10_15 = np.load(slope_extdir + 'crater_frac_diambin_10_15_newbool.npy')
    crater_frac_diambin_15_20 = np.load(slope_extdir + 'crater_frac_diambin_15_20_newbool.npy')
    crater_frac_diambin_20_25 = np.load(slope_extdir + 'crater_frac_diambin_20_25_newbool.npy')
    crater_frac_diambin_25_30 = np.load(slope_extdir + 'crater_frac_diambin_25_30_newbool.npy')
    crater_frac_diambin_30_35 = np.load(slope_extdir + 'crater_frac_diambin_30_35_newbool.npy')

    return crater_frac_diambin_1_1p25, crater_frac_diambin_1p25_1p5, crater_frac_diambin_1p5_1p75, crater_frac_diambin_1p75_2, \
    crater_frac_diambin_2_2p25, crater_frac_diambin_2p25_2p5, crater_frac_diambin_2p5_2p75, crater_frac_diambin_2p75_3, \
    crater_frac_diambin_3_3p25, crater_frac_diambin_3p25_3p5, crater_frac_diambin_3p5_3p75, crater_frac_diambin_3p75_4, \
    crater_frac_diambin_4_4p25, crater_frac_diambin_4p25_4p5, crater_frac_diambin_4p5_4p75, crater_frac_diambin_4p75_5, \
    crater_frac_diambin_5_6, crater_frac_diambin_6_7, \
    crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, \
    crater_frac_diambin_10_15, crater_frac_diambin_15_20, crater_frac_diambin_20_25, \
    crater_frac_diambin_25_30, crater_frac_diambin_30_35

def get_specific_diambin_and_Nvalue_arrays(slope, pix_frac, crater_ids_arr, crater_diam_m_arr, crater_id_in_pix_arr, start):

    # ----------------------- No Overlap stuff ----------------------- #
    # Create empty lists for each specific diameter bin
    density_diambin_1_1p25 = []
    density_diambin_1p25_1p5 = []
    density_diambin_1p5_1p75 = []
    density_diambin_1p75_2 = []

    density_diambin_2_2p25 = []
    density_diambin_2p25_2p5 = []
    density_diambin_2p5_2p75 = []
    density_diambin_2p75_3 = []

    density_diambin_3_3p25 = []
    density_diambin_3p25_3p5 = []
    density_diambin_3p5_3p75 = []
    density_diambin_3p75_4 = []

    density_diambin_4_4p25 = []
    density_diambin_4p25_4p5 = []
    density_diambin_4p5_4p75 = []
    density_diambin_4p75_5 = []

    # ---

    density_diambin_5_6   = []
    density_diambin_6_7   = []
    density_diambin_7_8   = []
    density_diambin_8_9   = []
    density_diambin_9_10  = []
    density_diambin_10_15 = []
    density_diambin_15_20 = []
    density_diambin_20_25 = []
    density_diambin_25_30 = []
    density_diambin_30_35 = []

    # ----------------------
    slope_diambin_1_1p25 = []
    slope_diambin_1p25_1p5 = []
    slope_diambin_1p5_1p75 = []
    slope_diambin_1p75_2 = []

    slope_diambin_2_2p25 = []
    slope_diambin_2p25_2p5 = []
    slope_diambin_2p5_2p75 = []
    slope_diambin_2p75_3 = []

    slope_diambin_3_3p25 = []
    slope_diambin_3p25_3p5 = []
    slope_diambin_3p5_3p75 = []
    slope_diambin_3p75_4 = []

    slope_diambin_4_4p25 = []
    slope_diambin_4p25_4p5 = []
    slope_diambin_4p5_4p75 = []
    slope_diambin_4p75_5 = []
    # ---

    slope_diambin_5_6   = []
    slope_diambin_6_7   = []
    slope_diambin_7_8   = []
    slope_diambin_8_9   = []
    slope_diambin_9_10  = []
    slope_diambin_10_15 = []
    slope_diambin_15_20 = []
    slope_diambin_20_25 = []
    slope_diambin_25_30 = []
    slope_diambin_30_35 = []

    # ----------------------- N value stuff ----------------------- #
    # create empty lists for slope and density for each Nvalue bin
    density_diambin_1 = []
    density_diambin_1p25 = []
    density_diambin_1p5 = []
    density_diambin_1p75 = []

    density_diambin_2 = []
    density_diambin_2p25 = []
    density_diambin_2p5 = []
    density_diambin_2p75 = []

    density_diambin_3 = []
    density_diambin_3p25 = []
    density_diambin_3p5 = []
    density_diambin_3p75 = []

    density_diambin_4 = []
    density_diambin_4p25 = []
    density_diambin_4p5 = []
    density_diambin_4p75 = []

    # ---

    density_diambin_5 = []
    density_diambin_6 = []
    density_diambin_7 = []
    density_diambin_8 = []
    density_diambin_9 = []
    density_diambin_10 = []
    density_diambin_15 = []
    density_diambin_20 = []
    density_diambin_25 = []
    density_diambin_30 = []

    # ----------------------
    slope_diambin_1 = []
    slope_diambin_1p25 = []
    slope_diambin_1p5 = []
    slope_diambin_1p75 = []

    slope_diambin_2 = []
    slope_diambin_2p25 = []
    slope_diambin_2p5 = []
    slope_diambin_2p75 = []

    slope_diambin_3 = []
    slope_diambin_3p25 = []
    slope_diambin_3p5 = []
    slope_diambin_3p75 = []

    slope_diambin_4 = []
    slope_diambin_4p25 = []
    slope_diambin_4p5 = []
    slope_diambin_4p75 = []

    # ---

    slope_diambin_5 = []
    slope_diambin_6 = []
    slope_diambin_7 = []
    slope_diambin_8 = []
    slope_diambin_9 = []
    slope_diambin_10 = []
    slope_diambin_15 = []
    slope_diambin_20 = []
    slope_diambin_25 = []
    slope_diambin_30 = []

    # read in arrays for crater contributions to specific bins
    crater_frac_diambin_1_1p25, crater_frac_diambin_1p25_1p5, crater_frac_diambin_1p5_1p75, crater_frac_diambin_1p75_2, \
    crater_frac_diambin_2_2p25, crater_frac_diambin_2p25_2p5, crater_frac_diambin_2p5_2p75, crater_frac_diambin_2p75_3, \
    crater_frac_diambin_3_3p25, crater_frac_diambin_3p25_3p5, crater_frac_diambin_3p5_3p75, crater_frac_diambin_3p75_4, \
    crater_frac_diambin_4_4p25, crater_frac_diambin_4p25_4p5, crater_frac_diambin_4p5_4p75, crater_frac_diambin_4p75_5, \
    crater_frac_diambin_5_6, crater_frac_diambin_6_7, \
    crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, \
    crater_frac_diambin_10_15, crater_frac_diambin_15_20, crater_frac_diambin_20_25, \
    crater_frac_diambin_25_30, crater_frac_diambin_30_35 = read_indiv_diambin_crater_frac()

    # create a list as a container for all crater fraction arrays
    allarr = [crater_frac_diambin_1_1p25, crater_frac_diambin_1p25_1p5, crater_frac_diambin_1p5_1p75, crater_frac_diambin_1p75_2, \
    crater_frac_diambin_2_2p25, crater_frac_diambin_2p25_2p5, crater_frac_diambin_2p5_2p75, crater_frac_diambin_2p75_3, \
    crater_frac_diambin_3_3p25, crater_frac_diambin_3p25_3p5, crater_frac_diambin_3p5_3p75, crater_frac_diambin_3p75_4, \
    crater_frac_diambin_4_4p25, crater_frac_diambin_4p25_4p5, crater_frac_diambin_4p5_4p75, crater_frac_diambin_4p75_5, \
    crater_frac_diambin_5_6, crater_frac_diambin_6_7, \
    crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, \
    crater_frac_diambin_10_15, crater_frac_diambin_15_20, crater_frac_diambin_20_25, \
    crater_frac_diambin_25_30, crater_frac_diambin_30_35]

    # create staggered sum arrays (see notes for explanation)
    crater_frac_diambin_1 = sum(allarr)
    crater_frac_diambin_1p25 = sum(allarr[1:])
    crater_frac_diambin_1p5 = sum(allarr[2:])
    crater_frac_diambin_1p75 = sum(allarr[3:])

    crater_frac_diambin_2 = sum(allarr[4:])
    crater_frac_diambin_2p25 = sum(allarr[5:])
    crater_frac_diambin_2p5 = sum(allarr[6:])
    crater_frac_diambin_2p75 = sum(allarr[7:])

    crater_frac_diambin_3 = sum(allarr[8:])
    crater_frac_diambin_3p25 = sum(allarr[9:])
    crater_frac_diambin_3p5 = sum(allarr[10:])
    crater_frac_diambin_3p75 = sum(allarr[11:])

    crater_frac_diambin_4 = sum(allarr[12:])
    crater_frac_diambin_4p25 = sum(allarr[13:])
    crater_frac_diambin_4p5 = sum(allarr[14:])
    crater_frac_diambin_4p75 = sum(allarr[15:])

    # ---

    crater_frac_diambin_5 = sum(allarr[16:])
    crater_frac_diambin_6 = sum(allarr[17:])
    crater_frac_diambin_7 = sum(allarr[18:])
    crater_frac_diambin_8 = sum(allarr[19:])
    crater_frac_diambin_9 = sum(allarr[20:])
    crater_frac_diambin_10 = sum(allarr[21:])
    crater_frac_diambin_15 = sum(allarr[22:])
    crater_frac_diambin_20 = sum(allarr[23:])
    crater_frac_diambin_25 = sum(allarr[24:])
    crater_frac_diambin_30 = sum(allarr[25:])

    # loop over all pixels
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
            current_diam = get_diam(crater_ids_arr, crater_diam_m_arr, current_id)

            if (current_diam > 35):
                continue

            else:
                append_to_density_slope_diambin_lists(current_diam, i, pix_frac, slope, \
                density_diambin_1_1p25, density_diambin_1p25_1p5, density_diambin_1p5_1p75, density_diambin_1p75_2, \
                density_diambin_2_2p25, density_diambin_2p25_2p5, density_diambin_2p5_2p75, density_diambin_2p75_3, \
                density_diambin_3_3p25, density_diambin_3p25_3p5, density_diambin_3p5_3p75, density_diambin_3p75_4, \
                density_diambin_4_4p25, density_diambin_4p25_4p5, density_diambin_4p5_4p75, density_diambin_4p75_5, \
                density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9,\
                density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, density_diambin_25_30,\
                density_diambin_30_35, \
                slope_diambin_1_1p25, slope_diambin_1p25_1p5, slope_diambin_1p5_1p75, slope_diambin_1p75_2, \
                slope_diambin_2_2p25, slope_diambin_2p25_2p5, slope_diambin_2p5_2p75, slope_diambin_2p75_3, \
                slope_diambin_3_3p25, slope_diambin_3p25_3p5, slope_diambin_3p5_3p75, slope_diambin_3p75_4, \
                slope_diambin_4_4p25, slope_diambin_4p25_4p5, slope_diambin_4p5_4p75, slope_diambin_4p75_5, \
                slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, slope_diambin_8_9, \
                slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20,\
                slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35, \
                crater_frac_diambin_1_1p25, crater_frac_diambin_1p25_1p5, crater_frac_diambin_1p5_1p75, crater_frac_diambin_1p75_2, \
                crater_frac_diambin_2_2p25, crater_frac_diambin_2p25_2p5, crater_frac_diambin_2p5_2p75, crater_frac_diambin_2p75_3, \
                crater_frac_diambin_3_3p25, crater_frac_diambin_3p25_3p5, crater_frac_diambin_3p5_3p75, crater_frac_diambin_3p75_4, \
                crater_frac_diambin_4_4p25, crater_frac_diambin_4p25_4p5, crater_frac_diambin_4p5_4p75, crater_frac_diambin_4p75_5, \
                crater_frac_diambin_5_6, crater_frac_diambin_6_7,\
                crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, crater_frac_diambin_10_15,\
                crater_frac_diambin_15_20, crater_frac_diambin_20_25, crater_frac_diambin_25_30, crater_frac_diambin_30_35)

                append_to_density_slope_Nvalue_lists(current_diam, i, pix_frac, slope, \
                density_diambin_1, density_diambin_1p25, density_diambin_1p5, density_diambin_1p75, \
                density_diambin_2, density_diambin_2p25, density_diambin_2p5, density_diambin_2p75, \
                density_diambin_3, density_diambin_3p25, density_diambin_3p5, density_diambin_3p75, \
                density_diambin_4, density_diambin_4p25, density_diambin_4p5, density_diambin_4p75, \
                density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8,\
                density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, density_diambin_25,\
                density_diambin_30, \
                slope_diambin_1, slope_diambin_1p25, slope_diambin_1p5, slope_diambin_1p75, \
                slope_diambin_2, slope_diambin_2p25, slope_diambin_2p5, slope_diambin_2p75, \
                slope_diambin_3, slope_diambin_3p25, slope_diambin_3p5, slope_diambin_3p75, \
                slope_diambin_4, slope_diambin_4p25, slope_diambin_4p5, slope_diambin_4p75, \
                slope_diambin_5, slope_diambin_6, slope_diambin_7, slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15,\
                slope_diambin_20, slope_diambin_25, slope_diambin_30, \
                crater_frac_diambin_1, crater_frac_diambin_1p25, crater_frac_diambin_1p5, crater_frac_diambin_1p75, \
                crater_frac_diambin_2, crater_frac_diambin_2p25, crater_frac_diambin_2p5, crater_frac_diambin_2p75, \
                crater_frac_diambin_3, crater_frac_diambin_3p25, crater_frac_diambin_3p5, crater_frac_diambin_3p75, \
                crater_frac_diambin_4, crater_frac_diambin_4p25, crater_frac_diambin_4p5, crater_frac_diambin_4p75, \
                crater_frac_diambin_5, crater_frac_diambin_6, crater_frac_diambin_7, \
                crater_frac_diambin_8, crater_frac_diambin_9, crater_frac_diambin_10, \
                crater_frac_diambin_15, crater_frac_diambin_20, crater_frac_diambin_25, crater_frac_diambin_30)

        elif len(current_crater_ids) > 1:
            for j in range(len(current_crater_ids)):
                current_id = current_crater_ids[j]
                current_diam = get_diam(crater_ids_arr, crater_diam_m_arr, current_id)

                if (current_diam > 35):
                    continue

                else:
                    append_to_density_slope_diambin_lists(current_diam, i, pix_frac, slope, \
                    density_diambin_1_1p25, density_diambin_1p25_1p5, density_diambin_1p5_1p75, density_diambin_1p75_2, \
                    density_diambin_2_2p25, density_diambin_2p25_2p5, density_diambin_2p5_2p75, density_diambin_2p75_3, \
                    density_diambin_3_3p25, density_diambin_3p25_3p5, density_diambin_3p5_3p75, density_diambin_3p75_4, \
                    density_diambin_4_4p25, density_diambin_4p25_4p5, density_diambin_4p5_4p75, density_diambin_4p75_5, \
                    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9,\
                    density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, density_diambin_25_30,\
                    density_diambin_30_35, \
                    slope_diambin_1_1p25, slope_diambin_1p25_1p5, slope_diambin_1p5_1p75, slope_diambin_1p75_2, \
                    slope_diambin_2_2p25, slope_diambin_2p25_2p5, slope_diambin_2p5_2p75, slope_diambin_2p75_3, \
                    slope_diambin_3_3p25, slope_diambin_3p25_3p5, slope_diambin_3p5_3p75, slope_diambin_3p75_4, \
                    slope_diambin_4_4p25, slope_diambin_4p25_4p5, slope_diambin_4p5_4p75, slope_diambin_4p75_5, \
                    slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, slope_diambin_8_9, \
                    slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20,\
                    slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35, \
                    crater_frac_diambin_1_1p25, crater_frac_diambin_1p25_1p5, crater_frac_diambin_1p5_1p75, crater_frac_diambin_1p75_2, \
                    crater_frac_diambin_2_2p25, crater_frac_diambin_2p25_2p5, crater_frac_diambin_2p5_2p75, crater_frac_diambin_2p75_3, \
                    crater_frac_diambin_3_3p25, crater_frac_diambin_3p25_3p5, crater_frac_diambin_3p5_3p75, crater_frac_diambin_3p75_4, \
                    crater_frac_diambin_4_4p25, crater_frac_diambin_4p25_4p5, crater_frac_diambin_4p5_4p75, crater_frac_diambin_4p75_5, \
                    crater_frac_diambin_5_6, crater_frac_diambin_6_7,\
                    crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, crater_frac_diambin_10_15,\
                    crater_frac_diambin_15_20, crater_frac_diambin_20_25, crater_frac_diambin_25_30, crater_frac_diambin_30_35)

                    append_to_density_slope_Nvalue_lists(current_diam, i, pix_frac, slope, \
                    density_diambin_1, density_diambin_1p25, density_diambin_1p5, density_diambin_1p75, \
                    density_diambin_2, density_diambin_2p25, density_diambin_2p5, density_diambin_2p75, \
                    density_diambin_3, density_diambin_3p25, density_diambin_3p5, density_diambin_3p75, \
                    density_diambin_4, density_diambin_4p25, density_diambin_4p5, density_diambin_4p75, \
                    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8,\
                    density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, density_diambin_25,\
                    density_diambin_30, \
                    slope_diambin_1, slope_diambin_1p25, slope_diambin_1p5, slope_diambin_1p75, \
                    slope_diambin_2, slope_diambin_2p25, slope_diambin_2p5, slope_diambin_2p75, \
                    slope_diambin_3, slope_diambin_3p25, slope_diambin_3p5, slope_diambin_3p75, \
                    slope_diambin_4, slope_diambin_4p25, slope_diambin_4p5, slope_diambin_4p75, \
                    slope_diambin_5, slope_diambin_6, slope_diambin_7, slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15,\
                    slope_diambin_20, slope_diambin_25, slope_diambin_30, \
                    crater_frac_diambin_1, crater_frac_diambin_1p25, crater_frac_diambin_1p5, crater_frac_diambin_1p75, \
                    crater_frac_diambin_2, crater_frac_diambin_2p25, crater_frac_diambin_2p5, crater_frac_diambin_2p75, \
                    crater_frac_diambin_3, crater_frac_diambin_3p25, crater_frac_diambin_3p5, crater_frac_diambin_3p75, \
                    crater_frac_diambin_4, crater_frac_diambin_4p25, crater_frac_diambin_4p5, crater_frac_diambin_4p75, \
                    crater_frac_diambin_5, crater_frac_diambin_6, crater_frac_diambin_7, \
                    crater_frac_diambin_8, crater_frac_diambin_9, crater_frac_diambin_10, \
                    crater_frac_diambin_15, crater_frac_diambin_20, crater_frac_diambin_25, crater_frac_diambin_30)

    # convert to numpy arrays so you can do array ops and save them
    # ----------------------- No Overlap stuff ----------------------- #
    density_diambin_1_1p25 = np.asarray(density_diambin_1_1p25)
    density_diambin_1p25_1p5 = np.asarray(density_diambin_1p25_1p5)
    density_diambin_1p5_1p75 = np.asarray(density_diambin_1p5_1p75)
    density_diambin_1p75_2 = np.asarray(density_diambin_1p75_2)

    density_diambin_2_2p25 = np.asarray(density_diambin_2_2p25)
    density_diambin_2p25_2p5 = np.asarray(density_diambin_2p25_2p5)
    density_diambin_2p5_2p75 = np.asarray(density_diambin_2p5_2p75)
    density_diambin_2p75_3 = np.asarray(density_diambin_2p75_3)

    density_diambin_3_3p25 = np.asarray(density_diambin_3_3p25)
    density_diambin_3p25_3p5 = np.asarray(density_diambin_3p25_3p5)
    density_diambin_3p5_3p75 = np.asarray(density_diambin_3p5_3p75)
    density_diambin_3p75_4 = np.asarray(density_diambin_3p75_4)

    density_diambin_4_4p25 = np.asarray(density_diambin_4_4p25)
    density_diambin_4p25_4p5 = np.asarray(density_diambin_4p25_4p5)
    density_diambin_4p5_4p75 = np.asarray(density_diambin_4p5_4p75)
    density_diambin_4p75_5 = np.asarray(density_diambin_4p75_5)

    # ---

    density_diambin_5_6 = np.asarray(density_diambin_5_6)
    density_diambin_6_7 = np.asarray(density_diambin_6_7)
    density_diambin_7_8 = np.asarray(density_diambin_7_8)
    density_diambin_8_9 = np.asarray(density_diambin_8_9)
    density_diambin_9_10 = np.asarray(density_diambin_9_10)
    density_diambin_10_15 = np.asarray(density_diambin_10_15)
    density_diambin_15_20 = np.asarray(density_diambin_15_20)
    density_diambin_20_25 = np.asarray(density_diambin_20_25)
    density_diambin_25_30 = np.asarray(density_diambin_25_30)
    density_diambin_30_35 = np.asarray(density_diambin_30_35)

    # -----------------------
    slope_diambin_1_1p25 = np.asarray(slope_diambin_1_1p25)
    slope_diambin_1p25_1p5 = np.asarray(slope_diambin_1p25_1p5)
    slope_diambin_1p5_1p75 = np.asarray(slope_diambin_1p5_1p75)
    slope_diambin_1p75_2 = np.asarray(slope_diambin_1p75_2)

    slope_diambin_2_2p25 = np.asarray(slope_diambin_2_2p25)
    slope_diambin_2p25_2p5 = np.asarray(slope_diambin_2p25_2p5)
    slope_diambin_2p5_2p75 = np.asarray(slope_diambin_2p5_2p75)
    slope_diambin_2p75_3 = np.asarray(slope_diambin_2p75_3)

    slope_diambin_3_3p25 = np.asarray(slope_diambin_3_3p25)
    slope_diambin_3p25_3p5 = np.asarray(slope_diambin_3p25_3p5)
    slope_diambin_3p5_3p75 = np.asarray(slope_diambin_3p5_3p75)
    slope_diambin_3p75_4 = np.asarray(slope_diambin_3p75_4)

    slope_diambin_4_4p25 = np.asarray(slope_diambin_4_4p25)
    slope_diambin_4p25_4p5 = np.asarray(slope_diambin_4p25_4p5)
    slope_diambin_4p5_4p75 = np.asarray(slope_diambin_4p5_4p75)
    slope_diambin_4p75_5 = np.asarray(slope_diambin_4p75_5)

    # ---

    slope_diambin_5_6 = np.asarray(slope_diambin_5_6)
    slope_diambin_6_7 = np.asarray(slope_diambin_6_7)
    slope_diambin_7_8 = np.asarray(slope_diambin_7_8)
    slope_diambin_8_9 = np.asarray(slope_diambin_8_9)
    slope_diambin_9_10 = np.asarray(slope_diambin_9_10)
    slope_diambin_10_15 = np.asarray(slope_diambin_10_15)
    slope_diambin_15_20 = np.asarray(slope_diambin_15_20)
    slope_diambin_20_25 = np.asarray(slope_diambin_20_25)
    slope_diambin_25_30 = np.asarray(slope_diambin_25_30)
    slope_diambin_30_35 = np.asarray(slope_diambin_30_35)

    # save the arrays
    np.save(slope_extdir + 'density_diambin_1_1p25.npy', density_diambin_1_1p25)
    np.save(slope_extdir + 'density_diambin_1p25_1p5.npy', density_diambin_1p25_1p5)
    np.save(slope_extdir + 'density_diambin_1p5_1p75.npy', density_diambin_1p5_1p75)
    np.save(slope_extdir + 'density_diambin_1p75_2.npy', density_diambin_1p75_2)

    np.save(slope_extdir + 'density_diambin_2_2p25.npy', density_diambin_2_2p25)
    np.save(slope_extdir + 'density_diambin_2p25_2p5.npy', density_diambin_2p25_2p5)
    np.save(slope_extdir + 'density_diambin_2p5_2p75.npy', density_diambin_2p5_2p75)
    np.save(slope_extdir + 'density_diambin_2p75_3.npy', density_diambin_2p75_3)

    np.save(slope_extdir + 'density_diambin_3_3p25.npy', density_diambin_3_3p25)
    np.save(slope_extdir + 'density_diambin_3p25_3p5.npy', density_diambin_3p25_3p5)
    np.save(slope_extdir + 'density_diambin_3p5_3p75.npy', density_diambin_3p5_3p75)
    np.save(slope_extdir + 'density_diambin_3p75_4.npy', density_diambin_3p75_4)

    np.save(slope_extdir + 'density_diambin_4_4p25.npy', density_diambin_4_4p25)
    np.save(slope_extdir + 'density_diambin_4p25_4p5.npy', density_diambin_4p25_4p5)
    np.save(slope_extdir + 'density_diambin_4p5_4p75.npy', density_diambin_4p5_4p75)
    np.save(slope_extdir + 'density_diambin_4p75_5.npy', density_diambin_4p75_5)

    # ---
    np.save(slope_extdir + 'density_diambin_5_6.npy', density_diambin_5_6)
    np.save(slope_extdir + 'density_diambin_6_7.npy', density_diambin_6_7)
    np.save(slope_extdir + 'density_diambin_7_8.npy', density_diambin_7_8)
    np.save(slope_extdir + 'density_diambin_8_9.npy', density_diambin_8_9)
    np.save(slope_extdir + 'density_diambin_9_10.npy', density_diambin_9_10)
    np.save(slope_extdir + 'density_diambin_10_15.npy', density_diambin_10_15)
    np.save(slope_extdir + 'density_diambin_15_20.npy', density_diambin_15_20)
    np.save(slope_extdir + 'density_diambin_20_25.npy', density_diambin_20_25)
    np.save(slope_extdir + 'density_diambin_25_30.npy', density_diambin_25_30)
    np.save(slope_extdir + 'density_diambin_30_35.npy', density_diambin_30_35)

    # ------------------
    np.save(slope_extdir + 'slope_diambin_1_1p25.npy', slope_diambin_1_1p25)
    np.save(slope_extdir + 'slope_diambin_1p25_1p5.npy', slope_diambin_1p25_1p5)
    np.save(slope_extdir + 'slope_diambin_1p5_1p75.npy', slope_diambin_1p5_1p75)
    np.save(slope_extdir + 'slope_diambin_1p75_2.npy', slope_diambin_1p75_2)

    np.save(slope_extdir + 'slope_diambin_2_2p25.npy', slope_diambin_2_2p25)
    np.save(slope_extdir + 'slope_diambin_2p25_2p5.npy', slope_diambin_2p25_2p5)
    np.save(slope_extdir + 'slope_diambin_2p5_2p75.npy', slope_diambin_2p5_2p75)
    np.save(slope_extdir + 'slope_diambin_2p75_3.npy', slope_diambin_2p75_3)

    np.save(slope_extdir + 'slope_diambin_3_3p25.npy', slope_diambin_3_3p25)
    np.save(slope_extdir + 'slope_diambin_3p25_3p5.npy', slope_diambin_3p25_3p5)
    np.save(slope_extdir + 'slope_diambin_3p5_3p75.npy', slope_diambin_3p5_3p75)
    np.save(slope_extdir + 'slope_diambin_3p75_4.npy', slope_diambin_3p75_4)

    np.save(slope_extdir + 'slope_diambin_4_4p25.npy', slope_diambin_4_4p25)
    np.save(slope_extdir + 'slope_diambin_4p25_4p5.npy', slope_diambin_4p25_4p5)
    np.save(slope_extdir + 'slope_diambin_4p5_4p75.npy', slope_diambin_4p5_4p75)
    np.save(slope_extdir + 'slope_diambin_4p75_5.npy', slope_diambin_4p75_5)

    # ---

    np.save(slope_extdir + 'slope_diambin_5_6.npy', slope_diambin_5_6)
    np.save(slope_extdir + 'slope_diambin_6_7.npy', slope_diambin_6_7)
    np.save(slope_extdir + 'slope_diambin_7_8.npy', slope_diambin_7_8)
    np.save(slope_extdir + 'slope_diambin_8_9.npy', slope_diambin_8_9)
    np.save(slope_extdir + 'slope_diambin_9_10.npy', slope_diambin_9_10)
    np.save(slope_extdir + 'slope_diambin_10_15.npy', slope_diambin_10_15)
    np.save(slope_extdir + 'slope_diambin_15_20.npy', slope_diambin_15_20)
    np.save(slope_extdir + 'slope_diambin_20_25.npy', slope_diambin_20_25)
    np.save(slope_extdir + 'slope_diambin_25_30.npy', slope_diambin_25_30)
    np.save(slope_extdir + 'slope_diambin_30_35.npy', slope_diambin_30_35)

    # ------------------------------ N value stuff ------------------------------ #
    # convert to numpy arrays so you can do array ops and save them
    density_diambin_1 = np.asarray(density_diambin_1)
    density_diambin_1p25 = np.asarray(density_diambin_1p25)
    density_diambin_1p5 = np.asarray(density_diambin_1p5)
    density_diambin_1p75 = np.asarray(density_diambin_1p75)

    density_diambin_2 = np.asarray(density_diambin_2)
    density_diambin_2p25 = np.asarray(density_diambin_2p25)
    density_diambin_2p5 = np.asarray(density_diambin_2p5)
    density_diambin_2p75 = np.asarray(density_diambin_2p75)

    density_diambin_3 = np.asarray(density_diambin_3)
    density_diambin_3p25 = np.asarray(density_diambin_3p25)
    density_diambin_3p5 = np.asarray(density_diambin_3p5)
    density_diambin_3p75 = np.asarray(density_diambin_3p75)

    density_diambin_4 = np.asarray(density_diambin_4)
    density_diambin_4p25 = np.asarray(density_diambin_4p25)
    density_diambin_4p5 = np.asarray(density_diambin_4p5)
    density_diambin_4p75 = np.asarray(density_diambin_4p75)

    # ---

    density_diambin_5 = np.asarray(density_diambin_5)
    density_diambin_6 = np.asarray(density_diambin_6)
    density_diambin_7 = np.asarray(density_diambin_7)
    density_diambin_8 = np.asarray(density_diambin_8)
    density_diambin_9 = np.asarray(density_diambin_9)
    density_diambin_10 = np.asarray(density_diambin_10)
    density_diambin_15 = np.asarray(density_diambin_15)
    density_diambin_20 = np.asarray(density_diambin_20)
    density_diambin_25 = np.asarray(density_diambin_25)
    density_diambin_30 = np.asarray(density_diambin_30)

    # --------------------------
    slope_diambin_1 = np.asarray(slope_diambin_1)
    slope_diambin_1p25 = np.asarray(slope_diambin_1p25)
    slope_diambin_1p5 = np.asarray(slope_diambin_1p5)
    slope_diambin_1p75 = np.asarray(slope_diambin_1p75)

    slope_diambin_2 = np.asarray(slope_diambin_2)
    slope_diambin_2p25 = np.asarray(slope_diambin_2p25)
    slope_diambin_2p5 = np.asarray(slope_diambin_2p5)
    slope_diambin_2p75 = np.asarray(slope_diambin_2p75)

    slope_diambin_3 = np.asarray(slope_diambin_3)
    slope_diambin_3p25 = np.asarray(slope_diambin_3p25)
    slope_diambin_3p5 = np.asarray(slope_diambin_3p5)
    slope_diambin_3p75 = np.asarray(slope_diambin_3p75)

    slope_diambin_4 = np.asarray(slope_diambin_4)
    slope_diambin_4p25 = np.asarray(slope_diambin_4p25)
    slope_diambin_4p5 = np.asarray(slope_diambin_4p5)
    slope_diambin_4p75 = np.asarray(slope_diambin_4p75)

    # ---

    slope_diambin_5 = np.asarray(slope_diambin_5)
    slope_diambin_6 = np.asarray(slope_diambin_6)
    slope_diambin_7 = np.asarray(slope_diambin_7)
    slope_diambin_8 = np.asarray(slope_diambin_8)
    slope_diambin_9 = np.asarray(slope_diambin_9)
    slope_diambin_10 = np.asarray(slope_diambin_10)
    slope_diambin_15 = np.asarray(slope_diambin_15)
    slope_diambin_20 = np.asarray(slope_diambin_20)
    slope_diambin_25 = np.asarray(slope_diambin_25)
    slope_diambin_30 = np.asarray(slope_diambin_30)

    # save the arrays
    np.save(slope_extdir + 'density_diambin_1.npy', density_diambin_1)
    np.save(slope_extdir + 'density_diambin_1p25.npy', density_diambin_1p25)
    np.save(slope_extdir + 'density_diambin_1p5.npy', density_diambin_1p5)
    np.save(slope_extdir + 'density_diambin_1p75.npy', density_diambin_1p75)

    np.save(slope_extdir + 'density_diambin_2.npy', density_diambin_2)
    np.save(slope_extdir + 'density_diambin_2p25.npy', density_diambin_2p25)
    np.save(slope_extdir + 'density_diambin_2p5.npy', density_diambin_2p5)
    np.save(slope_extdir + 'density_diambin_2p75.npy', density_diambin_2p75)

    np.save(slope_extdir + 'density_diambin_3.npy', density_diambin_3)
    np.save(slope_extdir + 'density_diambin_3p25.npy', density_diambin_3p25)
    np.save(slope_extdir + 'density_diambin_3p5.npy', density_diambin_3p5)
    np.save(slope_extdir + 'density_diambin_3p75.npy', density_diambin_3p75)

    np.save(slope_extdir + 'density_diambin_4.npy', density_diambin_4)
    np.save(slope_extdir + 'density_diambin_4p25.npy', density_diambin_4p25)
    np.save(slope_extdir + 'density_diambin_4p5.npy', density_diambin_4p5)
    np.save(slope_extdir + 'density_diambin_4p75.npy', density_diambin_4p75)

    # ---

    np.save(slope_extdir + 'density_diambin_5.npy', density_diambin_5)
    np.save(slope_extdir + 'density_diambin_6.npy', density_diambin_6)
    np.save(slope_extdir + 'density_diambin_7.npy', density_diambin_7)
    np.save(slope_extdir + 'density_diambin_8.npy', density_diambin_8)
    np.save(slope_extdir + 'density_diambin_9.npy', density_diambin_9)
    np.save(slope_extdir + 'density_diambin_10.npy', density_diambin_10)
    np.save(slope_extdir + 'density_diambin_15.npy', density_diambin_15)
    np.save(slope_extdir + 'density_diambin_20.npy', density_diambin_20)
    np.save(slope_extdir + 'density_diambin_25.npy', density_diambin_25)
    np.save(slope_extdir + 'density_diambin_30.npy', density_diambin_30)

    # ------------------------
    np.save(slope_extdir + 'slope_diambin_1.npy', slope_diambin_1)
    np.save(slope_extdir + 'slope_diambin_1p25.npy', slope_diambin_1p25)
    np.save(slope_extdir + 'slope_diambin_1p5.npy', slope_diambin_1p5)
    np.save(slope_extdir + 'slope_diambin_1p75.npy', slope_diambin_1p75)

    np.save(slope_extdir + 'slope_diambin_2.npy', slope_diambin_2)
    np.save(slope_extdir + 'slope_diambin_2p25.npy', slope_diambin_2p25)
    np.save(slope_extdir + 'slope_diambin_2p5.npy', slope_diambin_2p5)
    np.save(slope_extdir + 'slope_diambin_2p75.npy', slope_diambin_2p75)

    np.save(slope_extdir + 'slope_diambin_3.npy', slope_diambin_3)
    np.save(slope_extdir + 'slope_diambin_3p25.npy', slope_diambin_3p25)
    np.save(slope_extdir + 'slope_diambin_3p5.npy', slope_diambin_3p5)
    np.save(slope_extdir + 'slope_diambin_3p75.npy', slope_diambin_3p75)

    np.save(slope_extdir + 'slope_diambin_4.npy', slope_diambin_4)
    np.save(slope_extdir + 'slope_diambin_4p25.npy', slope_diambin_4p25)
    np.save(slope_extdir + 'slope_diambin_4p5.npy', slope_diambin_4p5)
    np.save(slope_extdir + 'slope_diambin_4p75.npy', slope_diambin_4p75)

    # ---

    np.save(slope_extdir + 'slope_diambin_5.npy', slope_diambin_5)
    np.save(slope_extdir + 'slope_diambin_6.npy', slope_diambin_6)
    np.save(slope_extdir + 'slope_diambin_7.npy', slope_diambin_7)
    np.save(slope_extdir + 'slope_diambin_8.npy', slope_diambin_8)
    np.save(slope_extdir + 'slope_diambin_9.npy', slope_diambin_9)
    np.save(slope_extdir + 'slope_diambin_10.npy', slope_diambin_10)
    np.save(slope_extdir + 'slope_diambin_15.npy', slope_diambin_15)
    np.save(slope_extdir + 'slope_diambin_20.npy', slope_diambin_20)
    np.save(slope_extdir + 'slope_diambin_25.npy', slope_diambin_25)
    np.save(slope_extdir + 'slope_diambin_30.npy', slope_diambin_30)

    print "Arrays for no overlap and N value crater fraction saved."

    return None

def append_to_density_slope_diambin_lists(current_diam, pix_idx, pix_frac, slope, \
    density_diambin_1_1p25, density_diambin_1p25_1p5, density_diambin_1p5_1p75, density_diambin_1p75_2, \
    density_diambin_2_2p25, density_diambin_2p25_2p5, density_diambin_2p5_2p75, density_diambin_2p75_3, \
    density_diambin_3_3p25, density_diambin_3p25_3p5, density_diambin_3p5_3p75, density_diambin_3p75_4, \
    density_diambin_4_4p25, density_diambin_4p25_4p5, density_diambin_4p5_4p75, density_diambin_4p75_5, \
    density_diambin_5_6, density_diambin_6_7, density_diambin_7_8, density_diambin_8_9,\
    density_diambin_9_10, density_diambin_10_15, density_diambin_15_20, density_diambin_20_25, density_diambin_25_30,\
    density_diambin_30_35, \
    slope_diambin_1_1p25, slope_diambin_1p25_1p5, slope_diambin_1p5_1p75, slope_diambin_1p75_2, \
    slope_diambin_2_2p25, slope_diambin_2p25_2p5, slope_diambin_2p5_2p75, slope_diambin_2p75_3, \
    slope_diambin_3_3p25, slope_diambin_3p25_3p5, slope_diambin_3p5_3p75, slope_diambin_3p75_4, \
    slope_diambin_4_4p25, slope_diambin_4p25_4p5, slope_diambin_4p5_4p75, slope_diambin_4p75_5, \
    slope_diambin_5_6, slope_diambin_6_7, slope_diambin_7_8, slope_diambin_8_9, \
    slope_diambin_9_10, slope_diambin_10_15, slope_diambin_15_20,\
    slope_diambin_20_25, slope_diambin_25_30, slope_diambin_30_35, \
    crater_frac_diambin_1_1p25, crater_frac_diambin_1p25_1p5, crater_frac_diambin_1p5_1p75, crater_frac_diambin_1p75_2, \
    crater_frac_diambin_2_2p25, crater_frac_diambin_2p25_2p5, crater_frac_diambin_2p5_2p75, crater_frac_diambin_2p75_3, \
    crater_frac_diambin_3_3p25, crater_frac_diambin_3p25_3p5, crater_frac_diambin_3p5_3p75, crater_frac_diambin_3p75_4, \
    crater_frac_diambin_4_4p25, crater_frac_diambin_4p25_4p5, crater_frac_diambin_4p5_4p75, crater_frac_diambin_4p75_5, \
    crater_frac_diambin_5_6, crater_frac_diambin_6_7,\
    crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, crater_frac_diambin_10_15,\
    crater_frac_diambin_15_20, crater_frac_diambin_20_25, crater_frac_diambin_25_30, crater_frac_diambin_30_35):

    if current_diam >= 1.0 and current_diam < 1.25:
        density_diambin_1_1p25.append(crater_frac_diambin_1_1p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1_1p25.append(slope[pix_idx])
    elif current_diam >= 1.25 and current_diam < 1.5:
        density_diambin_1p25_1p5.append(crater_frac_diambin_1p25_1p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1p25_1p5.append(slope[pix_idx])
    elif current_diam >= 1.5 and current_diam < 1.75:
        density_diambin_1p5_1p75.append(crater_frac_diambin_1p5_1p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1p5_1p75.append(slope[pix_idx])
    elif current_diam >= 1.75 and current_diam < 2.0:
        density_diambin_1p75_2.append(crater_frac_diambin_1p75_2[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1p75_2.append(slope[pix_idx])

    elif current_diam >= 2.0 and current_diam < 2.25:
        density_diambin_2_2p25.append(crater_frac_diambin_2_2p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2_2p25.append(slope[pix_idx])
    elif current_diam >= 2.25 and current_diam < 2.5:
        density_diambin_2p25_2p5.append(crater_frac_diambin_2p25_2p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2p25_2p5.append(slope[pix_idx])
    elif current_diam >= 2.5 and current_diam < 2.75:
        density_diambin_2p5_2p75.append(crater_frac_diambin_2p5_2p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2p5_2p75.append(slope[pix_idx])
    elif current_diam >= 2.75 and current_diam < 3.0:
        density_diambin_2p75_3.append(crater_frac_diambin_2p75_3[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2p75_3.append(slope[pix_idx])

    elif current_diam >= 3.0 and current_diam < 3.25:
        density_diambin_3_3p25.append(crater_frac_diambin_3_3p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3_3p25.append(slope[pix_idx])
    elif current_diam >= 3.25 and current_diam < 3.5:
        density_diambin_3p25_3p5.append(crater_frac_diambin_3p25_3p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3p25_3p5.append(slope[pix_idx])
    elif current_diam >= 3.5 and current_diam < 3.75:
        density_diambin_3p5_3p75.append(crater_frac_diambin_3p5_3p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3p5_3p75.append(slope[pix_idx])
    elif current_diam >= 3.75 and current_diam < 4.0:
        density_diambin_3p75_4.append(crater_frac_diambin_3p75_4[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3p75_4.append(slope[pix_idx])
    
    elif current_diam >= 4.0 and current_diam < 4.25:
        density_diambin_4_4p25.append(crater_frac_diambin_4_4p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4_4p25.append(slope[pix_idx])
    elif current_diam >= 4.25 and current_diam < 4.5:
        density_diambin_4p25_4p5.append(crater_frac_diambin_4p25_4p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4p25_4p5.append(slope[pix_idx])
    elif current_diam >= 4.5 and current_diam < 4.75:
        density_diambin_4p5_4p75.append(crater_frac_diambin_4p5_4p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4p5_4p75.append(slope[pix_idx])
    elif current_diam >= 4.75 and current_diam < 5.0:
        density_diambin_4p75_5.append(crater_frac_diambin_4p75_5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4p75_5.append(slope[pix_idx])

    # ---
    
    elif current_diam >= 5.0 and current_diam < 6.0:
        density_diambin_5_6.append(crater_frac_diambin_5_6[pix_idx] / pix_frac[pix_idx])
        slope_diambin_5_6.append(slope[pix_idx])
    
    elif current_diam >= 6.0 and current_diam < 7.0:
        density_diambin_6_7.append(crater_frac_diambin_6_7[pix_idx] / pix_frac[pix_idx])
        slope_diambin_6_7.append(slope[pix_idx])
    
    elif current_diam >= 7.0 and current_diam < 8.0:
        density_diambin_7_8.append(crater_frac_diambin_7_8[pix_idx] / pix_frac[pix_idx])
        slope_diambin_7_8.append(slope[pix_idx])
    
    elif current_diam >= 8.0 and current_diam < 9.0:
        density_diambin_8_9.append(crater_frac_diambin_8_9[pix_idx] / pix_frac[pix_idx])
        slope_diambin_8_9.append(slope[pix_idx])
    
    elif current_diam >= 9.0 and current_diam < 10.0:
        density_diambin_9_10.append(crater_frac_diambin_9_10[pix_idx] / pix_frac[pix_idx])
        slope_diambin_9_10.append(slope[pix_idx])

    elif current_diam >= 10.0 and current_diam < 15.0:
        density_diambin_10_15.append(crater_frac_diambin_10_15[pix_idx] / pix_frac[pix_idx])
        slope_diambin_10_15.append(slope[pix_idx])

    elif current_diam >= 15.0 and current_diam < 20.0:
        density_diambin_15_20.append(crater_frac_diambin_15_20[pix_idx] / pix_frac[pix_idx])
        slope_diambin_15_20.append(slope[pix_idx])

    elif current_diam >= 20.0 and current_diam < 25.0:
        density_diambin_20_25.append(crater_frac_diambin_20_25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_20_25.append(slope[pix_idx])
    
    elif current_diam >= 25.0 and current_diam < 30.0:
        density_diambin_25_30.append(crater_frac_diambin_25_30[pix_idx] / pix_frac[pix_idx])
        slope_diambin_25_30.append(slope[pix_idx])
    
    elif current_diam >= 30.0 and current_diam < 35.0:
        density_diambin_30_35.append(crater_frac_diambin_30_35[pix_idx] / pix_frac[pix_idx])
        slope_diambin_30_35.append(slope[pix_idx])

    return None

def plot_by_diam(slope_arr, pix_frac, start):

    crater_ids_arr, crater_diam_m_arr, crater_id_in_pix_arr = get_ids()

    # ------------ for cumulative numbers ------------ # 
    # get arrays where the crater diam has been identified by color
    #density_arr_color, slope_arr_color, color_arr = \
    #get_diam_ref_arrays(density, slope_arr, crater_vert_cat, crater_id_in_pix_arr, start)

    # ------------ to get crater contributions included only from within diambin ------------ # 
    get_specific_diambin_and_Nvalue_arrays(slope_arr, pix_frac, crater_ids_arr, crater_diam_m_arr, crater_id_in_pix_arr, start)
    return None

    # do the actual plotting
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 1, 2)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 2, 3)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 3, 4)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 4, 5)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 5, 6)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 6, 7)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 7, 8)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 8, 9)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 9, 10)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 10, 15)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 15, 20)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 20, 25)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 25, 30)
    make_plot_diam_bin(density_arr_color, slope_arr_color, color_arr, 30, 35)

    # Put all together
    # ----------------------------- GRID PLOT ----------------------------- #
    gs = gridspec.GridSpec(4,4)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.02, hspace=0.02)

    fig = plt.figure()
    ax_1to2   = fig.add_subplot(gs[0, 0])
    ax_2to3   = fig.add_subplot(gs[0, 1])
    ax_3to4   = fig.add_subplot(gs[0, 2])
    ax_4to5   = fig.add_subplot(gs[0, 3])
    ax_5to6   = fig.add_subplot(gs[1, 0])
    ax_6to7   = fig.add_subplot(gs[1, 1])
    ax_7to8   = fig.add_subplot(gs[1, 2])
    ax_8to9   = fig.add_subplot(gs[1, 3])
    ax_9to10  = fig.add_subplot(gs[2, 0])
    ax_10to15 = fig.add_subplot(gs[2, 1])
    ax_15to20 = fig.add_subplot(gs[2, 2])
    ax_20to25 = fig.add_subplot(gs[2, 3])
    ax_25to30 = fig.add_subplot(gs[3, 1])
    ax_30to35 = fig.add_subplot(gs[3, 2])

    all_axes = [ax_1to2, ax_2to3, ax_3to4, ax_4to5, ax_5to6, ax_6to7, ax_7to8, ax_8to9, ax_9to10, ax_10to15, ax_15to20, ax_20to25, ax_25to30, ax_30to35]
    diam_bins = ['1to2', '2to3', '3to4', '4to5', '5to6', '6to7', '7to8', '8to9', '9to10', '10to15', '15to20', '20to25', '25to30', '30to35']
    all_diam_idx = []
    min_val_list = []
    max_val_list = []
    for j in range(len(all_axes)):
        diam_bin_idx = get_diam_idx(color_arr, diam_bins[j])
        all_diam_idx.append(diam_bin_idx)

        diam_bin_min = int(diam_bins[j].split('to')[0])
        diam_bin_max = int(diam_bins[j].split('to')[1])
        min_val_list.append(4e6 / (np.pi * (diam_bin_max*1e3)**2))
        max_val_list.append(4e6 / (np.pi * (diam_bin_min*1e3)**2))

    ax_25to30.set_xlabel(r'$\mathrm{Slope}$', fontsize=15)
    ax_25to30.xaxis.set_label_coords(1.05, -0.25)

    ax_9to10.set_ylabel(r'$\mathrm{Density}$', fontsize=15)
    ax_9to10.yaxis.set_label_coords(-0.3, 1.05)

    for i in range(len(all_axes)):

        all_axes[i].scatter(slope_arr_color[all_diam_idx[i]], density_arr_color[all_diam_idx[i]], s=5, c=color_arr[all_diam_idx[i]], alpha=0.4, edgecolors='none')
        all_axes[i].set_yscale('log')
        all_axes[i].set_ylim(1e-8, 2.0)
        all_axes[i].set_xlim(0, 35)

        # plot contours for point density
        # make sure the arrays dont have NaNs
        slope_fin_idx = np.where(np.isfinite(slope_arr_color[all_diam_idx[i]]))[0]
        density_fin_idx = np.where(np.isfinite(density_arr_color[all_diam_idx[i]]))[0]
        fin_idx = np.intersect1d(slope_fin_idx, density_fin_idx)

        slope_arr_color_plot = slope_arr_color[all_diam_idx[i]][fin_idx]
        density_arr_color_plot = density_arr_color[all_diam_idx[i]][fin_idx]

        counts, xbins, ybins = np.histogram2d(slope_arr_color_plot, density_arr_color_plot, \
            bins=25, normed=False)
        levels_to_plot = [10, 50, 200, 500, 1e3, 2e3, 5e3]
        c = all_axes[i].contour(counts.transpose(), levels=levels_to_plot, \
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], \
            colors='lime', linestyles='solid', interpolation='None', zorder=10)
        all_axes[i].clabel(c, inline=True, colors='darkgoldenrod', inline_spacing=8, \
            fontsize=4, fontweight='black', fmt='%d', lw=2, ls='-')

        all_axes[i].axhline(y=min_val_list[i], ls='--', color='k', lw=1)  # min value of density from biggest crater in bin
        all_axes[i].axhline(y=max_val_list[i], ls='--', color='k', lw=1)  # max value of density from smallest crater in bin

        all_axes[i].minorticks_on()
        all_axes[i].tick_params('both', width=1, length=3, which='minor', labelsize=6)
        all_axes[i].tick_params('both', width=1, length=4.7, which='major', labelsize=6)

        # add text on figure to indicate diameter bin
        diam_bin_min = int(diam_bins[i].split('to')[0])
        diam_bin_max = int(diam_bins[i].split('to')[1])
        diambinbox = TextArea(str(diam_bin_min) + ' to ' + str(diam_bin_max) + ' km', textprops=dict(color='k', size=5))
        anc_diambinbox = AnchoredOffsetbox(loc=2, child=diambinbox, pad=0.0, frameon=False,\
                                             bbox_to_anchor=(0.65, 0.16),\
                                             bbox_transform=all_axes[i].transAxes, borderpad=0.0)
        all_axes[i].add_artist(anc_diambinbox)

        if i <= 11:
            all_axes[i].set_xticklabels([])

        if i == 1 or i == 2 or i == 3 or i == 5 or i == 6 or i == 7 or i == 9 or i == 10 or i == 11 or i == 13:
            all_axes[i].set_yticklabels([])

    ax_9to10.set_xticklabels(['0', '10', '20', ''])
    ax_20to25.set_xticklabels(['', '10', '20', '30'])
    ax_25to30.set_xticklabels(['0', '10', '20', '30'])
    ax_30to35.set_xticklabels(['0', '10', '20', '30'])

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

        ax.scatter(slope_arr_color[all_diam_idx[i]], density_arr_color[all_diam_idx[i]], s=1, c=color_arr[all_diam_idx[i]], alpha=0.4, edgecolors='none')

    ax.set_yscale('log')
    ax.set_ylim(1e-8, 2.0)
    ax.set_xlim(0, 35)
    
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

def append_to_density_slope_Nvalue_lists(current_diam, pix_idx, pix_frac, slope, \
    density_diambin_1, density_diambin_1p25, density_diambin_1p5, density_diambin_1p75, \
    density_diambin_2, density_diambin_2p25, density_diambin_2p5, density_diambin_2p75, \
    density_diambin_3, density_diambin_3p25, density_diambin_3p5, density_diambin_3p75, \
    density_diambin_4, density_diambin_4p25, density_diambin_4p5, density_diambin_4p75, \
    density_diambin_5, density_diambin_6, density_diambin_7, density_diambin_8,\
    density_diambin_9, density_diambin_10, density_diambin_15, density_diambin_20, density_diambin_25,\
    density_diambin_30, \
    slope_diambin_1, slope_diambin_1p25, slope_diambin_1p5, slope_diambin_1p75, \
    slope_diambin_2, slope_diambin_2p25, slope_diambin_2p5, slope_diambin_2p75, \
    slope_diambin_3, slope_diambin_3p25, slope_diambin_3p5, slope_diambin_3p75, \
    slope_diambin_4, slope_diambin_4p25, slope_diambin_4p5, slope_diambin_4p75, \
    slope_diambin_5, slope_diambin_6, slope_diambin_7, slope_diambin_8, slope_diambin_9, slope_diambin_10, slope_diambin_15,\
    slope_diambin_20, slope_diambin_25, slope_diambin_30, \
    crater_frac_diambin_1, crater_frac_diambin_1p25, crater_frac_diambin_1p5, crater_frac_diambin_1p75, \
    crater_frac_diambin_2, crater_frac_diambin_2p25, crater_frac_diambin_2p5, crater_frac_diambin_2p75, \
    crater_frac_diambin_3, crater_frac_diambin_3p25, crater_frac_diambin_3p5, crater_frac_diambin_3p75, \
    crater_frac_diambin_4, crater_frac_diambin_4p25, crater_frac_diambin_4p5, crater_frac_diambin_4p75, \
    crater_frac_diambin_5, crater_frac_diambin_6, crater_frac_diambin_7, \
    crater_frac_diambin_8, crater_frac_diambin_9, crater_frac_diambin_10, \
    crater_frac_diambin_15, crater_frac_diambin_20, crater_frac_diambin_25, crater_frac_diambin_30):

    if current_diam >= 1.0:
        density_diambin_1.append(crater_frac_diambin_1[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1.append(slope[pix_idx])
    if current_diam >= 1.25:
        density_diambin_1p25.append(crater_frac_diambin_1p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1p25.append(slope[pix_idx])
    if current_diam >= 1.5:
        density_diambin_1p5.append(crater_frac_diambin_1p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1p5.append(slope[pix_idx])
    if current_diam >= 1.75:
        density_diambin_1p75.append(crater_frac_diambin_1p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_1p75.append(slope[pix_idx])

    if current_diam >= 2.0:
        density_diambin_2.append(crater_frac_diambin_2[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2.append(slope[pix_idx])
    if current_diam >= 2.25:
        density_diambin_2p25.append(crater_frac_diambin_2p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2p25.append(slope[pix_idx])
    if current_diam >= 2.5:
        density_diambin_2p5.append(crater_frac_diambin_2p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2p5.append(slope[pix_idx])
    if current_diam >= 2.75:
        density_diambin_2p75.append(crater_frac_diambin_2p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_2p75.append(slope[pix_idx])

    if current_diam >= 3.0:
        density_diambin_3.append(crater_frac_diambin_3[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3.append(slope[pix_idx])
    if current_diam >= 3.25:
        density_diambin_3p25.append(crater_frac_diambin_3p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3p25.append(slope[pix_idx])
    if current_diam >= 3.5:
        density_diambin_3p5.append(crater_frac_diambin_3p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3p5.append(slope[pix_idx])
    if current_diam >= 3.75:
        density_diambin_3p75.append(crater_frac_diambin_3p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_3p75.append(slope[pix_idx])

    if current_diam >= 4.0:
        density_diambin_4.append(crater_frac_diambin_4[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4.append(slope[pix_idx])
    if current_diam >= 4.25:
        density_diambin_4p25.append(crater_frac_diambin_4p25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4p25.append(slope[pix_idx])
    if current_diam >= 4.5:
        density_diambin_4p5.append(crater_frac_diambin_4p5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4p5.append(slope[pix_idx])
    if current_diam >= 4.75:
        density_diambin_4p75.append(crater_frac_diambin_4p75[pix_idx] / pix_frac[pix_idx])
        slope_diambin_4p75.append(slope[pix_idx])

    # ---

    if current_diam >= 5.0:
        density_diambin_5.append(crater_frac_diambin_5[pix_idx] / pix_frac[pix_idx])
        slope_diambin_5.append(slope[pix_idx])
    
    if current_diam >= 6.0:
        density_diambin_6.append(crater_frac_diambin_6[pix_idx] / pix_frac[pix_idx])
        slope_diambin_6.append(slope[pix_idx])
    
    if current_diam >= 7.0:
        density_diambin_7.append(crater_frac_diambin_7[pix_idx] / pix_frac[pix_idx])
        slope_diambin_7.append(slope[pix_idx])
    
    if current_diam >= 8.0:
        density_diambin_8.append(crater_frac_diambin_8[pix_idx] / pix_frac[pix_idx])
        slope_diambin_8.append(slope[pix_idx])
    
    if current_diam >= 9.0:
        density_diambin_9.append(crater_frac_diambin_9[pix_idx] / pix_frac[pix_idx])
        slope_diambin_9.append(slope[pix_idx])

    if current_diam >= 10.0:
        density_diambin_10.append(crater_frac_diambin_10[pix_idx] / pix_frac[pix_idx])
        slope_diambin_10.append(slope[pix_idx])

    if current_diam >= 15.0:
        density_diambin_15.append(crater_frac_diambin_15[pix_idx] / pix_frac[pix_idx])
        slope_diambin_15.append(slope[pix_idx])

    if current_diam >= 20.0:
        density_diambin_20.append(crater_frac_diambin_20[pix_idx] / pix_frac[pix_idx])
        slope_diambin_20.append(slope[pix_idx])
    
    if current_diam >= 25.0:
        density_diambin_25.append(crater_frac_diambin_25[pix_idx] / pix_frac[pix_idx])
        slope_diambin_25.append(slope[pix_idx])
    
    if current_diam >= 30.0:
        density_diambin_30.append(crater_frac_diambin_30[pix_idx] / pix_frac[pix_idx])
        slope_diambin_30.append(slope[pix_idx])

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

    plot_by_diam(slope_arr, pix_frac, start)
    #plot_3d_hist(density, slope_arr)

    # total run time
    print '\n'
    print "Total time taken --", (time.time() - start)/60, "minutes."
    sys.exit(0)
