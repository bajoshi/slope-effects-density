from __future__ import division

import numpy as np

import sys
import os

home = os.getenv('HOME')
slopedir = home + '/Desktop/slope-effects-density/'
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

sys.path.append(slopedir)
import plot_maps as pm
import slope_utils as su
import slope_effects_density as se

if __name__ == '__main__':

    # read the crater frac files
    crater_frac_diambin_1_2, crater_frac_diambin_2_3, crater_frac_diambin_3_4, \
        crater_frac_diambin_4_5, crater_frac_diambin_5_6, crater_frac_diambin_6_7, \
        crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, \
        crater_frac_diambin_10_15, crater_frac_diambin_15_20, crater_frac_diambin_20_25, \
        crater_frac_diambin_25_30, crater_frac_diambin_30_35 = pm.read_indiv_diambin_crater_frac()

    allarr = [crater_frac_diambin_1_2, crater_frac_diambin_2_3, crater_frac_diambin_3_4, \
        crater_frac_diambin_4_5, crater_frac_diambin_5_6, crater_frac_diambin_6_7, \
        crater_frac_diambin_7_8, crater_frac_diambin_8_9, crater_frac_diambin_9_10, \
        crater_frac_diambin_10_15, crater_frac_diambin_15_20, crater_frac_diambin_20_25, \
        crater_frac_diambin_25_30, crater_frac_diambin_30_35]

    # create N value arrays
    crater_frac_diambin_1 = sum(allarr)

    # convert to fits so you can check with ds9
    pm.convert_npy_array_tofits(crater_frac_diambin_1, (2109,1949), slope_extdir, 'crater_frac_diambin_1')
    pm.convert_npy_array_tofits(crater_frac_diambin_1_2, (2109,1949), slope_extdir, 'crater_frac_diambin_1_2')
    pm.convert_npy_array_tofits(crater_frac_diambin_2_3, (2109,1949), slope_extdir, 'crater_frac_diambin_2_3')
    pm.convert_npy_array_tofits(crater_frac_diambin_3_4, (2109,1949), slope_extdir, 'crater_frac_diambin_3_4')
    pm.convert_npy_array_tofits(crater_frac_diambin_4_5, (2109,1949), slope_extdir, 'crater_frac_diambin_4_5')
    pm.convert_npy_array_tofits(crater_frac_diambin_5_6, (2109,1949), slope_extdir, 'crater_frac_diambin_5_6')

    # convert to ascii raster after check with ds9
    # get supplementary data first
    # read in pix frac array
    pix_frac = np.load(slope_extdir + 'pix_area_fraction_clipped.npy')
    pix_frac = pix_frac.ravel()

    # get pixel x and y arrays and rows and cols
    nrows = 2109
    ncols = 1949
    pix_x_cen_arr = np.arange(-3822217, -1874217 + 1000, 1000)
    pix_y_cen_arr = np.ones(ncols) * 380491
    for count in range(1,nrows):
        pix_x_cen_arr = np.append(pix_x_cen_arr, np.arange(-3822217, -1874217 + 1000, 1000))
        pix_y_cen_arr = np.append(pix_y_cen_arr, np.ones(ncols) * (380491 - 1000*count))

    rows, columns = se.get_rows_columns(pix_x_cen_arr, pix_y_cen_arr)

    # divide crater frac arrays by pixel frac to get densities
    density_2darray_diambin_1 = crater_frac_diambin_1 / pix_frac
    density_2darray_diambin_1_2 = crater_frac_diambin_1_2 / pix_frac
    density_2darray_diambin_2_3 = crater_frac_diambin_2_3 / pix_frac
    density_2darray_diambin_3_4 = crater_frac_diambin_3_4 / pix_frac
    density_2darray_diambin_4_5 = crater_frac_diambin_4_5 / pix_frac
    density_2darray_diambin_5_6 = crater_frac_diambin_5_6 / pix_frac
    
    # save N(1) crater fraction array and all density arrays
    np.save(slope_extdir + 'crater_frac_diambin_1.npy', crater_frac_diambin_1)
    np.save(slope_extdir + 'density_2darray_diambin_1.npy', density_2darray_diambin_1)
    np.save(slope_extdir + 'density_2darray_diambin_1_2.npy', density_2darray_diambin_1_2)
    np.save(slope_extdir + 'density_2darray_diambin_2_3.npy', density_2darray_diambin_2_3)
    np.save(slope_extdir + 'density_2darray_diambin_3_4.npy', density_2darray_diambin_3_4)
    np.save(slope_extdir + 'density_2darray_diambin_4_5.npy', density_2darray_diambin_4_5)
    np.save(slope_extdir + 'density_2darray_diambin_5_6.npy', density_2darray_diambin_5_6)    

    # now convert densities to ascii raster
    su.numpy_to_asciiraster(slope_extdir + 'density_2darray_diambin_1.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)
    su.numpy_to_asciiraster(slope_extdir + 'density_2darray_diambin_1_2.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)
    su.numpy_to_asciiraster(slope_extdir + 'density_2darray_diambin_2_3.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)
    su.numpy_to_asciiraster(slope_extdir + 'density_2darray_diambin_3_4.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)
    su.numpy_to_asciiraster(slope_extdir + 'density_2darray_diambin_4_5.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)
    su.numpy_to_asciiraster(slope_extdir + 'density_2darray_diambin_5_6.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)

    sys.exit(0)