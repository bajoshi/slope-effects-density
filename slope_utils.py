from __future__ import division  

import numpy as np

import sys
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def numpy_to_asciiraster(npy_array_path, final_shape, pix_x, pix_y, dir_for_raster='same'):

    # read in numpy array
    npy_array = np.load(npy_array_path)

    # decide where to save raster and create raster file
    if dir_for_raster == 'same':
        dir_for_raster = os.path.dirname(npy_array_path) + '/'

    npy_filename = os.path.splitext(os.path.basename(npy_array_path))[0]
    raster_file = open(dir_for_raster + '/' + npy_filename + '.asc', 'wa')

    # determine rows and columns and reshape numpy 1D array to 2D image
    rows = final_shape[0]
    cols = final_shape[1]

    # reshape numpy array to be a 2D array for the rest of the code to work correctly
    npy_array = npy_array.reshape((rows, cols))

    # get pixel values
    pix_x_ll = pix_x[0]
    pix_y_ll = pix_y[-1]  # these indices are the way they are because I know how the coord arrays are shaped
    cellsize = 1000  # pixel step in meters

    # write array to file
    raster_file.write('NCOLS ' + str(cols) + '\n')
    raster_file.write('NROWS ' + str(rows) + '\n')

    raster_file.write('XLLCENTER ' + str(pix_x_ll) + '\n')
    raster_file.write('YLLCENTER ' + str(pix_y_ll) + '\n')

    raster_file.write('CELLSIZE ' + str(cellsize) + '\n')
    raster_file.write('NODATA_VALUE ' + str(-9999.0) + '\n')

    # loop over all rows and write row data to raster file
    for i in range(rows):
        for j in range(cols):
        
            raster_file.write(str(npy_array[i,j]) + ' ')
        raster_file.write('\n')

    # close file to actually write it to disk and return
    raster_file.close()

    return None

def raster_to_numpy(raster_path, total_pixels, dir_for_array='same'):

    # read in raster file
    # opened only in read mode to stop from accidentally messing with it
    raster_file = open(raster_path, 'r')

    # decide where to save array
    if dir_for_array == 'same':
        dir_for_array = os.path.dirname(raster_path) + '/'

    array_filename = os.path.splitext(os.path.basename(raster_path))

    npy_array = np.empty(total_pixels)

    # convert data in raster to 1D numpy array
    for line in raster_file.readlines()[6:]:
        # MAKE SURE THAT THE RASTER FILE ACTUALLY DOES HAVE 6 
        # LINES IN ITS HEADER AND THAT NO ROWS ARE BEING SKIPPED.
        
        col_arr = np.array(float(line.split(' ')))
        np.append(npy_array, col_arr)

    # write file as numpy binary file
    np.save(dir_for_array + array_filename + '.npy', npy_array)

    # close file and return
    raster_file.close()

    return None

def blank_raster(final_shape, pix_x, pix_y, savedir):
    """
    This function will create a blank ascii 
    raster file with all ones in it.
    It needs the shape and the pixel coordinates.
    It also needs to know the directory path 
    for saving the raster.
    """

    # create raster file
    if savedir[-1] != '/':
        raster_file = open(savedir + '/' + 'blank_raster.asc', 'wa')
    elif savedir[-1] == '/':
        raster_file = open(savedir + 'blank_raster.asc', 'wa')

    # determine rows and columns
    rows = final_shape[0]
    cols = final_shape[1]

    # get pixel values
    pix_x_ll = pix_x[0]
    pix_y_ll = pix_y[-1]  # these indices are the way they are because I know how the coord arrays are shaped
    cellsize = 1000  # pixel step in meters

    # write array to file
    raster_file.write('NCOLS ' + str(cols) + '\n')
    raster_file.write('NROWS ' + str(rows) + '\n')

    raster_file.write('XLLCENTER ' + str(pix_x_ll) + '\n')
    raster_file.write('YLLCENTER ' + str(pix_y_ll) + '\n')

    raster_file.write('CELLSIZE ' + str(cellsize) + '\n')
    raster_file.write('NODATA_VALUE ' + str(-9999.0) + '\n')

    # loop over all rows and write row data to raster file
    for i in range(rows):
        for j in range(cols):
        
            raster_file.write(str(1.0) + ' ')
        raster_file.write('\n')

    # close file to actually write it to disk and return
    raster_file.close()

    return None

if __name__ == '__main__':
    
    print "\n"
    print "This module contains the following utility functions:"
    print "1. numpy_to_asciiraster: To convert 1D numpy arrays which are saved as "\
    "numpy binary files (i.e. with a .npy extension) to 2D ascii raster images."
    print "2. asciiraster_to_numpy: Inverse of the previous conversion."
    print "\nYou have to import this code as a module to use the functions."
    print "E.g. import slope_utils as su"
    print "Make sure that this module lies on your PYTHONPATH."
    print "Check documentation for individual functions for more details."
    print "\n"

    sys.exit(0)