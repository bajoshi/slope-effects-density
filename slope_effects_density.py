"""
This code will take arrays of veritces (i.e. pixel coords of vertices) defining the study area
and crater info (i.e. crater diam and crater center coords) and figure out the fraction of each 
pixel inside the study area and the summed fraction of crater area (i.e. summing over each crater 
that intersects a given pixel and giving a cumulative sum ) within each pixel.
"""

# To do list: 
# 0 - Give boolean option to user
# DONE# 1 - Convert and save csv
# DONE # 2 - Convert and save density numpy into ASCII Raster (or other raster type)
# DONE # 3 - Convert clipped ASCII Raster back to numpy
# 4 - Create 3D histogram from new clipped numpy
# DONE # 5 - Write function to read in crater vertices
# DONE # Probably easiest to put numpy2raster and raster2numpy in a separate utilities code and call here as needed.

from __future__ import division  # __future__ calls most recent version of a library or commands from python

import numpy as np
import Polygon as pg
from Polygon.Shapes import Circle

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import slope_utils as su

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

def check_point_polygon(x, y, poly):
    """
    This function, which is a very important one, came from a 
    solution online. I simply copy pasted it.
    It checks whether the supplied coordinates of a point are within the area
    enclosed by the supplied polygon vertices and returns True or False.
    The algorithm is called the Ray Casting algorithm. 
    """

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def part_of_old_main():
    """
    This is part of the main code that I was writing earlier,
    which was to be used with the get_pixel_area() function and others.
    I've kept it here so that I can come back to it should I ever need to.
    """

    for i in range(len(pix_centers)):

        # check if pixel center falls "well" inside the inner excluding rectangle
        if (min(in_rect_x) + 500 < pix_x_cen_arr[i]) and (pix_x_cen_arr[i] < max(in_rect_x) - 500) and \
        (min(in_rect_y) + 500 < pix_y_cen_arr[i]) and (pix_y_cen_arr[i] < max(in_rect_y) - 500):
            pix_area_arr[i] = 0
            continue

        # in any other case you'll have to define the corners and proceed
        tl_x = pix_centers[i][0] - 5e2
        tr_x = pix_centers[i][0] + 5e2
        bl_x = pix_centers[i][0] - 5e2
        br_x = pix_centers[i][0] + 5e2

        tl_y = pix_centers[i][1] + 5e2
        tr_y = pix_centers[i][1] + 5e2
        bl_y = pix_centers[i][1] - 5e2
        br_y = pix_centers[i][1] - 5e2

        tl = [tl_x, tl_y]
        tr = [tr_x, tr_y]
        bl = [bl_x, bl_y]
        br = [br_x, br_y]

        pixel_corners = [tl, tr, bl, br]  # top and bottom, left and right

        # if the pixel center is "well" inside the outer excluding rectangle then 
        # you only need to check the pixel corners with respect to the inner polygon
        if (min(out_rect_x) + 500 < pix_x_cen_arr[i]) and (pix_x_cen_arr[i] < max(out_rect_x) - 500) and \
        (min(out_rect_y) + 500 < pix_y_cen_arr[i]) and (pix_y_cen_arr[i] < max(out_rect_y) - 500):
            # so its inside the outer rect
            # check corners with inner polygon
            in_bool_tl = check_point_polygon(tl[0], tl[1], poly_inner)
            in_bool_tr = check_point_polygon(tr[0], tr[1], poly_inner)
            in_bool_bl = check_point_polygon(bl[0], bl[1], poly_inner)
            in_bool_br = check_point_polygon(br[0], br[1], poly_inner)

            # all cases are either pixels 
            # that intersect the edge or the pix in the annulus
            # the 5 cases are:
            # Case 1: All vertices True -- TTTT
            # Case 2: One False -- TTTF
            # Case 3: two False -- TTFF
            # Case 4: three False -- TFFF
            # Case 5: All veritces False -- FFFF

            # Case 1:
            if in_bool_tl and in_bool_tr and in_bool_bl and in_bool_br:
                pix_poly = get_pixel_intersect_shape()

                if pix_poly is None:
                    pix_area_arr[i] = 0
                    continue
                else:
                    continue

        # check corners 
        out_bool_tl = check_point_polygon(tl[0], tl[1], poly_outer)
        out_bool_tr = check_point_polygon(tr[0], tr[1], poly_outer)
        out_bool_bl = check_point_polygon(bl[0], bl[1], poly_outer)
        out_bool_br = check_point_polygon(br[0], br[1], poly_outer)
        

        # Case 1: All pixels True (i.e. TTTT)
        if out_bool_tl and out_bool_tr and out_bool_bl and out_bool_br:
            if in_bool_tl and in_bool_tr and in_bool_bl and in_bool_br:
                pix_area_arr[i] = 0

        pix_area_arr[i] = pixel_area(x_cen, y_cen, poly)

    return None

def get_pixel_area(polygon_x, polygon_y, pix_x_cen, pix_y_cen, case):
    """ 
    I'm not developing this function anymore. The problem of finding the 
    intersection of two polygons, which can be convex or concave, seems
    to be more complicated than I expected. A generic mathematical solution
    will take far too long to derive on my own. 
    I'm switching over to using the General Polygon Clipping (GPC) Library. 
    This was written by ALan Murta at U. Manchester. http://www.cs.man.ac.uk/~toby/gpc/
    Seems like it took him a couple years to write. I'm not sure what algorithms
    he uses to find polygon intersection but it appears to do what I want it 
    to do.
    # With respect to the problem at hand this function does not check if the 
    # supplied bounding polygon is the inner one or the outer one. 
    # The preceding code must make sure the correct polygon is given.
    # find number of bounding polygon points inside pixel
    # if there are none:
    # then find the closest two points on either side of the pixel
    # else if there are some points inside the pixel: 
    # then those are to be taken into account for finding the intersecting shape
    # again you need to find the two closest points on the bounding polygon 
    # that are closest to an edge
    """

    polygon_inside_idx_x = np.where((polygon_x > pix_x_cen - 500) & (polygon_x < pix_x_cen + 500))[0]
    polygon_inside_idx_y = np.where((polygon_y > pix_y_cen - 500) & (polygon_y < pix_y_cen + 500))[0]

    polygon_inside_idx = np.intersect1d(polygon_inside_idx_x, polygon_inside_idx_y)

    # define intersecting polygon as empty list
    intersect_poly = []

    if polygon_inside_idx.size:
        # ie. there are bounding polygon vertices inside the pixel
        # this can happen for all cases

        if (case == 'tttt') or (case == 'ffff'):

            if (0 in polygon_inside_idx) or ((len(polygon_x)-1) in polygon_inside_idx):
                # returning wrong value # will fix this later
                if case == 'tttt':
                    return 1
                elif case == 'ffff':
                    return 0 
            elif (max(polygon_inside_idx) - min(polygon_inside_idx)) == len(polygon_inside_idx)-1:
                start_idx = min(polygon_inside_idx)
                end_idx = max(polygon_inside_idx)

            # loop over all segments
            seg_count = 0
            total_segments = len() + 2
            for i in range(int(end_idx) - int(start_idx) + 2):

                if start_idx+i >= len(polygon_x):
                    # must check that end_idx + 1 does not lead to an index outside the array indices
                    # it should come back around to 0 if that is the case
                    x1 = polygon_x[start_idx+i-1 - len(polygon_x)]
                    y1 = polygon_y[start_idx+i-1 - len(polygon_x)]
                    x2 = polygon_x[start_idx+i - len(polygon_x)]
                    y2 = polygon_y[start_idx+i - len(polygon_x)]
                else:
                    x1 = polygon_x[start_idx+i-1]
                    y1 = polygon_y[start_idx+i-1]
                    x2 = polygon_x[start_idx+i]
                    y2 = polygon_y[start_idx+i] 

                # find if this segment can intersect any pixel edge
                seg_inter_begin = np.where((x1 >= pix_x_cen - 500) & (x1 <= pix_x_cen + 500) &\
                    (y1 > pix_y_cen - 500) & (y1 < pix_y_cen + 500))[0]

                seg_inter_finish = np.where((x2 >= pix_x_cen - 500) & (x2 <= pix_x_cen + 500) &\
                    (y2 >= pix_y_cen - 500) & (y2 <= pix_y_cen + 500))[0]

                # redundant checks
                # should never be trigerred if the program logic is correct
                if (seg_count == 0) and seg_inter_begin.size:
                    print "Something went wrong with assigning indices to bounding polygon vertices inside pixel. Exiting..."
                    sys.exit(0)
                elif (seg_count == 0) and not seg_inter_finish.size:
                    print "Something went wrong with assigning indices to bounding polygon vertices inside pixel. Exiting..."
                    sys.exit(0)

                # the beginning of the last segment by default must be part of the intersecting polygon 
                if seg_count == total_segments - 1:
                    intersect_poly.append((x1,y1))
                    seg_count += 1

                    continue

                # start checking which segment ends are inside pixel
                if seg_inter_begin.size or seg_inter_finish.size:
                    # i.e. the case where one or both of the ends of the segments is inside the pixel
                    # find the point of intersection

                    # if both segment ends are inside then 
                    # check that the previous 
                    if seg_inter_begin.size and seg_inter_finish.size:
                        intersect_poly.append((x2,y2))
                        # check that the previous one did actually a
                        seg_count += 1
                        continue

                    # if only one end is inside then find where the segment intersects the edge
                    for j in range(4):

                        x3 = pixel_corners[j-1][0]
                        y3 = pixel_corners[j-1][1]
                        x4 = pixel_corners[j][0]
                        y4 = pixel_corners[j][1]

                        m1 = (y2 - y1) / (x2 - x1)
                        m2 = (y4 - y3) / (x4 - x3)

                        c1 = y1 - m1 * x1
                        c2 = y3 - m2 * x3

                        xi = (c2 - c1) / (m1 - m2)
                        yi = (m1*c2 - m2*c1) / (m1 - m2)

                        # check that the intersection point found actually does lie on the edge and the segment
                        if (x1 <= xi <= x2) and (x3 <= xi <= x4) and (y1 <= yi <= y2) and (y3 <= yi <= y4):
                            intersect_poly.append((xi,yi))

                            # the end of the first segment by default must be part of the intersecting polygon 
                            if seg_count == 0:
                                intersect_poly.append((x2,y2))

                            edge_inter_count.append()
                            seg_count += 1
                            #if (in edge_inter_count) and (seg_count == total_segments - 1):
                            #    intersect_poly.append((pixel_corners[edge_count][0], pixel_corners[edge_count][1]))
                            break
                        else:
                            continue

                else:
                    # move to the next segment
                    seg_count += 1
                    continue

            return get_area(intersect_poly)

        elif case == 'tttf':
            # dummy return # to be fixed
            return None

    # check that there are no bounding polygon vertices on the pixel edge
    # this can only occur for the middle 3 cases with non-zero intersecting area
    if polygon_edge_idx.size:
        #ie. there are bounding polygon vertices on the pixel edge
        # dummy return # to be fixed
        return None

    else:
        # find the two closest points on either side
        # dummy return # to be fixed
        return None

    #else:
    #    # these are the couple cases where the pixel intersecting area is either 0 or 1
    #    # in this case there will be no polygon points that are inside the pixel
    #    if (case == 'tttt'):
    #        return 1

    #    elif (case == 'ffff'):
    #        return 0

    #    elif (case == 'tttf'):
    #        # dummy return # to be fixed
    #        return None

def return_unzipped_list(poly):
    """
    This function can take a polygon object from the Polygon module
    or it can take a numpy array.
    In both cases the polygon is a list of coordinate pairs.
    """

    if type(poly) is pg.cPolygon.Polygon:
        px, py = zip(*poly[0])
    elif type(poly) is np.ndarray:
        px, py = zip(*poly)

    return px, py

def polygon_plot_prep(poly):

    px, py = return_unzipped_list(poly)

    px = np.append(px, px[0])
    py = np.append(py, py[0])

    return px, py

def plot_polygon_intersection(poly1, poly2):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    px1, py1 = polygon_plot_prep(poly1)
    px2, py2 = polygon_plot_prep(poly2)

    ax.plot(px1, py1, '-', color='b', lw=2)
    ax.plot(px2, py2, '-', color='r', lw=2)

    poly_i = poly1 & poly2
    print poly_i
    px_i, py_i = polygon_plot_prep(poly_i)

    ax.plot(px_i, py_i, '-', color='k', lw=3)

    # shade intersection region
    poly_i_patch = Polygon(poly_i[0], closed=True, fill=True, color='gray', alpha=0.5)
    ax.add_patch(poly_i_patch)

    # get bounding boxes to set limits automatically
    xmin1, xmax1, ymin1, ymax1 = poly1.boundingBox()
    xmin2, xmax2, ymin2, ymax2 = poly2.boundingBox()

    xmin = min(xmin1, xmin2) - 3
    xmax = max(xmax1, xmax2) + 3
    ymin = min(ymin1, ymin2) - 3  #0.1 * min(ymin1, ymin2)
    ymax = max(ymax1, ymax2) + 3  #0.1 * max(ymax1, ymax2)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.show()

    return None

def plot_single_polygon(poly):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    px, py = polygon_plot_prep(poly)

    ax.plot(px, py, '-', color='k', lw=2)

    xmin, xmax, ymin, ymax = poly.boundingBox()

    xmin = xmin - 300
    xmax = xmax + 300
    ymin = ymin - 300
    ymax = ymax + 300

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.show()

    return None

def plot_region(vert_x, vert_y, vert_x_cen, vert_y_cen, eff_rad, valid_in, valid_out,\
    region_name='orientale', save=False, with_craters=False, show_rect=False):
    # plots the region of interest
    # showing the annulus with the inner and outer polygons in different colors

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # plot dividing circle
    circle = plt.Circle((vertices_x_center, vertices_y_center), eff_rad, color='black', fill=False, ls='--')
    ax.add_artist(circle)

    # plot vertices inside and outside dividing circle in different colors
    ax.plot(vertices_x[valid_out], vertices_y[valid_out], '.-', color='r', markersize=2, markeredgecolor='r')
    ax.plot(vertices_x[valid_in], vertices_y[valid_in], '.-', color='g', markersize=2, markeredgecolor='g')

    plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    #if with_craters:

    if show_rect:
        # Plot the rectangles within which the pixel corners are not checked
        in_rect_x = [-3.35e6, -2.55e6, -2.55e6, -3.35e6, -3.35e6]
        in_rect_y = [-9.5e5, -9.5e5, -3.3e5, -3.3e5, -9.5e5]

        out_rect_x = [-3.55e6, -2.32e6, -2.32e6, -3.59e6, -3.55e6]
        out_rect_y = [-1.47e6, -1.47e6, 5e4, 5e4, -1.47e6]

        ax.plot(in_rect_x, in_rect_y, '-', color='b')
        ax.plot(out_rect_x, out_rect_y, '-', color='b')

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    if save:
        fig.savefig(mdir + region_name + '.png', dpi=300)
    else:
        plt.show()

    return None

def get_rows_columns(pix_x_cen_arr, pix_y_cen_arr):

    orig = pix_x_cen_arr[0]
    for j in range(1, len(pix_x_cen_arr)):
        if pix_x_cen_arr[j] == orig:
            columns = j
            break

    orig = pix_y_cen_arr[0]
    rows = 1
    for j in range(0, len(pix_y_cen_arr)):
        if pix_y_cen_arr[j] != orig:
            rows += 1
            orig = pix_y_cen_arr[j]

    return rows, columns

def closest_pixel_indices(xp, yp, X, Y):

    x_dist_arr = np.abs(X-xp)
    y_dist_arr = np.abs(Y-yp)
    idx_arr = np.where((x_dist_arr == np.min(x_dist_arr)) & (y_dist_arr == np.min(y_dist_arr)))

    row_idx, col_idx = int(idx_arr[0]), int(idx_arr[1])

    return row_idx, col_idx

def get_pixels_in_bbox_old(bbox, pix_x_cen_arr, pix_y_cen_arr, mode='run'):
    """
    This is the old function I wrote which is too precise 
    (in getting the correct pixels in the bbox) and also
    uses far too much brute force to get to the answer. 
    The new function get_pixels_in_bbox() uses an analytical
    solution and is way faster.
    """

    # get limits from the bounding box
    xmin = bbox[0]
    xmax = bbox[1]
    ymin = bbox[2]
    ymax = bbox[3]

    # turn the limits into search area for pixels
    # i.e. be conservative and search an additional 500 units on each side
    # the _s is for search
    xmin_s = xmin - 500
    xmax_s = xmax + 500
    ymin_s = ymin - 500
    ymax_s = ymax + 500

    # now look for pixels within the search area
    # all you need to do is to find the pixel indices at the four corners 
    # of the bounding box and you can easily populate the rest of the array.
    pix_bbox_x = []
    pix_bbox_y = []

    # first, create a coordinate grid
    x1d_short = np.arange(np.min(pix_x_cen_arr), np.max(pix_x_cen_arr)+1000.0, 1000.0)
    y1d_short = np.arange(np.min(pix_y_cen_arr), np.max(pix_y_cen_arr)+1000.0, 1000.0)

    X, Y = np.meshgrid(x1d_short, y1d_short)

    # second, find the pixel coords (and their 
    # indices) that are closest to hte search corners
    len_x1d_arr = len(x1d_short)
    len_y1d_arr = len(y1d_short)
    
    bl_row_idx, bl_col_idx = closest_pixel_indices(xmin_s, ymin_s, X, Y)
    tr_row_idx, tr_col_idx = closest_pixel_indices(xmax_s, ymax_s, X, Y)
    tl_row_idx, tl_col_idx = closest_pixel_indices(xmin_s, ymax_s, X, Y)
    # the row and col indices in teh above lines are indices that will give you
    # the x and y values of the pixel center that is closest to the search corner
    # e.g. X[bl_row_idx, bl_col_idx] and Y[bl_row_idx, bl_col_idx] are the x and y
    # coords of the pixel closest to the bottom left corner of hte search area

    bl_xy_1d_idx = np.where((pix_x_cen_arr == X[bl_row_idx, bl_col_idx]) & (pix_y_cen_arr == Y[bl_row_idx, bl_col_idx]))[0]
    tr_xy_1d_idx = np.where((pix_x_cen_arr == X[tr_row_idx, tr_col_idx]) & (pix_y_cen_arr == Y[tr_row_idx, tr_col_idx]))[0]
    tl_xy_1d_idx = np.where((pix_x_cen_arr == X[tl_row_idx, tl_col_idx]) & (pix_y_cen_arr == Y[tl_row_idx, tl_col_idx]))[0]

    # will run the following lines in test mode to check
    # if the 1d array index assignment worked
    if mode == 'test':
        print '\n'
        print X[bl_row_idx, bl_col_idx], Y[bl_row_idx, bl_col_idx], pix_x_cen_arr[bl_xy_1d_idx], pix_y_cen_arr[bl_xy_1d_idx]
        print X[tr_row_idx, tr_col_idx], Y[tr_row_idx, tr_col_idx], pix_x_cen_arr[tr_xy_1d_idx], pix_y_cen_arr[tr_xy_1d_idx]
        print X[tl_row_idx, tl_col_idx], Y[tl_row_idx, tl_col_idx], pix_x_cen_arr[tl_xy_1d_idx], pix_y_cen_arr[tl_xy_1d_idx]

    # lastly, populate the pixel and corresponding index array
    # first populate the x and y arrays in bbox
    x_bbox_min = X[bl_row_idx, bl_col_idx]
    x_bbox_max = X[tr_row_idx, tr_col_idx]
    y_bbox_min = Y[bl_row_idx, bl_col_idx]
    y_bbox_max = Y[tr_row_idx, tr_col_idx]
    if mode == 'test':
        print "xmin and xmax values in search area:", x_bbox_min, x_bbox_max
        print "ymin and ymax values in search area:", y_bbox_min, y_bbox_max

    pix_bbox_x_short = np.arange(x_bbox_min, x_bbox_max+1000.0, 1000.0)
    pix_bbox_y_short = np.arange(y_bbox_min, y_bbox_max+1000.0, 1000.0)
    pix_bbox_x = pix_bbox_x_short
    pix_bbox_y = pix_bbox_y_short

    for u in range(len(pix_bbox_y_short) - 1):
        pix_bbox_x = np.append(pix_bbox_x, pix_bbox_x_short)

    for v in range(len(pix_bbox_x_short) - 1):
        pix_bbox_y = np.vstack((pix_bbox_y, pix_bbox_y_short))

    pix_bbox_y = pix_bbox_y.T.flatten()
    pix_bbox_y = pix_bbox_y[::-1]
    # I'm reversing the y array because the original pix_y_cen_arr is also reversed
    # i.e. it goes from max to min 
    # because the origin of the coordinates originally given in the slope file
    # is to the top left. i.e. the coordinates are mostly in the fourth quadrant
    if mode == 'test':
        print len(pix_bbox_x), len(pix_bbox_y)  # should be equal lengths

    # populate pixel index array
    pixel_indices = []

    rows, columns = get_rows_columns(pix_x_cen_arr, pix_y_cen_arr)
    rows_in_bbox, columns_in_bbox = get_rows_columns(pix_bbox_x, pix_bbox_y)
    # this will (almost?) always be square
    # the bounding box for hte circle will always be square 
    # but because I add 500 to the search area, that might 
    # cause the returned shape to be rectangular

    current_start = int(tl_xy_1d_idx)
    current_row_indices = np.arange(int(tl_xy_1d_idx), int(current_start + columns_in_bbox), 1)
    row_count = 0
    while 1:
        if row_count == rows_in_bbox:
            break

        for w in range(len(current_row_indices)):
            pixel_indices.append(current_row_indices[w])
        
        row_count += 1
        current_start += columns
        current_row_indices = np.arange(int(current_start), int(current_start + columns_in_bbox), 1)

    pixel_indices = np.asarray(pixel_indices)

    return pix_bbox_x, pix_bbox_y, pixel_indices

def get_pixels_in_bbox(bbox, pix_x_cen_arr, pix_y_cen_arr, rows, columns):

    # get limits from the bounding box
    xmin = bbox[0]
    xmax = bbox[1]
    ymin = bbox[2]
    ymax = bbox[3]

    # turn the limits into search area for pixels
    # i.e. be conservative and search an additional 
    # 1000 units (i.e. 1 pixel) on each side. 
    # the _s is for search
    xmin_s = xmin - 1000
    xmax_s = xmax + 1000
    ymin_s = ymin - 1000
    ymax_s = ymax + 1000

    # now get the coordinates of hte top left pixels 
    # in the full array and the bounding box.
    x_tl = pix_x_cen_arr[0]
    y_tl = pix_y_cen_arr[0]

    x_bbox_tl = xmin_s
    y_bbox_tl = ymax_s

    # now get the differences to find hte number
    # of rows and columns in between these two coords.
    delta_x = abs(x_bbox_tl - x_tl)
    delta_y = abs(y_bbox_tl - y_tl)

    delta_rows = int(delta_y / 1000)
    delta_cols = int(delta_x / 1000)

    # now populate the pixel values and indices arrays
    bbox_rows = int(((ymax_s - ymin_s) / 1000) + 1) + 2  # padding it by 2 pixels just so that I don't miss any
    bbox_cols = int(((xmax_s - xmin_s) / 1000) + 1) + 2

    pixel_indices = []
    for row_count in range(bbox_rows):
        tl_idx = columns*(delta_rows + row_count) + delta_cols
        tr_idx = tl_idx + bbox_cols
        pixel_indices.append(np.arange(tl_idx, tr_idx, 1))
    
    # make sure that pixel_indices is 1d 
    # and get pix values and return
    pixel_indices = np.asarray(pixel_indices)
    pixel_indices = pixel_indices.ravel()

    # remove all pixel indices that are greater than
    # the total length of the full array
    # this is a check to make sure that the indices 
    # supplied don't cause the next bit of code to 
    # look for pixel centers that don't exist.
    valid_idx = np.where(pixel_indices <= len(pix_x_cen_arr))[0]
    pixel_indices = pixel_indices[valid_idx]
    
    pix_bbox_x = pix_x_cen_arr[pixel_indices]
    pix_bbox_y = pix_y_cen_arr[pixel_indices]

    return pix_bbox_x, pix_bbox_y, pixel_indices

def crater_test(pix_x_cen_arr, pix_y_cen_arr):

    # define sample craters
    c1 = Circle(radius=2.5e5, center=(-3.4e6,0), points=128)
    c2 = Circle(radius=3.5e5, center=(-3.5e6,-0.1e6), points=128)
    c3 = Circle(radius=1e5, center=(-3.3e6,0), points=128)
    c4 = Circle(radius=2e5, center=(-2e6,-1.5e6), points=128)
    c5 = Circle(radius=1e5, center=(-2.2e6,-1.4e6), points=128)

    # plot all craters
    fig = plt.figure()
    ax = fig.add_subplot(111)

    c1x, c1y = polygon_plot_prep(c1)
    c2x, c2y = polygon_plot_prep(c2)
    c3x, c3y = polygon_plot_prep(c3)
    c4x, c4y = polygon_plot_prep(c4)
    c5x, c5y = polygon_plot_prep(c5)

    # do the crater calc
    craters_x = np.array([c1x, c2x, c3x, c4x, c5x])
    craters_y = np.array([c1y, c2y, c3y, c4y, c5y])
    all_craters = [c1,c2,c3,c4,c5]

    pix_crater_area = np.zeros(len(pix_x_cen_arr))

    for i in range(len(craters_x)):

        current_crater_x_cen = craters_x[i]
        current_crater_y_cen = craters_y[i]

        crater_poly = all_craters[i]

        pix_bbox_x, pix_bbox_y, pixel_indices =\
         get_pixels_in_bbox(crater_poly.boundingBox(), pix_x_cen_arr, pix_y_cen_arr, mode='test')

        # first check that the lengths of indices array and 
        # returned x and y array are equal and then
        # check if the x and y elements given by pixel indices
        # are indeed the elements in pix_bbox_x and pix_bbox_y
        print "Returned x, y, and index arrays are:"
        print pix_bbox_x
        print pix_bbox_y
        print pixel_indices

        if len(pix_bbox_x) == len(pix_bbox_y) == len(pixel_indices):
            print "Equal length. Now checking for equality of pixel coord values that are obtained via two different ways."
        else:
            print "Lengths:", len(pix_bbox_x), len(pix_bbox_y), len(pixel_indices)
            print "Returned arrays are not of equal length. Exiting."
            sys.exit(0)
        print np.array_equal(pix_bbox_x, pix_x_cen_arr[pixel_indices])
        print np.array_equal(pix_bbox_y, pix_y_cen_arr[pixel_indices])

        for j in range(len(pix_bbox_x)):

            current_pix_x_cen = pix_bbox_x[j]
            current_pix_y_cen = pix_bbox_y[j]

            # define a polygon using pixel corners in exactly the same way as done for the pixel fraction case
            tl_x = current_pix_x_cen - 5e2
            tr_x = current_pix_x_cen + 5e2
            bl_x = current_pix_x_cen - 5e2
            br_x = current_pix_x_cen + 5e2

            tl_y = current_pix_y_cen + 5e2
            tr_y = current_pix_y_cen + 5e2
            bl_y = current_pix_y_cen - 5e2
            br_y = current_pix_y_cen - 5e2

            tl = [tl_x, tl_y]
            tr = [tr_x, tr_y]
            bl = [bl_x, bl_y]
            br = [br_x, br_y]

            pixel_corners = [tl, tr, br, bl]  # top and bottom, left and right going clockwise

            pixel_corners = pg.Polygon(pixel_corners)

            # find the area of intersection between the pixel and crater
            inter_area = (pixel_corners & crater_poly).area()

            # find pixel index using pixel center to append to the correct array element
            pix_index = pixel_indices[j]
            pix_crater_area[pix_index] += inter_area

    # pix_crater_area /= 1e6 -- normalized to 1 sq km if needed (comment out if using fractions)

    print np.where(pix_crater_area != 0)
    rows, columns = get_rows_columns(pix_x_cen_arr, pix_y_cen_arr)
    im = ax.imshow(pix_crater_area.reshape(rows, columns), cmap='bone')
    plt.colorbar(im, ax=ax)

    plt.show()

    plt.clf()
    plt.cla()
    plt.close()

    return None

def get_crater_circle_poly(x_cen, y_cen, dia):

    rad = dia / 2

    circ_poly = Circle(radius=rad, center=(current_crater_x_cen,current_crater_y_cen), points=128)
    # the crater circle is approximated using a polygon of 128 vertices

    return circ_poly

def get_density(crater_frac, pix_frac, total_pix):

    density = np.zeros(total_pix)
    for i in range(total_pix):

        if pix_frac[i] != 0.0:
            density[i] = crater_frac[i] / pix_frac[i] 
            # Because pixel area is currently 1 sq. km, we do not need to multiply 
            # pix_frac by the area of a single pixel. 
            # Density = # of craters per sq. km.
            # If you ever have pixels that are not 1 sq. km in area then you should use
            # density[i] = crater_frac[i] / pix_area_arr_phys[i] (pix_area_arr_phys calculation currently commented out)
            # with area_single_pix (see code above where it saves the pixel area frac) 
            # set to the area of a pixel in the physical units you need.
        
        if (pix_frac[i] == 0.0) or (pix_frac[i] == -9999.0):
            density[i] = np.nan

    return density

def write_LofL_pickle(l, outname):
    """Function came from SO solution:
    see here -- https://stackoverflow.com/questions/25464295/how-to-pickle-a-list
    """

    import pickle

    with open(outname + ".pkl", "wb") as f:
        pickle.dump(l, f)

    return None

def plot_im(arr, outname):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(arr, cmap='bone')
    plt.colorbar(im, ax=ax)

    fig.savefig(slope_extdir + outname + '.png', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.cla()
    plt.close()

    return None

def save_csv(arr, arr_name, datatype, savepath, hdr):

    # save as csv
    data = np.array(zip(arr), dtype=[(arr_name, datatype)])
    # the string in the dtype here should match the array variable
    np.savetxt(savepath, data, fmt=['%.4f'], delimiter=',', header=hdr)

    return None

def convert_arr_to_bool_int(arr, returntype):

    if returntype == 'bool' or returntype == 'boolean':
        return arr.astype(bool)
    else:
        nonzero_idx = np.nonzero(arr)
        arr[nonzero_idx] = 1.0
        return arr

if __name__ == '__main__': 
    
    # Start time
    start = time.time()
    dt = datetime.datetime
    print "Starting at --", dt.now()

    # read in catalogs
    vertices_cat = np.genfromtxt(slope_extdir + 'HF_vertices_m.csv', dtype=None, names=True, delimiter=',')

    # Old code using crater centers
    #craters_cat = np.genfromtxt(slope_extdir + 'CRATER_FullHF_m.csv', dtype=None, names=True, delimiter=',')
    #craters_x = craters_cat['x_coord_m']
    #craters_y = craters_cat['y_coord_m']
    #craters_diam = craters_cat['Diameter_m']

    # create arrays for more convenient access
    vertices_x = vertices_cat['x_coord_m']
    vertices_y = vertices_cat['y_coord_m']

    # delete offending points -- those that cause the polygon to cross itself
    #argmin takes an array of difference and it gives us the argument that 
    #gave the minimum difference (so finding the closest point in an array).
    #Here, we store the x and y positions of the closest points and then delete them.
    off_x_idx1 = np.argmin(abs(vertices_x - -2.41165e6))
    off_y_idx1 = np.argmin(abs(vertices_y - 176717))

    off_x_idx2 = np.argmin(abs(vertices_x - -3.61074e6))
    off_y_idx2 = np.argmin(abs(vertices_y - 41808.2))

    off_x_idx3 = np.argmin(abs(vertices_x - -3.61526e6))
    off_y_idx3 = np.argmin(abs(vertices_y - 41295.4))

    off_x_idx = np.array([off_x_idx1, off_x_idx2, off_x_idx3])
    off_y_idx = np.array([off_y_idx1, off_y_idx2, off_y_idx3])
    vertices_x = np.delete(vertices_x, off_x_idx, axis=None)
    vertices_y = np.delete(vertices_y, off_y_idx, axis=None)

    # define radius and centers for vertices
    # eyeballed for now
    vertices_x_center = -2.94e6
    vertices_y_center = -5.87e5
    rad_vertices = np.sqrt((vertices_x - vertices_x_center)**2 + (vertices_y - vertices_y_center)**2)
    eff_rad = np.min(rad_vertices) + 25e4  # number put in by trial and error
    
    # define valid indices for vertices inside and outside dividing effective radius 
    #i.e. if the position is farther from the center of the study area than the exent 
    #of the radius of the outer circle, then it's outside and vice versa for the radius 
    #of the inner circle. np.where gives the array indices that satisfy a given condition.
    valid_out = np.where(rad_vertices > eff_rad)[0]
    valid_in = np.where(rad_vertices < eff_rad)[0]

    #plot_region(vertices_x, vertices_y, vertices_x_center, vertices_y_center, eff_rad, valid_in, valid_out,\
    #region_name='orientale', save=True, with_craters=False, show_rect=False)

    # define inner and outer polygons
    poly_outer = zip(vertices_x[valid_out], vertices_y[valid_out])
    poly_inner = zip(vertices_x[valid_in], vertices_y[valid_in])

    poly_outer = pg.Polygon(poly_outer)
    poly_inner = pg.Polygon(poly_inner)

    # ----------------  measure and populate pixel area fraction array  ---------------- # 
    # read in pixel slope info
    # this file also gives the x and y centers of pixels
    do_pix_frac = False
    if do_pix_frac:
        print "Will compute pixel fractions now."
    else:
        print "Will skip computing pixel fractions. Moving to crater fractions now."

    # Old file for slope values and pixel centers that
    # was in a bigger area than the clipped area which 
    # is used for the newer arrays with secondary craters
    # stamped out.
    #slope_arr = np.load(slope_extdir + '3km_slope_points.npy')
    #pix_x_cen_arr = slope_arr['pix_x_cen']
    #pix_y_cen_arr = slope_arr['pix_y_cen']
    #slope = slope_arr['slope_val']

    # New arrays that have compatible shapes with the 
    # clipped files. 
    # The limits have been hardcoded in for now.
    # Will have to be changed if the clipped area changes.
    nrows = 2109
    ncols = 1949
    pix_x_cen_arr = np.arange(-3822217, -1874217 + 1000, 1000)
    pix_y_cen_arr = np.ones(ncols) * 380491
    for count in range(1,nrows):
        pix_x_cen_arr = np.append(pix_x_cen_arr, np.arange(-3822217, -1874217 + 1000, 1000))
        pix_y_cen_arr = np.append(pix_y_cen_arr, np.ones(ncols) * (380491 - 1000*count))

    rows, columns = get_rows_columns(pix_x_cen_arr, pix_y_cen_arr)

    # zip combines the given arrays.
    # In this case, the pixel centers were taken from the slope raster, 
    # which was converted to a numpy binary file to minimize computation time. 
    pix_centers = zip(pix_x_cen_arr, pix_y_cen_arr)
    pix_area_arr = np.zeros(len(pix_x_cen_arr))

    # define rectangles within which pixels can be skipped, i.e. well within 
    # the inner rectangle or well beyond the outer rectangle. You need 5 points 
    # to define a rectangle in order to make sure the polygon closes. However, 
    # you only really need the min/max extent of the area shapefile.
    inner_rect_x = [-3.35e6, -2.55e6, -2.55e6, -3.35e6, -3.35e6]
    inner_rect_y = [-9.5e5, -9.5e5, -3.3e5, -3.3e5, -9.5e5]

    outer_rect_x = [-3.55e6, -2.32e6, -2.32e6, -3.59e6, -3.55e6]
    outer_rect_y = [-1.47e6, -1.47e6, 5e4, 5e4, -1.47e6]
    
    # loop over all pixels, range just designates the iterable -- in this case, 
    # the pix_centers array.
    if do_pix_frac:
        for i in range(len(pix_centers)):

            if (i % 100000) == 0.0:
                print '\r',
                print "At pixel number:",'{0:.2e}'.format(i),\
                "; time taken up to now:",'{0:.2f}'.format((time.time() - start)/60),"minutes.",
                sys.stdout.flush()

            # check if pixel center falls "well" inside the inner excluding rectangle
            if (min(inner_rect_x) + 500 < pix_x_cen_arr[i]) and (pix_x_cen_arr[i] < max(inner_rect_x) - 500) and \
            (min(inner_rect_y) + 500 < pix_y_cen_arr[i]) and (pix_y_cen_arr[i] < max(inner_rect_y) - 500):
                pix_area_arr[i] = 0.0
                continue
            # pix_area_arr defines the starting array for the fractional area of pixels within the study area.
            # in any other case you'll have to define the corners and proceed.
            tl_x = pix_centers[i][0] - 5e2
            tr_x = pix_centers[i][0] + 5e2
            bl_x = pix_centers[i][0] - 5e2
            br_x = pix_centers[i][0] + 5e2

            tl_y = pix_centers[i][1] + 5e2
            tr_y = pix_centers[i][1] + 5e2
            bl_y = pix_centers[i][1] - 5e2
            br_y = pix_centers[i][1] - 5e2

            tl = [tl_x, tl_y]
            tr = [tr_x, tr_y]
            bl = [bl_x, bl_y]
            br = [br_x, br_y]

            pixel_corners = [tl, tr, br, bl]  # top and bottom, left and right going clockwise

            pixel_corners = pg.Polygon(pixel_corners) # creates a polygon for each pixel as it iterates

            # The Polygon module is capable of finding the area of intersection between two polygons, which is what we've implemented below.
            # case 1: check if the pixel is completely inside both polygons
            if ((pixel_corners & poly_inner).area() == 1e6) and ((pixel_corners & poly_outer).area() == 1e6):
                # if it is completely inside the inner polygon then it is not part of the
                # annulus of interest, so it gets assigned zero.
                pix_area_arr[i] = 0.0
                continue

            # case 2: check if the pixel is completely outside both polygons, also assigned zero.
            if ((pixel_corners & poly_inner).area() == 0.0) and ((pixel_corners & poly_outer).area() == 0.0):
                pix_area_arr[i] = 0.0
                continue

            # case 3: check if the pixel is completely outside the inner polygon but completely inside the outer polygon
            if ((pixel_corners & poly_inner).area() == 0.0) and ((pixel_corners & poly_outer).area() == 1e6):
                # if it is outside the inner polygon but inside the outer one (i.e. completely within the annulus)
                pix_area_arr[i] = 1.0
                continue

            # case 4: check if the pixel is completely inside the outer polygon but intersects the inner polygon
            if ((pixel_corners & poly_inner).area() < 1e6) and ((pixel_corners & poly_inner).area() != 0.0) and\
             ((pixel_corners & poly_outer).area() == 1e6):
                pix_area_arr[i] = 1.0 - (pixel_corners & poly_inner).area() / 1e6  # stores the fraction of the pixel area that is within the annulus
                continue

            # case 5: check if the pixel is completely outside the inner polygon but intersects the outer polygon
            if ((pixel_corners & poly_outer).area() < 1e6) and ((pixel_corners & poly_outer).area() != 0.0) and\
             ((pixel_corners & poly_inner).area() == 0.0):
                pix_area_arr[i] = (pixel_corners & poly_outer).area() / 1e6  # stores the fraction of the pixel area that is within the annulus
                continue

        # write all zeros as -9999.0 which is the NODATA_VALUE (for the ascii raster and numpy array)
        invalid_idx = np.where(pix_area_arr == 0.0)[0]
        pix_area_arr[invalid_idx] = -9999.0

        # This array will give pixel area in physical units
        # Uncomment the line for saving physical area arr if the area is not unity 
        # check this part of the code again and make sure that it is okay with the NODATA_VALUE
        #area_single_pix = 1.0  # in square km
        #pix_area_arr_phys = pix_area_arr * area_single_pix

        # save as numpy binary
        np.save(slope_extdir + 'pix_area_fraction.npy', pix_area_arr)
        #np.save(slope_extdir + 'pix_area_km.npy', pix_area_arr_phys)

        # save as csv
        data = np.array(zip(pix_area_arr), dtype=[('pixel_area_frac', float)])
        # the string in the dtype here should match the array variable
        np.savetxt(slope_extdir + 'pixel_area_fraction.csv', data, fmt=['%.4f'], delimiter=',', \
            header='pixel_area_fraction')

        # save as ascii raster
        su.numpy_to_asciiraster(slope_extdir + 'pix_area_fraction.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)

        print "\n","Pixel fractional area computation done and saved."
        print "Moving to craters now.", '\n'

        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(pix_area_arr.reshape(rows, columns), cmap='bone')
        plt.colorbar(im, ax=ax)

        fig.savefig(slope_extdir + 'pix_area_frac.png', dpi=300, bbox_inches='tight')
        plt.clf()
        plt.cla()
        plt.close()

    # ----------------  measure and populate crater pixel fraction array  ---------------- # 
    # read crater vertices file
    crater_vert = np.genfromtxt(slope_extdir + 'CRATER_FullHF_Vertices_coords.txt', \
        dtype=None, names=True, delimiter=',')
    crater_ids = np.unique(crater_vert['ORIG_FID'])

    # loop over all craters
    # for each crater get its bounding box and 
    # loop over all pixels in that bounding box
    # find intersecting area for each pixel and keep a running sum
    pix_crater_area = np.zeros(len(pix_x_cen_arr))

    # also create a blank array for associating 
    # crater ids with each pixel. I need a blank 
    # array of lists so that I can append to each element.
    pix_crater_id = np.zeros(len(pix_x_cen_arr))
    pix_crater_id = pix_crater_id.tolist()

    for w in range(len(pix_crater_id)):
        pix_crater_id[w] = []

    for i in range(len(crater_ids)):

        if (i % 1000) == 0.0:
            print '\r',
            print 'Analyzing crater', i+1, "of", len(crater_ids),
            sys.stdout.flush()

        # this line was originally used to create a circular polygon for each crater
        # based on its x,y center and diameter
        #get_crater_circle_poly(craters_x[i], craters_y[i], craters_diam[i])

        # Now, using the explicit x and y crater vertices to create a polygon
        current_crater_vert_idx = np.where(crater_vert['ORIG_FID'] == crater_ids[i])
        current_x_vert = crater_vert['x_coord_m'][current_crater_vert_idx]
        current_y_vert = crater_vert['y_coord_m'][current_crater_vert_idx]

        crater_poly = pg.Polygon(zip(current_x_vert, current_y_vert))

        # Check the diameter of the crater and only proceed if 
        # the diameter is equal to or greater than 1 KM.
        # The area will be in sq. meters. 
        current_diam = np.sqrt(crater_poly.area() * 4 / np.pi)
        if current_diam < 1000:
            continue

        # get all pixels within the crater's bounding box
        pix_bbox_x, pix_bbox_y, pixel_indices = \
        get_pixels_in_bbox(crater_poly.boundingBox(), pix_x_cen_arr, pix_y_cen_arr, rows, columns)

        #pix_bbox_x_old, pix_bbox_y_old, pixel_indices_old = \
        #get_pixels_in_bbox_old(crater_poly.boundingBox(), pix_x_cen_arr, pix_y_cen_arr, mode='run')

        #print '\n'
        #if pixel_indices_old.tolist() in pixel_indices.tolist(): 
        #    print i, 'True'
        #else:
        #    print i, 'False'
        #    print len(pix_bbox_x), len(pix_bbox_x_old)
        #    print len(pix_bbox_y), len(pix_bbox_y_old)
        #    print len(pixel_indices), len(pixel_indices_old)
        #    print pixel_indices
        #    print pixel_indices_old
        #continue

        # loop over all pixels within the crater's bounding box and assign crater area fraction to each
        for j in range(len(pix_bbox_x)):

            current_pix_x_cen = pix_bbox_x[j]
            current_pix_y_cen = pix_bbox_y[j]

            # define a polygon using pixel corners in exactly the same way as done for the pixel fraction case
            tl_x = current_pix_x_cen - 5e2
            tr_x = current_pix_x_cen + 5e2
            bl_x = current_pix_x_cen - 5e2
            br_x = current_pix_x_cen + 5e2

            tl_y = current_pix_y_cen + 5e2
            tr_y = current_pix_y_cen + 5e2
            bl_y = current_pix_y_cen - 5e2
            br_y = current_pix_y_cen - 5e2

            tl = [tl_x, tl_y]
            tr = [tr_x, tr_y]
            bl = [bl_x, bl_y]
            br = [br_x, br_y]

            pixel_corners = [tl, tr, br, bl]  # top and bottom, left and right going clockwise

            pixel_corners = pg.Polygon(pixel_corners)

            # find the area of intersection between the pixel and crater and
            # the fraction of original crater that area amounts to
            inter_area = (pixel_corners & crater_poly).area()
            inter_area_crater_frac = inter_area / crater_poly.area() # store the fraction of the crater occupying that pixel

            # find pixel index using pixel center to append to the correct array element
            pix_index = pixel_indices[j]
            pix_crater_area[pix_index] += inter_area_crater_frac #for each pixel, keep a running sum of the fractions of craters within it
            if inter_area_crater_frac != 0.0:
                pix_crater_id[pix_index].append(crater_ids[i])

    # pix_crater_area /= 1e6 -- normalized to 1 sq km if needed (comment out if using fractions)

    """
    write all zeros as -9999.0 which is the NODATA_VALUE (for the ascii raster and numpy array)
    For the crater area fraction array, I want to save the pixels outside the study area with the
    NODATA_VALUE and those inside the study area but NOT intersecting a crater to 0.0 which should
    be the case after the previous for loop is done. So I only need to change the pixels outside
    study area to the NODATA_VALUE using the same invalid_idx as before.
    x_area_arr = np.load(slope_extdir + 'pix_area_fraction.npy')  
    """
    #invalid_idx = np.where(pix_area_arr == 0.0)[0]
    # these lines is here just in case you're running the code only for the 
    # craters and the invalid_idx definition would be commented out otherwise
    #pix_crater_area[invalid_idx] = -9999.0

    # save as numpy binary array and csv
    np.save(slope_extdir + 'crater_area_frac_in_pix_fastcomp.npy', pix_crater_area)
    save_csv(pix_crater_area, 'pix_crater_area', float, slope_extdir + 'crater_area_fraction_in_pixel_fastcomp.csv', 'crater_area_fraction_in_pixel')

    # save the list of lists as csv
    write_LofL_pickle(pix_crater_id, slope_extdir + 'pix_crater_id_fastcomp')

    # save as ascii raster
    su.numpy_to_asciiraster(slope_extdir + 'crater_area_frac_in_pix_fastcomp.npy', (rows, columns), pix_x_cen_arr, pix_y_cen_arr)

    print "\n","Crater fractional area in each pixel computation done and saved."

    sys.exit(0)

    # ----------------  loop over all pixels and measure crater density  ---------------- # 
    # read in pix area and crater area files 
    pix_frac = np.load(slope_extdir + 'pix_area_fraction.npy')
    crater_frac = np.load(slope_extdir + 'crater_area_frac_in_pix.npy')

    density = get_density(crater_frac, pix_frac, len(pix_centers))

    np.save(slope_extdir + 'density.npy', density.reshape(rows,columns))
    np.save(slope_extdir + 'slope_val.npy', slope.reshape(rows,columns))

    # plot pixel crater fraction with slope overplotted
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(np.log10(density), slope, 'o', color='k', markersize=1)

    fig.savefig(slope_extdir + 'density_slope.png', dpi=300, bbox_inches='tight')
    plt.show()

    # total run time
    print '\n'
    print "Total time taken --", (time.time() - start)/60, "minutes."
    sys.exit(0)
