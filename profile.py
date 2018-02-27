#!/usr/bin/env python
# encoding: utf-8
# filename: profile.py

#import pstats, cProfile
import line_profiler

import numpy as np
import os

import cython_util_funcs

home = os.getenv('HOME')
slope_extdir = home + '/Documents/plots_codes_for_heather/slope_effects_files/'

if __name__ == '__main__':
	
    # Read in arrays 
    crater_ids_arr = np.load(slope_extdir + 'crater_ids_arr.npy')
    crater_x_arr = np.load(slope_extdir + 'crater_x_arr.npy')
    crater_y_arr = np.load(slope_extdir + 'crater_y_arr.npy')
    crater_diam_m_arr = np.load(slope_extdir + 'crater_diam_m_arr.npy')

    # get unique ids
    crater_ids = np.unique(crater_ids_arr)
    total_craters = len(crater_ids)

    # Run profiling
    #cProfile.runctx("cython_util_funcs.get_crater_diams(crater_diam_m_arr, crater_ids, crater_ids_arr, total_craters)", \
    #    globals(), locals(), "Profile.prof")
    #s = pstats.Stats("Profile.prof")
    #s.strip_dirs().sort_stats("time").print_stats()

    profile = line_profiler.LineProfiler(cython_util_funcs.get_crater_diams)
    profile.runcall(cython_util_funcs.get_crater_diams, crater_diam_m_arr, crater_ids, crater_ids_arr, total_craters)
    profile.print_stats()
