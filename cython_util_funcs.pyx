from __future__ import division

import numpy as np
cimport numpy as np

DTYPE1 = np.float64
DTYPE2 = np.int
ctypedef np.float64_t DTYPE_t1
ctypedef np.int_t DTYPE_t2

def get_crater_diams(np.ndarray[DTYPE_t1, ndim=1] crater_diam_m_arr, np.ndarray[DTYPE_t2, ndim=1] crater_ids, \
	np.ndarray[DTYPE_t2, ndim=1] crater_ids_arr, int total_craters):

    cdef np.ndarray[DTYPE_t1, ndim=1] crater_diams = np.zeros((total_craters), dtype=DTYPE1)

    cdef int i
    cdef int current_id
    cdef int current_id_idx

    for i in xrange(total_craters):
        current_id = crater_ids[i]
        current_id_idx = np.where(crater_ids_arr == current_id)[0][0]

        crater_diams[i] = crater_diam_m_arr[current_id_idx]

    return crater_diams