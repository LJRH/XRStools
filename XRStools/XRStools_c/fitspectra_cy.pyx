# -*- coding: utf-8 -*-
###################################################################################
# Fits, Elastic constants fits : Alessandro Mirone
# European Synchrotron Radiation Facility
###################################################################################
#distutils: extra_compile_args = -fopenmp
# distutils: language = c++


import cython
from cython.parallel cimport prange
from cpython cimport bool
cimport numpy
from numpy cimport ndarray
import math
from libc.stdlib cimport free
from libc.string cimport memcpy
from libcpp.vector cimport vector

from libc.stdio cimport printf



cdef extern from "math.h":
    double  fabs(float)nogil
import numpy

cdef extern from "fitspectra.h" :
    void spectra_roi_by_line_c( int Nsp,
                              double *spettro_byline,
                              double* error_sum,
                              double* frequencies_sum,
                              int Nenes,
                              float* enes_data,
                              float* denominator ,
                              int ny,
                              int nx,
                              float* mms,
                              float* MASK,
                              float mine,		 
                              float de,	
                              float discard_threshold ,
                              float threshold_fraction,
                              float hline,
                              float slopeline,
                              float DHoverDI,
                              float deltaEref,
                              int useref,
                              int weight_by_response,
                              float Xintercept ,
                              float CRX ,
                              float fNmiddle,
                              float Xslope ,
                              int Nresp,
                              double* response_line_intensity
                          ) 
     
 
def spectra_roi_by_line(
        ndarray[numpy.float64_t, ndim = 1] spettro_byline,
        ndarray[numpy.float64_t, ndim = 1] error_sum ,
        ndarray[numpy.float64_t, ndim = 1] frequencies_sum ,
        ndarray[numpy.float32_t, ndim = 1] enes_data ,
        ndarray[numpy.float32_t, ndim = 1] denominator ,
        ndarray[numpy.float32_t, ndim = 3] mms ,
        ndarray[numpy.float32_t, ndim = 2] MASK ,

	float mine,		 
        float de,	
        float discard_threshold,
        float threshold_fraction,
        float hline,
        float slopeline,
        float DHoverDI,
        float deltaEref,
        int useref,
        int weight_by_response,
        float Xintercept ,
        float CRX ,
        float fNmiddle,
        float Xslope ,
        ndarray[numpy.float64_t, ndim = 1] response_line_intensity):
    
    
    cdef  Nenes =  enes_data.shape[0]
    cdef  ny  = mms.shape[1]
    cdef  nx  = mms.shape[2]
    cdef  Nsp = spettro_byline.shape[0]
    cdef Nresp = response_line_intensity.shape[0]
    
    assert Nenes == mms.shape[0]
    assert Nenes == denominator.shape[0]
    assert MASK.shape[0] == ny
    assert MASK.shape[1] == nx
    assert Nsp == error_sum.shape[0]
    assert Nsp == frequencies_sum.shape[0]
    
    assert   spettro_byline.flags["C_CONTIGUOUS"]
    assert   error_sum.flags["C_CONTIGUOUS"]
    assert   frequencies_sum.flags["C_CONTIGUOUS"]
    assert   enes_data.flags["C_CONTIGUOUS"]
    assert   mms.flags["C_CONTIGUOUS"]
    assert   denominator.flags["C_CONTIGUOUS"]
    assert   MASK.flags["C_CONTIGUOUS"]
    assert   response_line_intensity.flags["C_CONTIGUOUS"]

    print(" CHIAMO SPPCCC ")
    spectra_roi_by_line_c( Nsp,
                         &(spettro_byline[0]),
                         &(error_sum[0]),
                         &(frequencies_sum[0]),                       
                         Nenes,
                         &(enes_data[0]),
                         &(denominator [0]),
                         ny,
                         nx,
                         &(mms[0,0,0]),
                         &(MASK[0,0]),
			 mine,
                           de,
                         discard_threshold ,
                         threshold_fraction,
                         hline,
                         slopeline,
                         DHoverDI,
                         deltaEref,
                         useref,
                         weight_by_response,
                         Xintercept ,
                         CRX ,
                         fNmiddle,
                         Xslope ,
                         Nresp,
                         &(response_line_intensity[0]) 
                    ) 
    print(" CHIAMO SPPCCC ")

    return None

