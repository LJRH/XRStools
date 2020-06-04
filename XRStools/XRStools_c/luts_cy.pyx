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

cdef extern from "luts.h" :

    void lutprod(
        int n_1, float *lut1,
        int n_2, float *lut2,
        int na2,
        int nb2,
        int dim1,
        int dim2,
        float *reponse_pixel,
        int     &n_result,
        float * &result,
        float * rois
    )

    void lutprod4reponse(
        int n_1, float *lut1,
        int n_2, float *lut2,
        int na1,
        int na2,
        int nb2,
        int dim1,
        int dim2,
        float *reponse_pixel,
        int dim_sol ,
        float *solution,
        int     &n_result,
        float * &result,
        int simmetrizza,
        float * rois
    )

  

def get_product4reponse(
        ndarray[numpy.float32_t, ndim = 2] lut_1,
        ndarray[numpy.float32_t, ndim = 2] lut_2,
        int na1,
        int na2,
        int nb2 ,
        ndarray[numpy.float32_t, ndim = 2]reponse_pixel,
        ndarray[numpy.float32_t, ndim = 2]solution,
        int simmetrizza,
        ndarray[numpy.float32_t, ndim = 2] rois):
    
    cdef  int n_1 =  lut_1.shape[0]
    cdef  int n_2 =  lut_2.shape[0]
    
    cdef  int dim1 =  reponse_pixel.shape[0]
    cdef  int dim2 =  reponse_pixel.shape[1]

    cdef int dim_sol = solution.shape[1]  
    
    cdef int n_result = -1
    cdef float *result = NULL

    assert   lut_1.flags["C_CONTIGUOUS"]
    assert   lut_2.flags["C_CONTIGUOUS"]
    assert   reponse_pixel.flags["C_CONTIGUOUS"]
    assert   solution.flags["C_CONTIGUOUS"]

     
    lutprod4reponse(
        n_1 , &lut_1[0,0],n_2 , &lut_2[0,0],
        na1, na2, nb2,
        dim1, dim2,
        &reponse_pixel[0,0],
        dim_sol,
        &solution[0,0],
        n_result,
        result,
        simmetrizza,
        &rois[0,0]
    )

    cdef  float[:,:]  result_py = numpy.zeros(  [n_result,3] , dtype=numpy.float32)	
    memcpy( &(result_py[0,0]), result , n_result*3*sizeof(float)   )
    free(result)

    return result_py



def get_product(
        ndarray[numpy.float32_t, ndim = 2] lut_1,
        ndarray[numpy.float32_t, ndim = 2] lut_2,
        int na2,
        int nb2 ,
        ndarray[numpy.float32_t, ndim = 2]reponse_pixel,
        ndarray[numpy.float32_t, ndim = 2] rois):

        cdef  int n_1 =  lut_1.shape[0]
        cdef  int n_2 =  lut_2.shape[0]

        cdef  int dim1 =  reponse_pixel.shape[0]
        cdef  int dim2 =  reponse_pixel.shape[1]

        cdef int n_result = -1
        cdef float *result = NULL
        assert   lut_1.flags["C_CONTIGUOUS"]
        assert   lut_2.flags["C_CONTIGUOUS"]
        assert   reponse_pixel.flags["C_CONTIGUOUS"]
        lutprod(
            n_1 , &lut_1[0,0],n_2 , &lut_2[0,0],
            na2, nb2,
            dim1, dim2,
            &reponse_pixel[0,0],
            n_result,
            result,
            &rois[0,0]
        )
        cdef  float[:,:]  result_py = numpy.zeros(  [n_result,3] , dtype=numpy.float32)	

        memcpy( &(result_py[0,0]), result , n_result*3*sizeof(float)   )
        free(result)
        return result_py


        
