from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import h5py
import math
from six.moves import range
myrank=0
try:
    import skimage.restoration
except:
    print( " ATTENTION : SKIMAGE RESORATION NOT LOADED ")




def mdiv(grad):
    res = np.zeros(grad.shape[1:])
    res[ :-1  , :, :  ] += this_grad[0, :-1  , :,:]
    res[ 1:-1 , :, :  ] -= this_grad[0, :-2  , :,:]
    res[ -1    , :, : ] -= this_grad[0,  -2  , :,:]
    res[ :, :-1  , :  ] += this_grad[1, :, :-1  ,:]
    res[ :, 1:-1 , :  ] -= this_grad[1, :, :-2  ,:]
    res[ :, -1   , : ]  -= this_grad[1, :,  -2  ,:]
    res[ :, : , :-1   ] += this_grad[1, : ,:, :-1 ]
    res[ :, : , 1:-1  ] -= this_grad[1, : ,:, :-2 ]
    res[ :, : , -1   ]  -= this_grad[1, : ,:,  -2 ]
    return res
def mygradient(img):
    shape = [3 ] + list(img.shape)
    gradient = np.zeros(shape, dtype=img.dtype)
    gradient[0,:,:,: ] = np.diff(img, axis=0)
    gradient[1,:,:,: ] = np.diff(img, axis=1)
    gradient[2,:,:,: ] = np.diff(img, axis=2)
    return gradient
def v_project(v,weight ):
    norms = np.minimum( weight, np.sqrt(  v[0]*v[0] + v[1]*v[1]     ))
    return v/ norms
def my_denoise_tv_chambolle_positive(image, weight=0.1,  n_iter_max=200):
    ndim = image.ndim
    g = np.zeros_like(p)
    x       =   np.zeros_like(image)
    tmpxa   =   np.zeros_like(image)
    v = np.zeros((image.ndim, ) + image.shape, dtype=image.dtype)
    i = 0
    sigma = 1.0/math.sqrt(8.0)
    tau   = 1.0/math.sqrt(8.0)
    while i < n_iter_max:
        tmpxa[:] = x + sigma * (  ( image-x)       + mydiv(  v   ) )
        tmpxa[:] = np.maximum (tmpxa)
        tmpxa[:] = tmpxa-x
        x[:]     =  x +  tmpxa
        tmpxa[:] =  x +  tmpxa
        v[:]  = v + tau * mygrad(tmpxa)
        v = v_project(v,weight )
    return x
        

def _denoise_tv_chambolle_nd(image, weight=0.1, eps=2.e-4, n_iter_max=200,
                             positivity=False):
    """Perform total-variation denoising on n-dimensional images.

    Parameters
    ----------
    image : ndarray
        n-D input data to be denoised.
    weight : float, optional
        Denoising weight. The greater `weight`, the more denoising (at
        the expense of fidelity to `input`).
    eps : float, optional
        Relative difference of the value of the cost function that determines
        the stop criterion. The algorithm stops when:

            (E_(n-1) - E_n) < eps * E_0

    n_iter_max : int, optional
        Maximal number of iterations used for the optimization.

    positivity : bool, optional
        Adds positivity constraint

    Returns
    -------
    out : ndarray
        Denoised array of floats.

    Notes
    -----
    Rudin, Osher and Fatemi algorithm.

    LICENCE
    -------

                    Copyright (C) 2011, the scikit-image team
                    All rights reserved.
                    
                    Redistribution and use in source and binary forms, with or without
                    modification, are permitted provided that the following conditions are
                    met:
                    
                     1. Redistributions of source code must retain the above copyright
                        notice, this list of conditions and the following disclaimer.
                     2. Redistributions in binary form must reproduce the above copyright
                        notice, this list of conditions and the following disclaimer in
                        the documentation and/or other materials provided with the
                        distribution.
                     3. Neither the name of skimage nor the names of its contributors may be
                        used to endorse or promote products derived from this software without
                        specific prior written permission.
                    
                    this software is provided by the author ``as is'' and any express or
                    implied warranties, including, but not limited to, the implied
                    warranties of merchantability and fitness for a particular purpose are
                    disclaimed. in no event shall the author be liable for any direct,
                    indirect, incidental, special, exemplary, or consequential damages
                    (including, but not limited to, procurement of substitute goods or
                    services; loss of use, data, or profits; or business interruption)
                    however caused and on any theory of liability, whether in contract,
                    strict liability, or tort (including negligence or otherwise) arising
                    in any way out of the use of this software, even if advised of the
                    possibility of such damage.



    """

    ndim = image.ndim
    p = np.zeros((image.ndim, ) + image.shape, dtype=image.dtype)
    g = np.zeros_like(p)
    d = np.zeros_like(image)
    i = 0
    while i < n_iter_max:
        if i > 0:
            # d will be the (negative) divergence of p
            d = -p.sum(0)
            slices_d = [slice(None), ] * ndim
            slices_p = [slice(None), ] * (ndim + 1)
            for ax in range(ndim):
                slices_d[ax] = slice(1, None)
                slices_p[ax+1] = slice(0, -1)
                slices_p[0] = ax
                d[tuple(slices_d)] += p[tuple(slices_p)]
                slices_d[ax] = slice(None)
                slices_p[ax+1] = slice(None)
            out_nopos = image + d
        else:
            out_nopos = image

        if not positivity:
            out = out_nopos
        else:
            out = np.maximum(0, out_nopos)
            removed = np.minimum(out_nopos, 0) 
            d = d-removed

        E = (d ** 2).sum()

        # g stores the gradients of out along each axis
        # e.g. g[0] is the first order finite difference along axis 0
        slices_g = [slice(None), ] * (ndim + 1)
        for ax in range(ndim):
            slices_g[ax+1] = slice(0, -1)
            slices_g[0] = ax
            g[tuple(slices_g)] = np.diff(out, axis=ax)
            slices_g[ax+1] = slice(None)

        norm = np.sqrt((g ** 2).sum(axis=0))[np.newaxis, ...]
        E += weight * norm.sum()
        tau = 1. / (2.*ndim)
        norm *= tau / weight
        norm += 1.
        p -= tau * g
        p /= norm
        E /= float(image.size)
        if i == 0:
            E_init = E
            E_previous = E
        else:
            if np.abs(E_previous - E) < eps * E_init:
                break
            else:
                E_previous = E
        i += 1
    return out
    
def superr( scalDD, scalDS, scalSS, niter=15, beta=1.0e-8):
    """ 
    -    scalDS  which is  an array  [ZDIM,YDIM,XDIM]  , type "d" .
    -    scalDD  which is the total sum of the squared datas.
    -    scalSS  which is an array [XDIM,XDIM]  , type "d" .
    """
    ZDIM,YDIM,XDIM = scalDS.shape
    Volume = np.zeros( [ZDIM,YDIM,XDIM]  ,"f" )
    assert( scalSS.shape == (XDIM,XDIM))


    scalSS = scalSS.astype("f")
    scalDS = scalDS.astype("f")


    print( "  SHAPES ", scalDS.shape, scalSS.shape )

    
    Fista (  scalDD, scalDS, scalSS, Volume,niter=niter, beta=beta)
    
    return Volume

def calculate_grad( scalDD, scalDS , scalSS,   solution, grad) :
    grad [:]  = np.tensordot(  solution, scalSS, axes=[[-1],[-1]])
    err  = (grad*solution).sum()
    if scalDS is not None:
        err -= (scalDS*solution).sum()*2
        err +=  scalDD
        grad [:] -= scalDS
    return err/2
    
def    Fista( scalDD, scalDS , scalSS,  solution      , niter=500, beta=0.1 ):


    grad   = np.zeros_like(solution)
    grad2  = np.zeros_like(solution)
    x_old  = np.zeros_like(solution)
    y      = np.zeros_like(solution)

    err = 0.0
    err=calculate_grad( scalDD, scalDS , scalSS,   solution, grad) 
    for i in range(20):
        
        calculate_grad(None, None , scalSS,   grad, grad2)
        Lip = math.sqrt( np.linalg.norm(grad2/100000.0)   )*100000
        grad[:]  = grad2/ Lip
        if myrank==0:
            print( "LIP ", Lip)
    Lip = Lip*1.2

    t=1.0
    y[:] = solution
    x_old[:] = solution
    for iter in range(abs(niter)):
        err = calculate_grad(scalDD, scalDS , scalSS,   y, grad) 

        solution[:] =  y - grad/Lip
        
        # solution[:]=skimage.restoration.denoise_tv_chambolle(solution, weight=beta, eps=0.000002) 
        solution[:]=_denoise_tv_chambolle_nd(solution, weight=beta, eps=0.000002, positivity=True)


        ## solution[:] = np.maximum(solution, 0)

        tnew = ( 1+math.sqrt(1.0+4*t*t) )/2


        y[:] = solution +(t-1)/tnew *( solution - x_old )
        t = tnew
        if niter<0:
            t=1
        x_old[:] = solution
        if myrank==0:
            print( " errore est %e  mod_grad est  %e\n" % (  err, grad.std()) )










