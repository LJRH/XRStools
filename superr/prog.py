from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import h5py
import math
from six.moves import filter
from six.moves import range
myrank=0
import skimage.restoration

def main2():
    file = "demo.hdf5"
    group_delta  = "ROI_B_FIT8/scanXX/scans/Scan273/"
    group_sample = "ROI_B_FIT8/images/scans/"
    
    h5f = h5py.File(file,"r")
    h5  = h5f[group_delta]

    print( list(h5.keys()))
    sonde = {}

    width = 5
    XDIM = None
    roi_keys  = list(h5.keys()) 

    roi_keys = [60, 64, 35, 69, 34, 24, 5, 6, 71, 70, 39, 58, 56, 33]
    roi_keys = [60,64,71,70,69]
    roi_keys = [str(t) for t in roi_keys]

    tot_dizio={}

    for t in roi_keys:
        m = np.array(h5[t+"/matrix"][:],"d")
        mm = m.sum(axis=0)
        tot_dizio[t]=[mm,0]

        m=m[:(m.shape[0]//width)*width].reshape(-1, width, m.shape[1],m.shape[2]).sum(axis=1)/width
        sonde [t]  =  m
        if XDIM is None:
            XDIM = m.shape[0]
        else:
            assert(XDIM==m.shape[0])

    h5  = h5f[group_sample]

    zscan_keys =  sorted(  list(h5.keys())     , key = lambda x:  int(list(filter(str.isdigit, str(x) ))) ) 
    ZDIM = len(zscan_keys)
    
    m = h5[zscan_keys[0]][roi_keys[0]]["matrix"][:]
    YDIM = m.shape[0]
    
    scalDS = np.zeros( [ZDIM,YDIM,XDIM]  ,"f" )
    scalDD = 0.0
    scalSS = np.zeros( [XDIM,XDIM]  ,"f" )
    Volume = np.zeros( [ZDIM,YDIM,XDIM]  ,"f" )


    for rk in roi_keys:
        print( rk)
        probes      = sonde [rk]
        scalSS[:]  += np.tensordot( probes, probes, axes = [  [1,2], [1,2] ] )

    for iz in range(ZDIM):
        print( iz)
        zkey = zscan_keys[iz]
        for rk in roi_keys:
            m = np.array(h5[ zkey ][ rk ]["matrix"][:],"d")
            tot_dizio[rk][1]=tot_dizio[rk][1]+m

            probes = sonde [rk]
            assert( probes.shape[1:] == m.shape[1:])
            assert( XDIM == probes.shape[0] )
            assert( YDIM == m.shape[0]      )
            
            plane_contrib  = np.tensordot( m, probes, axes = [  [1,2], [1,2] ] ) 
            scalDS[iz] += plane_contrib

            scalDD     += (m*m).sum() 
    h5f.close()


    # import pickle
    # f=open("mats.pic","w")
    # pickle.dump(tot_dizio, f)
    # f.close()
# for n in d.keys():
#     print n
#     B=d[n][1]
#     A=d[n][0]
#     B=B.sum(axis=0)
#     pesiA = A.sum(axis=0)
#     pesiB = B.sum(axis=0)
#     medieA = (arange(A.shape[0])[:,None]*A).sum(axis=0)/pesiA
#     medieB = (arange(B.shape[0])[:,None]*B).sum(axis=0)/pesiB
#     plot(medieA)
#     plot(medieB)
#     show()
# 60 64 35 69 34 24 5 6 71 70 39 58 56 33
# 67 68  47 12 

# 20 21 44 45 40 1 3 4 59 19 53
# 28 2 15 31 52 55


    Fista (  scalDD, scalDS, scalSS, Volume,niter=15, beta=1.0e-8)
    
    h5f = h5py.File("Volume.h5","r+")
    h5f["Volume_onedectonly_15iters_beta1.0em8"] =  Volume
    h5f.close()

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
        solution[:]=skimage.restoration.denoise_tv_chambolle(solution, weight=beta, eps=0.000002) 


        ## solution[:] = np.maximum(solution, 0)

        tnew = ( 1+math.sqrt(1.0+4*t*t) )/2


        y[:] = solution +(t-1)/tnew *( solution - x_old )
        t = tnew
        if niter<0:
            t=1
        x_old[:] = solution
        if myrank==0:
            print( " errore est %e  mod_grad est  %e\n" % (  err, grad.std()) )





main()







