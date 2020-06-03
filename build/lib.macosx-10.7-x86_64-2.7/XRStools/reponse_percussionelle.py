from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import math
import scipy
import scipy.stats
import scipy.ndimage.filters as filt
import sys
import pickle
import os

from scipy.optimize import minimize
from six.moves import filter
from six.moves import range
from six.moves import zip
scipymin=1

from . import minimiser

import h5py

if(sys.argv[0][-12:]!="sphinx-build"):
    from . import luts_cy

try:
    from mpi4py import MPI
    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()
    procnm = MPI.Get_processor_name()
    comm = MPI.COMM_WORLD
    print( "MPI LOADED , nprocs = ", nprocs)
except:
    nprocs=1
    myrank = 0
    print( "MPI NOT LOADED ")


global indent 
indent = ""

def filterRoiList(l):
    return [t for t in l if t not in [ "motorDict"] ]



def get_LUT_1d( na ,nb, cp, cb , nref):
    # ecco un caso dove bisogna considerare che il pixel va da -0.5 a 0.5
    res=[]
    # print na, nb, cp, cb , nref
    for i in range(na):   ## qui sotto si pensa che lo zero e' all' inizio del pixel
        X0 = (cb+0.5) +(i-(cp+0.5))*nref
        X1 = (cb+0.5) +(i+1-(cp+0.5))*nref
        I0 = int(math.floor(X0))
        I1 = int(math.ceil(X1))
        for j in range(I0,I1):
            x0 = float(max(j,X0))
            x1 = float(min(j+1,X1))
            if j>=0 and j<nb:
                targetx0   =  ( x0-(cb+0.5) )/nref +(cp+0.5) -i  ## la posizione esatta nel pixel su cui cade
                targetx1   =  ( x1-(cb+0.5) )/nref +(cp+0.5) -i  ##  per maschera, sono compresi fra i e i+1
              
                res.append([i,j,(x1-x0)/nref ,  targetx0  ,targetx1 ])
    return res
            
def get_product(lut_1,lut_2, na2, nb2 , reponse_pixel):
    res=[]
    dim1,dim2 = reponse_pixel.shape
    for i1,j1,f1, y0,y1 in lut_1:

        Y0 = dim1*y0
        Y1 = dim1*y1
        iY0  = int(math.ceil(Y0))
        iY1  = int(math.floor(Y1))
        facts = (reponse_pixel[ iY0:iY1]).sum(axis=0)

        fiY0 = min(  iY0, Y1 )
        fiY1 = max(  iY1, Y0 )
        
        if  fiY0-Y0>1.0e-10 :
            facts[:] += (fiY0-Y0)  *  reponse_pixel[int(math.floor(Y0))]
        if(iY1>=iY0):
            if(Y1-fiY1>1.0e-10):
                facts[:] += (Y1-fiY1)  *  reponse_pixel[iY1]
            
        facts[:] /= (Y1-Y0) # media. Il fattore area c' e' gia in f1, f2
        
       
        for i2,j2,f2, x0,x1 in lut_2:
            # print j2
            X0 = dim2*x0
            X1 = dim2*x1
            iX0  = int(math. ceil(X0))
            iX1  = int(math.floor(X1))
            
            Fatt    = (facts[ iX0:iX1]).sum(axis=0)

            
            fiX0 = min(  iX0, X1 )
            fiX1 = max(  iX1, X0 )

            
            if ( fiX0-X0>1.0e-10):
                Fatt   += (fiX0-X0)  *  facts[int(math.floor(X0))]
            if(iX1>=iX0):
                if(X1-fiX1>1.0e-10):
                    Fatt  += (X1-fiX1)  *  facts[iX1]
                
            Fatt /= (X1-X0) # media. Il fattore area c' e' gia in f1, f2
            if Fatt>1.1:
                print( Fatt, " Fatt ")
                print( facts, Y0, Y1, iY0, iY1)
                print( "  X0,X1 " , X0,X1)
                print( " iX0, iX1 " ,iX0, iX1)
                print( "(facts[ iX0:iX1]).sum(axis=0)  " , (facts[ iX0:iX1]).sum(axis=0))
                print( " (fiX0-X0)  *  facts[int(math.floor(X0))] " ,(fiX0-X0)  *  facts[int(math.floor(X0))])
                print( " (X1-fiX1)  *  facts[iX1] " ,(X1-fiX1)  *  facts[iX1])

                
                raise
            newel =    [  i1*na2+i2   ,  j1*nb2+j2  , f1*f2* Fatt ] 
            res.append(newel)
    return res


def get_product4reponse(lut_1,lut_2, na2, nb2 , reponse_pixel, solution):
    res=[]
    dim1,dim2 = reponse_pixel.shape
    repnseq = np.arange(reponse_pixel.size)
    repnseq.shape = reponse_pixel.shape
    for i1,j1,f1, y0,y1 in lut_1:

        Y0 = dim1*y0
        Y1 = dim1*y1
        iY0  = int(math.ceil(Y0))
        iY1  = int(math.floor(Y1))
        facts = (reponse_pixel[ iY0:iY1]).sum(axis=0)
        
        recolte = np.array(repnseq[ iY0:iY1] )
        factrec = np.ones( repnseq[ iY0:iY1].shape   )

        fiY0 = min(  iY0, Y1 )
        fiY1 = max(  iY1, Y0 )
        
        
        if fiY0-Y0>1.0e-10 :
            where = int(math.floor(Y0))
            recolte = np.concatenate(  [recolte ,   [repnseq[where]]   ]   )            
            factrec = np.concatenate(  [factrec , [ [ (fiY0-Y0) ]  *  repnseq.shape[1]]   ]  )

        if(iY1 > Y0):            
            if( Y1-fiY1 > 1.0e-10 ):
                recolte = np.concatenate( [ [repnseq[iY1]],  recolte ]    )
                factrec = np.concatenate(  [ [[(Y1-fiY1)]*repnseq.shape[1]],     factrec  ]  )
            
            
        ### facts[:]   /= (Y1-Y0) # media. Il fattore area c' e' gia in f1, f2
        factrec[:] /= (Y1-Y0) # *dim1

        # print " =========== "
        # TERM = factrec.sum(axis=0)
        # print factrec.sum(axis=0)
        
        for i2,j2,f2, x0,x1 in lut_2:

            X0 = dim2*x0
            X1 = dim2*x1
            iX0  = int(math. ceil(X0))
            iX1  = int(math.floor(X1))
            
            recolteX = np.array(recolte  [ :, iX0:iX1] )
            factrecX = np.array(factrec  [ :, iX0:iX1] ) 

            # print " PEZZZO " , TERM[ iX0:iX1]

            
            fiX0 = min(  iX0, X1 )
            fiX1 = max(  iX1, X0 )
        
            if(fiX0-X0>1.0e-10):
                where = int(math.floor(X0))
                recolteX  = np.concatenate(  [ recolte[:, where:where+1]               ,  recolteX]    , axis=1 )
                factrecX  = np.concatenate(  [   (fiX0-X0) * factrec [:, where:where+1] ,  factrecX]    , axis=1  )
                # print " PEZZZO " , (fiX0-X0) * TERM[ where:where+1]
           
            if(iX1>=iX0):
                if(X1-fiX1>1.0e-10):
                    recolteX = np.concatenate(  [  recolteX  , recolte[:, iX1:iX1+1]                 ] , axis=1   )
                    factrecX = np.concatenate(  [  factrecX  ,  (X1-fiX1) * factrec[ :, iX1:iX1+1 ]   ] , axis=1 )
                    # print " PEZZZO " ,  (X1-fiX1) * TERM[ iX1:iX1+1]
                
                
            factrecX     /= (X1-X0) # *dim2  # media. Il fattore area c' e' gia in f1, f2


            # print " ---> " , factrecX.sum()
            if  (   factrecX.sum()>1.1   ) :
                print(    X0,X1,Y0,Y1)
                print(    iX0,iX1,iY0,iY1)
                print(    fiX0,fiX1,fiY0,fiY1)
                raise
            ## factrecX[:] = factrecX/factrecX.sum()

            
            for J,F in zip( recolteX.flatten()   ,   factrecX.flatten()  ):
                # print j1, nb2, j2
                newel =    [  i1*na2+i2   ,  J  , f1*f2* solution[j1,j2]*F      ] 
            # newel =    [  i1*na2+i2   ,  j1*nb2+j2  , f1*f2, Fatt ] 
                res.append(newel)
    return res
            
                


def get_LUT( mat_a,mat_b, center_pic, nref, reponse_pixel, doproduct = 1, soluzione=None, ROI = None):

    na1, na2  = mat_a.shape
    nb1, nb2  = mat_b.shape

    if ROI is None:
        ROI = np.ones_like(mat_a)
    ROI = ROI.astype("f")
    
    center_b = np.array(    [ (nb1-1)/2.0, (nb2-1)/2.0 ]   )

    lut_1 = get_LUT_1d( na1 ,nb1 , center_pic[0] ,center_b[0]  , nref)
    lut_2 = get_LUT_1d( na2 ,nb2 , center_pic[1] ,center_b[1]  , nref)
    if doproduct:
        if doproduct ==1:
            # LUT = get_product(lut_1,lut_2, na2, nb2 , reponse_pixel)
            
            # print " cy_product ", np.array(LUT).sum()
            #print np.array(lut_1,"f").shape
            #print np.array(lut_2,"f").shape
            # print " SHAPE "
            # print np.array(reponse_pixel,"f").shape, np.array(lut_1,"f").shape, np.array(lut_2,"f").shape, na2, nb2 , len( np.array(lut_1,"f")  ), len( np.array(lut_2,"f") ) 


            if len(lut_1)==0 or len(lut_2)==0:
                return None
            tmp = luts_cy.get_product( np.array(lut_1,"f"),
                                                np.array(lut_2,"f"),
                                                na2, nb2 ,
                                                np.array(reponse_pixel,"f"), ROI)
            LUT = np.array(tmp)
            
            
            #print LUT.sum()
            
        else :
            #print soluzione.shape
            #print mat_a.shape
            #print " ------------------ " 
            #print np.array(lut_2)[:,1].max()
            #print np.array(lut_1)[:,1].max()
            #print soluzione.shape
            if 0:
                LUT= get_product4reponse(  lut_1,
                                           lut_2,
                                           na2, nb2 ,
                                           reponse_pixel,
                                           soluzione)
            else:
                global SYMM_RESPO
                tmp = luts_cy.get_product4reponse(np.array(lut_1,"f"),np.array(lut_2,"f"), na1,
                                                           na2, nb2 , np.array(reponse_pixel,"f"), np.array(soluzione,"f"),
                                                           SYMM_RESPO, ROI)
                LUT = np.array(tmp)
                

            
        return LUT
    else:
        return lut_1, lut_2

def calculate_grad(grad,data ,solution , s2d, d2s, solution_shape , parallel = 0 , beta=0.1 ) :
    nsol, ndata = d2s.shape
    
    proj =   s2d.dot(solution)
    if data is not None:
        err  =   data - proj
    else:
        err  =    -proj
    fid  =   np.dot(err,err)
    grad[:] =   d2s.dot(err)

    #print solution_shape

    if parallel and nprocs>1:
        
        grad_res = np.zeros_like(grad)
        comm.Allreduce([ grad , MPI.FLOAT  ], [grad_res, MPI.FLOAT  ], op=MPI.SUM)
        grad[:]= grad_res

        a    = np.array( [fid]    ,  "f"  )
        a_r  = np.array( [0]    ,  "f"  )
        comm.Allreduce([ a , MPI.FLOAT  ], [a_r, MPI.FLOAT  ], op=MPI.SUM)
        fid = a_r[0]        
        
    sol_on_surf = np.reshape( solution, solution_shape )
    laplacian    = filt.laplace(sol_on_surf , mode='reflect')
    result =  fid/2.0 - beta*        np.dot( laplacian.flatten(), sol_on_surf.flatten()  ) /2
    grad[:] +=   beta * laplacian.flatten()


    # if myrank==0:
    #    print " fid ----<" , result, parallel
    return result


def    cg( data , solution   ,   s2d, d2s, solution_shape  ):
    maxdim = max(  s2d.shape[0], s2d.shape[1]   )
    auxs = [   np.zeros( [maxdim]  ,"f")  for i in range(5)                            ]

    nsol, ndata = d2s.shape
    
    grad=auxs[0][:nsol]
    grad_old=auxs[1][:nsol]
    p=auxs[2][:nsol]

    err = 0.0
    err=calculate_grad(grad,data ,solution , s2d, d2s, solution_shape  ) 
    # calculate_grad(grad, Volume, donnees, dim, n_row,n_col,  nnz, Bp, Bj, Bx, beta,auxs+3 ) ; 
      
    p[:] = grad[:]
    
    rrold=0.0
    rr   = 0.0
    
    rrold +=  np.dot( grad, grad  )
    rr=rrold;

    for iter in range(500):
        grad_old[:] = grad[:]
        calculate_grad(grad,None ,p , s2d, d2s, solution_shape  )
        pap = np.dot(  p, grad )

        solution[:] -= p*(rr/pap)

        err=calculate_grad(grad , data ,solution , s2d, d2s, solution_shape  ) 

        rrold=rr;
        
        rr = np.dot(  grad,   (grad-grad_old)  )

        
        beta = rr/rrold
        if beta<0 :
            beta = 0

        p[:] =  grad[:]+  p[:]*beta

        if myrank==0:
            if iter%1 ==0:
                sys.stdout.write(" "*60+"\r"+indent+( "iter %d errore est %e  mod_grad est  %e" % ( iter,  err,rr )    )  )
                sys.stdout.flush()
    if myrank==0:
        print( " ")


def    Fista( data , solution   ,   s2d, d2s, solution_shape , parallel = 0 , niter=500, beta=0.1 ):
    global indent
    maxdim = max(  s2d.shape[0], s2d.shape[1]   )
    auxs = [   np.zeros( [maxdim]  ,"f")  for i in range(5)                            ]

    nsol, ndata = d2s.shape
    
    grad =auxs[0][:nsol]
    grad2=auxs[1][:nsol]
    x_old=auxs[2][:nsol]
    y    =auxs[3][:nsol]

    err = 0.0
    err=calculate_grad(grad,data ,solution , s2d, d2s, solution_shape , parallel = parallel , beta=beta) 
    # calculate_grad(grad, Volume, donnees, dim, n_row,n_col,  nnz, Bp, Bj, Bx, beta,auxs+3 ) ; 


    if myrank==0:
        print( indent+"CALCULATING LIPSCHITZ FACTOR ")
    for i in range(100):
        calculate_grad(grad2,None ,grad , s2d, d2s, solution_shape , parallel = parallel , beta=beta)
        Lip = math.sqrt( np.linalg.norm(grad2/100000)   )*100000
        grad[:]  = grad2/ Lip
        if myrank==0:
            sys.stdout.write( " "*60+"\r"+indent+"LIP  %e"% Lip    ) 
            sys.stdout.flush()
    #if myrank==0:
        #print ""
    

    t=1.0
    y[:] = solution
    x_old[:] = solution
    for iter in range(abs(niter)):
        err = calculate_grad(grad,data ,y , s2d, d2s, solution_shape , parallel = parallel, beta=beta )
        solution[:] =  y + grad/Lip
        solution[:] = np.maximum(solution, 0)
        tnew = ( 1+math.sqrt(1.0+4*t*t) )/2
        y[:] = solution +(t-1)/tnew *( solution - x_old )
        t = tnew
        if niter<0:
            t=1
        x_old[:] = solution
        if myrank==0:
            if iter%1 ==0:
                # sys.stdout.write(" "*60+"\r"+indent+("FISTA iter %d  errore est %e  mod_grad est  %e" % ( iter,  err, grad.std()) ))
                sys.stdout.write(indent+("FISTA iter %d  errore est %e  mod_grad est  %e\n" % ( iter,  err, grad.std()) ))
                sys.stdout.flush()

    if myrank==0:
        print( " " )


        
def   trajectory_error( XYXY    ,  O_spots, opticalPSF,  nref, reponse_pixel  , retrieve_spots, suggerimento=None, ROI = None):
    if myrank==0: print( XYXY)

    Nspots = O_spots.shape[0]

    if(len(XYXY)==4):
        X1,Y1,X2,Y2 = XYXY
    else:
        X1,X2 = XYXY
        Y1    = suggerimento[0] + X1* suggerimento[1]
        Y2    = suggerimento[0] + X2* suggerimento[1]


    trajectory = Trajectory()
    trajectory.set_extrema(X1,Y1,X2,Y2, Nspots  )

    err = 0.0
    dd = 0.0
    ss = 0.0
    ds = 0.0

    total_I = []
    total_J = []
    total_F = []
    total_data = []

    for n, data in enumerate( O_spots  ):
        # print n
        center_pic = [   trajectory.Y.intercept + n* trajectory.Y.slope  ,      trajectory.X.intercept + n*  trajectory.X.slope         ]
        # lut1, lut2 =  get_LUT(data  ,  opticalPSF , center_pic, nref, reponse_pixel, doproduct = 0)
        LUT =  get_LUT(data  ,  opticalPSF , center_pic, nref, reponse_pixel, doproduct = 1, ROI=ROI)
        if LUT is not  None:        
            I,J,F = np.array(LUT).T
            assert(J.max()<    opticalPSF.size )
            total_I.append( n*data.size + I  )
            total_J.append( J  )
            total_F.append(F)
        total_data.append( data.flatten()  )

        # tmp = np.zeros(  [ data.shape[0] , opticalPSF.shape[1]   ],"f"  )
        # tmp2= np.zeros(  [ data.shape[0] ,       data.shape[1]   ],"f"  )
        # for i,j,f, dum1, dum2 in lut1:
        #     tmp [i,:] += opticalPSF[j,: ] *f
        # for i,j,f , dum1, dum2 in lut2:
        #     tmp2[:,i] += tmp[:,j]  *f 
        # ss += (tmp2*tmp2).sum()
        # dd += (data*data).sum()
        # ds += (tmp2*data).sum()

    #print " CONCATENO " 
    I = np.concatenate(total_I)
    J = np.concatenate(total_J)
    F = np.concatenate(total_F)
    data = np.concatenate(total_data )
    
    s2d_coo = scipy.sparse.coo_matrix( (F,(I,J)) , shape = [   O_spots.size   ,  opticalPSF.size   ])
    s2d =  s2d_coo.tocsr()
    sim = s2d.dot( opticalPSF.flatten())

    if retrieve_spots:
        return sim

    dd = (data*data).sum()
    ss=(sim*sim).sum()
    ds=(data*sim).sum()
    
    error = dd - ds*ds/ss
    
    ### f2 * ss + dd -2*ds*f    ==>  2f*ss -2 ds *f = 0
    ###  f = ds/ss        error =  ds*ds/ss +dd -2ds*ds/ss = dd - ds*ds/ss
    
    if myrank==0:
        sys.stdout.write( " "*60+"\r"+indent+"trajectory error %e" % error )
        sys.stdout.flush()
    return error

# class Tobject:
#     def __init__(self,  variables ,  O_spots, opticalPSF,  nref , reponse_pixel   ) :
#         self.variables      =  variables
#         self.O_spots        =  O_spots
#         self.opticalPSF     =  opticalPSF
#         self.nref           =  nref
#         self.reponse_pixel  =  reponse_pixel
        
#     def error(self):
#         X1 = self.variables[0].value
#         X2 = self.variables[1].value
#         Y1 = self.variables[2].value
#         Y2 = self.variables[3].value
        
#         res = trajectory_error(np.array([X1,Y1,X2,Y2]),  self.O_spots, self.opticalPSF,  self.nref , self.reponse_pixel,0  )
#         if myrank==0:
#             print " error = " , res
#         return res


def refine_trajectory(O_spots, opticalPSF, trajectory   , nref  , reponse_pixel , retrieve_spots = 0 , suggerimento = None , ROI = None):
    global indent

    Nspots = O_spots.shape[0]



    X1,Y1 = trajectory.get_coords(0)
    X2,Y2 = trajectory.get_coords(Nspots-1)

    if retrieve_spots:
        res = trajectory_error(np.array([X1,Y1,X2,Y2]),  O_spots, opticalPSF,  nref , reponse_pixel ,   1 , None, ROI)
        return res

    
    # if not scipymin or 0:
    #     variables = [ minimiser.Variable( X1,X1-2, X1+2),
    #                   minimiser.Variable( X2,X2-2, X2+2),
    #                   minimiser.Variable( Y1,Y1-2, Y1+2),
    #                   minimiser.Variable( Y2,Y2-2, Y2+2)
    #     ]
    #     tominimise=Tobject(variables ,  O_spots, opticalPSF,  nref , reponse_pixel  )
        
    #     miniob=minimiser.minimiser(tominimise,variables)
    #     miniob.amoeba(0.0001)  
        
    #     X1 = variables[0].value
    #     X2 = variables[1].value
    #     Y1 = variables[2].value
    #     Y2 = variables[3].value

    
    #     trajectory.set_extrema( X1,Y1,X2,Y2  , Nspots)
    #     return trajectory    

    if suggerimento is not None:

        if myrank==0:
            print( indent +"IMPROVING TRAJECTORY ")

        res = minimize(trajectory_error,np.array([X1,X2]),( O_spots, opticalPSF,  nref , reponse_pixel ,0 ,suggerimento , None, ROI),     # method='Powell',
                       method='Nelder-Mead',
                       options={'disp': False, 'maxiter': 40,
                                'return_all': False,
                                'maxfev': None, 'xtol': 0.001,
                       'ftol': 0.001 }  )
        if myrank==0:
            print( " ...  MINIMO IN ",    res.x,)
            print( " INIZIALE  ",    [X1,X2])
           #  print " FINALE   ",    res.x[0], res.x[1]

        trajectory.set_extrema_suggestion( res.x[0], res.x[1],   Nspots  , suggerimento)
       

    else:
        dodebug = False
        if myrank==0:
            print( indent + "IMPROVING TRAJECTORY ")
            dodebug = True
        res = minimize(trajectory_error,np.array([X1,Y1,X2,Y2]),( O_spots, opticalPSF,  nref , reponse_pixel ,0, None, ROI ),     # method='Powell',
                       method='Nelder-Mead',
                       options={'disp': dodebug, 'maxiter': 40, 'return_all': False,  'maxfev': None,  'ftol': 0.001})
                       # options={'disp': dodebug, 'maxiter': 40, 'return_all': False,  'maxfev': None, 'xtol': None, 'ftol': 0.001})
        # res = minimize(trajectory_error,np.array([X1,Y1,X2,Y2]),( O_spots, opticalPSF,  nref ))
        if myrank==0:
            print( " ...  MINIMO IN ",    res.x,)
            print( " INIZIALE  ",    [X1,Y1,X2,Y2])
            # print " FINALE   ",    res.x[0], res.x[1], res.x[2],  res.x[3]
            
        trajectory.set_extrema( res.x[0], res.x[1], res.x[2],  res.x[3],   Nspots)


    return trajectory    


def fit_reponse(trajectory, O_spots , nref, reponse_pixel, beta=0.1, niter=500, ROI = None )    :

    solution = np.zeros( [O_spots[0].shape[0]*nref, O_spots[0].shape[1]*nref    ] , "f" )
    total_I = []
    total_J = []
    total_F = []
    total_data = []

    for n, data in enumerate( O_spots  ):
        # if myrank==0:
        #    print n
        center_pic = [   trajectory.Y.intercept + n* trajectory.Y.slope ,      trajectory.X.intercept + n*  trajectory.X.slope         ]
        LUT =  get_LUT(data  ,solution , center_pic, nref, reponse_pixel, ROI=ROI)
        if LUT is not None:
            I,J,F = np.array(LUT).T
            total_I.append( n*data.size + I  )
            total_J.append( J  )
            total_F.append(F)
        total_data.append( data.flatten()  )

    # if myrank==0:
    #     print " CONCATENO "
    
    I = np.concatenate(total_I)
    J = np.concatenate(total_J)
    F = np.concatenate(total_F)
    data = np.concatenate(total_data )

    # if myrank==0:
    #     print " CREO MATRICE "
    #     print data.shape
    #     print I.shape
    #     print J.shape
        
    s2d_coo = scipy.sparse.coo_matrix( (F,(I,J)) , shape = [   O_spots.size   ,  solution.size   ]    )
    d2s_coo = scipy.sparse.coo_matrix( (F,(J,I)) , shape = [   solution.size  ,  O_spots.size  ]    )


    # if myrank==0:
    #     print " CAMBIO FORMATO "
    
    s2d =  s2d_coo.tocsr()
    d2s =  d2s_coo.tocsr()

    # if myrank==0:
    #     print " FISTA  "
        
    b=solution[:]
    solution_shape = solution.shape
    solution=np.reshape(solution,[-1])
    # cg( data , solution    ,   s2d, d2s, solution_shape  )
    Fista( data , solution    ,   s2d, d2s, solution_shape, niter=niter, beta=beta )
    solution.shape = solution_shape
    return solution



def fit_reponse_pixel(trajectory_list, O_spots_list , nref, reponse_pixel, solution_list, niter=500, beta=100 , rois = None):


    if type(trajectory_list)!=type([]):
        trajectory_list = [trajectory_list]
        O_spots_list = [O_spots_list]
        solution_list = [solution_list]
        
    total_I = []
    total_J = []
    total_F = []
    total_data = []

    tot_size_data = 0
    for trajectory , O_spots, solution, ROI  in zip(trajectory_list,O_spots_list, solution_list, rois ) :

        if trajectory is None:
            continue

        assert( solution.shape == (O_spots[0].shape[0]*nref, O_spots[0].shape[1]*nref ) )


        for n, data in enumerate( O_spots  ):
            # if myrank==0:
        #     print n
            center_pic = [   trajectory.Y.intercept + n* trajectory.Y.slope ,      trajectory.X.intercept + n*  trajectory.X.slope         ]
            LUT =  get_LUT(data  ,solution , center_pic, nref, reponse_pixel, soluzione = solution, doproduct = 2, ROI = ROI )
            if LUT is not None:
                I,J,F = np.array(LUT).T
                total_I.append( tot_size_data+ n*data.size + I  )
                total_J.append( J  )
                total_F.append(F)
            total_data.append( data.flatten()  )
            
        tot_size_data +=O_spots .size
        
    if myrank==0:
        print( " CONCATENO ", )
    I = np.concatenate(total_I)
    J = np.concatenate(total_J)
    F = np.concatenate(total_F)
    data = np.concatenate(total_data )

    # if myrank==0:
    #     print " CREO MATRICE "
    #     print data.shape
    #     print I.shape
    #     print J.shape

    s2d_coo = scipy.sparse.coo_matrix( (F,(I,J)) , shape = [   tot_size_data   ,  reponse_pixel.size   ]    )
    d2s_coo = scipy.sparse.coo_matrix( (F,(J,I)) , shape = [   reponse_pixel.size  ,  tot_size_data  ]    )


    # if myrank==0:
    #     print " CAMBIO FORMATO "
        
    s2d =  s2d_coo.tocsr()
    d2s =  d2s_coo.tocsr()

    datasim = s2d.dot( reponse_pixel.flatten()   )
    datasim = datasim - data

    if myrank==0:
        print( " FIDELITY ", (datasim*datasim).sum(),)
        print( " FISTA  ",)
        
    reponse_pixel_shape = reponse_pixel.shape
    reponse_pixel=np.reshape(reponse_pixel,[-1])
    
    Fista( data , reponse_pixel    ,   s2d, d2s, reponse_pixel_shape , parallel = 1 , niter=niter, beta=beta )
    reponse_pixel.shape = reponse_pixel_shape
    return reponse_pixel



def get_spots_list(  filename  , groupname , filter_rois =1 ):
    res=[]
    nomi_scan = []
    stats   = []
    origini = {}
    
    h5 = h5py.File(filename,"r")

    xscales = {}
    enescan = None
    
    for sn in filterRoiList(list(h5[groupname].keys())):
        
        print( groupname+"/"+sn+"/matrix")
        m  = h5[groupname+"/"+sn+"/matrix"][:]
        print( m.shape)

        if groupname+"/"+sn+"/xscale" in h5:
            xscales[sn] = h5[groupname+"/"+sn+"/xscale"][:]
        else:
            xscales[sn] = np.arange( m.shape[0]).astype("f") 


        
        if groupname+"/motorDict/energy" in h5:
            enescan = h5[groupname+"/motorDict/energy"].value

        
        if m.shape!=(0,):
            stat = m.sum(axis=-1).sum(axis=-1)
            
            if (not filter_rois) or (stat.min()>stat.max()/2.0 and stat.min() >250.0):
                res.append(m)
                nomi_scan.append(  sn   )
                stats.append(stat)
                origini[sn] =  h5[groupname+"/"+sn+"/cornerpos"][:]
                
    h5.close()
    rois = []
                                    ## METTERE D UFFICIO ANCHE LE ROI QUANDO SI SCRIVE UNO SCAN E CAMBIARE QUI IL ALLINDIETRO
    for i in range(1):
        pos = groupname.rfind("/")
        groupname=groupname[:pos]

    # groupname=groupname+"/rois_definition/rois_dict/"

    print( " NOMI SCAN ", nomi_scan)
    
    h5 = h5py.File(filename,"r")
    for sn in nomi_scan:

        print( " DATANAME ", groupname+"/ROI%02d/mask" %int(sn), filename)
        m  = h5[groupname+"/ROI%02d/mask" %int(sn)][:]
        rois.append(m)
    h5.close()

    print( res)
    for m in res:
        print( m.sum())
    print( nomi_scan)
    print( rois)
    return res, nomi_scan, stats, origini, rois, xscales, enescan

   
class pippo:
    pass

class Trajectory :
    def __init__(self):
        self.X=pippo()
        self.Y=pippo()

    def get_coords(self,n):
        X1 = self.X.intercept + self.X.slope*n
        Y1 = self.Y.intercept + self.Y.slope*n
        return X1, Y1
    def set_extrema(self, X1, Y1, X2, Y2,N):
        self.N = N
        self.X.intercept = X1
        self.Y.intercept = Y1
        self.X.slope     = (X2-X1)/(N-1)
        self.Y.slope     = (Y2-Y1)/(N-1)
        
    def set_extrema_suggestion( self, X1, X2,   Nspots  , suggerimento):
        Nspots = self.N
        Y1    = suggerimento[0] + X1* suggerimento[1]
        Y2    = suggerimento[0] + X2* suggerimento[1]            
        self. set_extrema( X1, Y1, X2, Y2,Nspots)
        

    def  projectOnHint( self, suggerimento  ):
        Nspots = self.N

        X1,Y1 = self.get_coords(0)
        X2,Y2 = self.get_coords(Nspots-1)
        Y1    = suggerimento[0] + X1* suggerimento[1]
        Y2    = suggerimento[0] + X2* suggerimento[1]
            
        self. set_extrema( X1, Y1, X2, Y2,Nspots)
    def __repr__(self):
        s=""
        s=s+" X.intercept = %e " % self.X.intercept
        s=s+" X.slope     = %e " % self.X.slope
        s=s+" Y.intercept = %e " % self.Y.intercept
        s=s+" Y.slope     = %e " % self.Y.slope
        s=s+" Nspots     = %d " % self.N
        
        return s


class NoSignal(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
def get_Trajectory_byregression( O_spots    ):
    Xs=[]
    Ys=[]
    Ws=[]
    for slice in O_spots:
        mmax  = slice.max()
        thres = mmax / 20
        slice = np.array(slice)
        slice[ slice<thres  ] = 0

        m0 = slice.sum()
        if m0>0:
            X  = (  slice * np.arange(slice.shape[1] )  ).sum()/m0
            Y  = (  slice.T * np.arange(slice.shape[0] )  ).sum()/m0
            Xs.append(X)
            Ys.append(Y)
            Ws.append(m0)
    if  len(Xs)<=1:
        raise NoSignal(" All images were zero  ")
    
    ret_val = Trajectory()

    ret_val.X.slope, ret_val.X.intercept = np.polyfit( np.arange(len(Xs)), Xs ,1,  w = Ws)
    ret_val.Y.slope, ret_val.Y.intercept = np.polyfit( np.arange(len(Ys)), Ys ,1,w= Ws )

    
    # ret_val.X.slope, ret_val.X.intercept, Xr_value, Xp_value, Xstd_err = scipy.stats.linregress(np.arange(len(Xs)),Xs , )
    # ret_val.Y.slope, ret_val.Y.intercept, Yr_value, Yp_value, Ystd_err = scipy.stats.linregress(np.arange(len(Ys)),Ys , )
    ret_val.N = len(O_spots)
    
    return ret_val


def DOFIT(filename=None, groupname=None, nref=5, niter_optical=500, beta_optical=0.1 ,
          beta_pixel=1000.0, niter_pixel = -20,
          niter_global  = 50, pixel_dim=6, simmetrizza=1, do_refine_trajectory=1, target_file="responses.h5",target_groupname="FIT",
          trajectory_reference_scansequence_filename = None, trajectory_reference_scansequence_groupname = None ,  trajectory_threshold = None,
          trajectory_file=None  , filter_rois=1 , fit_lines = False):

    ###### filename = "../nonregressions/demo_imaging.hdf5"
    ###### groupname = "ROI_B/foil_scanXX/scans/Scan273/"
    # print filename,  groupname

    O_spots_list, nomi_rois, stats, origini_foil, rois, xscales, energy = get_spots_list( filename,  groupname,  filter_rois= filter_rois )
    stats = np.array(stats)
    
    if do_refine_trajectory ==2:

        h5f = h5py.File(trajectory_reference_scansequence_filename,"r")
        h5  = h5f[trajectory_reference_scansequence_groupname]

        zscan_keys =  sorted(  list(filterRoiList(list(h5.keys())))     , key = lambda x:  int(list(filter(str.isdigit, str(x) ))) )
        ZDIM = len(zscan_keys)
        myrange = np.array_split( np.arange( ZDIM ), nprocs   )[myrank]
        myZDIM =   len(myrange)
        
        roi_keys  = nomi_rois
        integrated     = {}
        origini_ref    = {}
        for rk in roi_keys:
            zkey = zscan_keys[ 0 ]
            m  = h5[zkey][rk]["matrix"][:]
            
            mm = m.sum(axis=0)
            integrated[rk ] = np.zeros_like( mm ) .astype("d")
            origini_ref   [rk ] = h5[zkey][rk]["cornerpos"][:]

        for iz in range(myZDIM):        
            zkey = zscan_keys[myrange[iz] ]
            for rk in roi_keys:
                m = h5[ zkey ][ rk ]["matrix"][:]
                msum = m.sum(axis=0)
                integrated[rk] = integrated[rk]+msum
                
        if nprocs>1:    
            for n in  integrated.keys():
                    comm.Allreduce(  [np.array(integrated[n]), MPI.DOUBLE], [integrated[n], MPI.DOUBLE],  op=MPI.SUM)

        trajectory_hints={}
        for n in  integrated.keys():
            C     = integrated[n]
            pesiC = C.sum(axis=0)
            medieC = (np.arange(C.shape[0])[:,None]*C).sum(axis=0)/pesiC
            coords  = np.arange(len( medieC  )) 

            maxP=pesiC.max()            
            whereBigger = np.where( pesiC > maxP*trajectory_threshold   )[0]
            
            inizio = whereBigger.min()
            fine   = whereBigger.max()
            pesifit = np.zeros( len(pesiC),"f")
            pesifit[inizio:fine+1]   =   ( pesiC[inizio:fine+1]> maxP*trajectory_threshold  )*1
            pfit=np.polynomial.polynomial.polyfit(coords, medieC, 1,  w=pesifit)
            tmp_intercept = pfit[0]
            tmp_slope     = pfit[1]
            
            DY = origini_foil[n][0] -origini_ref[n][0]
            DX = origini_foil[n][1] -origini_ref[n][1]

            hint_slope = tmp_slope
            hint_intercept  = tmp_intercept - DY  +  hint_slope*DX
            
            trajectory_hints[n] = hint_intercept,hint_slope
    else:
        trajectory_hints={}

    global SYMM_RESPO
    SYMM_RESPO = simmetrizza

    reponse_pixel = np.ones([ pixel_dim , pixel_dim  ],"f")
    
    #### nref = 5
    
    trajectory_list = []
    solution_list   = []
    trajectories_from_file = None
    if trajectory_file is not None:
        ### file = open(trajectory_file,"r")
        print( "RELOADING TRAJECTORIES ")
        trajectories_from_file = reload_trajectories(trajectory_file,  nomi_rois)

    
    for iterm,  (O_spots, name, ROI)  in enumerate(zip(O_spots_list, nomi_rois, rois)):
        # print "MYRANK %d fa "%myrank, iterm
        if (iterm)%nprocs == myrank:

            if trajectories_from_file is not None:
                " REUSING TRAJECTORIES "
                trajectory = trajectories_from_file[iterm]
                trajectory.N = O_spots.shape[0]
            else:
                try:
                    trajectory = get_Trajectory_byregression( O_spots    )
                except NoSignal as inst:
                    msg = " Was not able to determine trajactory because not enough signal in ROI %s"%name
                    print(msg)
                    inst.value = msg
                    raise

            print( " TRAJECTORY ", trajectory)

            if  do_refine_trajectory  ==2 :      
  
                PRIMA = trajectory.Y.intercept, trajectory.Y.slope
                print( "============= projectOnHint============= ", )


                trajectory.projectOnHint(  trajectory_hints[  name  ]  )
                
  
                print( "Intercept : ", PRIMA[0] , "===>", trajectory.Y.intercept , "Slope : ", PRIMA[1] , "===>", trajectory.Y.slope)

            if myrank ==0:
                print( "Process 0 FITTING OPTICAL RESPONSE for scan ", name)

            solution =  fit_reponse( trajectory, O_spots , nref  , reponse_pixel,niter=niter_optical, beta= beta_optical, ROI=ROI)   ## niter 500    beta_optical 0.1

        else:
            solution = None
            trajectory = None
       
        solution_list.append(solution)
        trajectory_list.append(trajectory)

    if myrank ==0:
        print( "NOW ALL PROCESS FITTING  PIXEL RESPONSE")

    if pixel_dim>1:
        newreponse = fit_reponse_pixel(trajectory_list , O_spots_list , nref, reponse_pixel, solution_list ,beta = beta_pixel  , niter  = niter_pixel , rois=rois)    
    global indent 
    indent = ""

    for iter in range(niter_global):
        indent = "      "
        if myrank ==0:
            print( "ITERATIVE REFINEMENT, cycle No", iter)


        print(  nomi_rois)
        print(  O_spots_list)
        for iterm, (O_spots, name, ROI )  in enumerate(zip(O_spots_list, nomi_rois, rois)):
            if (iterm)%nprocs == myrank:
                # print " MYRANK %d fa "%myrank, iterm
                trajectory =  trajectory_list[iterm]

                if myrank ==0:
                    print( "  Process 0 FITTING OPTICAL RESPONSE for roi ", name)


                solution =  fit_reponse( trajectory, O_spots , nref  , reponse_pixel,niter=niter_optical, beta=beta_optical, ROI=ROI)

                if do_refine_trajectory:
                    suggerimento = None
                    if do_refine_trajectory  ==2 :  
                        suggerimento = trajectory_hints[name]

                    if myrank ==0:
                        print( "  Process %d REFINING TRAJECTORY for roi "%myrank, name)

                    trajectory = refine_trajectory(O_spots, solution , trajectory   , nref  , reponse_pixel , suggerimento = suggerimento, ROI=ROI )
                    #####solution =  fit_reponse( trajectory, O_spots , nref  , reponse_pixel)
            else:
                solution = None
                trajectory = None

            solution_list[iterm]=  solution
            trajectory_list[iterm] = trajectory

        if nprocs>1:
            # print " process ", myrank, " aspetta " 
            comm.Barrier()
            # print " TUTTI I PROCESSI SONO ARRIVATI QUI "
            #for i in range(nprocs):
            #    if i==myrank:
            #        print (" %d " %i )*10 
            #        for t in trajectory_list:
            #            if t is not None:
            #                print t
            #    comm.Barrier()

        if myrank ==0:
            print( "REFITTING PIXEL RESPONSE ",)
            
        # import pickle
        # file = open("tra%d"%myrank,"r")
        # # pickle.dump( trajectory_list, file  )
        # trajectory_list = pickle.load(  file  )
        # file.close()


        if pixel_dim>1:
            newreponse = fit_reponse_pixel(trajectory_list, O_spots_list , nref, reponse_pixel, solution_list ,beta = beta_pixel, niter  =  niter_pixel , rois = rois)
        # reponse_pixel=newreponse
    
    if myrank ==0:
        print( "FINISHED. NOW WRITING ")

    trajectory_dic = {}
    for i,t in enumerate(trajectory_list):
        if t is not None:
            trajectory_dic[i] = t
    tosend = pickle.dumps(trajectory_dic)
        
    for iproc in range(nprocs):
        if iproc:
            if myrank==iproc:
                comm.send(  tosend, dest=0, tag=11)
            elif myrank==0:
                data = comm.recv(source=iproc, tag=11)
                tdic = pickle.loads(data)
                trajectory_dic.update(tdic)



        if fit_lines :
            lines = []
            for solution in solution_list :
                if solution is None:
                    lines.append(None)
                    continue
                
                weigths = np.array(solution.flat)
                indici  = np.indices(solution.shape)
                                     
                ind_y   = np.array(indici[0].flat)
                ind_x   = np.array(indici[1].flat)
                print( " SHAPES ", ind_x.shape, ind_y.shape, weigths.shape)

                pfit = np.polynomial.polynomial.polyfit(ind_x.astype("d"),ind_y.astype("d"),1,w= weigths.astype("d")  )

                lines.append( [pfit[0], pfit[1] ]    )
                
                                     
                

                
        if nprocs>1:
            comm.Barrier()
        if iproc==myrank:
            if myrank==0:
                h5f = h5py.File(target_file,"a")
                
                h5  = h5f.require_group(target_groupname)
                if "response" in  h5:
                    del  h5["response"]
                h5["response"] = reponse_pixel

                if energy is not None:
                    h5.require_group("motorDict")
                    if "energy" in h5["motorDict"]:
                        del h5["motorDict"]["energy"]
                    h5["motorDict/energy"] = energy
                
            else:
                os.system("touch %s" %target_file )
                h5f = h5py.File(target_file,"a")
                h5  = h5f[target_groupname]

            for iterm, (O_spots, name, stat)  in enumerate(zip(O_spots_list, nomi_rois, stats)):
                if (iterm)%nprocs == myrank:
                    trajectory =  trajectory_list[iterm]
                    solution =    solution_list[iterm]

                    if name in h5:
                        del h5[name]

                    h5group  =  h5.require_group(name )
                    README=""
                    h5group[ "data"   ] = solution
                    README +="data : The reponse which fits all the respones if translated as below\n"

                    h5group["Xintercept"] = trajectory.X.intercept 
                    README +="Xintercept : The X-shift of the fitted response for response numero 0\n"

                    h5group["Xslope"    ] = trajectory.X.slope 
                    README +="Xslope : The step of X-shift of the fitted response going from response n + n+1\n"

                    h5group["Yintercept"] = trajectory.Y.intercept
                    README +="Yintercept : same as for Xintercept but for Y \n"
                    
                    h5group["Yslope"    ] = trajectory.Y.slope
                    README +="Yslope : same as for Xslope but for Y \n"

                    h5group["stat"    ] = stat
                    README +="stat : this relates to original responses. it is the total energy contained in reference n for all the n's\n"

                    h5group["nref"    ] = nref
                    README +="nref : how many pixel of the fitted response we have for one pixel of the original response. this rescale the response but not all the other coordinates which are referred to the detector,\n"
                    
                    h5group["xscale"    ]   = xscales[name]
                    README +="xscale : the scan varying variable values for the range over which the fit has been done\n"
                    
                    if fit_lines:
                        h5group["line"    ] = lines[iterm]
                        README +="nref : how many pixel of the fitted response we have for one pixel of the original response\n"
                        h5group["python_plot_imageAndLine"]=("fig,ax=pylab.plt.subplots()\n"
                                                             "ax.imshow(self.data)\n"
                                                             "xs= pylab.arange(self.data.shape[1])\n"
                                                             "ys= xs*self.line[1]+self.line[0]\n"
                                                             "ax.plot(xs,ys)\n"
                                                             "pylab.plt.show()\n")
                        README +="python_plot_imageAndLine : double click on it to execute the code and have a plot\n"
                                                             
                    h5group["README"    ] = README

                    
            h5f.flush()        
            h5f.close()
            h5f = None
            
    # if myrank == 0:
    #     tl = [trajectory_dic[i]  for i in range( len(trajectory_list) )]
    #     file = open( "trajectories.pickle"  ,"w")
    #     pickle.dump(   tl   , file   )
    #     file.close()

        
def reload_trajectories(trajectory_file, nomi_rois) : ### , trajectory_file_group, nomi_rois):
    
    h5 =  h5py.File( trajectory_file   ,"r")
    ### h5group = h5[trajectory_file_group]
    trajectory_list = []
    for iterm,  name  in enumerate( nomi_rois[:]):
        h5group = h5[name]
        trajectory = Trajectory()
        trajectory.X.intercept  = h5group["Xintercept"].value
        trajectory.X.slope      = h5group["Xslope"].value
        trajectory.Y.intercept  = h5group["Yintercept"].value
        trajectory.Y.slope      = h5group["Yslope"].value
        # trajectory.N            = O_spots.shape[0]
        trajectory_list.append( trajectory ) 
    return trajectory_list


    


def DOROIS(filename = "../nonregressions/demo_imaging.hdf5" , groupname = "ROI_B/foil_scanXX/scans/Scan273/",
           roisgroupname = "ROI_B/",
           target_filename="newscan.h5", 
           roisgroupname_target = "ROI_B_FIT8/",
           newscanstarget    = "scanXX/scans/Scan273",
           responsefilename =  "responses.h5",
           responsepath =  "fit",
           nex = 3,
           filter_rois=1,
           recenterings= None,
           filterMask = None):
    
    O_spots_list, nomi_rois, stats , origini, rois, xscales, energy= get_spots_list( filename,  groupname ,  filter_rois= filter_rois  )

    trajectory_list = []
    solution_list   = []
    h5f =  h5py.File( responsefilename   ,"r")
    print( responsefilename)
    h5 = h5f[responsepath]
    print( responsepath)
    
    for iterm,  (O_spots, name)  in enumerate(zip(O_spots_list[:], nomi_rois[:])):


        if (iterm)%nprocs == myrank:

            if myrank ==0 :
                print( "Process 0 reading responses and trajectories for roi ",  name)

            h5group = h5[name]
            print( name, "data")
            solution =  h5group["data"][:]
            trajectory = Trajectory()
            trajectory.X.intercept  = h5group["Xintercept"].value
            trajectory.X.slope      = h5group["Xslope"].value
            trajectory.Y.intercept  = h5group["Yintercept"].value
            trajectory.Y.slope      = h5group["Yslope"].value
            nref                    = h5group["nref"].value
            trajectory.N            = O_spots.shape[0]

            if recenterings is not None:
                rec = recenterings[int(name)]
                recy,recx = rec
                irecy, irecx = int(round(recy)), int(round(recx))
                recy,recx  = recy-irecy, recx-irecx
            else:
                recy,recx   = 0,0
                irecy,irecx = 0,0

            trajectory.X.intercept += -recx
            trajectory.Y.intercept += -recy

            if "line" in h5group:
                trajectory.line = h5group["line"].value
                lh,lslope =trajectory.line
                lh = lh-recy + recx*lslope
                trajectory.line = np.array([lh,lslope])
            
            trajectory.recy = recy
            trajectory.recx = recx
            trajectory.irecy = irecy
            trajectory.irecx = irecx

        else:
            solution = None
            trajectory = None
        solution_list.append(solution)
        trajectory_list.append(trajectory)

    reponse_pixel =  h5["response" ][:]
    h5f.close()
    
    del solution
    del trajectory
    
    masklist = []
    for iterm, (O_spots, name, ROI)  in enumerate(zip(O_spots_list[:], nomi_rois[:], rois)):

        if (iterm)%nprocs == myrank:
            # print " MYRANK fa ",  myrank, iterm, solution_list[iterm]

            if myrank ==0 :
                print( "Process 0 expanding roi ",  name)

            assert(    (np.array(solution_list[iterm].shape) -    np.array(O_spots[0].shape)*nref).sum() == 0           )
            
            trajectory =  trajectory_list[iterm]
            
            # trajectory = refine_trajectory(O_spots, solution_list[iterm] , trajectory   , nref  , reponse_pixel  )
            # solution =  fit_reponse( trajectory, O_spots , nref  , reponse_pixel,niter=500, beta=0.001)
            ## solution_list[iterm] = solution


            # print( filename)
            # print( roisgroupname+ "/rois_definition/image")
            
            h5f = h5py.File(filename,"r")
            
            # imagealldect =  h5f [roisgroupname+ "/rois_definition/image"][:]
            imagealldect =  h5f [roisgroupname+ "/image"][:]

            # h5 = h5f [roisgroupname+  ("/rois_definition/rois_dict/ROI%02d"%int(name)) ]
            h5 = h5f [roisgroupname+  ("/ROI%02d"%int(name)) ]
            mask   = h5["mask"  ][:]
            origin = h5["origin"][:]
            h5f.close()

            if nex==0:
                newmask = mask
                neworigin=origin
                
            h,w = mask.shape
            del mask
            nspots = trajectory.N
            sezione = [int(nspots*(nex)) ,int(nspots*(nex))+nspots ]
            if nex:

                newmask = np.zeros([ (1+4*nex)*h,(1+4*nex)*w ]   )

                # print " CI SONO nspots" ,  nspots


                for i in range(  - int(nspots*(nex+0.5)),  int(nspots*(1+nex+0.5))):
                    X1,Y1 = trajectory.get_coords(i)
                    DX = int( trajectory.X.slope  ) 
                    X1shift = X1 + 2*nex*w
                    Y1shift = Y1 + 2*nex*h
                    newmask[ max(0,int(Y1shift)-h//2) :  min(int(Y1shift)+h//2+1, newmask.shape[0]) , max(0,int(X1shift-1+DX)):min(int(X1shift+1+DX)+1 , newmask.shape[1])  ] =1 
                yinds, xinds = np.where(newmask)

                ymin,ymax = yinds.min(), yinds.max()
                xmin,xmax = xinds.min(), xinds.max()

                # print " VA DA " , ymin,ymax, xmin,xmax 

                newmask = newmask[  ymin:ymax+1,   xmin:xmax+1  ]
                SHIFTORIGIN = np.array(  [ (ymin-2*nex*h), (xmin-2*nex*w) ] )
                neworigin = origin + SHIFTORIGIN

                if neworigin[0] <0:
                    newmask = newmask[-neworigin[0]:]
                    neworigin[0]=0
                    
                if neworigin[1] <0:
                    newmask = newmask[ : , :-neworigin[1] ]
                    neworigin[1]=0
                
                
                corry  =  neworigin[0] + newmask.shape[0] -   imagealldect.shape[0]
                corrx  =  neworigin[1] + newmask.shape[1] -   imagealldect.shape[1]
                if corry>0:
                    newmask = newmask[:-corry]
                if corrx>0:
                    newmask = newmask[:,:-corrx]

                new_Ospots  = np.zeros( [  nspots*(1+2*nex)  ,   newmask.shape[0], newmask.shape[1]   ] ,"f"     )
                trajectory.X.intercept =  trajectory.X.intercept - nex*nspots*trajectory.X.slope  - SHIFTORIGIN[1]
                trajectory.Y.intercept =  trajectory.Y.intercept - nex*nspots*trajectory.Y.slope  - SHIFTORIGIN[0]


                if filterMask is not None:
                    newmask[:] =  newmask * filterMask[ neworigin[0]:neworigin[0]+newmask.shape[0],  neworigin[1]:neworigin[1]+newmask.shape[1]   ]

                
                if nex>0:
                    ROI=None
            else:
                new_Ospots = O_spots
                    
            expandedSpots = refine_trajectory(new_Ospots, solution_list[iterm] , trajectory   , nref  , reponse_pixel , retrieve_spots = 1, ROI=ROI )
            expandedSpots.shape = new_Ospots.shape
            
        else:
            newmask   = None
            neworigin = None
            expandedSpots = None
            sezione = None
            
        masklist .append(   [newmask,neworigin, expandedSpots , sezione , trajectory   ]       )

    if target_filename!=filename:
        Oh5f = h5py.File(filename,"r")
        
    for iterm, (O_spots, name)  in enumerate(zip(O_spots_list[:], nomi_rois[:])):
        if nprocs>1:
            comm.Barrier()
        if (iterm)%nprocs == myrank:
            newmask,neworigin, newspots, sezione, trajectory = masklist[iterm]
            solution =    solution_list[iterm] 
            if(iterm==0):  h5f = h5py.File(target_filename,"a")
            else         :  h5f = h5py.File(target_filename,"a")
            if target_filename==filename:
                Oh5f = h5f

            if iterm ==0 and (roisgroupname_target in h5f):
                del h5f[roisgroupname_target]
                
            h5f.require_group( roisgroupname_target )
            print( " h5g " , target_filename , roisgroupname_target)
            h5g = h5f [roisgroupname_target]
            h5g.require_group( "rois_definition"  )
            h5g.require_group( "rois_definition/rois_dict/"  )
            h5g.require_group( "rois_definition/rois_dict/ROI%02d"%int(name)  )

            if iterm==0: 
                h5g [ "rois_definition/image"] = imagealldect
            
            h5 = h5g [ "rois_definition/rois_dict/ROI%02d"%int(name)   ]
            h5["mask"]   = newmask
            h5["origin"] = neworigin - np.array( [trajectory.irecy, trajectory.irecx] )
            h5["sezioneold"] = sezione
            print( "newscanstarget  " , newscanstarget)

            h5g.require_group( newscanstarget  )
            h5 = h5g [newscanstarget   ]

            
            h5.require_group( str(name)  )
            h5 = h5 [str(name)   ]
            h5["matrix"] = newspots
            h5["Xintercept"] =trajectory.X.intercept
            h5["Yintercept"] =trajectory.Y.intercept
            h5["Xslope"] =trajectory.X.slope
            h5["Yslope"] =trajectory.Y.slope
            h5["monitor"] = np.ones([newspots[0].shape[0]],"d")
            h5["monitor_divider"] = 1.0




            if hasattr(trajectory,"line"):
                h5["line"] =trajectory.line
                h5[ "optical_response" ] = solution
                h5["nref"] = nref

            Oh5  = Oh5f[groupname]
            Oh5g = Oh5[name]

            for key in Oh5g.keys():
                if key=="cornerpos":
                    h5[key] = neworigin - np.array( [trajectory.irecy, trajectory.irecx] )
                elif key!="matrix" and key not in h5 :
                    print( " aggiungo " , key)
                    h5.copy(  Oh5g[key],  key)

            if iterm == 0 :
                if "motorDict" in Oh5:
                    print( newscanstarget)
                    print( list(h5g [newscanstarget ].keys()))
                    print( list(Oh5.keys()))
                    
                    h5g [newscanstarget ].copy( Oh5 ["motorDict"],  h5g [newscanstarget]    )
                    
            h5f.close()
        if nprocs>1:
            comm.Barrier()

    if target_filename!=filename:
        Oh5f.close()


if __name__ == "__main__":

    # todo="fit"
    todo="rois"

    if todo=="fit":

        filename = "../nonregressions/demo_imaging.hdf5"
        groupname = "ROI_B/foil_scanXX/scans/Scan273/"
        nref=5
        niter_optical=500
        beta_optical=0.1
        beta_pixel=1000.0
        niter_pixel = -20
        niter_global  = 10
        pixel_dim=6
        simmetrizza=1
        do_refine_trajectory=0
        target_file="responses.h5"

        DOFIT(filename=filename, groupname=groupname, nref=nref, niter_optical=niter_optical, beta_optical=beta_optical ,
              beta_pixel=beta_pixel  , niter_pixel = niter_pixel,
              niter_global  = niter_global, pixel_dim=pixel_dim , simmetrizza=simmetrizza , do_refine_trajectory=do_refine_trajectory,
              target_file=target_file
        )

    if todo=="rois":

        filename = "../nonregressions/demo_imaging.hdf5"  # dove c' e lo scan originale
        groupname = "ROI_B/foil_scanXX/scans/Scan273/"    # Dove si trova esattament lo scan
        roisgroupname = "ROI_B/"                          # dove c' e' la roi originale


        target_filename = "newscan.h5"
        roisgroupname_target = "ROI_B_FIT8/"
        newscanstarget    = "scanXX/scans/Scan273"

        responsefilename =  "responses.h5"
        nex = 3

        DOROIS(filename = filename , groupname =  groupname,
               roisgroupname = roisgroupname ,
               target_filename = target_filename,
               roisgroupname_target =  roisgroupname_target ,
               newscanstarget    = newscanstarget,
               responsefilename =  "responses.h5",
               nex = nex)








