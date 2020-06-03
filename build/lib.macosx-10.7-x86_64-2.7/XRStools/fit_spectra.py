from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import scipy
import math
import sys
from six.moves import range
from six.moves import zip
import pickle
if(sys.argv[0][-12:]!="sphinx-build"):
    from XRStools import fitspectra_cy


def Fista(solution   ,  problem,   niter, niterLip):
    # print " SHAPE " , solution.shape
    dim = solution.shape[0]
    err = 0.0
    beta=problem.beta
    grad , err=problem. calculate_grad(    solution ) 
    Lip = math.sqrt( np.linalg.norm(grad)   )
    grad   = grad/ Lip

    print( "CALCULATING LIPSCHITZ FACTOR ")
    
    for i in range(niterLip):
        grad2,err2 = problem.calculate_grad(grad,  quadratic_only=1)
        Lip = math.sqrt( np.linalg.norm(grad2)   )
        grad   = grad2/ Lip
        print( "LIP ", Lip)

    Lip = Lip*1.05
    
    t=1.0
    y      = solution
    x_old  = solution

    errList=[]
    for i in range(1000):
        grad, err = problem.calculate_grad(y)
        print( "err " , err)
        errList.append(err)
        y = y-grad/Lip
        y=np.maximum(y-beta/Lip,0)
    return y,errList

    errList=[]
    for iter in range(abs(niter)):
        grad, err = problem.calculate_grad(y)
        errList.append(err)
        solution =  y + grad/Lip
        solution = np.maximum(solution-beta/Lip, 0)
        tnew = ( 1+math.sqrt(1.0+4*t*t) )/2
        y[:] = solution +(t-1)/tnew *( solution - x_old )
        t = tnew
        if niter<0:
            t=1
        x_old[:] = solution
        if iter%1 ==0:
            # sys.stdout.write(" "*60+"\r"+("FISTA iter %d  errore est %e  mod_grad est  %e" % ( iter,  err, grad.std()) ))
            sys.stdout.write(("FISTA iter %d  errore est %e  mod_grad est  %e\n" % ( iter,  err, grad.std()) ))
            sys.stdout.flush()
        print( " ")
    return solution, errList


def     calculate_grad(problem, spettro,  quadratic_only=0):

    err = problem.DD_scal/2
    grad = -problem.i2f.dot( problem.SD_scal_flat    )
    err = err +np.dot( spettro, grad   )

    tmp = problem.f2i.dot(spettro)
    tmp = np.reshape(tmp,[-1, problem.SS_scal.shape[0]]    )
    tmp = np.tensordot(  tmp, problem.SS_scal, axes=[[-1],[0]]     )
    tmp = np.reshape(tmp,[-1]    )
    add = problem.i2f.dot(tmp)
    err = err + np.dot(add, spettro    )/2
    if not quadratic_only:
        grad[:]+=add
    else:
        grad = add

    err=err + problem.beta*( np.abs(spettro) ).sum()
    return grad, err
    

def  fitta( DD_scal,  SS_scal,  SD_scal,    energie_spettro,  energie_for_SD     , spettro   , beta,   niter, niterLip  ):
    
    de = energie_spettro[1]-energie_spettro[0]
    e_min = energie_spettro[0]
    e_max = energie_spettro[-1]
    
    energie_for_SD_flat = np.array(energie_for_SD.flat)
    SD_scal_flat        = np.array(SD_scal.flat)
    
    sparse_elements = []
    for j in range(len(energie_for_SD_flat)):
        eneSD =   energie_for_SD_flat[j]
        assert(eneSD>e_min )
        assert(eneSD<e_max )

        i = int((eneSD-e_min)/de)
        f = (eneSD-e_min)/de -i
        
        sparse_elements.append([ j,i  , (1-f)   ])
        sparse_elements.append([ j,i+1,    f    ])


    J,I,F = np.array(sparse_elements).T
        
    f2i_coo =    scipy.sparse.coo_matrix( (F,(J,I)) , shape = [   energie_for_SD_flat.size   ,   spettro.size  ])   
    i2f_coo =    scipy.sparse.coo_matrix( (F,(I,J)) , shape = [   spettro.size  , energie_for_SD_flat.size     ])
    
    f2i = f2i_coo.tocsr() 
    i2f = i2f_coo.tocsr() 

    problem  = type('MyObject', (object,), {"calculate_grad":calculate_grad  ,  "f2i":f2i,"i2f":i2f, "SS_scal":SS_scal, "SD_scal":SD_scal,  "SD_scal_flat":SD_scal_flat, "DD_scal":DD_scal, "beta":beta })()


    solution, errList = Fista(spettro   ,  problem,   niter, niterLip )

    ii = f2i.dot( solution  )
  
    return ii, solution, errList
    

def fit_spectra_main( references  , sample_s , DE , beta,   niter, niterLip , slopeInfos  ,  discard_threshold = 0 , threshold_fraction = 0    ):
    
    chiavi = list(references.keys())
    res = {}
    for k in chiavi :
        ref = references[k]
        SInfo  = slopeInfos[k]
        data = []
        for scan in sample_s.keys():
            data.append(    sample_s[scan][k] )
          
        (sintesi, energie_spettro, spettro_byline, errors,
         solution,   errList      )=  fit_spectra_roi( ref,  data, DE  , beta,   niter, niterLip  , SInfo ,  discard_threshold = discard_threshold , threshold_fraction = threshold_fraction )


        res[k] = {"energies" : energie_spettro,
                  "spectraByLine" : spettro_byline,
                  "errors" : errors, 
                  "spectraByFit": solution,
                  "fit_errList": errList,
                  "sintesi":sintesi 
                  ##"spectra_byscalprod" : spectra_byscalprod
                  
        }
    return res



def do_spettro_byscal(energie_spettro, spettro_byscal, SD_scal, energie_for_SD,  S_L1  ):
    
    byscal_sum = np.zeros_like(energie_spettro)


    Slong_L1 = np.zeros_like(SD_scal)
    Slong_L1[:,:] = Slong_L1[:,:]+S_L1[:,None]


    El = energie_spettro[0]
    Eh = energie_spettro[-1] + energie_spettro[1]-energie_spettro[0]

    NS = len(energie_spettro)
    
    SD_scal = np.array(SD_scal.flat)
    energie_for_SD = np.array(energie_for_SD.flat)
    Slong_L1 = np.array(Slong_L1.flat)
    
    for E,Scal,Sl1 in zip(  energie_for_SD  , SD_scal,  Slong_L1    ) :
        fpos = NS*((  E-El)/(Eh-El))
        if fpos>0 and fpos<NS-1:
            ipos = int(fpos)
            ipos1= ipos+1
            f = fpos-ipos
            
            byscal_sum[ipos] +=   (1-f)*Sl1
            byscal_sum[ipos1] +=   (f)*Sl1
            
            spettro_byscal[ipos] +=   (1-f)*Scal
            spettro_byscal[ipos1] +=   (f)*Scal
    
    for i in range(len(spettro_byscal)) :
        if byscal_sum[i]>0:
            spettro_byscal[i] = spettro_byscal[i]/byscal_sum[i]

def fit_spectra_roi( ref,  datas, user_de  , beta,   niter, niterLip  , SInfo ,  discard_threshold = 0 , threshold_fraction = 0  ):
    ##  reference scan energy ( the analyser energies, in an array_
    if ref is not None:
        enes_ref = ref.zscale
        DE_ref = (enes_ref[-1]-enes_ref[0])
        ME_ref = (enes_ref[-1]+enes_ref[0])/2

        if   hasattr( ref,  "incidentE"  ):
            incidentE = ref.incidentE
            if incidentE is not None:
                SHIFT = incidentE*1000 - ME_ref
                print(  " SHIFT " , SHIFT)
        else:
            SHIFT = 0

        deltaEref = (enes_ref[-1]-enes_ref[0])/(len(enes_ref)-1)
        enes_ref = ME_ref+SHIFT - enes_ref  # da aggiunger  # (ME_ref+SHIFT) is nothing more than incident energy for reference scan

        MASK = ref.mask
        
        CRY = (MASK.shape[0]-1)/2.0  ## The intercepts refer the positions of the center of the ROI 
        CRX = (MASK.shape[1]-1)/2.0  ## Here the center coordinates

        fNmiddle =  ( len(enes_ref) - 1 )/2.0 + SHIFT/deltaEref  ##   this is the reference step's  number  for which the analyser energy = incidentE 

        
    else:
        MASK = SInfo["mask"]
        CRY = (MASK.shape[0]-1)/2.0  ## The intercepts refer the positions of the center of the ROI
        CRX = (MASK.shape[1]-1)/2.0  ## The intercepts refer the positions of the center of the ROI
        # fNmiddle = ( len(enes_ref)-1.0 )/2.0
        fNmiddle = None
        SHIFT=0
        
    count=0
    for data in datas:
        
        count+=1
        enes_data = data.zscale     
        # print " ENES DATA ", enes_data
        
        mymine = enes_data.min()
        mymaxe = enes_data.max()
        myde   = (mymaxe-mymine)/(len(enes_data)-1)
        
        if count==1:
            mine = mymine
            maxe = mymaxe
            de   = myde
        else:
            mine = min(mine, mymine  )
            maxe = max(maxe, mymaxe  )
            de   = min(de, myde )
            
    if user_de !=0:
        de = max(user_de,de)
    
    if ref is not None:
        
        mine = mine+enes_ref.min()-de
        maxe = maxe+enes_ref.max()+de
        
    else:
        
        mine = mine+  (-CRY/ abs( SInfo["zrate"]))*SInfo["estep"]          - de
        maxe = maxe+  ( CRY/ abs( SInfo["zrate"]))*SInfo["estep"]          + de
    
    nsteps = int(   round((maxe-mine)/de) )
    maxe = mine + (nsteps)*de
    # print maxe,mine,de, (maxe-mine)/de

    energie_spettro = np.linspace(mine -mine, maxe-mine, num=nsteps+1  ,  endpoint=True)
    spettro = np.zeros(nsteps+1,"d")

    ###############################################################
    
    spettro_byline = np.zeros(nsteps+1,"d") ## where intensity will be accumulated
    error_sum = np.zeros(nsteps+1,"d")      ## where the N of occurencies will be accumulated
                                            ##  needed to do the average 
    frequencies_sum = np.zeros(nsteps+1,"d") ## to keep statisthics



    if ref is not None:

        hline =ref.line_infos.line [0]             ## the height of the line
        slopeline = ref.line_infos.line [1]        ##  the slope of the line
        response_line_intensity = ref.line_infos.optical_response.sum(axis=0) # This is used to do some kind of weigthing
                                                                              # The reponse might be stronger at the middle of the line  and weaker
                                                                              #  at the left and right extrema
                                                                              # The response function is a 2D image , with this instruction we
                                                                              # project it over the line abscissa

        ## Yintercept and Xintercept  are always subtracted by CRY and CRX respectively because in reponsepercussionelle.py
        ## they are the coordinate of the moving center of the response function
                                                                              
        hline += (ref.line_infos.Yintercept-CRY) + fNmiddle*ref.line_infos.Yslope - ((ref.line_infos.Xintercept-CRX) + fNmiddle*ref.line_infos.Xslope ) *slopeline
                                                                                                                              ##  This isthe line height 
                                                                                                                              ## in the middle of the reference scan
                                                                                                                              ## We count energy from this middle
        DHoverDI = ref.line_infos.Yslope - ref.line_infos.Xslope *slopeline           ## at the begginiing of the scan.
                                                                                      ## The line start from a shiftx = ref.line_infos.Xintercept and this
                                                                                      ## cross-talk with the slopeline
                                                                                      ## DHoverDI is the height variation at each energy step, it keeps into account also
                                                                                      ## the slope
    else:
        DHoverDI =  SInfo["zrate"]  ## when given in input zrate must correspond to:   Yslope - Xslope *slopeline
        slopeline = SInfo["slope"]
        deltaEref = SInfo["estep"]
        hline = CRY #### - CRX*slopeline
        
    for data in datas:
        ## here data is a stack of 2D images for a given ROI
        
        enes_data = data.zscale   ## this is supposed to be the analyser_energy

        denominator = data.denominator  ## This is the factor applied for renormalisation
                                        ## it will be used to properly account for statisthical error

        mms        = data.mm

        # P_enes=[]
        # P_mms=[]
        # P_sp_e=[]
        # P_sp_s=[]
        # P_MASK=MASK



        usecython=1
        if usecython:


            
            if ref is not None:
                useref=1
                weight_by_response = ref.line_infos.weight_by_response
                Xintercept = ref.line_infos.Xintercept
                # CRX
                # fNmiddle
                Xslope =ref.line_infos.Xslope
                response_line_intensity=response_line_intensity
            else:
                useref=0
                weight_by_response =0
                Xintercept =0.0
                CRX =0.0
                fNmiddle=0.0
                Xslope =0.0
                response_line_intensity=np.array([0.0,0.0],"d")

            print( " CHIAMO FITSP")
            print (" enes_data  " , enes_data.dtype)
            print (" denominator  " , denominator.dtype)
            print ("  mms " , mms.dtype)
            print ("  MASK " , MASK.dtype)

            #  spettro_byline = spettro_byline.astype("d")  ne pas faire ca enligne de l'argument
            #  car spettro_byline est valeur de retour, ainsi que
            
            frequencies_sum = frequencies_sum.astype("d")
            error_sum       = error_sum.astype("d")


            
            fitspectra_cy.spectra_roi_by_line(spettro_byline,
                                              error_sum,
                                              frequencies_sum,
                                              enes_data.astype("f"),
                                              denominator.astype("f") ,
                                              mms.astype("f"),
                                              MASK.astype("f") ,
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
                                              response_line_intensity.astype("d")
                                          )
            print( "  in uscita da  spectra_roi_by_line "  ,  spettro_byline.sum() ) 
        else:


            mask_npix = MASK.sum()

            for iE , (ene , mm, deno) in enumerate(zip(enes_data,mms, denominator)):
        
                
                if discard_threshold :
                    mm_npix=( np.less(   discard_threshold ,mm )  ).sum()
                    print (" DISCARD ", discard_threshold, threshold_fraction , mask_npix, mm_npix)
                    if mm_npix > mask_npix*threshold_fraction:
                        continue



                ## one energy, one 2D image from the stack, one value for the denominator

                # P_enes.append(ene)
                # P_mms.append(mm)
                # P_sp_e.append(   ene -(np.arange(mm.shape[0]) -hline)/DHoverDI  * deltaEref   )
                # P_sp_s.append(  mm[:,int(CRX)]   )

                for iy in range(mm.shape[0]):
                    for ix in range(mm.shape[1]):
                        # slow loop in python. You know why it is orribly slow.

                        if not MASK[iy,ix]:
                            ## Good. Mask is already taken into account
                            ## But needs to be passed . Issue to be  followed
                            continue

                        ## Assigning an energy to the pixel
                        H0 = iy- ( hline + ix * slopeline)
                        ii = H0 /DHoverDI
                        E = ene - ii *deltaEref    #   + DE_ref*0

                        ## calculating the longitudinal position along the line ( to weight with the projected response_line_intensity)
                        if ref is not None   and  ref.line_infos.weight_by_response:
                            posx_inref =   int( round(   (ix-( (ref.line_infos.Xintercept-CRX+ fNmiddle*ref.line_infos.Xslope ) + ii * ref.line_infos.Xslope  )   )))
                            if posx_inref>=0 and posx_inref< len(response_line_intensity):
                                freq = response_line_intensity[posx_inref] # a weigthing factor
                            else:
                                freq = None
                        else:
                            freq = 1.0

                        if   freq is not  None:                        
                            fipos = (E-mine)/de                        # the position in pixel units of the contribution to the spectra array (which starts from mine)
                            ipos  = int(fipos)                         # the integer part of fipos
                            f = fipos - ipos                           # the fractional residu of fipos

                            if ipos>0 and ipos < nsteps:              # If I am withing the range of the spectra
                                frequencies_sum[ipos] += (1-f)* freq       # I distribute the contribution : 100% if f=0, to ipos , with the weigth given by response
                                spettro_byline [ipos] += (1-f)* mm[iy,ix]  # Same thing : distributing intensity to spectro_..

                                frequencies_sum[ipos+1] +=  f*freq         # Same thing as above, just 100% if f=1 because we are distributing to the upper pixel
                                spettro_byline [ipos+1] +=  f* mm[iy,ix]

                                ## Calculating the error by hoping that the final result  be  gaussian
                                error_sum[ipos] += (1-f)*(1-f)* mm[iy,ix] /deno   # 
                                error_sum[ipos+1] += f*f* mm[iy,ix] /deno   # 

            # ff = open("/tmp/sp.p","wb")
            # todump =  [P_enes, P_mms, P_sp_e, P_sp_s, P_MASK]
            # pickle.dump(   todump  , ff)

        
    for i in range(len(spettro_byline)):
        if spettro_byline [i] >0:
            spettro_byline [i] = spettro_byline [i]  /frequencies_sum[i]
            error_sum[ i ]     = math.sqrt(error_sum[ i ])/frequencies_sum[i]
    ####################################################################
    if niter:

        len_data = 0

        for data in datas:
            len_data +=   len(data.zscale)

        len_ref = len(enes_ref)
        SS_scal = np.zeros([len_ref, len_ref],"d")

        S_L1    = np.zeros([len_ref],"d")

        SD_scal = np.zeros([len_ref, len_data],"d")

        DD_scal = 0.0
        for data in datas:
            DD_scal += (  data.mm* data.mm  ).sum()

        mm = ref.mm
        for i in range(len_ref):
            S_L1[i] = np.abs(mm[i]).sum()

            for j in range(len_ref):
                SS_scal[i,j] = ( mm[i] * mm[j]  ).sum()

        count=0
        for data in datas:
            for image in data.mm:
                for i in range(len_ref):
                    SD_scal[ i , count] = ( mm[i] * image  ).sum()
                count+=1
        energies_sample =  np.zeros( len_data , "d")
        len_data = 0
        for data in datas:
            energies_sample[len_data:len_data+len(data.zscale)]=data.zscale
            len_data +=   len(data.zscale)

        energie_for_SD   =  enes_ref[:,None]         +  (energies_sample-mine)

        # spettro_byscal = np.zeros_like(energie_spettro)
        # byscal_sum = np.zeros_like(energie_spettro)
        #do_spettro_byscal(energie_spettro, spettro_byscal, SD_scal, energie_for_SD,  S_L1  )
        ##solution, errList = fitta( DD_scal,  SS_scal,  SD_scal,    energie_spettro,  energie_for_SD     , spettro    , beta,   niter, niterLip   )


        ii, solution, errList = fitta( DD_scal,  SS_scal,  SD_scal,    energie_spettro,  energie_for_SD     , spettro    , beta,   niter, niterLip   )

        ii.shape = energie_for_SD.shape
        # print mm.shape
        # print ii.T.shape
        sintesi  =  np.tensordot(ii.T , mm , axes = [ (-1),(0)   ] )
        len_data = 0

    else:
        sintesi , solution, errList = None, None, None
        
    # return sintesi, energie_spettro, spettro_byline, solution,   spettro_byscal, errList
    return sintesi, energie_spettro+mine, spettro_byline, error_sum, solution,    errList








