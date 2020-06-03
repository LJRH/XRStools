from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
import scipy.misc
from scipy import signal
import scipy

from scipy.ndimage import maximum_filter
import scipy.ndimage.morphology as morph
import scipy.ndimage.measurements as meas
## import skimage
## import skimage.transform


from numpy import gradient
from six.moves import range
from six.moves import zip

LENMIN=10

def gradX(x):
    res=np.zeros_like(x)
    res[:,  :-1]=x[:,:-1]-x[:,1:]
    return res
def gradXm(x):
    res=np.zeros_like(x)
    res[:,1: ]=-x[:,:-1]+x[:,1:]
    return res
def gradY(x):
    res=np.zeros_like(x)
    res[  :-1]=x[:-1]-x[1:]
    return res
def gradYm(x):
    res=np.zeros_like(x)
    res[1: ]=-x[:-1]+x[1:]
    return res
  
  
def lap(x):
  res = -(    gradX(gradXm(x))+gradXm(gradX(x))+    gradY(gradYm(x)) + gradYm(gradY(x))       )/2.0
  return res


def Diff(x ,coeffs):
  print( coeffs.shape)
  print( x.shape)

  DY,DX = gradient(x)

  yy , xy  = gradient(DY)

  dum , xx = gradient(DX)

  # xx = gradient(gradient(x,axis=1),axis=1)
  # yy = gradient(gradient(x,axis=0),axis=1)
  # xy = gradient(gradient(x,axis=0),axis=0)

  res = xx*coeffs[2] + yy*coeffs[0]+ +2*coeffs[1]*xy

  return res

def shrink(tmp, binning):
  n_cols =  (tmp.shape[0] // binning)
  n_rows =  (tmp.shape[1] // binning)
  tmp = tmp[:n_rows * binning,:n_cols * binning]
  tmp.shape = [ n_rows, binning, n_cols, binning  ]
  shrinked  =  tmp.max(axis=3).max(axis=1)
  return shrinked

def   deshrink(   large   ,  shrinked, binning    ):
  n_cols =  (large.shape[0] // binning)
  n_rows =  (large.shape[1] // binning) 

  tmp = large[:n_rows * binning,:n_cols * binning]  
  tmp.shape = [ n_rows, binning, n_cols, binning  ]
  shrinked=shrinked[:]
  shrinked.shape=[ n_rows,1 ,   n_cols, 1  ]
  tmp[:] =  shrinked
 

def ReadFile(name):
  data=scipy.misc.imread(name)*1.0
  data=np.max(data, axis=-1)
  return data


def CercaAnelli(A, lines=False):
  A=A.astype("d")
  # Nmarge=5
  # A[:Nmarge ,:]=1
  # A[:,:Nmarge ]=1
  # A[-Nmarge:,:]=1
  # A[:,-Nmarge:]=1


  Ax = np.concatenate( [A[1:], A[:1] ]  )- np.concatenate( [A[-1:], A[:-1] ]  )
  Ax = 2*Ax + np.concatenate( [Ax[:, 1:], Ax[:, :1] ], axis=1  ) + np.concatenate( [Ax[:,-1:], Ax[:,:-1] ] , axis=1  )


  Ay = np.concatenate( [A[:, 1:], A[:, :1] ] , axis=1  ) - np.concatenate( [A[:, -1:], A[:, :-1] ]  , axis=1   ) 
  Ay = 2*Ay + np.concatenate( [Ay[ 1:], Ay[ :1] ], axis=0  ) + np.concatenate( [Ay[-1:], Ay[:-1] ] , axis=0  )

  if lines:
      res=Canny_lines(Ax,Ay)
  else:
      res=Canny(Ax,Ay)
  # print len(res)
  
  return res


def IsMaximum( i,j,slope):
  totry = [[i-1,j-1],[i,j-1],[i+1,j-1],[i-1,j],[i+1,j],[i-1,j+1],[i,j+1],[i+1,j+1], ]
  for punto in totry:
   if( slope[punto[0], punto[1]  ]> slope[i,j]):
     return 0
  return 1

def Canny(Ax,Ay):
  Angles = np.arctan2(    Ax, -Ay )

  Angles=Angles*180/np.pi

  Angles = np.floor( Angles/45  + 0.5)
   
  slope = Ax*Ax + Ay*Ay
  
  Nmarge=6
  MaximumRatio=100.0
  

  # maxfits = maximum_filter(slope, size=3)
  # indici = np.where(slope==maxfits)
  # ListLocalMaxima=[ [i,j] for i,j in zip(indici[0], indici[1])   ]
  
  Ni=Ax.shape[0]
  Nj=Ax.shape[1]
  ListLocalMaxima=[]
  for i in range( Nmarge, Ni-Nmarge):
   for j in range( Nmarge, Nj-Nmarge):
     if( IsMaximum( i,j,slope) ):
        ListLocalMaxima.append( [i,j] )

  # print len(ListLocalMaxima)
  # raise

  steps={-5:[-1,1], -4:[-1,0],-3:[-1,-1], -2:[0,-1],-1:[1,-1],0:[1,0],1:[1,1],2:[0,1],3:[-1,1],4:[-1,0], 5:[-1,-1]}
  
  EdgeList=[]

  endpoints=[]



  EdgePoints=np.zeros(Ax.shape  )
  EdgePoints[:Nmarge ,:]=-1
  EdgePoints[:,:Nmarge ]=-1
  EdgePoints[-Nmarge:,:]=-1
  EdgePoints[:,-Nmarge:]=-1
  StartPoints = ListLocalMaxima
  
  for starting in StartPoints:
    EdgePointsTmp = np.zeros(Ax.shape  )  
    s0=starting[0]
    s1=starting[1]
    pentevalue= slope[s0, s1 ]
    if(EdgePoints[ s0, s1  ]==0):
      Edge=[]
      Edge.append( [s0,s1] ) 
      prendi=0
      traccia=0
      while(1):        
        if( slope[s0, s1 ]/pentevalue < 1.0/MaximumRatio):
          prendi=0
          traccia=1
          break  
        EdgePointsTmp[ s0, s1  ]=1
        direction = Angles[s0,s1]
        ss =  [ steps[direction-1],  steps[direction],  steps[direction+1]]
        pttrs = [  [s0+ss[0][0],s1+ss[0][1]], [s0+ss[1][0],s1+ss[1][1]], [s0+ss[2][0],s1+ss[2][1]]  ]
        values =  [ [slope[ pttrs[i][0], pttrs[i][1]  ],i]  for i in range(3)  ]
        imax=    max(values)[1]           
        s0=pttrs[imax][0]
        s1=pttrs[imax][1]                
        if(EdgePoints[ s0, s1  ]):
          prendi=0
          traccia=1
          break        
        if( EdgePointsTmp[ s0, s1  ] ):
          prendi=1
          traccia=1
          Edge.reverse()
          newedge = [ ]
          for p in Edge :
            if tuple(p)== (s0,s1):
              break
            newedge.append(p)
          Edge = [[s0,s1]]+newedge
          break
        else:
          Edge.append( [s0,s1] ) 
      if(prendi):
        endpoints.append([s0,s1])
        EdgeList.append(Edge)
      if traccia:
        for p in Edge:
          EdgePoints[p[0],p[1]]=1
  

  # values =  [ [len(edge),edge]  for edge  in  EdgeList]
  # edge =    max(values)[1]           

  values =  [ edge   for edge  in  EdgeList   if   len(edge)>10   ]
  # edge =    max(values)[1]           

  # return [edge]  
  return  values


def Canny_lines(Ax,Ay):
  Angles = np.arctan2(    Ax, -Ay )

  Angles=Angles*180/np.pi

  Angles = np.floor( Angles/45  + 0.5)
   
  slope = Ax*Ax + Ay*Ay
  
  Nmarge=1
  MaximumRatio=10000.0
  

  # maxfits = maximum_filter(slope, size=3)
  # indici = np.where(slope==maxfits)
  # ListLocalMaxima=[ [i,j] for i,j in zip(indici[0], indici[1])   ]
  
  Ni=Ax.shape[0]
  Nj=Ax.shape[1]
  ListLocalMaxima=[]
  for i in range( Nmarge, Ni-Nmarge):
   for j in range( Nmarge, Nj-Nmarge):
     if( IsMaximum( i,j,slope) ):
        ListLocalMaxima.append( [i,j] )

  # print len(ListLocalMaxima)
  # raise

  steps={-5:np.array([-1,1]), -4:np.array([-1,0]),-3:np.array([-1,-1]), -2:np.array([0,-1]),-1:np.array([1,-1]),
          0:np.array([1,0]),1:np.array([1,1]),2:np.array([0,1]),3:np.array([-1,1]),4:np.array([-1,0]), 5:np.array([-1,-1])}
  
  EdgeList=[]

  endpoints=[]

  EdgePoints=np.zeros(Ax.shape  )
  EdgePoints[:Nmarge ,:]=-1
  EdgePoints[:,:Nmarge ]=-1
  EdgePoints[-Nmarge:,:]=-1
  EdgePoints[:,-Nmarge:]=-1
  StartPoints = ListLocalMaxima
  

  stack = []
  for starting in StartPoints:
    s0=starting[0]
    s1=starting[1]
    pentevalue= slope[s0, s1 ]
    stack.append( [  (s0,  s1), pentevalue, 1    ] ) 
  

  while len(stack):
    starting, pentevalue, direction_fact = stack[-1]
    stack=stack[:-1]
    
    print( " inizio da ", starting)
    EdgePointsTmp = np.zeros(Ax.shape  )  
    s0=starting[0]
    s1=starting[1]
    pentevalue= slope[s0, s1 ]
    if(EdgePoints[ s0, s1  ]==0):
      Edge=[]
      Edge.append( [s0,s1] ) 
      prendi=0
      traccia=0
      while(1):   
        print( s0,s1)
        if( slope[s0, s1 ]/pentevalue < 1.0/MaximumRatio):

          if direction_fact==1  :
            prendi=0
            traccia=0
            if  len(Edge)>4:
              stack.append( [  Edge[-4], pentevalue, -1    ]) 
            print( " troppo debole Inverto")
          else:
            prendi= len(Edge)>5
            traccia=1
            print( " troppo debole finisco")
          break  
        EdgePointsTmp[ s0, s1  ]=1
        direction = Angles[s0,s1]
        ss =  np.array([ steps[direction-1],  steps[direction],  steps[direction+1]])*direction_fact
        pttrs = [  [s0+ss[0][0],s1+ss[0][1]], [s0+ss[1][0],s1+ss[1][1]], [s0+ss[2][0],s1+ss[2][1]]  ]
        values =  [ [slope[ pttrs[i][0], pttrs[i][1]  ],i]  for i in range(3)  ]
        imax=    max(values)[1]           
        s0=pttrs[imax][0]
        s1=pttrs[imax][1]                
        if(EdgePoints[ s0, s1  ]):
          print( " scontro vecchio in ", s0, s1)
          if direction_fact==-1:
            Edge=Edge[:-10]
            prendi=1
            traccia=1
            break
          else:
            EdgePointsTmp[:]=0
            prendi=0
            traccia=0
            if  len(Edge)>10:
              stack.append( [  Edge[-4], pentevalue, -1    ]) 
            print( " troppo debole Inverto")
            break


            
        if( EdgePointsTmp[ s0, s1  ] ):
          print( " scontro nuovo ")
          prendi=1
          traccia=1
          Edge.reverse()
          newedge = [ ]
          for p in Edge :
            if tuple(p)== (s0,s1):
              break
            newedge.append(p)
          Edge = [[s0,s1]]+newedge
          break
        else:
          print( " continup ")
          Edge.append( [s0,s1] ) 
      if(prendi):
        endpoints.append([s0,s1])
        EdgeList.append(Edge)
      if traccia:
        for p in Edge:
          EdgePoints[p[0],p[1]]=1
  

  # values =  [ [len(edge),edge]  for edge  in  EdgeList]
  # edge =    max(values)[1]           

  values =  [ edge   for edge  in  EdgeList   if   len(edge)>10   ]
  # edge =    max(values)[1]           

  # return [edge]  
  return  values




def divergenza( Ax,Ay, edge):
  if( edge[-1] != edge[0]):
      edge.append(edge[0])
  res=0
  for i in range(len(edge)-1):
     p1=edge[i]
     p2=edge[i+1]
     Ax1=Ax[p1[0],p1[1] ]
     Ax2=Ax[p2[0],p2[1] ]
     Ay1=Ay[p1[0],p1[1] ]
     Ay2=Ay[p2[0],p2[1] ]
     Dx = p2[0]-p1[0]
     Dy = p2[1]-p1[1]
    
     res=res+ ( Dx*(Ay1+Ay2) -Dy*(Ax1+Ax2) )*0.5
  return res 



def threshold_mask(mask, image,  value  ) :
  n=mask.max()
  npix=3

  image = signal.medfilt2d(image, kernel_size=3)

  for i in range(1,n+1):
  
    icount=0
    tmp_mask_old = 0
    while(icount<100):
      tmp_mask = np.equal( i,  mask  )
      
      if not np.any(tmp_mask):
        continue

      mask[:]=mask*(1-tmp_mask)
      
      massimo = (image*tmp_mask).max()
      print( " MASSIMO " , massimo)
      tmp_mask_2 =   morph.grey_dilation(tmp_mask,  footprint=np.ones([npix,npix]),
                                         structure=np.zeros([npix,npix]))

      print( " Npunti , value " ,   tmp_mask.sum(), tmp_mask_2.sum() , value)

      tmp_mask = np.less( value*massimo,  image  )*tmp_mask_2

      print( " Npunti  " ,   tmp_mask.sum())
     
      
      mask[:]+=tmp_mask*i
     
      if ( tmp_mask^tmp_mask_old).sum()==0:
        break
      tmp_mask_old = tmp_mask


      print( " ================ ")
      print( mask.sum())

      icount+=1

  return mask


def grow_mask(input, npix  ) :
  output =  morph.grey_dilation(input,  footprint=np.ones([npix,npix]),
                                structure=np.zeros([npix,npix]))
  return output
def shrink_mask(input, npix  ) :
  output =  morph.grey_erosion(input,  footprint=np.ones([npix,npix]),
                               structure=np.zeros([npix,npix]))
  return output


def get_spots_mask( A, rrA, median_size=5, nofroi=12,  give_borders=False, tval=-1, Hough=False ) :
  if not Hough:
    res=get_spots_mask_Normal( A, rrA, median_size=median_size, nofroi=nofroi,  give_borders=False, tval=-1 )
    return res
  else:
    res=get_spots_mask_Lines( A,  rrA, median_size=median_size, nofroi=nofroi,  give_borders=False, tval=-1 )
    return res
    
def get_spots_mask_Normal( A, rrA, median_size=None, nofroi=12,  give_borders=False, tval=-1, Hough=False ) :
    

      print( " qui Hough " , Hough)
      A = A*(1-rrA)
      A = signal.medfilt2d(A, kernel_size=median_size)
      mmax = A.max()
      A[A<tval*mmax]=0
  
      A=scipy.ndimage.filters.gaussian_filter(A, 3.0) # order=0, output=None, mode='reflect', cval=0.0, truncate=4.0)[source]
      # if Hough:
      #   print( " CHIAMO " )
      #   hspace, angles, dists = skimage.transform.hough_line(A)
      #   print( angles.shape)
      #   print( dists.shape)
      #   hspace, angles, dists = skimage.transform.hough_line_peaks(A, angles, dists, min_distance=4, min_angle=10,)
      #   print( angles , dists)


      if True or not Hough:      
        mask=np.zeros(A.shape,"i")
        cerchi = CercaAnelli(A)
        for c in cerchi[:]:
          if len(c)<LENMIN: continue
          for y,x in c:
            mask[y,x] = 1

        filled=morph.binary_fill_holes(mask)

        newmask = relabelise(filled,A, nofroi)

        if give_borders:
          return newmask, mask
        else:
          return newmask
      else:
        pass


def relabelise(filled, A , nofroi):
      labels, nlabs = meas.label(filled, structure=np.ones([3,3]))
      print( " nlabs ", nlabs)
      aves = meas.mean(A ,labels=labels, index = list(range(1,nlabs+1)) )
      # if type(aves)!=type([]):
      #   aves=[aves]
      avlabs = list(zip(aves, np.arange(1,nlabs+1)))  
      avlabs.sort()
      intes = [ a for a,l in avlabs[-nofroi:]]
      amed = np.median(intes)
      labs = [ l for a,l in avlabs[-nofroi:] if a>amed/100]
      newmask = np.zeros(filled.shape,"i")
      i=1
      for l in labs:
        newmask[np.equal(labels,l)]=i
        i+=1
      return newmask


def get_spots_mask_Lines( A,rrA,  median_size=None, nofroi=12,  give_borders=False, tval=-1 ):
  A=A*(1-rrA)  
  thematrix = np.zeros(    [A.shape[0]+4, A.shape[1]+4 ]  , A.dtype)
  thematrix[2:-2   ,2:-2 ] = A
  thematrix=scipy.ndimage.filters.gaussian_filter(thematrix, 3.0) 

  Nrows , Ncols =  thematrix.shape
  
  submatrix = np.array(thematrix)
  mask, diff_coeffs = get_spots_mask_Lines_slave( submatrix  , nofroi=nofroi, give_diff_coeffs=True , tval=tval) 



#  ######################
#  nspots = mask.max()
#  if nspots!=12:
#      print "WARNING: LESS SPOTS WERE FOUND :  " ,  nspots
#  else:
#      print "GOOD!:  FOUND :  " ,  nspots, "  SPOTS " 

#  for l in range(1,nspots+1) :
#      maskzone = np.equal(mask,l)
#      submatrix[maskzone] = submatrix.max()
        
#  import pylab
#  f1=pylab.figure()
#  from matplotlib.colors import LogNorm
#  pylab.imshow(submatrix, norm=LogNorm(vmin=0.01))
#  pylab.show()
#  #################################


  
  ##############################################################
  submatrix = np.array(thematrix)
  for iter in range( (median_size*2)**2):
    submatrix  =  submatrix + Diff(submatrix,diff_coeffs )*0.5*0.2
    newsub=np.zeros_like(submatrix)
    newsub[2:-2   , 2:-2 ] = submatrix [2:-2   , 2:-2 ]
    submatrix=newsub

          
  mask = get_spots_mask_Lines_slave( submatrix    ,  nofroi=nofroi, give_diff_coeffs=False , tval=tval ) 



#  ######################
#  nspots = mask.max()
#  if nspots!=12:
#      print "WARNING: LESS SPOTS WERE FOUND :  " ,  nspots
#  else:
#      print "GOOD!:  FOUND :  " ,  nspots, "  SPOTS " 

#  for l in range(1,nspots+1) :
#      maskzone = np.equal(mask,l)
#      submatrix[maskzone] = submatrix.max()
        
#  import pylab
#  f1=pylab.figure()
#  from matplotlib.colors import LogNorm
#  pylab.imshow(submatrix, norm=LogNorm(vmin=0.01))
#  pylab.show()
#  #################################



  
  nspots = mask.max()
  if nspots!= nofroi:
    print( "WARNING: LESS SPOTS WERE FOUND :  " ,  nspots)
  else:
    print( "GOOD!:  FOUND :  " ,  nspots, "  SPOTS " )

  return mask[2:-2,2:-2]

def get_spots_mask_Lines_slave( A, median_size=None, nofroi=120,  give_borders=False, give_diff_coeffs=False, tval=-1) :
      mmax = A.max()
      A[A<tval*mmax]=0
  
      mask=np.zeros(A.shape,"i")
      cerchi = CercaAnelli(A, lines=True)
      maxlen = 0
      for c in cerchi[:]:
        maxlen = max(maxlen,len(c))

      momenti = []
      for c in cerchi[:]:
        c=np.array(c)
        
        m = np.array([0.0,0.0])
        m2 = np.zeros([2,2],"d")
        
        if len(c)<LENMIN  or len(c)<maxlen/5: continue
        for cel in c:
          mask[cel[0],cel[1]] = 1

          m[:]+= cel
          m2[:]+= cel[:,None]*cel
          
        N = len(c)
        
        m=m/N
        m2 = m2/N -m[:,None]*m
        m2=m2/(m2[0,0]+m2[1,1])
        print( " m2 " , m2)
      
        momenti.append( [m,   np.array([ m2[0,0], m2[0,1], m2[1,1]])  ])

      if give_diff_coeffs:
        coefficients_fissi = np.zeros([3]+list(A.shape)      ,"f")
        coefficients_large = np.zeros([3]+list(A.shape)      ,"f")
        for icoeff in range(3):
          print( " icoeff " , icoeff)
          for (m,m2), c in zip(momenti, cerchi):
            c=np.array(c)
            coefficients_fissi[icoeff][c[:,0], c[:,1]] = m2[icoeff]

        for icoeff in [2,1,0]:
          fisso =  coefficients_fissi[icoeff]
          mask_large =  np.less(0,np.abs(fisso ) )
          coefficients_large[icoeff][ mask_large  ] = fisso[ mask_large ] 

          for binning,niters in zip([12,8,4,2,1] , [2000,500, 400,50,50]  ):

            fissi_shrinked = shrink( coefficients_fissi[icoeff] ,binning   )
            coeff_shrinked = shrink( coefficients_large[icoeff] ,binning   )
            mask_shrink = np.less(0,np.abs(fissi_shrinked) )
            coeff_shrinked[ mask_shrink  ] =   fissi_shrinked[ mask_shrink  ]

            for iter in range(niters):
              if iter%50==0:
                print( iter)
              coeff_shrinked += 0.2* lap(  coeff_shrinked )
              coeff_shrinked[ mask_shrink  ] =   fissi_shrinked[ mask_shrink  ]

            deshrink(   coefficients_large[icoeff]   ,  coeff_shrinked, binning    )

        # import pylab
        # pylab.imshow(coefficients_large[icoeff]); pylab.show()
        # pylab.plot(  coefficients_large[icoeff][:,50])
        # pylab.plot(  coefficients_large[icoeff][:,100])
        # pylab.plot(  coefficients_large[icoeff][:,150])
        # pylab.plot(  coefficients_large[icoeff][:,200])
        # pylab.show()
          
      ##############################################################
      # for iter in range(500):
      #   A  =  A + Diff(A,coefficients)*0.5*0.2

      newmask = mask
      print( " NEWMASK ", newmask)
      
      if not give_diff_coeffs:

#          maskzone = np.less(0,mask)
#          A[maskzone] = A.max()
        
#          import pylab
#          f1=pylab.figure()
#          from matplotlib.colors import LogNorm
#          pylab.imshow(A, norm=LogNorm(vmin=0.01))
#          pylab.show()
#          raise



          
          filled=morph.binary_fill_holes(mask)
          newmask = filled
          newmask = relabelise(filled,A, nofroi)
          pass
      else:
          newmask = mask

      if give_diff_coeffs:
          return newmask, coefficients_large
      else:
          return newmask






from PyMca5 import EdfFile
# from PyMca5 import EdfFile


if __name__=="__main__":

  thematrix=EdfFile.EdfFile("test.edf").GetData(0)

  # A=A[:256,:256]

  Nrows , Ncols =  thematrix.shape

  # thematrix=thematrix[Nrows-256:, Ncols-256:  ]
  # Nrows , Ncols =  thematrix.shape

  for i  in range(0,Nrows,256):
    for j  in range(0,Ncols,256):

      submatrix = thematrix[i:i+256,j:j+256]

      A = signal.medfilt2d(submatrix, kernel_size=5)
      submatrix[:]=A

      mask,borders  = get_spots_mask( submatrix, submatrix*0 ,   median_size=5, give_borders=True) 

      nspots = mask.max()

      if nspots!=12:
        print( "WARNING: LESS SPOTS WERE FOUND :  " ,  nspots)
      else:
        print( "GOOD!:  FOUND :  " ,  nspots, "  SPOTS " )



      for l in range(1,nspots+1) :
        maskzone = np.equal(mask,l)*np.equal(borders,1)
        submatrix[maskzone] = submatrix.max()


  EdfFile.EdfFile("res.edf","w+").WriteImage({},thematrix)  









