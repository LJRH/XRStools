from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy  as np
from six.moves import range
from six.moves import zip


def match(distorted,refgrid, coors,  targets, possibili):
    for target in set(tuple(coors))-set(tuple(targets)):
        new_targets=targets+[tuple(target)]
        if validate( distorted,refgrid, coors,  new_targets):
            if len(new_targets)==len(distorted):
                print( new_targets)
                possibili.append(   new_targets  )
            else:
                match(distorted,refgrid, coors,  new_targets, possibili)
        

 
def validate(distorted,refgrid, coors,  new_targets)  :
    C=new_targets[-1]
    targets= new_targets[:-1]
    nt=len(targets)

    tC=distorted[nt]

    for i in range(nt):
        A=targets[i]
        tA=distorted[i]
        for j in range(i+1,nt):
            B=targets[j]
            tB=distorted[j]
            v1 = [ C[0]-A[0] , C[1]-A[1]   ] 
            v2 = [ C[0]-B[0] , C[1]-B[1]   ] 
            vv=abs(np.array( [ C[0]-A[0] , C[1]-A[1], C[0]-B[0] , C[1]-B[1]]))
            Area = v1[0]*v2[1]-v1[1]*v2[0]
            dmax= vv.max()
            if abs(Area)>= dmax*dmax/2 :
                v1 = [ tC[0]-tA[0] , tC[1]-tA[1]   ] 
                v2 = [ tC[0]-tB[0] , tC[1]-tB[1]   ] 
                tArea = v1[0]*v2[1]-v1[1]*v2[0]           
                if Area*tArea<0:
                    return False
    return True

def merit( distorted,refgrid, coors, choice):
    result=0.0
    dic={}
    minx=1.0e20
    maxx=-1.0e20
    for (k,l),b in zip(choice,distorted):        
        dic[(int(k), int(l))]=b
        minx=min(minx, b[0])
        maxx=max(maxx, b[0])

    D = maxx-minx
    print( dic)

    NX=np.max(np.array(coors)[:,0])+1
    NY=np.max(np.array(coors)[:,1])+1
    for i in range(NX-1):
        for j in range(NY-1):
            try:
                a1=dic[(i,j)]
                presenti=1
            except:
                presenti=0
            if presenti:
                result=result - D/NX * i*a1[0]*2
                result=result - D/NX * j*a1[1]*2


    for i in range(NX-1):
        for j in range(NY-1):
            try:
                a1=dic[(i,j)]
                a2=dic[(i+1,j+1)]
                b1=dic[(i+1,j)]
                b2=dic[(i,j+1)]                
                presenti=1
            except:
                presenti=0

            
            if  presenti:

                d1=a1-a2
                d2=b1-b2

                c1=(a1+a2)/2.0
                c2=(b1+b2)/2.0
                
                D1=np.sqrt((d1*d1).sum())
                D2=np.sqrt((d2*d2).sum())
                
                DC =  ((c1-c2)*(c1-c2)).sum()   
                DD=(D1-D2)*(D1-D2)

                result=result+ (DC+DD)-((d1*d1).sum()+(d1*d1).sum())/2.0

    for i in range(NX-1):
        for j in range(NY-1):
            try:
                a1=dic[(i,j)]
                b1=dic[(i+1,j)]
                b2=dic[(i-1,j)]
                presenti=1
            except:
                presenti=0
            if  presenti:

                v1=b2-a1
                v2=b1-a1
                Area = v1[0]*v2[1]-v1[1]*v2[0]
                result=result+ abs(Area)


    for i in range(NX-1):
        for j in range(NY-1):
            try:
                a1=dic[(i,j)]
                b1=dic[(i,j+1)]
                b2=dic[(i,j-1)]
                presenti=1
            except:
                presenti=0
            if  presenti:

                v1=b2-a1
                v2=b1-a1
                Area = v1[0]*v2[1]-v1[1]*v2[0]
                result=result+ abs(Area)

    return result


def register(distorted):

    refgrid=np.zeros([4,3,2],"f")
    refgrid[:,:,1] = np.arange(3)
    refgrid[:,:,0] = np.reshape(np.arange(4),[4,1])
    coors=np.reshape(refgrid,[-1,2])
    targets=[]
    possibili=[]
    match(distorted,refgrid, [tuple(tok) for tok in coors], targets, possibili)
    print( possibili)
    best=1.0e30
    coors=coors.astype("i")
    for choice in possibili:
        m=merit( distorted , refgrid  ,[tuple(tok.tolist()) for tok in  coors],  choice)
        if m<best:
            best=m
            choosen=choice
    return choosen


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Qt4Agg" )
    import pylab


    refgrid=np.zeros([4,3,2],"f")
    refgrid[:,:,1] = np.arange(3)
    refgrid[:,:,0] = np.reshape(np.arange(4),[4,1])
    distorted = np.reshape(refgrid+0.6*np.random.random(refgrid.shape),[12,2])*1.2
    distorted=np.array(np.random.permutation(distorted))
    distorted=distorted[:-1]
    choosen = register( distorted  ) 
    # pylab.plot(coors[:,0],coors[:,1],"o")
    pylab.plot(distorted[:,0],distorted[:,1],"o")
    pylab.xlim( -0.5,5.5   )
    pylab.ylim( -0.5,5.5   )
    
    labels=[ "(%d %d)"%tuple(tok)  for tok in  choosen]
    
    
    for label, x, y in zip(labels, distorted[:, 0], distorted[:, 1]):
        pylab.annotate(
            label, 
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))




    pylab.show()







