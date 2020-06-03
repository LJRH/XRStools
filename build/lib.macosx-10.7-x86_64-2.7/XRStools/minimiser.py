from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


#/*##########################################################################
# Copyright (C) 2004-2012 European Synchrotron Radiation Facility
#
# PPM  : Alessandro Mirone.
# GPU   Cuda     ( OCL is in progress.. )    Dimitris Karkoulis
#   Qt : Interface : Vorobeva Anastasiya   Vorobyeva Veronika
#           nvorobeva@hotmail.fr  vorobyevav@yahoo.com 
#                     and Alessandro Mirone
#  European Synchrotron Radiation Facility, Grenoble,France
#
#
# PPM is  developed at
# the ESRF by the SciSoft  group.
# PPM CUDA is developed by  Dimitris Karkoulis, financed by:
#        LinkSCEEM-2 (INFRA-2010-1.2.3) Work Package 12 project 
#           (grant number RI-261600)
#
# This toolkit is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) 
# any later version.
#
# PPM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# PPM; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA 02111-1307, USA.
#
# PPM follows the dual licensing model of Trolltech's Qt and Riverbank's PyQt
# and cannot be used as a free plugin for a non-free program. 
#
# Please contact the ESRF industrial unit (industry@esrf.fr) if this license 
# is a problem for you.
#############################################################################*/



############################################################
# Alessandro Mirone
# April 2001
# ESRF
from math import exp
import math
from six.moves import map
from six.moves import range

def difforzero(a):
     if(a<0):
          return 0
     else:
          return a




ISPARALLEL=0
ISMASTER=1
try:
     import mpi
     ISPARALLEL=1
     if(mpi.rank==0):
       ISMASTER=1
     else:
       ISMASTER=0
except:
     pass


from numpy  import *
import random
import os
import glob
import re
from PyMca5 import Gefit

colpo=[10]

try:
  import beansgui
  beans_OK=1
except:
  beans_OK=0


def update( var, x):
     var.value=(var.max+var.min)*0.5 + (var.max-var.min)*math.sin(x)*0.5

def getx( var):
     x=( 2*var.value -var.min - var.max)/(var.max-var.min)
     return math.asin(x)


# def update( var, x):
#      var.value=(   (1-tanh(x))*var.min+(1+tanh(x))*var.max)/2.0

# def getx( var):
#      x=( 2*var.value -var.min - var.max)/(var.max-var.min)
#      # print "x ", x
#      if (x==1): return 40
#      if(x==-1): return -40
#      if(x>1 or x<-1):
#           print var.value," ", var.min," ", var.max
#           raise " out of range"
     
#      x=0.5*log(   (1+x) /(1-x)  )
#      return x

def generateverticesNormal(X,dx):
     pointlist=[array(X)]
     for i in range(0,len(X)):
        nX=array(X)

        nX[i]=nX[i]+dx
        pointlist.append(nX)
     return array(pointlist) 

def generateverticesRandomNormal(X,dx, centershiftFact=0):
     pointlist=[array(X)]
     add=0
     centershift=0
     if( ISPARALLEL==0 or mpi.rank==0):
       centershift=centershiftFact*dx*array([ (-0.5+random.random()) for i in range(0, len(X)) ]     )
     if( ISPARALLEL ):
       centershift =mpi.bcast(centershift ,0   )
       

       
     for   i in range(0,len(X)):
       if( ISPARALLEL==0 or mpi.rank==0):
         add=dx[i]*(-0.5+random.random())+centershift[i]
       if( ISPARALLEL ):
         add =mpi.bcast(add ,0   )
       pointlist[0][i]+=add
     for i in range(0,len(X)):
        nX=array(X)
        if( ISPARALLEL==0 or mpi.rank==0):
          add=dx[i]*(-0.5+random.random())+centershift[i]
        if( ISPARALLEL ):
          add =mpi.bcast(add ,0   )
        nX[i]=nX[i]+add
        pointlist.append(nX)
     return array(pointlist) 
   

class minimiser:
  extrapotfactor=10.0
  stoppa=[0,0]      #  emergency stop, delayed stop
  history=[]
  def __init__(self, EcartObject, VariableList):
    self.EcartObject=EcartObject
    self.VariableList=[]

    for tok in VariableList:
      if(isinstance(tok,list)):
        self.VariableList = self.VariableList + tok
      else:
           if tok.min<tok.max:
                self.VariableList.append(tok)

    minimiser.history=[]

  def calculateerror(self, x):
     colpo[0]=colpo[0]+1
     for i in range (0, len(self.VariableList)):
         update( self.VariableList[i], x[i])
     extrapot=0.0
     limit=3.0
     for i in range (0, len(self.VariableList)):
        # print x[i]
        if(abs(x[i])>limit):
           extrapot=extrapot+(exp(abs(x[i]))-exp(limit))
     # print " EXTRAPOT__ " , extrapot*self.extrapotfactor

     if( beans_OK):
       if( beansgui.BLock is not None):
          beansgui.BLock.acquire()

     error = self.EcartObject.error()

     if( beans_OK):
       if( beansgui.BLock is not None):
          beansgui.BLock.release()



     # print " error ", error

     return error +extrapot*self.extrapotfactor

  
    
  def amoeba(self, ftol, arret=None):
    X=[]
    for var in self.VariableList:
      X.append(getx(var) )
    x=array(X)
    p=generateverticesNormal(X, 0.1)
    y=list(map(self.calculateerror,p))
    return self.amoebaspecial(ftol,p,y, arret=arret)

  
  def amoebaAnnealed(self, ftol, temperature_fn, max_refusedcount=100,
                     centershiftFact=0.0,
                     arret=None,
                     max_isthesame=10,
                     expandby_if_isthesame=3.0,
                     shrinkby_if_accepted =0.7,
                     triggerby_if_accepted=1,
                     proportionality_whenlower=2.5,
                     maxium_ratio_max_min_whenlower=33.33,
                     callbackIT=None,
                     callbackLocal=None,
                     maxiters=100000
                     ):
    count=0

    minimiser.history=[]
    history=self.history
    ###############################################
    # If EcartObject has the method setitercounter
    # then set it
    if(hasattr(self.EcartObject, 'setitercounter' ) ):
        self.EcartObject.setitercounter(0)

    ###############################################
    #  initialization of x variable.
    #  x maps [-inf,+inf] to the limited interval
    #  that has been assigned to the given variable
    #  X contains the x for every variable.
    #
    self.history.append(str(list(map(par,self.VariableList))))
    X=[]
    for var in self.VariableList:
      X.append(getx(var) )
    p=generateverticesNormal(X, 0.1)
    y=list(map(self.calculateerror,p))
    self.history[-1]=self.history[-1 ]+ str([" ERROR ", y[0] ])


    if( ISPARALLEL==0 or mpi.rank==0):
         fperformance=open("performances","a")
         fperformance.write("At ITERATION %d, ERROR %e \n" % (0, y[0] ))
         fperformance.close()

         if( hasattr(self.EcartObject,"writeminimum")):
                      self.EcartObject.writeminimum( 0 ,  y[0]  )





    res=self.amoebaspecial(ftol,p,y, arret=arret, callbackIT=callbackIT)
    if( callbackLocal is not None):
            callbackLocal()



    X=res[0]
    Y=res[1]   
    Yold=Y

    Xmin=X
    Ymin=res[1]
    lower=1
    
    Dpar=ones(len(X))*0.1    
    refusedcount=0
    isthesamecount=0

    while(refusedcount<max_refusedcount  and  isthesamecount<max_isthesame
                  and self.stoppa[0]==0  and self.stoppa[1]==0 and count < maxiters
          ) :
       count=count+1
       if(lower ):
         if( count > 1):
            Dpar=proportionality_whenlower*abs(res[0]-Xold)
            Dpar=Dpar+max(Dpar)/maxium_ratio_max_min_whenlower
       else:
          if(isthesame):
             Dpar*=expandby_if_isthesame
          else:
             mask=less(Dpar,triggerby_if_accepted)
             Dpar=Dpar*mask+(1-mask)*triggerby_if_accepted
             Dpar*=shrinkby_if_accepted
       Xold=X*1.0
       mask=less(Dpar,10.0)
       Dpar=Dpar*mask+(1-mask)*10

       # ---print " X " , X
       # ---print " Xold " , Xold


       p=generateverticesRandomNormal(X, Dpar, centershiftFact)
       # ---print p
       y=list(map(self.calculateerror,p))



       res=self.amoebaspecial(ftol,p,y, arret=arret, callbackIT=callbackIT)
       if( callbackLocal is not None):
         callbackLocal()
    
       self.history.append(str(list(map(par,self.VariableList))))
       self.history[-1]= self.history[-1]+ str([" ERROR ", res[1] ])


       ###############################################
       # If EcartObject has the method getitercounter
       # then use it getitercounter
       if(hasattr(self.EcartObject,'getitercounter') ):
            if( ISPARALLEL==0 or mpi.rank==0):
                 fperformance=open("performances","a")
                 fperformance.write("At ITERATION %d, ERROR %e \n" % (self.EcartObject.getitercounter(), res[1] ))
                 fperformance.close()
                 if( hasattr(self.EcartObject,"writeminimum")):
                      self.EcartObject.writeminimum( self.EcartObject.getitercounter(), res[1]   )

       # ---print "********************************************"
       # ---print " newY = ", res[1]
       # ---print "    Y = ", Y

       # ---print " Ymin= ", Ymin
       if(abs(res[1]-Y)/abs(res[1]+Y)<4*ftol):
         isthesame=1
         isthesamecount=isthesamecount+1
         accepted=0
         lower=0
         X=res[0]*1
         Y=res[1]
       else:
         isthesame=0
         if(res[1]<Y):
            Y=res[1]
            lower=1
            accepted=1
            isthesamecount=0
            if(Y<Ymin):
               Ymin=Y
               Xmin=res[0]*1
         else:
            lower=0
            espo=(res[1]-Y)/temperature_fn(count)
            if(espo>30):
              accepted = 0
            prob=exp(-espo)

            if( ISPARALLEL==0 or mpi.rank==0):
              probcomp=random.random()
            if( ISPARALLEL ):
              probcomp =mpi.bcast( probcomp,0   )


            
            if(probcomp<prob):
              isthesamecount=0
              accepted=1
            else:
              accepted=0
         if(accepted):
            X=res[0]*1
            Y=res[1]
            refusedcount=0
         else:
            refusedcount+=1
       # ---print " accepted = ", accepted
       # ---print " refusedcount = ", refusedcount
       # ---print " TEMP = ",  temperature_fn(count)
       # ---print " count = " , count
       # ---print " isthesame ", isthesame
       # ---print " isthesamecount ", isthesamecount
       # ---print " lower = ", lower
       # ---print "Dpar = ", Dpar
       if(arret != None and Ymin<arret):
          print( " mi fermo " , Ymin, " " , Y)
          break
    for i in range (0, len(self.VariableList)):
     update( self.VariableList[i], Xmin[i])
         





  #########################################################3
  # run an amoeba minimisation taking as initial simplex p
  # with initial values given by y.
  # It stops when either 
  # 2*abs(y[ihighest]-y[ilowest])/(abs(y[ihighest])+abs(y[ilowest]))<ftol
  #
  # or  y[ilowest]< arret if arret!=None
  #
  def amoebaspecial(self, ftol,p,y, arret=None, callbackIT=None):
    ndim = len(y)-1   

    iter=0
    maxiter=10000
    while(1):
      ilowest=ihighest=isecondhighest=0
      for i in range(0, len(p)):
         if( y[i]>y[ihighest]):
                isecondhighest=ihighest
                ihighest=i
         elif(y[i]>y[isecondhighest]):
                isecondhighest=i
         if( y[i]<y[ilowest] ):
                ilowest=i
      center= (sum(p, axis=0)-p[ihighest])/ndim
      # print "-----------------------"
      # print p
      # print y
      # print "------------------------"

      rtol=2*abs(y[ihighest]-y[ilowest])/(abs(y[ihighest])+abs(y[ilowest]))
      if(rtol<ftol): break
      # print dir()
      if (  self.stoppa[0] ): break

      
      if(arret!=None and arret>y[ilowest]): break
      iter=iter+1
      if(iter>maxiter):
         return (p[ilowest]*1.0,y[ilowest])

      ptry= 2*center- p[ihighest]
      # print ptry
      ytry=self.calculateerror(ptry)
      if( callbackIT is not None):
        callbackIT()
      if(ytry< y[ilowest]):
          ptry2= 2*ptry-center 

          # print ptry2
          ytry2= self.calculateerror(ptry2)
          if( callbackIT is not None):
            callbackIT()
          if(ytry2 < y[ilowest] ):
              p[ihighest]=ptry2
              y[ihighest]=ytry2
          else:
              p[ihighest]=ptry
              y[ihighest]=ytry
      elif( ytry>=y[isecondhighest]):
          if(ytry<y[ihighest]):
             p[ihighest]=ptry
             y[ihighest]=ytry
          ptry=0.5*p[ihighest]+0.5*center
          # print ptry
          ytry=self.calculateerror(ptry)
          if( callbackIT is not None):
            callbackIT()
          if(ytry<y[ihighest]):
             y[ihighest]=ytry
             p[ihighest]=ptry
          else:
            for i in range(0,ndim+1):
               p[i]=0.5*(p[i]+p[ilowest])
               # print p[i]
               y[i]=self.calculateerror(p[i])
      else:
          p[ihighest]=ptry
          y[ihighest]=ytry
    return (p[ilowest]*1.0,y[ilowest])


class VariableShifted:
     def getvalue(self):
          return self.shift+self.varia.getvalue()
     def __init__(self, varia, shift):
          self.varia=varia
          self.shift=shift


class Variable:
     def getvalue(self):
        return self.value

     def __init__(self, value,min,max):
        self.shift=0.0
        self.value=value*1.0
        self.max=max*1.0
        self.min=min*1.0
        # self.getvalue=self.getvalue

     def setrandom(self):
       probcomp=0
       if( ISPARALLEL==0 or mpi.rank==0):
         probcomp=random.random()
       if( ISPARALLEL ):
         probcomp =mpi.bcast( probcomp,0   )
       self.value=self.min+(self.max-self.min)*probcomp



     def __repr__(self):
         res="Variable(%e,%e,%e)"%(self.value, self.min,self.max)
         return res

 
     def __add__(self, other):
         return VariableShifted(self, other)

def CreateVariableArray(Np=None, Xmin=None,Xmax=None, value=1.0, min=0.7, max=1.4):
    res=[]
    for i in range(Np):
         X=(Xmin*(Np-i-1.0) +Xmax*(i+0.0)  )/(Np-1)
         ob = Variable(value, min,max)
         res.append(ob)
         ob.X=X
    return res

GDIC={}

class DependentVariable:
     depth=0.0
     def getvalue(self):
        command="res="+self.formula
        exec(command)
        return res

     def VarCounter( self):
         if( self.vardepthstart>=self.depth):
            self.vardepthstart = self.depth
            self.varcounter=0
            self.varlastdepth=-1

         result = self.varcounter
         if(result==0):
             self.vardepthstart=self.depth

         if( self.depth >  self.varlastdepth):
             self.varcounter=self.varcounter+1

         self.varlastdepth=self.depth

    # print " for ", self, " DependentVariable.VarCounter has been called, result is ", result
    # print " depth is ", self.depth
         return result

     def __init__(self,stringa , dictio={}):

        self.varcounter=0
        self.vardepthstart=0
        self.varlastdepth=0


##################################################################################3
        posold=0
        while(1):
             pos = stringa[posold:].find("par(")
             if pos==-1: break
             pos2= stringa[posold+pos:].find(")")
             pos2=pos2+pos
             key  = stringa[posold+pos+4:posold+pos2]
             if key in dictio:
                  stringa1=stringa[:posold+pos+4]+"self."+key+")"
                  stringa=stringa1+stringa[posold+pos2+1:]
                  posold=len(stringa1)
                  setattr(self,key,dictio[key])
             else:
                  msg= " KEY "+  key+ "\n STRINGA " +stringa+"\n error processing dependent variable string " 
                  
                  raise Exception(msg)
  
        # for key in dictio.keys():
        #     while(1):
        #         res=re.search("par *\( *%s *\)"%key, stringa  )
        #         if(res is not None):
        #             stringa=stringa[0:res.start()]+"par(self.%s)"%key + stringa[res.end():]
        #             setattr(self,key,dictio[key])
        #         else:
        #             break


        self.getvalue=self.getvalue
        self.formula=stringa

#     def __repr__(self):
#         res="DependentVariable ->getvalue() = %e"% self.getvalue()
#         return res



def par(a):
    if( hasattr(a,'getvalue' )):
          return a.getvalue()
    else:
          return a      



class wrapperForFit:
     
     def __init__(self, fit, variable_list):
          self.fit=fit
          self.variable_list=variable_list

     def get_angles(self):
          res=[]
          for v in self.variable_list:
               ss=(v.max+v.min)*0.5
               sd=(v.max-v.min)*0.5
               if sd==0:
                    angle=0
               else:
                    angle = arcsin((v.value-ss)/sd)
               res.append(angle)
          return res

     
     def apply_angles(self,angles):
          for i in range(len(self.variable_list)):
               v=self.variable_list[i]
               x=angles[i]
               ss=(v.max+v.min)*0.5
               sd=(v.max-v.min)*0.5
               v.value =  ss+sd*sin(x)
               
               
     def theory(self, paras , xs=None, ):
          self.apply_angles(paras)
          colpo[0]=colpo[0]+1
          self.fit.error()
          return self.fit.GetResForGeFit()

     def getData(self):
           return self.fit.GetDataForGeFit()


if __name__=='__main__':




  class scarto:
    def __init__(self, a,b,c):
       self.a=a
       self.b=b
       self.c=c
    def error(self):
      # print par(self.a)
      # print dir(self.a)
      res= 1+par(self.a)*par(self.a)+par(self.b)*par(self.b)+par(self.c)*par(self.c)
      # print "error=", res
      return res


  a=Variable(2.,5.,-6.)
  b=Variable(7.,9.,-10.)
  c=Variable(7.,10.,-11.)

  ecob=scarto(a,b,c)

  miniob=minimiser(ecob,[a,b,c])

  (p,y)=miniob.amoeba(0.000000001)  

  miniob.amoebaAnnealed( 0.001, lambda x: 1000.0/(1.0+x*0.3), 100,
      arret=1.0e-12, maxiters=10) 

  # ---print p
  # ---print y


  class scartoF2:
    def __init__(self, a,b):
       self.a=a
       self.b=b
    def error(self):
      x1=par(self.a)
      x2=par(self.b)
      a=x2-x1*x1
      res= 100*a*a + (1-x1)*(1-x1)
      # ---print "error=", res, " itc= ", self.itc
      self.itc=self.itc+1
      return res
    def setitercounter(self, n):
      self.itc=n
    def getitercounter(self):
      return self.itc


  for ini in [ (1001,1001), (1001,-999),  (-999,-999),    (-999,1001),
                 (1443,1), (1,1443),(1.2,1)]   :

    a=Variable(ini[0],-2000.048,2000.048)
    b=Variable(ini[1],-2000.048,2000.048)

    tominimise=scartoF2(a,b)

    miniob=minimiser(tominimise,[a,b])

    os.system("rm performances")
    
    miniob.amoebaAnnealed( 0.001, lambda x: 1000.0/(1.0+x*0.3), 100,
    arret=1.0e-12) 
    
    s=open("performances","r").read()
    x1=par(a)
    x2=par(b)
    dist=(x1-1)*(x1-1)+(x2-1)*(x2-1)
    dist=sqrt(dist)

    open("resocontoRosF2d","a").write( "##### ini=(%f %f)\n%s// finale= %e   x1= %e x2=%e dist=%e\n" %
      ( ini[0],ini[1], s, tominimise.error(), x1,x2,dist) )









