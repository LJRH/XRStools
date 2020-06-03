from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: superresolution.py

from .helpers import *

import os
import numpy as np
import pylab
from scipy import io
from itertools import groupby
from scipy.interpolate import Rbf, RectBivariateSpline
from scipy.optimize import leastsq, fmin

from  .xrs_imaging import *
from six.moves import range

class imageset:
    """
    class to make SR-images from list of LR-images
    """
    def __init__(self):
        self.list_of_images = []
        self.xshifts  = []
        self.yshifts  = []
        self.shifts   = []
        self.srimage  = []
        self.srxscale = []
        self.sryscale = []
        self.refimagenum = []

    def estimate_xshifts(self,whichimage=None):
        if not whichimage:
            ind = 0 # first image in list_of_images is the reference image
        else:
            ind = whichimage
        origx   = self.list_of_images[ind].xscale
        origy   = self.list_of_images[ind].yscale
        origim  = self.list_of_images[ind].matrix
        xshifts = []        
        for image in self.list_of_images:
            newx  = image.xscale
            newy  = image.yscale
            newim = image.matrix
            xshifts.append(estimate_xshift(origx,origy,origim,newx,newy,newim))
        self.refimagenum = ind
        self.xshifts = xshifts

    def estimate_yshifts(self,whichimage=None):
        if not whichimage:
            ind = 0 # first image in list_of_images is the reference image
        else:
            ind = whichimage
        origx   = self.list_of_images[ind].xscale
        origy   = self.list_of_images[ind].yscale
        origim  = self.list_of_images[ind].matrix
        yshifts = []        
        for image in self.list_of_images:
            newx  = image.xscale
            newy  = image.yscale
            newim = image.matrix
            yshifts.append(estimate_yshift(origx,origy,origim,newx,newy,newim))
        self.refimagenum = ind
        self.yshifts = yshifts

    def estimate_shifts(self,whichimage=None):
        if not whichimage:
            ind = 0 # first image in list_of_images is the reference image
        else:
            ind = whichimage
        origx   = self.list_of_images[ind].xscale
        origy   = self.list_of_images[ind].yscale
        origim  = self.list_of_images[ind].matrix
        shifts = []        
        for image in self.list_of_images:
            newx  = image.xscale
            newy  = image.yscale
            newim = image.matrix
            shifts.append(estimate_shift(origx,origy,origim,newx,newy,newim))
        self.refimagenum = ind
        self.shifts = shifts

    def interpolate_xshift_images(self,scaling,whichimages=None):
        if not whichimages:
            inds = list(range(len(self.list_of_images)))
        elif not isinstance(whichimages,list):
            inds = []
            inds.append(whichimages)
        else:
            inds = whichimages

        newim = np.zeros((len(inds),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling,np.shape(self.list_of_images[inds[0]].matrix)[1]))
        newx  = np.linspace(self.list_of_images[inds[0]].xscale[0]-self.xshifts[inds[0]],self.list_of_images[inds[0]].xscale[-1]-self.xshifts[inds[0]],len(self.list_of_images[inds[0]].xscale)*scaling)
        
        newy  = self.list_of_images[inds[0]].yscale

        for n in range(len(inds)):
            print( self.xshifts[inds[n]])
            oldim = self.list_of_images[inds[n]].matrix
            oldx  = self.list_of_images[inds[n]].xscale-self.xshifts[inds[n]]
            oldy  = self.list_of_images[inds[n]].yscale
            newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

        self.srimage = np.sum(newim,axis=0)
        self.srxscale = newx
        self.sryscale = newy

    def interpolate_yshift_images(self,scaling,whichimages=None):
        if not whichimages:
            inds = list(range(len(self.list_of_images)))
        elif not isinstance(whichimages,list):
            inds = []
            inds.append(whichimages)
        else:
            inds = whichimages

        newim = np.zeros((len(inds),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling,np.shape(self.list_of_images[inds[0]].matrix)[1]))
        newx  = self.list_of_images[0].xscale
        newy  = np.linspace(self.list_of_images[inds[0]].yscale[0]-self.yshifts[inds[0]],self.list_of_images[inds[0]].yscale[-1]-self.yshifts[inds[0]],len(self.list_of_images[inds[0]].yscale)*scaling)
        for n in range(len(inds)):
            oldim = self.list_of_images[inds[n]].matrix
            oldx  = self.list_of_images[inds[n]].xscale
            oldy  = self.list_of_images[inds[n]].yscale-self.yshifts[inds[n]]
            newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

        self.srimage = np.sum(newim,axis=0)
        self.srxscale = newx
        self.sryscale = newy

    def interpolate_shift_images(self,scaling,whichimages=None):
        if not whichimages:
            inds = list(range(len(self.list_of_images)))
        elif not isinstance(whichimages,list):
            inds = []
            inds.append(whichimages)
        else:
            inds = whichimages

        if len(scaling)<2:
            scaling = [scaling, scaling]
        print( inds, self.list_of_images[inds[0]].xscale[0], self.shifts[inds[0]], self.list_of_images[inds[0]].xscale[-1])
        newim = np.zeros((len(self.list_of_images),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling[0],np.shape(self.list_of_images[inds[0]].matrix)[1]*scaling[1]))
        newx  = np.linspace(self.list_of_images[inds[0]].xscale[0]-self.shifts[inds[0]][0],self.list_of_images[inds[0]].xscale[-1]-self.shifts[inds[0]][0],len(self.list_of_images[inds[0]].xscale)*scaling[0])
        newy  = np.linspace(self.list_of_images[inds[0]].yscale[0]-self.shifts[inds[0]][1],self.list_of_images[inds[0]].yscale[-1]-self.shifts[inds[0]][1],len(self.list_of_images[inds[0]].yscale)*scaling[1])

        for n in range(len(inds)):
            oldim = self.list_of_images[inds[n]].matrix
            oldx  = self.list_of_images[inds[n]].xscale-self.shifts[inds[n]][0]
            oldy  = self.list_of_images[inds[n]].yscale-self.shifts[inds[n]][1]
            newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

        self.srimage = np.sum(newim,axis=0)
        self.srxscale = newx
        self.sryscale = newy

    def plotSR(self):
        X,Y = pylab.meshgrid(self.srxscale,self.sryscale)
        pylab.pcolor(X,Y,np.transpose(self.srimage))
        pylab.show(block=False)

    def plotLR(self,whichimage):
        if not isinstance(whichimage,list):
            inds = []
            inds.append(whichimage)
        else:
            inds = list(whichimage)

        for ind in inds:
            X,Y = pylab.meshgrid(self.list_of_images[ind].xscale,self.list_of_images[ind].yscale)
            pylab.figure(ind)
            pylab.pcolor(X,Y,np.transpose(self.list_of_images[ind].matrix))
            pylab.show(block=False)

    def save(self):
        pass

    def load(self):
        pass

    def loadkimberlite(self,matfilename):
        data_dict   = io.loadmat(matfilename)
        sorted_keys = sorted(data_dict.keys())
        sy    = data_dict['sy'][0]
        allsx = []
        for key in sorted_keys[3:12]:
            allsx.append(data_dict[key][0])

        allmats = []
        for key in sorted_keys[13:]:
            allmats.append(data_dict[key])

        alllengths = []
        for sx in allsx:
            alllengths.append(len(sx))
        
        ind = np.where(alllengths == np.max(alllengths))[0][0]
        # spline everything onto longest sx-scale
        for n in range(len(allmats)):
            print( np.shape(allsx[n]), np.shape(sy), np.shape(allmats[n]))
            ip         = RectBivariateSpline(allsx[n],sy,allmats[n])
            allmats[n] = ip(allsx[ind],sy)
            allsx[n]   = allsx[ind]

        allimages = []
        for n in range(len(allmats)):
            allimages.append(LRimage(allmats[n],allsx[n],sy))
        self.list_of_images = allimages

    def loadhe3070(self,matfilename):
        data_dict  = io.loadmat(matfilename)
        sy    = data_dict['det'][0][0]['sy'][0][0][0]
        allsx = []
        allmats = []
        for n in range(9):
            allsx.append(np.reshape(data_dict['det'][0][n]['sx'][0][0]-data_dict['det'][0][n]['sx'][0][0][0],len(data_dict['det'][0][n]['sx'][0][0],)))
            allmats.append(data_dict['det'][0][n]['img'][0][0])

        alllengths = []
        for sx in allsx:
            alllengths.append(len(sx))

        ind = np.where(alllengths == np.max(alllengths))[0][0]
        for n in range(len(allmats)):
            print( np.shape(allsx[n]), np.shape(sy), np.shape(np.transpose(allmats[n])))
            ip         = RectBivariateSpline(allsx[n],sy,np.transpose(allmats[n]))
            allmats[n] = ip(allsx[ind],sy)
            allsx[n]   = allsx[ind]

        allimages = []
        for n in range(len(allmats)):
            allimages.append(LRimage(allmats[n],allsx[n],sy))
        self.list_of_images = allimages

    def loadimage(self,filename):
        f = open(filename,'rb')
        self.list_of_images.append(pickle.load(f))
        f.close()

    def correct_scattering_angle(self,whichimage,tth):
        pass









