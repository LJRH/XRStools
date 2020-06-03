from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import range
from six.moves import input
#!/usr/bin/python
# Filename: extraction.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

import os
import math

import numpy as np
import array as arr
import matplotlib.pyplot as plt
import pickle

from matplotlib.widgets import Cursor
from itertools import groupby
from scipy.integrate import trapz
from scipy import interpolate, signal, integrate, constants, optimize
from re import findall
from scipy.ndimage import measurements
from scipy.optimize import leastsq, fmin
from scipy.interpolate import Rbf, RectBivariateSpline

from scipy.integrate import odeint


#fcomp = 1/(-i*lex)*(-2*((abb0*(abb8 + abb7*sgbeta*t) + abb1) + i*y0) *(y(1) + i*y(2)) + c1*(1 + (y(1) + i*y(2)).^2));





    

#print mpr_compds([1,2,3],['SiO2','CO2'],[0.5,0.5],9.86,[1 1])



# xrs_read


class scan:
    """
    this is a container class for inelastic x-ray scattering scans with 2D detectors    
    each scan is an instance of this class, the scans of an experiment can then be 
    grouped and ordered by this class's attributes scantype, starting energy, etc. 
    """
    def __init__(self,mats,num,en,monoa,moni,counts,mots,data,scantype='generic'):
        # rawdata
        self.edfmats  = np.array(mats)
        self.number   = num
        self.scantype = scantype
        self.energy   = np.array(en)
        self.monoangle= np.array(monoa)
        self.monitor  = np.array(moni)
        # some things maybe not imediately necessary
        self.counters = np.array(counts)
        self.motors   = mots
        self.specdata = np.array(data)
        # data (to be filled after defining rois)
        self.eloss    = []
        self.signals  = []
        self.errors   = []
        self.signals_orig = [] # keep a copy of uninterpolated data
        # would like to keep uninterpolated signals/eloss/errors, too
        
    def applyrois(self,rois,scaling=None):
        """
        sums up each 2D matrix of a scan over the indices given in rois,
        i.e. turns the 3D matrix of a scan (stack of 2D detector images)
        into a matrix of size (len(energy),number of rois)
        rois are a list of tuples
        """
        data = np.zeros((len(self.edfmats),len(rois)))
        for n in range(len(rois)): # each roi (default is 9)
            for m in range(len(self.edfmats)): # each energy point along the scan
                for l in range(len(rois[n])): # each pixel on the detector
                    data[m,n] += self.edfmats[m,rois[n][l][1],rois[n][l][0]]
        self.signals = np.array(data)
        self.signals_orig = np.array(data)
        if np.any(scaling):
            assert len(scaling) == len(rois) # make sure, there is one scaling factor for each roi
            for ii in range(len(rois)):
                self.signals[:,ii] *= scaling[ii]
                self.signals_orig[:,ii] *= scaling[ii]

    def applyrois_old(self,rois):
        """
        sums up each 2D matrix of a scan over the indices given in rois,
        i.e. turns the 3D matrix of a scan (stack of 2D detector images)
        into a matrix of size (len(energy),number of rois)
        rois are a list of tuples
        this seems a bit faster than the old version
        """
        data     = np.zeros((len(self.edfmats),len(rois)))
        roixinds = []
        roiyinds = []
        for r in range(len(rois)):
            for n in range(len(rois[r])):
                roixinds.append(rois[r][n][0])
                roiyinds.append(rois[r][n][1])
            data[:,r] = np.sum(np.sum(self.edfmats[:,np.amin(roiyinds):np.amax(roiyinds)+1,np.amin(roixinds):np.amax(roixinds)+1],axis=1),axis=1)
        self.signals = data

    def gettype(self):
        return self.scantype

    def getnumber(self):
        return self.number

    def getshape(self):
        if not np.any(self.signals):
            print( 'please apply the ROIs first.')
            return
        else:
            return np.shape(self.signals)

    def getnumofrois(self):
        if not self.signals.any():
            print( 'please apply the ROIs first.')
            return
        else:
            return np.shape(self.signals)[1]

class scangroup:
    """
    container class for scans with the same 'scantype' attribute
    each group of scans will summed into an instance of this class
    different goups of scans will then be stitched together based on 
    their type, starting energy, etc. 
    """
    def __init__(self,energy,signals,errors,grouptype='generic'):
        self.energy    = energy
        self.eloss     = []
        self.signals   = signals
        self.errors    = errors
        self.grouptype = grouptype
        self.signals_orig = signals # keep a copy of uninterpolated data

    def gettype(self):
        return self.grouptype

    def getmeanenergy(self):
        return np.mean(self.energy)

    def getestart(self):
        return self.energy[0]

    def geteend(self):
        return self.energy[-1]

    def getmeanegridspacing(self):
        return np.mean(np.diff(self.energy))

    def getmaxediff(self):
        return (self.energy[-1]-self.energy[0])

# superresolution/imaging


#################################
# ROIs
#################################

class rois:
    """
    a class to define ROIs 
    ROIs are saved in a container (self.rois)
    """
    def __init__(self,scans,scannumbers):
        self.scandata = scans # dictionary with instances of the onescan class
        if not isinstance(scannumbers,list):
            thescannumbers = []
            thescannumbers.append(scannumbers)
        else:
            thescannumbers = scannumbers  
        self.scannums   = thescannumbers
        self.DET_PIXEL_NUM = 256
        self.rois       = container() # single variable that holds all relevant info about the defined rois
        self.rois.inds  = [] # list of lists of touples (for each ROI)
        self.rois.xinds = [] # list of lists of x-indices (for each ROI)
        self.rois.yinds = [] # list of lists of y-indices (for each ROI)
        self.rois.roiNumber = 0 # number of ROIs defined
        self.rois.kind    = [] # what kind of ROI (auto, zoom, linear)
        self.rois.height  = [] # list of hights (in pixels, for each ROI)
        self.rois.width   = [] # list of widths (in pixels, for each ROI)

    def preparemats(self):
        """
        sums and squeezes all edf-files of a scan/all scans into one matrix 
        """
        # take shape of the first edf-matrices
        dim = np.shape(self.scandata['Scan%03d' % self.scannums[0]].edfmats)
        edfmatrices = np.zeros((1,dim[1],dim[2]))
        for number in self.scannums:
            scanname = 'Scan%03d' % number
            edfmatrices = np.append(edfmatrices,self.scandata[scanname].edfmats,axis=0)
        return sumx(edfmatrices)
    
    def prepareedgemats(self,index):
        """
        difference between two summed and squeezed edf-files of a scan from below and above energyval 
        """
        # take shape of the first edf-matrices 
        dim = np.shape(self.scandata['Scan%03d' % self.scannums[0]].edfmats)
        abovemats = np.zeros((1,dim[1],dim[2]))
        belowmats = np.zeros((1,dim[1],dim[2]))
        for number in self.scannums:
            scanname = 'Scan%03d' % number
            abovemats = np.append(abovemats,self.scandata[scanname].edfmats[index:,:,:],axis=0)
            belowmats = np.append(belowmats,self.scandata[scanname].edfmats[:index,:,:],axis=0)
        return np.absolute(np.squeeze(np.sum(abovemats,axis=0))-np.squeeze(np.sum(belowmats,axis=0)))

    def getlinrois(self,numrois,logscaling=True,height=5,colormap='jet'):
        """
        define ROIs by clicking two points
        """
        # clear all existing rois
        self.deleterois()

        plt.clf()

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, axisbg='#FFFFCC')

        thematrix = self.preparemats()
        plt.title('Click start and endpoints of the linear ROIs you whish.',fontsize=14)        
        if logscaling:
            theimage = plt.imshow(np.log(np.absolute(thematrix)))
        else:
            theimage = plt.imshow(thematrix)
        # Choose a color palette
        theimage.set_cmap(colormap)
        plt.colorbar()

        therois   = [] # list of lists of tuples with (x,y) coordinates
        theroisx  = [] # list of lists of x-indices
        theroisy  = [] # list of lists of y-indices

        cursor = Cursor(ax, useblit=True, color='red', linewidth=1 )

        for m in range(numrois):
            plt.clf()
            if m>0:
                if therois[m-1]:
                    for index in therois[m-1]:
                        thematrix[index[1],index[0]] = np.amax(thematrix)
            print ("chose two points as endpoints for linear ROI No. %s" % (m+1))

            if logscaling:
                theimage = plt.imshow(np.log(np.absolute(thematrix)))
            else:
                theimage = plt.imshow(thematrix)    
            theimage.set_cmap(colormap)
            plt.draw()    

            endpoints = []
            endpts    = plt.ginput(2)
            for point in endpts:
                point = [round(element) for element in point]
                endpoints.append(np.array(point))
            roix = np.arange(endpoints[0][0],endpoints[1][0])
            roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]

            theheight = np.arange(-height,height)
            
            eachroi  = []
            eachroix = []
            eachroiy = []
            for n in range(len(roix)):
                for ind in theheight:
                    eachroi.append((roix[n],roiy[n]+ind))
                    eachroix.append(roix[n])
                    eachroiy.append(roiy[n]+ind)
            therois.append(eachroi)
            theroisx.append(eachroix)
            theroisy.append(eachroiy)

            self.rois.height.append(len(theheight)) # save the hight of each roi (e.g. for imaging)
            self.rois.width.append(endpoints[1][0]-endpoints[0][0]) # save the width of each roi
            del endpoints
            del eachroi
            del eachroix
            del eachroiy
            
        self.rois.kind    = 'linear'
        self.rois.inds  = therois
        self.rois.xinds = theroisx
        self.rois.yinds = theroisy
        self.rois.roiNumber = numrois

    def getzoomrois(self,numrois,logscaling=True,colormap='jet'):
        """
        define ROIs by zooming in on plot
        """
        # clear all existing rois
        self.deleterois()

        plt.clf()
        thematrix = self.preparemats()
        limmax    = np.shape(thematrix)[1]
        therois   = [] # list of lists of tuples with (x,y) coordinates
        theroisx  = [] # list of lists of x-indices
        theroisy  = [] # list of lists of y-indices
        plt.title('Zoom around ROI, change to shell, and press Return')
        if logscaling:
            theimage = plt.imshow(np.log(np.absolute(thematrix)))
        else:
            theimage = plt.imshow(thematrix)    
        theimage.set_cmap(colormap)
        plt.ion()    
        for m in range(numrois):
            plt.clf()
            if m>0:
                if therois[m-1]:
                    for index in therois[m-1]:
                        thematrix[index[0],index[1]] = np.amax(thematrix)
            plt.title('Zoom around ROI, change to shell, and press Return')
            if logscaling:
                theimage = plt.imshow(np.log(np.absolute(thematrix)))
            else:
                theimage = plt.imshow(thematrix)    
            theimage.set_cmap(colormap)
            plt.draw()    
                        
            thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
            wait = input(thestring)
                                
            limits = np.floor(plt.axis())
            
            inds = limits < 1
            limits[inds] = 1
            inds = limits > limmax
            limits[inds] = limmax
            print( 'selected limits are: ', limits)

            T = np.zeros(np.shape(thematrix))
            T[limits[3]:limits[2],:] = 1
            t = np.zeros(np.shape(thematrix))
            t[:,limits[0]:limits[1]] = 1
            indsx,indsy = np.where(t*T == 1)
            eachroi  = []
            eachroix = []
            eachroiy = []
            for n in range(len(indsx)):
                eachroi.append((indsy[n],indsx[n]))
                eachroix.append(indsx[n])
                eachroiy.append(indsy[n])
            therois.append(eachroi)
            theroisx.append(eachroix)
            theroisy.append(eachroiy)
        plt.show()
        self.rois.kind      = 'zoom'
        self.rois.inds      = therois
        self.rois.xinds     = theroisx
        self.rois.yinds     = theroisy
        self.rois.roiNumber = numrois

    def getzoomrois_frommatrix(self,matrix,numrois,logscaling=True,colormap='jet'):
        """
        define ROIs by zooming in on plot
        """
        # clear all existing rois
        self.deleterois()

        plt.clf()
        thematrix = matrix
        limmax = np.shape(thematrix)[1]
        therois = []
        plt.title('Zoom around ROI, change to shell, and press Return')
        if logscaling:
            theimage = plt.imshow(np.log(np.absolute(thematrix)))
        else:
            theimage = plt.imshow(thematrix)    
        theimage.set_cmap(colormap)
        plt.ion()
        for m in range(numrois):
            
            plt.title('Zoom around ROI, change to shell, and press Return')
            if logscaling:
                theimage = plt.imshow(np.log(np.absolute(thematrix)))
            else:
                theimage = plt.imshow(thematrix)    
            theimage.set_cmap(colormap)
            plt.ion()

            thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
            wait = input(thestring)

            limits = np.floor(plt.axis())
            plt.axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0])

            inds = limits < 1
            limits[inds] = 1
            inds = limits > limmax
            limits[inds] = limmax
            print( 'selected limits are: ', limits)

            T = np.zeros(np.shape(thematrix))
            T[limits[3]:limits[2],:] = 1
            t = np.zeros(np.shape(thematrix))
            t[:,limits[0]:limits[1]] = 1
            indsx,indsy = np.where(t*T == 1)
            eachroi = []
            for n in range(len(indsx)):
                eachroi.append((indsy[n],indsx[n]))
            therois.append(eachroi)

        plt.show()
        self.rois.kind      = 'zoom'
        self.rois.inds      = therois
        self.rois.roiNumber = numrois

    def getlinedgerois(self,index,numrois=9,logscaling=True,height=5,colormap='jet'):
        """
        define ROIs by clicking two points in difference picture
        """
        # clear all existing rois
        self.deleterois()

        plt.clf()
        thematrix = self.prepareedgemats(index) # change this to use an index
        if logscaling:
            theimage = plt.imshow(np.log(np.absolute(thematrix)))
        else:
            theimage = plt.imshow(thematrix)

        # Choose a color palette
        theimage.set_cmap(colormap)
        plt.colorbar()

        therois = []
        for m in range(numrois):
            print ("chose two points as endpoints for linear ROI No. %s" % (m+1))
            endpoints = []
            endpts    = plt.ginput(2)
            for point in endpts:
                point = [round(element) for element in point]
                endpoints.append(np.array(point))
            roix = np.arange(endpoints[0][0],endpoints[1][0])
            roiy = [round(num) for num in np.polyval(np.polyfit([endpoints[0][0],endpoints[1][0]],[endpoints[0][1],endpoints[1][1]],1),roix)]
            del endpoints
            theheight = np.arange(-height,height)
            eachroi = []
            for n in range(len(roix)):
                for ind in theheight:
                    eachroi.append((roix[n],roiy[n]+ind))
            therois.append(eachroi)
            del eachroi            
        plt.show()
        self.rois.kind      = 'linear'
        self.rois.inds      = therois
        self.rois.roiNumber = numrois

    def getzoomedgerois(self,index,numrois=9,logscaling=True,colormap='jet'):
        """
        define ROIs by zooming in on plot
        """
        # clear all existing rois
        self.deleterois()

        plt.clf()
        thematrix = self.prepareedgemats(index)
        limmax = np.shape(thematrix)[1]
        therois = []
        title('Zoom around ROI, change to shell, and press Return')
        if logscaling:
            theimage = plt.imshow(np.log(np.absolute(thematrix)))
        else:
            theimage = plt.imshow(thematrix)
        theimage.set_cmap(colormap)
        plt.ion()
        for m in range(numrois):

            plt.title('Zoom around ROI, change to shell, and press Return')
            if logscaling:
                theimage = plt.imshow(np.log(np.absolute(thematrix)))
            else:
                theimage = plt.imshow(thematrix)
            theimage.set_cmap(colormap)
            plt.ion()
                        
            thestring = 'Zoom around ROI No. %s and press enter to continue.' % (m+1)
            wait = input(thestring)
                                
            limits = np.floor(axis())
            plt.axis([0.0,self.DET_PIXEL_NUM,self.DET_PIXEL_NUM,0.0])

            inds = limits < 1
            limits[inds] = 1
            inds = limits > limmax
            limits[inds] = limmax
            print( 'selected limits are: ', limits)

            T = np.zeros(np.shape(thematrix))
            T[limits[3]:limits[2],:] = 1
            t = np.zeros(np.shape(thematrix))
            t[:,limits[0]:limits[1]] = 1
            indsx,indsy = np.where(t*T == 1)
            eachroi = []
            for n in range(len(indsx)):
                eachroi.append((indsy[n],indsx[n]))
            therois.append(eachroi)
        plt.show()
        self.rois.kind      = 'zoom'
        self.rois.inds      = therois
        self.rois.roiNumber = numrois

    def getautorois(self,kernel_size=5,colormap='jet',thresholdfrac=100):
        """
        define ROIs by choosing a threshold through a slider bar on the plot window
        """
        # clear all existing rois
        self.deleterois()

        plt.clf()
        thematrix = self.preparemats() # the starting matrix to plot
        plt.title('Crank the threshold and close the plotting window, when you are satisfied.',fontsize=14)        
        ax = plt.subplot(111) 
        plt.subplots_adjust(left=0.05, bottom=0.2)
        thres0 = 0 # initial threshold value
        theimage = plt.imshow(np.log(np.absolute(thematrix)))
        theimage.set_cmap(colormap)
        plt.colorbar()
        thresxcolor = 'lightgoldenrodyellow'
        thresxamp  = plt.axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
        sthres = plt.Slider(thresxamp, 'Threshold', 0.0, np.floor(np.amax(thematrix)/thresholdfrac), valinit=thres0)

        def update(val):
            thres = sthres.val
            newmatrix = signal.medfilt2d(thematrix, kernel_size=kernel_size)
            belowthres_indices = newmatrix < thres
            newmatrix[belowthres_indices] = 0            
            labelarray,numfoundrois = measurements.label(newmatrix)
            print( str(numfoundrois) + ' ROIs found!'    )
            # organize rois
            therois = []
            for n in range(numfoundrois):
                rawindices = np.nonzero(labelarray == n+1)
                eachroi = []                
                for m in range(len(rawindices[0])):
                    eachroi.append((rawindices[1][m],rawindices[0][m]))
                therois.append(eachroi)
            self.rois.inds = therois
            self.rois.roiNumber = numfoundrois
            theimage.set_data(newmatrix)
            plt.draw()
        sthres.on_changed(update)
        plt.show()
        self.rois.kind = 'auto'

    def getautorois_eachdet(self,kernel_size=5,colormap='jet',thresholdfrac=100):
        """
        autoroi, detector for detector for ID20 (6 detector setup) so that different 
        thresholds can be choosen for the different detectors
        define ROIs by choosing a threshold through a slider bar on the plot window
        """
        # clear all existing rois
        self.deleterois()

        wholematrix = np.array(self.preparemats()) # full 6-detector image
        imagelims = [[0,256,0,256],[0,256,256,512],[0,256,512,768],[256,512,0,256],[256,512,256,512],[256,512,512,768]] 

        for n in range(len(imagelims)):
            plt.clf()
            # thematrix now is a single detector image (256x256 pixels)
            thematrix = wholematrix[imagelims[n][0]:imagelims[n][1],imagelims[n][2]:imagelims[n][3]]
            plt.title('Crank the threshold and close the plotting window, when you are satisfied.',fontsize=14)        
            ax = plt.subplot(111) 
            plt.subplots_adjust(left=0.05, bottom=0.2)
            thres0 = 0 # initial threshold value
            theimage = plt.imshow(np.log(np.absolute(thematrix)))
            theimage.set_cmap(colormap)
            plt.colorbar()
            thresxcolor = 'lightgoldenrodyellow'
            thresxamp  = plt.axes([0.2, 0.10, 0.55, 0.03], axisbg=thresxcolor)
            sthres = plt.Slider(thresxamp, 'Threshold', 0.0, np.floor(np.amax(thematrix)/thresholdfrac), valinit=thres0)

            # using the "container" class to pass the results of the rois from the nested function
            roi_result = container()

            def update(val):
                thres = sthres.val
                newmatrix = signal.medfilt2d(thematrix, kernel_size=kernel_size)
                belowthres_indices = newmatrix < thres
                newmatrix[belowthres_indices] = 0            
                labelarray,numfoundrois = measurements.label(newmatrix)
                print( str(numfoundrois) + ' ROIs found!'    )
                # organize rois
                therois = []
                for n in range(numfoundrois):
                    rawindices = np.nonzero(labelarray == n+1)
                    eachroi = []
                    for m in range(len(rawindices[0])):
                        eachroi.append((rawindices[1][m],rawindices[0][m]))
                    therois.append(eachroi)
                theimage.set_data(newmatrix)
                plt.draw()
                roi_result.therois =  therois
                roi_result.numfoundrois = numfoundrois

            sthres.on_changed(update)
            plt.show()
            thestring = 'Press enter to continue.'
            wait = input(thestring)
            self.rois.inds.extend(roi_result.therois)
            self.rois.roiNumber += roi_result.numfoundrois

        self.rois.kind = 'auto'

    def saverois_old(self,filename):
        """
        save the ROIs in file with name filename
        """    
        f = open(filename, 'wb')
        theobject = [self.rois, self.roinums, self.roikind]
        pickle.dump(theobject, f,protocol=-1)
        f.close()

    def loadrois_old(self,filename):
        """
        load ROIs from file with name filename
        """
        f = open(filename,'rb')
        theobject = pickle.load(f)
        f.close()
        self.rois    = theobject[0]
        self.roinums = theobject[1]
        self.roikind = theobject[2]

    def saverois(self,filename):
        """
        save the ROIs in file with name filename
        """    
        f = open(filename, 'wb')
        pickle.dump(self.rois, f,protocol=-1)
        f.close()

    def loadrois(self,filename):
        """
        load ROIs from file with name filename
        """
        f = open(filename,'rb')
        self.rois = pickle.load(f)
        f.close()
        
    def deleterois_old(self):
        """
        delete the existing ROIs
        """
        self.rois    = []
        self.roinums = []
        self.roikind = []

    def deleterois(self):
        """
        delete existing ROIs (same as at initialization)
        """
        self.rois       = container() # single variable that holds all relevant info about the defined rois
        self.rois.inds  = [] # list of lists of touples (for each ROI)
        self.rois.xinds = [] # list of lists of x-indices (for each ROI)
        self.rois.yinds = [] # list of lists of y-indices (for each ROI)
        self.rois.roiNumber = 0 # number of ROIs defined
        self.rois.kind    = [] # what kind of ROI (auto, zoom, linear)
        self.rois.height  = [] # list of hights (in pixels, for each ROI)
        self.rois.width   = [] # list of widths (in pixels, for each ROI)

    def plotrois(self):
        """
        returns a plot of the ROI shapes
        """
        plt.figure()
        thematrix = np.zeros((self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
        for roi in self.rois:
            for index in roi:
                thematrix[index[1],index[0]] = 1
        theimage = plt.imshow(thematrix)
        plt.show()

class container:
    """
    random container class to hold values
    """
    def __init__(self):
        pass



#########################################################################
#
# trash
#

#    def createrois(self,scannumbers):
#        """
#        Creates a Rois object from information from this class and scannumbers.
#        INPUT:
#        scannumbers = integer or list of integers of scannumbers from which should be used to define the ROIs. 
#        """
#        rois_object = xrs_rois.rois(self.scans,scannumbers)
#        return rois_object
#
#    def getautorois(self,scannumbers,kernel_size=5,colormap='jet',threshfrac=100):
#        """
#        Define ROIs automatically using median filtering and a variable threshold.
#       INPUT:
#        scannumbers = 
#        kernel_size =
#        colormap    =
#        threshfrac  = 
#        (This function should be deprecated once the new style of ROI definitions work properly)
#        """
#        rois_object = self.createrois(scannumbers)
#        rois_object.getautorois(kernel_size,colormap,thresholdfrac=threshfrac)
#        self.rois = rois_object.rois
#        # add the new xrs_rois's functions and classes
#        new_obj = xrs_rois.roi_object()
#        new_obj.indices = self.rois.inds
#        #print np.shape(self.rois.inds), np.shape(new_obj.indices)
#        self.roi_obj = new_obj
#
#    def getautorois_eachdet(self,scannumbers,kernel_size=5,colormap='jet',threshfrac=100):
#        """
#        define ROIs automatically one detector at a time using median filtering and a variable threshold
#        (This function should be deprecated once the new style of ROI definitions work properly)
#        """
#        rois_object = self.createrois(scannumbers)
#        rois_object.getautorois_eachdet(kernel_size,colormap,thresholdfrac=threshfrac)
#        self.rois = rois_object.rois        
#        # add the new xrs_rois's functions and classes
#        new_obj = xrs_rois.roi_object()
#        new_obj.indices = rois_object.rois.inds
#        self.roi_obj = new_obj
#                
#    def getlinrois(self,scannumbers,numrois=72,logscaling=True,height=5,colormap='jet'):
#        """
#        define ROIs by clicking two points
#        (This function should be deprecated once the new style of ROI definitions work properly)
#        """
#        rois_object = self.createrois(scannumbers)
#        rois_object.getlinrois(numrois,logscaling,height,colormap)
#        self.rois = rois_object.rois
#        # add the new xrs_rois's functions and classes
#        new_obj = xrs_rois.roi_object()
#        new_obj.indices = self.rois.inds
#        self.roi_obj = new_obj
        
#    def getzoomrois(self,scannumbers,numrois=72,logscaling=True,colormap='jet'):
#        """
#        define ROIs by zooming in on a region of interest
#        (This function should be deprecated once the new style of ROI definitions work properly)
#        """
#        rois_object = self.createrois(scannumbers)
#        rois_object.getzoomrois(numrois,logscaling,colormap)
#        self.rois = rois_object.rois
#        # add the new xrs_rois's functions and classes
#        new_obj = xrs_rois.roi_object()
#        new_obj.indices = rois_object.rois.inds
#        self.roi_obj = new_obj
#
#    def getlinedgerois(self,scannumbers,energyval,numrois=72,logscaling=True,height=5,colormap='jet'):
#        """
#        define ROIs by clicking two points on a matrix, which is a difference of edf-files above and 
#        below an energy value
#        (This function should be deprecated once the new style of ROI definitions work properly)
#        """
#        if not isinstance(scannumbers,list):
#            scannums = []
#            scannums.append(scannumbers)
#        else:
#            scannums = scannumbers
#        # search for correct index in the first of the given scans
#        scanname = 'Scan%03d' % scannums[0]
#        index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
#        rois_object = self.createrois(scannumbers)
#        rois_object.getlinedgerois(index,numrois,logscaling,height,colormap)
#        self.rois = rois_object.rois
#        # add the new xrs_rois's functions and classes
#        new_obj = xrs_rois.roi_object()
#        new_obj.indices = rois_object.rois.inds
#        self.roi_obj = new_obj
#
#    def getzoomedgerois(self,scannumbers,energyval,numrois=72,logscaling=True,colormap='jet'):
#        """
#        define ROIs by zooming in on a matrix, which is a difference of edf-files above and 
#        below an energy value
#        (This function should be deprecated once the new style of ROI definitions work properly)
#        """
#        if not isinstance(scannumbers,list):
#            scannums = []
#            scannums.append(scannumbers)
#        else:
#            scannums = scannumbers
#        # search for correct index in the first of the given scans
#        scanname = 'Scan%03d' % scannums[0]
#        index = self.scans[scanname].energy.flat[np.abs(self.scans[scanname].energy - energyval).argmin()]
#        rois_object = self.createrois(scannumbers)
#        rois_object.getlinedgerois(index,numrois,logscaling,colormap)
#        self.rois = rois_object.rois
#        # add the new xrs_rois's functions and classes
#        new_obj = xrs_rois.roi_object()
#        new_obj.indices = rois_object.rois.inds
#        self.roi_obj = new_obj
#
#    def saverois(self,filename):
#        """
#        saves the rois into an file using pickle
#        filename = absolute path to file and filename
#        """
#        import pickle
#        f = open(filename, 'wb')
#        pickle.dump(self.rois, f,protocol=-1)
#        f.close()








