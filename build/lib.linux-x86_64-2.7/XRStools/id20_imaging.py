from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: id20_imaging.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import io
from itertools import groupby
from scipy.interpolate import Rbf, RectBivariateSpline
from scipy.optimize import leastsq, fmin

from helpers import *
from .xrs_read import rois
from six.moves import range

class imaging:
    """
    a class to make images from id20 data, it will use some of the methods the read_id20 class has also,
    maybe it would make sense to use read_id20 as a superclass for this one and just add the imaging
    functionality
    """
    def __init__(self,absfilename,energycolumn='energy_cc',monitorcolumn='monitor'):
        self.scans         = {} # dictionary of scans that are loaded
        if not os.path.isfile(absfilename):
            raise Exception('IOError! No such file, please check filename.')
        self.path          = os.path.split(absfilename)[0] + '/'
        self.filename      = os.path.split(absfilename)[1]
        self.scannumbers   = []
        self.EDF_PREFIX    = 'edf/'
        self.EDF_POSTFIX   = '.edf'
        self.DET_PIXEL_NUM = 256
        self.TTH_OFFSETS   = np.array([-13.0,-6.5,-6.5,0.0,0.0,6.5,6.5,13.0,13.0])
        # which column in the SPEC file to be used for the energy and monitor
        self.encolumn      = energycolumn.lower()
        self.monicolumn    = monitorcolumn.lower()
        # ROIs: only rectangular rois should be possible (line or zoom ROIs)
        self.rois     = []
        # the stored images
        self.images   = {} # dictionary of images, which will be stored under their scan name

    def readscan(self,scannumber):
        """
        returns the data, motors, counter-names, and edf-files from the calss's specfile of "scannumber"
        should use PyMca's routines to read the Spec- and edf-files
        """
        fn = self.path + self.filename
        # loading spec file
        print ("parsing edf- and SPEC-files of scan No. %s" % scannumber)
        data, motors, counters = specread(fn,scannumber)
        # loading according edf files
        edfmats = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUM,self.DET_PIXEL_NUM)))
        for m in range(len(counters['ccdno'])):
            ccdnumber = counters['ccdno'][m]+1
            edfname   = self.path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
            edfmats[m,:,:] = edfread(edfname)
        self.scannumbers.extend([scannumber]) # keep track of the scnanumbers
        return data, motors, counters, edfmats

    def loadscan(self,scannumbers,scantype='generic'):
        """
        loads the files belonging to scan No "scannumber" and puts it into an instance
        of the container-class scan. the default scantype is 'generic', later the scans
        will be grouped (then added) based on the scantype
        """
        # make sure scannumbers are iterable (list)
        if not isinstance(scannumbers,list):
            scannums = []
            scannums.append(scannumbers)
        else:
            scannums = scannumbers 
        for number in scannums:
            scanname = 'Scan%03d' % number
            data, motors, counters, edfmats = self.readscan(number)
            # can assign some things here already (even if maybe redundant)
            monitor   = counters[self.monicolumn]
            # energy related columns cannot be assigned now, since scans here could also be position scans (like z-scans)
            monoangle = None #counters['pmonoa']
            energy    = None #counters[self.encolumn]
            # create an instance of the "scan" class for every scan
            onescan = scan(edfmats,number,energy,monoangle,monitor,counters,motors,data,scantype)
            # assign one dictionary entry to each scan 
            self.scans[scanname] = onescan

    def createrois(self,scannumbers):
        """
        create rois object from this class and scannumbers
        """
        rois_object = rois(self.scans,scannumbers)
        return rois_object

    def getlinrois(self,scannumbers,numrois=9,logscaling=True,height=5,colormap='jet'):
        """
        define ROIs by clicking two points
        """
        rois_object = self.createrois(scannumbers)
        rois_object.getlinrois(numrois,logscaling,height,colormap)
        # asign the entire rois object to self.rois to also have access to the roikind attribute e.g.
        self.rois = rois_object

    def getzoomrois(self,scannumbers,numrois=9,logscaling=True,colormap='jet'):
        """
        define ROIs by zooming in on a region of interest
        """
        rois_object = self.createrois(scannumbers)
        rois_object.getzoomrois(numrois,logscaling,colormap)
        # asign the entire rois object to self.rois to also have access to the roikind attribute e.g.
        self.rois = rois_object

    def image_from_line(self,scannumber,scanmotor,whichroi=[0]):
        """
        create images from a point focus, either along a sample position (scanmotor = 'sx', 'sz', ...)
        or along energy (scanmotor = 'energy_cc', 'nrj', ...) 
        """
        scanname = 'Scan%03d' % scannumber

        # make mask out of rois (this can be deleted once the ROIs are always handled as masks of zeros and ones)
        maskrois = []
        for roi in self.rois.rois:
            mask = np.zeros((self.DET_PIXEL_NUM,self.DET_PIXEL_NUM))
            for pixel in roi:
                mask[pixel[1],pixel[0]] = 1
            maskrois.append(mask)

        roisalongscan = [] # list of rectangular matricies (shaped like the rois) for each roi (len(scan)*width*height matrix)
        lineimages    = [] # list of images from lines (len(scan)*width)
        for n in whichroi: # for each roi
            inds = np.nonzero(maskrois[n])
            roialongscan = []
            lineimage     = []
            for edfmat in self.scans[scanname].edfmats: # each edf file of a given scan
                roialongscan.append(np.reshape(edfmat[inds],(self.rois.roiheight[n],self.rois.roiwidth[n]))) # the whole roi
                lineimage.append(np.sum( np.reshape(edfmat[inds],(self.rois.roiheight[n],self.rois.roiwidth[n])),axis=0))  # summed over height
            roisalongscan.append(roialongscan)
            lineimages.append(np.squeeze(lineimage))
        print( np.shape(roisalongscan), np.shape(lineimages))
        plt.imshow(lineimages[0])


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
            allimages.append(image(allmats[n],allsx[n],sy))
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
            allimages.append(image(allmats[n],allsx[n],sy))
        self.list_of_images = allimages

class image:
    """
    container class to hold info of a single LR-image to be put togther in a SR-image by the imageset class
    """
    def __init__(self,matrix,xscale,yscale):
        self.name    = []
        self.matrix    = matrix
        self.xscale    = xscale
        self.yscale    = yscale
        self.tth    = []

    def plotimage(self):
        pass

    def shiftx(self):
        pass

    def shifty(self):
        pass

    def save(self):
        pass

    def load(self):
        pass


def interpolate_image(oldx,oldy,oldIM,newx,newy):
    """
    2d interpolation
    """
    interp = RectBivariateSpline(oldx,oldy,oldIM)
    return interp(newx,newy)

def estimate_xshift(x1,y1,im1,x2,y2,im2):
    """
    estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
    """
    funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2+a,y2))**2.0)
    res = leastsq(funct,0.0)
    return res[0]

def estimate_yshift(x1,y1,im1,x2,y2,im2):
    """
    estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
    """
    funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2,y2+a))**2.0)
    res = leastsq(funct,0.0)
    return res[0]

def estimate_shift(x1,y1,im1,x2,y2,im2):
    """
    estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
    """
    funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2+a[0],y2+a[1]))**2.0)
    res = fmin(funct,[0.0,0.0],disp=0)
    return res


class LRimage:
    """
    container class to hold info of a single LR-image to be put togther in a SR-image by the imageset class
    """
    def __init__(self,matrix,xscale,yscale):
        self.name    = []
        self.matrix    = matrix
        self.xscale    = xscale
        self.yscale    = yscale
        self.tth    = []

    def plotimage(self):
        pass

    def shiftx(self):
        pass

    def shifty(self):
        pass

    def save(self):
        pass

    def load(self):
        pass








