from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import range
from six.moves import zip
from six.moves import input
#!/usr/bin/python
# Filename: rixs_read.py

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

from . import xrs_rois, xrs_scans, xrs_utilities, math_functions, xrs_fileIO
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# try to import the fast PyMCA file parsers

try:
    import PyMca5.PyMcaIO.EdfFile as EdfIO
    import PyMca5.PyMcaIO.specfilewrapper as SpecIO
    use_PyMca = True
except:
    use_PyMca = False

if not use_PyMca:
    
    try:
        import PyMca.EdfFile as EdfIO
        import PyMca.specfilewrapper as SpecIO
        use_PyMca = True
    except:
        use_PyMca = False


print( " >>>>>>>>  use_PyMca " , use_PyMca)

__metaclass__ = type # new style classes

class read_id20:
    """
    Main class for handling raw data from XES experiments on ESRF's ID20. This class
    is used to read scans from SPEC files and the according EDF-files, it provides access
    to all tools from the xrs_rois module for defining ROIs, it can be used to integrate
    scans, sum them up, stitch them together, and define the energy loss scale.
    INPUT:
    absFilename   = path and filename of the SPEC-file
    energyColumn  = name (string) of the counter for the energy as defined in the SPEC session (counter mnemonic)
    monitorColumn = name (string) of the counter for the monitor signals as defined in the SPEC session (counter mnemonic)
    edfName       = name/prefix (string) of the EDF-files (default is the same as the SPEC-file name)
    single_image  = boolean switch, 'True' (default) if all 6 detectors are merged in a single image,
                    'False' if two detector images per point exist.
    """
    def __init__(self,absFilename,energyColumn='Anal Energy',monitorColumn='kaprixs',edfName=None, edf_shape=(255, 1295), EinCoor=[0,0] ):

        self.scans         = {} # dictionary of scans

        try:
            self.path     = os.path.split(absFilename)[0] + '/'
            self.filename = os.path.split(absFilename)[1]
        except IOError:
            print('file does not exist.')

        if not edfName:
            self.edfName = os.path.split(absFilename)[1]
        else:
            self.edfName = edfName
            
        self.scannumbers   = []
        self.EDF_PREFIX    = 'edf/'
        self.EDF_POSTFIX   = '.edf'
        self.DET_PIXEL_NUMx = edf_shape[1]
        self.DET_PIXEL_NUMy = edf_shape[0]
        self.DET_PIXEL_NUM  = edf_shape[0]

        # which column in the SPEC file to be used for the energy and monitor
        self.encolumn      = energyColumn.lower()  # name of energy motor scanned
        self.monicolumn    = monitorColumn.lower() # name of the monitor counter
        self.EinCoor       = EinCoor               # incident energy coordinate

        # compensation factor
        self.comp_factor   = 0.0
        self.pixel_size    = 0.055 # pixel size in mm

        # here are the attributes of the old rawdata class
        self.energy   = []    # common energy scale for all analyzers
        self.energy2  = []  # analyzer energy 
        self.signals  = []    # signals for all analyzers
        self.errors   = []    # poisson errors

        self.groups   = {}  # dictionary of groups (such as 2 'elastic', or 5 'edge1', etc.)
        self.tth      = []  # list of scattering angles (one for each ROI)
        self.resolution   = []  # list of FWHM of the elastic lines for each analyzer
        self.signals_orig = []  # signals for all analyzers before interpolation

        # ROI object
        self.roi_obj = [] # an instance of the roi_object class from the xrs_rois module (new)

    def set_roiObj(self,roiobj):
        self.roi_obj = roiobj

    def readscan(self,scannumber):
        """
        Returns the data, motors, counter-names, and edf-files from the SPEC file defined when
        the xrs_read object was initiated.
        There should be an alternative that uses the PyMca module if installed.
        INPUT:
        scannumber = number of the scan to be loaded
        fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
        """
        # load SPEC-file
        print( 'Parsing EDF- and SPEC-files of scan No. %s' % scannumber)
        fn = self.path + self.filename
        if use_PyMca == True:
            data, motors, counters = xrs_fileIO.PyMcaSpecRead(fn,scannumber)
        else:
            data, motors, counters = xrs_fileIO.SpecRead(fn,scannumber)
        # load EDF-files
        edfmats = xrs_fileIO.ReadEdfImages(counters['ccdno'], self.DET_PIXEL_NUMx, self.DET_PIXEL_NUMy, self.path, self.EDF_PREFIX, self.edfName, self.EDF_POSTFIX)
        # add the scannumber to self.scannumbers, if not already present
        if not scannumber in self.scannumbers:
            self.scannumbers.extend([scannumber])
        return data, motors, counters, edfmats 

























#---------



    def SumDirect(self,scannumbers):
        Sum=None
        for number in scannumbers:
            data, motors, counters, edfmats = readscan(number, self.path, self.filename, self.DET_PIXEL_NUMx, self.DET_PIXEL_NUMy, self.EDF_PREFIX, self.EDF_POSTFIX, self.edfName)
            if Sum is None:
                Sum = np.zeros(edfmats[0].shape ,"f") 
            Sum[:] += edfmats.sum(axis=0) 
        return Sum

    def loadscan(self, scannumbers, scantype='generic', direct=True):
        """
        Loads the files belonging to scan No. "scannumber" and puts it into an instance
        of the xrs_scan-class 'scan'. The default scantype is 'generic', later the scans
        will be grouped (and added) based on the scantype.
        INPUT:
        scannumbers = integer or list of scannumbers that should be loaded
        scantype    = string describing the scan to be loaded (e.g. 'edge1' or 'K-edge') 
        """
        # make sure scannumbers are iterable
        if not isinstance(scannumbers,list):
            scannums = []
            scannums.append(scannumbers)
        else:
            scannums = scannumbers 

        # load scan/scans
        for number in scannums:
            scanname = 'Scan%03d' % number
            data, motors, counters, edfmats = readscan(number, self.path, self.filename, self.DET_PIXEL_NUMx, self.DET_PIXEL_NUMy, self.EDF_PREFIX, self.EDF_POSTFIX, self.edfName)
            # create an instance of "scan" class for every scan
            onescan = xrs_scans.scan(edfmats,number,counters[self.encolumn],counters[self.monicolumn],counters,motors,data,scantype)
            # assign one dictionary entry to each scan 
            self.scans[scanname] = onescan
            if not number in self.scannumbers:
                self.scannumbers.extend([number])
            # save the incident energy
            self.scans[scanname].Ein = motors[self.EinCoor[0]][self.EinCoor[1]]

            if direct:
                energy = self.scans[scanname].energy * 1e3 # energy in eV
                self.scans[scanname].signals = np.zeros((len(energy),self.roi_obj.number_of_rois))
                self.scans[scanname].errors  = np.zeros((len(energy),self.roi_obj.number_of_rois))
                for key,col in zip(self.roi_obj.red_rois, list(range(self.roi_obj.number_of_rois)) ):
                    (pos,M) = self.roi_obj.red_rois[key]
                    S = M.shape
                    x = np.arange(S[0])*self.pixel_size
                    inset = (slice(pos[0]  , pos[0]+(S[0])   ), slice(  pos[1]  , pos[1]+(S[1]) ) )
                    A =  self.scans[scanname].edfmats[:,inset[0],inset[1]]
                    y = np.squeeze(np.sum(A,2))
                    meanmon = np.mean(self.scans[scanname].monitor)
                    yn = np.zeros_like(y)
                    for jj in range(len(self.scans[scanname].monitor)):
                        yn[jj,:]=y[jj,:]/self.scans[scanname].monitor[jj]*meanmon
                    meanii = len(list(range(yn.shape[1])))//2
                    yc = np.zeros_like(y)
                    for ii in range(yn.shape[1]):
                        sort  = np.argsort(energy)
                        yc[:,ii] = np.interp(energy[sort]+(ii-meanii)*self.comp_factor*self.pixel_size, energy[sort],yn[sort,ii],left=float('nan'),right=float('nan'))
                    self.scans[scanname].signals[:,col] = xrs_utilities.nonzeroavg(yc)
                print('Deleting EDF-files of scan No. %03d'%number)
                self.scans[scanname].edfmats = np.array([])

    def deletescan(self,scannumbers):
        """
        Deletes scans from the class.
        INPUT:
        scannumbers = integer or list of integers (SPEC scan numbers) to delete
        """
        # make sure scannumbers are iterable
        numbers = []
        if not isinstance(scannumbers,list):
            numbers.append(scannumbers)
        else:
            numbers = scannumbers
        # delete scans
        for number in numbers:
            scanname = 'Scan%03d' % number
            del(self.scans[scanname])
            self.scannumbers.remove(number)


    def SumDirect(self,scannumbers):
        Sum=None
        for number in scannumbers:
            data, motors, counters, edfmats = self.readscan(number)
                        
            print( edfmats.shape)
                        
            if Sum is None:
                Sum = np.zeros(edfmats[0].shape ,"f") 
            Sum[:] += edfmats.sum(axis=0) 
        return Sum

    def getCompensationFactor(self,scannumber,roiNumber):
        if not self.roi_obj:
            print('Please set a ROI object first.')
            return
        else:
            data, motors, counters, edfmats = readscan(scannumber, self.path, self.filename, self.DET_PIXEL_NUMx, self.DET_PIXEL_NUMy, self.EDF_PREFIX, self.EDF_POSTFIX, self.edfName)
            plt.cla()
            plt.ion()
            el_positions = np.zeros(len(edfmats))
            for ii in range(len(edfmats)):
                key = list(self.roi_obj.red_rois.keys())[roiNumber]
                (pos,M) = self.roi_obj.red_rois[key]
                S = M.shape
                x = np.arange(S[0])*self.pixel_size+self.pixel_size
                inset = (slice(pos[0]  , pos[0]+(S[0])   ), slice(  pos[1]  , pos[1]+(S[1]) ) )
                y = np.sum( edfmats[ii,inset[0],inset[1]] ,axis=1)
                try:
                    popt = optimize.curve_fit(math_functions.gauss_forcurvefit, x, y)[0]
                    g = math_functions.gauss_forcurvefit(x,popt[0],popt[1],popt[2])
                    el_positions[ii] = (popt[1])
                except:
                    g = np.zeros_like(x)
                    el_positions[ii] = 0.0
                plt.cla()
                plt.plot(x,y,x,g)
                plt.draw()

            #print S, len(el_positions), len(edfmats), y.shape[0]
            #center = np.interp(range(len(edfmats)),range(y.shape[0]),el_positions)

            # fit the slope
            plt.cla()
            sort =  np.argsort(np.array(counters[self.encolumn])*1e3)
            energy = (np.array(counters[self.encolumn])*1e3)[sort]
            center = el_positions[sort]
            plt.plot(center,energy)
            plt.draw()
            input('Zoom in and press enter to continue')
            limits = plt.axis()
            inds = np.where(np.logical_and(center >= limits[0], center <= limits[1]))[0] 
            fact = np.polyfit( center[inds], energy[inds],1)
            plt.plot(center,energy,center,np.polyval(fact,center))
            self.comp_factor = comp_factor = fact[0]




def readscan(scannumber, path, filename, DET_PIXEL_NUMx, DET_PIXEL_NUMy, EDF_PREFIX, EDF_POSTFIX, edfName):
    """
    Returns the data, motors, counter-names, and edf-files from the SPEC file defined when
    the xrs_read object was initiated.
    There should be an alternative that uses the PyMca module if installed.
    INPUT:
    scannumber = number of the scan to be loaded
    fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
    """
    # load SPEC-file
    print( 'Parsing EDF- and SPEC-files of scan No. %s' % scannumber)
    fn = path + filename
    #if use_PyMca == True:
    #    data, motors, counters = xrs_fileIO.PyMcaSpecRead(fn,scannumber)
    #else:
    data, motors, counters = xrs_fileIO.SpecRead(fn,scannumber)
    # load EDF-files
    edfmats = xrs_fileIO.ReadEdfImages(counters['ccdno'], DET_PIXEL_NUMx, DET_PIXEL_NUMy, path, EDF_PREFIX, edfName, EDF_POSTFIX)
    return data, motors, counters, edfmats 









