from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: xes_read.py

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

# try to import the fast PyMCA file parsers
try:
    import PyMca5.PyMcaIO.EdfFile as EdfIO
    import PyMca5.PyMcaIO.specfilewrapper as SpecIO
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
    def __init__(self,absFilename,energyColumn='Anal Energy',monitorColumn='kaprixs',edfName=None):
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
        self.DET_PIXEL_NUMx = 1296
        self.DET_PIXEL_NUMy = 256
        self.DET_PIXEL_NUM  = 256

        # which column in the SPEC file to be used for the energy and monitor
        self.encolumn      = energyColumn.lower()
        self.monicolumn    = monitorColumn.lower()

        # here are the attributes of the old rawdata class
        self.energy   = []    # common energy scale for all analyzers
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

    def loadscan(self,scannumbers,scantype='generic'):
        """
        Loads the files belonging to scan No. "scannumber" and puts it into an instance
        of the xrs_scan-class 'scan'. The default scantype is 'generic', later the scans
        will be grouped (and added) based on the scantype.
        INPUT:
        scannumbers = integer or list of scannumbers that should be loaded
        scantype    = string describing the scan to be loaded (e.g. 'edge1' or 'K-edge') 
        fromtofile  = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
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
            monoangle = 1 # have to check for this later ... !!! counters['pmonoa']
            energy    = counters[self.encolumn]
            # create an instance of "scan" class for every scan
            onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
            # assign one dictionary entry to each scan 
            self.scans[scanname] = onescan
    
    def loadscandirect(self,scannumbers,scantype='generic',scaling=None):
        """
        Loads a scan without saving the edf files in matrices.
        scannumbers = integer or list of integers defining the scannumbers from the SPEC file
        scantype    = string describing the scan to be loaded (e.g. 'edge1' or 'K-edge') 
        fromtofile  = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
        scaling     = list of scaling factors to be applied, one for each ROI defined
        """
        # make sure scannumbers are iterable (list)
        if not isinstance(scannumbers,list):
            scannums = []
            scannums.append(scannumbers)
        else:
            scannums = scannumbers 
        # check if there are ROIs defined        
        if not self.roi_obj:
            print( 'Please define some ROIs first')
            return
        for number in scannums:
            scanname = 'Scan%03d' % number
            data, motors, counters, edfmats = self.readscan(number)
            # can assign some things here already (even if maybe redundant)
            monitor   = counters[self.monicolumn]
            monoangle = 1 # counters['pmonoa'] # this still needs checking
            energy    = counters[self.encolumn]
            # create an instance of "scan" class for every scan

            onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)

            onescan.applyrois(self.roi_obj.indices,scaling=scaling)

            print( 'Deleting EDF-files of Scan No. %03d' % number)
            onescan.edfmats = [] # delete the edfmats
            self.scans[scanname] = onescan

    def deletescan(self,scannumbers):
        """
        Deletes scans from the class.
        INPUT:
        scannumbers = integer or list of integers (SPEC scan numbers) to delete
        """
        numbers = []
        if not isinstance(scannumbers,list):
            numbers.append(scannumbers)
        else:
            numbers = scannumbers
        for number in numbers:
            scanname = 'Scan%03d' % number
            del(self.scans[scanname])
            self.scannumbers.remove(number)

    def SumDirect(self,scannumbers):
        Sum=None
        for number in scannumbers:
            data, motors, counters, edfmats = self.readscan(number)
            if Sum is None:
                Sum = np.zeros(edfmats[0].shape ,"f") 
            Sum[:] += edfmats.sum(axis=0) 
        return Sum

    def getXESspectrum(self):
        """
        Groups the instances of the scan class by their scantype attribute, 
        adds equal scans (each group of equal scans) and appends them.
        INPUT:
        include_elastic = boolean flag, skips the elastic line if set to 'False' (default)
        """
        # find the groups 
        allgroups = xrs_scans.findgroups(self.scans)
        for group in allgroups:
            # self.makegroup(group)
            onegroup = xrs_scans.makegroup_nointerp(group)
            self.groups[onegroup.get_type()] = onegroup
        self.energy,self.signals,self.errors = xrs_scans.appendXESScans(self.groups)

    def dump_data(self,filename):
        data = np.zeros((len(self.eloss),3))
        data[:,0] = self.eloss
        data[:,1] = self.signals
        data[:,2] = self.errors
        np.savetxt(filename,data)










