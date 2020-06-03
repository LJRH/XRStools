from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import range
from collections import Iterable



#!/usr/bin/python
# Filename: ixs_offDiagonal.py

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


#from helpers import *
from . import xrs_rois, xrs_scans, xrs_utilities, math_functions, xrs_fileIO, roifinder_and_gui
import h5py
from numpy import array
import scipy.io
import traceback
import sys
import os
import numpy as np
import array as arr
import pickle
from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from pylab import *
from scipy import signal
from scipy.ndimage import measurements
from itertools import groupby

import matplotlib.pyplot as plt
import warnings

# try to import the fast PyMCA parsers
try:
    import PyMca5.PyMcaIO.EdfFile as EdfIO
    import PyMca5.PyMcaIO.specfilewrapper as SpecIO
    use_PyMca = True
except:
    use_PyMca = False

print( " >>>>>>>>  use_PyMca " , use_PyMca)
__metaclass__ = type # new style classes


def print_citation_message():
	"""Prints plea for citing the XRStools article when using this software.

	"""
	print ('                                                                                ')
	print (' ############################# Welcome to XRStools #############################')
	print (' # If you are using this software, please cite the following work:             #')
	print (' # Ch.J. Sahle, A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari: #') 
	print (' # "Planning, performing, and analyzing X-ray Raman scattering experiments."   #')
	print (' # Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.                 #')
	print (' ###############################################################################')
	print ('                                                                                ')


class offDiagonal:
    """ **offDiagonal**
    Class for reading scans from off-diagonal IXS experiments on the high-resolution setup at ID20.

    Arguments:
    ----------
    absFilename (string): Absolute path and filename of the SPEC-file.
    scanMotor (string): Mnemonic of the motor that is scanned.
    monitorName (string): Mnemonic of the counter used for normalization.
    edfName (string): EDF-file base file name (default is None, i.e. same as SPEC-file).
    armLength (float): Legth (in m) of the spectrometer arm used (either 1.0 or 2.0).
    """
    def __init__( self, path, SPECfname='fourc', EDFprefix='/edf/', EDFname='fourc_', \
                        EDFpostfix='.edf', en_column='sry', moni_column='izero' ):

        self.path        = path
        self.SPECfname   = SPECfname

        if not os.path.isfile(os.path.join(path, SPECfname)):
            raise Exception('IOError! No such file or directory.')

        self.EDFprefix   = EDFprefix
        self.EDFname     = EDFname
        self.EDFpostfix  = EDFpostfix
        self.en_column   = en_column.lower()
        self.moni_column = moni_column.lower()

        self.scans        = {}
        self.scan_numbers = []

        self.eloss    = np.array([])
        self.energy   = np.array([])
        self.signals  = np.array([])
        self.errors   = np.array([])
        self.q_values = []
        self.groups   = {}
        self.cenom    = []
        self.E0       = []
        self.tth      = []
        self.resolution  = []
        self.comp_factor = None
        self.cenom_dict  = {}
        self.raw_signals = {}
        self.raw_errors  = {} 

        self.TTH_OFFSETS1 = np.array([5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0])
        self.TTH_OFFSETS2 = np.array([-9.71, -9.75, -9.71, -3.24, -3.25, -3.24, 3.24, 3.25, 3.24, 9.71, 9.75, 9.71]) 
        self.PIXEL_SIZE   = 0.055 # pixel size in mm

        self.roi_obj = None

        # specific to off-diagonal experiments
        self.scanMatrix     = np.array([])
        self.offDiaDataSets = []

        print_citation_message()

    def SumDirect(self, scan_numbers):
        """ **SumDirect**

        Creates a summed 2D image of a given scan or list of scans.

        Args:
            scan_numbers (int or list): Scan number or list of scan numbers to be added up.

        Returns:
            A 2D np.array of the same size as the detector with the summed image.

        """

        # make sure scannumbers are iterable (list)
        numbers = []
        if not isinstance(scan_numbers,list):
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers

        im_sum    = None
        en_column = None # uses first column in SPEC file as scanned motor
        for number in numbers:
            scan = xrs_scans.Scan()
            print(" SONO IN SumDirect   number " , number, " en_column  " , en_column,
                  " moni_column " , self.moni_column, " scan_type  " , None, " scaling  " , None
            ) 
            scan.load(self.path, self.SPECfname, self.EDFprefix, self.EDFname, self.EDFpostfix, number, \
                      direct=False, roi_obj=None, scaling=None, scan_type='generic', \
                      en_column=en_column, moni_column=self.moni_column)

            print( " IN SumDirect  la shape est   ",  scan.edfmats.shape  ) 
            
            if im_sum is None:
                im_sum = np.zeros(scan.edfmats[0].shape ,"f") 
            im_sum[:] += scan.edfmats.sum(axis=0) 
        return im_sum


    def set_roiObj(self,roiobj):
        """ **set_roiObj**
        Assigns an object of the roi_obj class to this class.
        """
        self.roi_obj = roiobj

    def load_scan( self, scan_numbers, scan_type='generic', direct=True, scaling=None, method='sum' ):
        """ **load_scan**

        Load a single or multiple scans.

        Note:
            When composing a spectrum later, scans of the same 'scan_type' will
            be averaged over. Scans with scan type 'elastic' or 'long' in their
            names are recognized and will be treated specially.

        Args:
            scan_numbers (int or list): Integer or iterable of scan numbers to be loaded.
            scan_type            (str): String describing the scan to be loaded (e.g. 'edge1' or 'K-edge').
            direct           (boolean): Flag, 'True' if the EDF-files should be deleted after loading/integrating the scan.

        """

        # make sure scannumbers are iterable (list)
        numbers = []
        if not isinstance(scan_numbers,list):
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers

        # make sure there is a cenom_dict available if direct=True AND method='sum'
        if direct and method=='pixel' and not self.cenom_dict:
            print('Please run the get_compensation_factor method first for pixel-wise compensation.')
            return

        # go throught list of scan_numbers and load scans		
        for number in numbers:

            # create a name for each scan
            scan_name = 'Scan%03d' % number

            # create an instance of the Scan class
            scan = xrs_scans.Scan()

            # load the scan
            scan.load( self.path, self.SPECfname, self.EDFprefix, self.EDFname, self.EDFpostfix, number, \
                direct=direct, roi_obj=self.roi_obj, scaling=scaling, scan_type=scan_type, \
                en_column=self.en_column, moni_column=self.moni_column, method=method, cenom_dict=self.cenom_dict )

            # add it to the scans dict
            self.scans[scan_name] = scan

            # add the number to the list of scan numbers
            if not number in self.scan_numbers:
                self.scan_numbers.extend([number])

    def save_state_hdf5( self,  file_name, group_name, comment="" , overwriteFile = True,  overwriteGroup=True):
        if overwriteFile:
            h5 = h5py.File(file_name,"w")
        else:
            h5 = h5py.File(file_name,"a")
        if overwriteGroup:
            if group_name in h5:
                del       h5[group_name]
        h5group =  h5.require_group(group_name)
        h5group["scan_numbers"] = self.scan_numbers
        for scan_name, scan in self.scans.items():
            h5group_scan  =  h5.require_group(group_name+"/"+scan_name)
            scan.save_hdf5( h5group_scan)
            h5group_scan["offdia_energy"]  = scan.offdia_energy
            h5group_scan["RCmonitor"]      = scan.RCmonitor
                    
        h5group["comment"]  = comment
        h5.flush()
        h5.close()

    def load_state_hdf5( self,  file_name, group_name ):
        h5 = h5py.File(file_name,"r")
        
        h5group =  h5.require_group(group_name)

        self.scan_numbers = h5group["scan_numbers"][()]

        self.scans() = {}

        for key in h5group:
            if str(key)[:4] == "Scan":
                scan_name = str(key)
                
                # for scan_name, scan in self.scans.items():
            
                h5group_scan  =  h5group[h5group_scan]

                scan xrs_scans.Scan()
                
                scan.load_hdf5( h5group_scan) 

                
                scan.offdia_energy = h5group_scan["offdia_energy"]  [()]
            
                scan.RCmonitor  = h5group_scan["RCmonitor"]   [()]   

                self.scans[scan_name] = scan
        h5.flush()
        h5.close()

                
    def loadRockingCurve(self, scan_numbers, en_colum='energy', RCmoni='alirixs', direct=False, scaling=None, method='sum', scan_type='RC', storeInsets = False):
        """ **loadRockingCurve**
        Load one or more rocking curves.

        Arguments:
        ----------
        scanNumbers (int or list of ints): Scan number or list of scan numbers of rocking curve scans to be loaded.
        scanType (string): String describing the scan for later automatic stitching/interpolation. Few special types exist: elastic, long.
        energyCoor (list): Indices to find the energy during the scan based on the SPEC-file header.
        direct (boolean): Keyword if EDF-files should be deleted or kept (default).
        """
        # make sure scanNumbers are iterable (list)
        numbers = []
           
        if not isinstance( scan_numbers , Iterable) :
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers

        for number in numbers:
            # create a name for each scan
            scan_name = 'Scan%03d' % number

            # create an instance of the Scan class
            scan = xrs_scans.Scan()

            # load the scan
            scan.load( self.path, self.SPECfname, self.EDFprefix, self.EDFname, self.EDFpostfix, number, \
                direct=direct, roi_obj=self.roi_obj, scaling=scaling, scan_type=scan_type, \
                       en_column=self.en_column, moni_column=self.moni_column, method=method, cenom_dict=self.cenom_dict , storeInsets = storeInsets)

            scan.offdia_energy = scan.motors[en_colum]
            scan.RCmonitor     = scan.counters[RCmoni]

            self.scans[scan_name] = scan
            if not number in self.scan_numbers:
                self.scan_numbers.extend([number])

    def stitchRockingCurves(self,RCmoni='kaprixs',I0moni='izero',addColumns = 0):
        """ **stitchRockingCurves**
        Go through all rocking curves and stitch them together to a 3D matrix.
        """
        RcScans = xrs_scans.findRCscans(self.scans)

        sorted_RcScans = sorted(RcScans,key=lambda x:x.offdia_energy)
        energy_points  = sorted(list(set([scan.offdia_energy for scan in sorted_RcScans])))

        dim1 = len(energy_points)
        dim2 = int(sum(list(set([scan.signals.shape[0] for scan in sorted_RcScans])))+addColumns) # dirty fix for 2 scans of the same length

        for ii in range(len(self.roi_obj.red_rois)):
            dataset = xrs_scans.offDiaDataSet()
            dataset.ROIno  = ii
            dataset.energy = energy_points
            moniMatrix   = np.zeros((dim1,dim2))
            motorMatrix  = np.zeros((dim1,dim2))
            signalMatrix = np.zeros((dim1,dim2))
            I0Matrix     = np.zeros((dim1,dim2))
            errorMatrix  = np.zeros((dim1,dim2))
            for jj in range(len(energy_points)):
                moniCol   = np.array([])
                motorCol  = np.array([])
                signalCol = np.array([])
                I0Col     = np.array([])<
                for scan in sorted_RcScans:
                    if scan.offdia_energy == energy_points[jj]:
                        moniCol   = np.insert(moniCol,np.searchsorted(moniCol,scan.counters[RCmoni]),scan.counters[RCmoni])
                        motorCol  = np.insert(motorCol,np.searchsorted(motorCol,scan.energy),scan.energy)
                        signalCol = np.insert(signalCol,np.searchsorted(signalCol,scan.signals[:,ii]),scan.signals[:,ii])
                        I0Col     = np.insert(I0Col,np.searchsorted(I0Col,scan.counters[I0moni]),scan.counters[I0moni])
                moniMatrix[jj,:]   = moniCol
                signalMatrix[jj,:] = signalCol
                motorMatrix[jj,:]  = motorCol
                I0Matrix[jj,:]     = I0Col
                errorMatrix[jj,:]  = np.sqrt(signalCol)

            dataset.signalMatrix = signalMatrix
            dataset.motorMatrix  = motorMatrix
            dataset.RCmonitor    = moniMatrix
            dataset.I0Matrix     = I0Matrix
            dataset.errorMatrix  = errorMatrix
            self.offDiaDataSets.append(dataset)

    def getrawdata(self):
        """ **getrawdata**
        Iterates through all instances of the scan class and calls it's applyrois method
        to sum up over all rois.
        """
        if not np.any(self.roi_obj.indices):
            print( 'Please define some ROIs first.')
            return
        for scan in self.scans:
            if len(self.scans[scan].edfmats):
                print ("integrating "+scan)
                self.scans[scan].applyrois(self.roi_obj.indices)

    def getrawdata_pixelwise(self):
        """
        Goes through all instances of the scan class and calls it's applyrois_pw method
        to extract intensities for all pixels in each ROI.
        """
        if not np.any(self.roi_obj.indices):
            print( 'Please define some ROIs first.')
            return
        for scan in self.scans:
            if len(self.scans[scan].edfmats):
                print ("integrating pixelwise "+scan)
                self.scans[scan].applyrois_pw(self.roi_obj.indices)

        def SumDirect(self,scannumbers):
            """ **SumDirect**
            """
            sum = None
            for number in scannumbers:
                data, motors, counters, edfmats = self.readscan(number)
                if sum is None:
                    sum = np.zeros(edfmats[0].shape ,"f") 
                    sum[:] += edfmats.sum(axis=0)
            return sum

    def deletescan(self,scannumbers):
        """ **deletescan**
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

    def save_raw_data(self,filename):
        data = np.zeros((len(self.eloss),len(self.signals[0,:])))
        data[:,0]   = self.eloss
        data[:,1::] = self.signals
        np.savetxt(filename,data)








