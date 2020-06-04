from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import range
from six.moves import zip
from six.moves import input
#!/usr/bin/python
# Filename: xrs_read.py

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
__author__    = "Christoph J. Sahle - ESRF"
__contact__   = "christoph.sahle@esrf.fr"
__license__   = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

from . import xrs_rois, xrs_scans, xrs_utilities, math_functions, xrs_fileIO, roifinder_and_gui

import h5py
import scipy.io
import traceback
import sys
import os, re
import numpy as np
import array as arr
import pickle
import matplotlib.pyplot as plt

from numpy import array
from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d, Rbf
from scipy import signal, optimize
from scipy.ndimage import measurements

# try to import the fast PyMCA parsers
try:
    import PyMca5.PyMcaIO.EdfFile as EdfIO
    import PyMca5.PyMcaIO.specfilewrapper as SpecIO
    use_PyMca = True
except:
    use_PyMca = False

print( " >>>>>>>>  use_PyMca " , use_PyMca)
__metaclass__ = type # new style classes

# These values are used in read_Lerix class but may be useful elsewhere? LJRH
TINY = 1.e-7
MAX_FILESIZE = 100*1024*1024  # 100 Mb limit
COMMENTCHARS = '#;%*!$'
NAME_MATCH = re.compile(r"[a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)*$").match
VALID_SNAME_CHARS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'
VALID_NAME_CHARS = '.%s' % VALID_SNAME_CHARS
RESERVED_WORDS = ('and', 'as', 'assert', 'break', 'class', 'continue',
                'def', 'del', 'elif', 'else', 'eval', 'except', 'exec',
                'execfile', 'finally', 'for', 'from', 'global', 'if',
                'import', 'in', 'is', 'lambda', 'not', 'or', 'pass',
                'print', 'raise', 'return', 'try', 'while', 'with',
                'group', 'end', 'endwhile', 'endif', 'endfor', 'endtry',
                'enddef', 'True', 'False', 'None')

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

print_citation_message()

class Hydra:
    """Main class for handling XRS data from ID20's multi-analyzer spectrometer 'Hydra'.

    This class is intended to read SPEC- and according EDF-files and generate spectra from
    multiple individual energy loss scans.

    Note:
        Hydra is the name of the multi-analyzer x-ray Raman scattering spectrometer at ESRF's
        ID20 beamline. This class has been adopted specifically for this spectrometer.

        If you are using this program, please cite the following work:

        Sahle, Ch J., A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari.
        "Planning, performing and analyzing X-ray Raman scattering experiments."
        Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.

    Args:
          * path        (str): Absolute path to directory holding the data.
          * SPECfname   (str): Name of the SPEC-file ('hydra' is the default).
          * EDFprefix   (str): Prefix for the EDF-files ('/edf/' is the default).
          * EDFname     (str): Filename of the EDF-files ('hydra_' is the default).
          * EDFpostfix  (str): Postfix for the EDF-files ('.edf' is the default).
          * en_column   (str): Counter mnemonic for the energy motor ('energy' is the default).
          * moni_column (str): Mnemonic for the monitor counter ('izero' is the default).

    Attributes:
          * path        (str): Absolute path to directory holding the data.
          * SPECfname   (str): Name of the SPEC-file ('hydra' is the default).
          * EDFprefix   (str): Prefix for the EDF-files ('/edf/' is the default).
          * EDFname     (str): Filename of the EDF-files ('hydra_' is the default).
          * EDFpostfix  (str): Postfix for the EDF-files ('.edf' is the default).
          * en_column   (str): Counter mnemonic for the energy motor ('energy' is the default).
          * moni_column (str): Mnemonic for the monitor counter ('izero' is the default).
          * scans        (dict): Dictionary holding all loaded scans.
          * scan_numbers (list): List with scan number of all loaded scans.
          * eloss    (np.array): Array with the common energy loss scale for all analyzers.
          * energy   (np.array): Array with the common energy scale for all analyzers.
          * signals  (np.array): Array with the signals for all analyzers (one column per anayzer).
          * errors   (np.array): Array with the poisson errors for all analyzers (one column per anayzer).
          * qvalues      (list): List of momentum transfer values for all analyzers.
          * groups       (dict): Dictionary of groups of scans (instances of the 'scangroup' class,
            such as 2 'elastic', or 5 'edge1', etc.).
          * cenom        (list): List of center of masses of the elastic lines.
          * E0            (int): Elastic line energy value, mean value of all center of masses.
          * tth          (list): List of all scattering angles (one value for each ROI).
          * resolution   (list): List of FWHM of the elastic lines (one for each analyzer).
          * comp_factor (float): Compensation factor for line-by-line energy dispersion compensation.
          * cenom_dict   (dict): Dictionary holding center-of-masses for of the elastic line.
          * raw_signals  (dict): Dictionary holding pixel- or line-wise signals.
          * raw_errors   (dict): Dictionary holding pixel- or line-wise Poisson errors.
          * TTH_OFFSETS1 (np.array): Two-Theta offsets between individual analyzers inside each analyzer
            module in one direction (horizontal for V-boxes, vertical for H-boxes).
          * TTH_OFFSETS2 (np.array): Two-Theta offsets between individual analyzers inside each analyzer
            module in one direction (horizontal for H-boxes, vertical for V-boxes).
          * roi_obj (instance): Instance of the roi_object class from the xrs_rois module defining all
            ROIs for the current dataset (default is 'None').

    """

    def __init__( self, path, SPECfname='hydra', EDFprefix='/edf/', EDFname='hydra_', \
                        EDFpostfix='.edf', en_column='energy', moni_column='izero' ):

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

        print_citation_message()

    def save_state_hdf5( self, file_name, group_name, comment="" ):
        """ **save_state_hdf5**

        Save the status of the current instance in an HDF5 file.

        Args:
            file_name  (str): Path and file name for the HDF5-file to be created.
            group_name (str): Group name under which to store status in the HDF5-file.
            comment    (str): Optional comment (no comment is default).

        """
        h5 = h5py.File(file_name,"a")

        h5.require_group(group_name)
        h5group =  h5[group_name]

        for key in self.__dict__.keys():
            if key in h5group.keys():
                raise Exception( 'Data \'' + key + '\' already present in  ' + file_name + ':' + group_name )
            else:
                h5group[key] = getattr( self, key )

        h5group["comment"]  = comment
        h5.flush()
        h5.close()

    def load_state_hdf5( self, file_name, group_name ):
        """ **load_state_hdf5**

        Load the status of an instance from an HDF5 file.

        Args:
            file_name  (str): Path and filename for the HDF5-file to be created.
            group_name (str): Group name under which to store status in the HDF5-file.

        """
        h5 = h5py.File( file_name,"r" )

        h5group =  h5[group_name]
        keys    = {"eloss":array, "energy":array, "signals":array, "errors":array, "q_values":array,
                    "cenom":array, "E0":float, "tth":array, "resolution":array }

        for key in keys:
            setattr(self, key, keys[key](array(h5group[key])))

        h5.flush()
        h5.close()

    def set_roiObj( self,roiobj ):
        """ **set_roiObj**

        Assign an instance of the 'roi_object' class to the current data set.

        Args:
            roiobj (instance): Instance of the 'roi_object' class holding all
                information about the definition of the ROIs.

        """
        self.roi_obj = roiobj

    def load_scan( self, scan_numbers, scan_type='generic', direct=True, scaling=None, method='sum' ):
        """**load_scan**

        Load a single or multiple scans.
        Note:
            When composing a spectrum later, scans of the same scan_type will
            be averaged over. Scans with scan type 'elastic' or long in their
            names are recognized and will be treated specially.

        Args:
              * scan_numbers (int or list): Integer or iterable of scan numbers to be loaded.
              * scan_type            (str): String describing the scan to be loaded (e.g. 'edge1' or 'K-edge').
              * direct           (boolean): Flag, 'True' if the EDF-files should be deleted after loading/integrating the scan.


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
                en_column=self.en_column, moni_column=self.moni_column, method=method, \
                cenom_dict=self.cenom_dict, comp_factor=self.comp_factor )

            # add it to the scans dict
            self.scans[scan_name] = scan

            # add the number to the list of scan numbers
            if not number in self.scan_numbers:
                self.scan_numbers.extend([number])

    def load_loop( self, beg_nums, num_of_regions, direct=True, method='sum' ):
        """ **load_loop**

        Loads a whole loop of scans based on their starting numbers and
        the number of single scans in the loop.

        Args:
        beg_nums      (list): List of scan numbers of the first scans in each loop.
        num_of_regions (int): Number of scans in each loop.

        """

        type_names = []
        for n in range(num_of_regions):
            type_names.append('edge'+str(n+1))

        numbers = []
        for n in range(len(beg_nums)):
            for m in range(num_of_regions):
                numbers.append((beg_nums[n]+m))

        type_names = type_names*len(beg_nums)

        for n in range(len(type_names)):
            number = []
            number.append(numbers[n])
            self.load_scan( number, type_names[n], direct=True, method=method )

    def delete_scan( self, scan_numbers ):
        """ **delete_scan**

        Deletes scans from the dictionary of scans.

        Args:
            scan_numbers (int or list): Integer or list of integers (SPEC scan numbers) to be deleted.

        """

        # make sure scannumbers are iterable (list)
        numbers = []
        if not isinstance(scan_numbers,list):
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers

        # delete the scan
        for number in numbers:
            scan_name = 'Scan%03d' % number
            del(self.scans[scan_name])
            self.scan_numbers.remove(number)

    def get_raw_data( self, method='sum', scaling=None ):
        """ **get_raw_data**

        Applies the ROIs to the EDF-files.

        This extracts the raw-data from the EDF-files subject to three
        different methods: 'sum' will simply sum up all pixels inside a
        ROI, 'row' will sum up over the dispersive direction, and 'pixel'
        will not sum at all.

        Args:
            method (str): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', or 'pixel'. Default is 'sum'.

        """

        if not self.roi_obj:
            print ( 'Did not find a ROI object, please set one first.' )
            return

        for scan in self.scans:
            if len(self.scans[scan].edfmats):
                print ( "Integrating " + scan + " using method \'" + method + '\'.')
                self.scans[scan].get_raw_signals( self.roi_obj , method=method, scaling=scaling )


    def get_data_new(self, method='sum', scaling=None):
        """ **get_data_new**

        Applies the ROIs to the EDF-files.

        This extracts the raw-data from the EDF-files subject to three
        different methods: 'sum' will simply sum up all pixels inside a
        ROI, 'row' will sum up over the dispersive direction, and 'pixel'
        will not sum at all.

        Args:
            method (str): Keyword specifying the selected choice of data treatment:
                can be 'sum', 'row', or 'pixel'. Default is 'sum'.

        """

        if not self.cenom_dict:
            print ( 'No compensation factor/center-of-mass info found, please provide first.' )
            return

        for scan in self.scans:
            self.scans[scan].get_signals( method=method, cenom_dict=self.cenom_dict, comp_factor=self.comp_factor, scaling=scaling )

    def get_data(self):
        """ **get_data**

        Applies the ROIs and sums up intensities.

        Returns:
            'None', if no ROI object is available.

        """

        if not self.roi_obj:
            print ( 'Did not find a ROI object, please set one first.' )
            return

        for scan in self.scans:
            if len(self.scans[scan].edfmats):
                print ( "Integrating " + scan )
                self.scans[scan].apply_rois( self.roi_obj )

    def get_data_pw(self):
        """ **get_data_pw**

        Extracts intensities for each pixel in each ROI.

        Returns:
            'None', if no ROI object is available.

        """

        if not np.any(self.roi_obj.indices):
            print ( 'Did not find a ROI object, please set one first.' )
            return

        for scan in self.scans:
            if len(self.scans[scan].edfmats):
                print ( "Integrating pixelwise " + scan )
                self.scans[scan].apply_rois_pw(self.roi_obj)

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
            scan.load(self.path, self.SPECfname, self.EDFprefix, self.EDFname, self.EDFpostfix, number, \
                    direct=False, roi_obj=None, scaling=None, scan_type='generic', \
                    en_column=en_column, moni_column=self.moni_column)

            if im_sum is None:
                im_sum = np.zeros(scan.edfmats[0].shape ,"f")
            im_sum[:] += scan.edfmats.sum(axis=0)
        return im_sum

    def get_eloss_new(self, method='sum'):
        """ **get_eloss_new**

        Defines the energy loss scale for all ROIs and applies dispersion
        compensation if applicable.

        Args:
            method (str): Keyword describing which dispersion compensation
                method to use. Possible choices are 'sum' (no compensation),
                'pixel' (pixel-by-pixel compensation), or 'row' (line-by-line)
                compensation.

        """

        # make sure an elastic line was loaded
        if not 'elastic' in self.groups:
            print( 'Please load/integrate at least one elastic scan first!' )
            return

        # make sure there is data
        elif not np.any(self.groups['elastic'].raw_signals):
            self.get_raw_data( method=method )

        # make sure there is center of masses available
        elif not self.cenom_dict and  method != 'row':
            for scan in self.scans:
                if self.scans[scan].get_type() == 'elastic':
                    print('GETIING COMP FACTOR!!!')
                    self.get_compensation_factor( self.scans[scan].scan_number, method=method )

        first_key = list(self.raw_signals.keys())[0]
        # 'sum'
        if method == 'sum':
            # master eloss scale in eV is the one of the first ROI
            self.signals = np.zeros((len(self.energy),len(self.cenom_dict)))
            self.errors  = np.zeros((len(self.energy),len(self.cenom_dict)))
            master_eloss = (self.energy - np.median([self.cenom_dict[key] for key in  self.cenom_dict]))*1.0e3
            self.E0      = np.median([self.cenom_dict[key] for key in  self.cenom_dict])
            for key,ii in zip(sorted(self.cenom_dict), range(len(self.cenom_dict))):
                # signals
                x = ( self.energy - self.cenom_dict[key] )*1.0e3
                y = self.raw_signals[key][:]
                #try:
                #	rbfi = Rbf( x, y, function='linear' )
                #	self.signals[:,ii] = rbfi( master_eloss )
                #except:
                self.signals[:,ii] = np.interp( master_eloss, x, y )
                # errors
                y = self.raw_errors[key][:]
                #try:
                #	rbfi = Rbf( x, y, function='linear' )
                #	self.errors[:,ii] = rbfi( master_eloss )
                #except:
                self.errors[:,ii] = np.interp( master_eloss, x, y )
            self.eloss = master_eloss

        # 'pixel'
        elif method ==  'pixel':
            # master eloss scale in eV is the one of central pixel in first ROI
            self.signals = np.zeros((len(self.energy),len(self.cenom_dict)))
            self.errors  = np.zeros((len(self.energy),len(self.cenom_dict)))
            master_eloss = ( self.energy - np.median([np.median(self.cenom_dict[key]) for key in self.cenom_dict]) )*1.0e3
            self.E0      = np.mean(self.cenom_dict[first_key][self.cenom_dict[first_key] > 0.0])
            for key,ii in zip(sorted(self.cenom_dict), range(len(self.cenom_dict))):
                print ('Pixel-by-pixel compensation for ' + key +'.')
                signal = np.zeros(len(master_eloss))
                error  = np.zeros(len(master_eloss))
                for dim1 in range(self.cenom_dict[key].shape[0]):
                    for dim2 in range(self.cenom_dict[key].shape[1]):
                        x = ( self.energy - self.cenom_dict[key][dim1, dim2] )*1.0e3
                        # signals
                        y = self.raw_signals[key][:, dim1, dim2]
                        #print "Y AMAX", signal.max()
                        #rbfi    = Rbf( x, y, function='linear' )
                        rbfi = interp1d(x, y,bounds_error=False, fill_value=0.0)
                        #print "rbf AMAX", Rbf( x, y, function='linear' )
                        #signal += rbfi( master_eloss )
                        signal += rbfi( master_eloss )
                        #print "SIGNAL AMAX", signal.max()
                        # errors
                        y = self.raw_errors[key][:, dim1, dim2]
                        rbfi    = Rbf( x, y, function='linear' )
                        error  += rbfi( master_eloss )**2
                self.signals[:,ii] = signal
                self.errors[:,ii]  = np.sqrt(error)
            self.eloss = master_eloss

        # 'row'
        elif method == 'row':
            self.signals = np.zeros((len(self.energy),len(self.roi_obj.red_rois)))
            self.errors  = np.zeros((len(self.energy),len(self.roi_obj.red_rois)))
            energy = self.energy * 1e3 # energy in eV

            for key,ii in zip(sorted(self.raw_signals), range(len(self.raw_signals))):
                y  = self.raw_signals[key]
                meanii = len(range(y.shape[1]))/2
                yc = np.zeros_like(y)
                for jj in range(y.shape[1]):
                    sort  = np.argsort(energy)
                    yc[:,jj] = np.interp( energy[sort] + (jj-meanii)*self.comp_factor*self.PIXEL_SIZE, energy[sort], y[sort,jj],
                                        left=float('nan'), right=float('nan') )
                self.signals[:,ii] = xrs_utilities.nonzeroavg(yc)
                self.errors[:,ii]  = np.sqrt(self.signals[:,ii])

        else:
            print('Method \''+method+'\' not supported, use either \'sum\', \'pixel\', or \'row\'.')
            return

    def get_eloss( self ):
        """ **get_eloss**

        Finds the energy loss scale for all ROIs by calculating the center
        of mass (COM) for each ROI's elastic line. Calculates the resolution
        function (FWHM) of the elastic lines.

        """

        # make sure an elastic line was loaded
        if not 'elastic' in self.groups:
            print( 'Please load/integrate at least one elastic scan first!' )
            return

        # make sure there is data
        elif not self.groups['elastic'].signals.shape[1] == len(self.roi_obj.indices):
            self.get_data()

        else:
            # reset values, in case this is run several times
            self.cenom      = []
            self.resolution = []
            Origin          = None
            valid_cenoms    = []

            for n in range(self.groups['elastic'].signals.shape[1]):
                # find the center of mass for each ROI
                cofm = xrs_utilities.find_center_of_mass(self.groups['elastic'].energy,self.groups['elastic'].signals[:,n])
                self.cenom.append(cofm)
                if self.there_is_a_valid_roi_at( n ):
                    valid_cenoms.append(cofm)
                    if Origin is None:
                        Origin = cofm

                # find the FWHM/resolution for each ROI
                en_scale  = (self.groups['elastic'].energy - self.cenom[n])*1e3
                intensity = self.groups['elastic'].signals_orig[:,n]
                FWHM, x0 = xrs_utilities.fwhm(en_scale,intensity)

            # find master E0
            self.E0 = np.mean(valid_cenoms)

            # define last eloss scale as the 'master' scale for all ROIs
            self.eloss = (self.energy - cofm)*1e3 # energy loss in eV

            # define eloss-scale for each ROI and interpolate onto the 'master' eloss-scale
            for n in range(self.signals.shape[1]):
                # inserting zeros at beginning and end of the vectors to avoid interpolation errors
                x = (self.energy-self.cenom[n])*1e3
                x = np.insert(x,0,-1e10)
                x = np.insert(x,-1,1e10)
                y = self.signals[:,n]
                y = np.insert(y,0,0)
                y = np.insert(y,-1,0)
                f = interp1d(x,y, bounds_error=False,fill_value=0.0)
                self.signals[:,n] = f(self.eloss)

            # do the same for the errors
            for n in range(self.signals.shape[1]):
                # inserting zeros at beginning and end of the vectors to avoid interpolation errors
                x = (self.energy-self.cenom[n])*1e3
                x = np.insert(x,0,-1e10)
                x = np.insert(x,-1,1e10)
                y = self.errors[:,n]
                y = np.insert(y,0,0)
                y = np.insert(y,-1,0)
                f = interp1d(x,y, bounds_error=False,fill_value=0.0)
                self.errors[:,n] = f(self.eloss)

    def get_spectrum( self, include_elastic=False, abs_counts=False ):
        """ **get_spectrum**

        Constructs a spectrum based on the scans loaded so far. Defines the energy
        loss scale based on the elastic lines.

        Args:
            include_elastic (boolean): Boolean flag, does not include the elastic line if
                set to 'False' (this is the default).
            abs_counts      (boolean): Boolean flag, constructs the spectrum in absolute
                counts if set to 'True' (default is 'False')

        """

        # find the groups
        all_groups = xrs_scans.findgroups(self.scans)

        # append scans
        for group in all_groups:
            one_group = xrs_scans.makegroup_nointerp(group, absCounts=abs_counts)
            self.groups[one_group.get_type()] = one_group

        self.energy, self.signals, self.errors = xrs_scans.appendScans(self.groups,include_elastic)

        # define the energy loss scale
        self.get_eloss()

    def get_spectrum_new( self, method='sum', include_elastic=False, abs_counts=False, interpolation=False ):
        """ **get_spectrum_new**

        Constructs a spectrum from the scans loaded so far based on the chosen method:
        detector pixels can either be summed up (method='sum'), pixel-by-pixel
        compensation can be applied (method='pixel'), or a line-by-line compensation
        scheme can be applied (method='row').

        The energy loss scale will be defined in the process and the data will be
        interpolated onto this scale using norm-conserving wavelet interpolation.

        Args:
            method              (str): Keyword describing the kind of integration scheme
                to be used (possible values are 'sum', 'pixel', or 'row'), default is 'sum'.
            include_elastic (boolean): Boolean flag, does not include the elastic line if
                set to 'False' (this is the default).
            abs_counts      (boolean): Boolean flag, constructs the spectrum in absolute
                counts if set to 'True' (default is 'False')
            interpolation   (boolean): Boolean flag, if True, signals are interpolated
                onto energy grid of the first scan in each group of scans.

        """

        if not method in ['pixel', 'sum', 'row']:
            print('Unknown integration method. Use either \'pixel\', \'sum\', or \'row\'')
            return

        # make sure there is an elastic line available
        elastic_number = None
        for key in self.scans:
            if self.scans[key].scan_type == 'elastic':
                elastic_number = self.scans[key].scan_number
                print ('Using scan No. %d for CENOMs.'%elastic_number)
                break
            else:
                pass

        if not elastic_number:
            print( 'Please load/integrate at least one elastic scan first!' )
            return

        # get compensation factor/CENOM for an elastic line
        if method in [ 'sum' , 'pixel']  and not self.cenom_dict:
            self.get_compensation_factor( elastic_number, method=method )

        if method == 'row' and not self.comp_factor:

            self.get_compensation_factor( elastic_number, method=method )

        # get raw data
        for key in self.scans:
            if not self.scans[key].raw_signals:
                self.scans[key].get_raw_signals( self.roi_obj, method=method )

        # sum up similar scans
        # find all groups of scans
        all_groups = xrs_scans.findgroups(self.scans)

        # initiate groups
        self.groups = {}

        # sum up similar scans
        for group in all_groups:
            self.groups[group[0].get_type()] = xrs_scans.sum_scans_to_group( group, method=method, interp=interpolation )

        # stitch groups together into a spectrum
        spectrum = xrs_scans.stitch_groups_to_spectrum( self.groups, method=method, include_elastic=include_elastic )
        self.energy      = spectrum.energy
        self.signals     = spectrum.signals
        self.errors      = spectrum.errors
        self.raw_signals = spectrum.raw_signals
        self.raw_errors  = spectrum.raw_errors

        # define energy loss scale and apply compensation if applicable
        self.get_eloss_new( method=method )

    def get_q_values( self, inv_angstr=False, energy_loss=None ):
        """ **get_q_values**

        Calculates the momentum transfer for each analyzer.

        Args:
            inv_angstr (boolean): Boolean flag, if 'True' momentum transfers are calculated in
                inverse Angstroms.
            energy_loss  (float): Energy loss value at which the momentum transfer is to be
                calculated. If 'None' is given, the momentum transfer is calculated for every
                energy loss point of the spectrum.

        Returns:
            If an energy loss value is passed, the function returns the momentum transfers
            at this energy loss value for each analyzer crystal.

        """

        # one q-value per analyzer and energy loss point
        q_vals = np.zeros_like(self.signals)
        if inv_angstr:
            for n in range( self.signals.shape[1] ):
                q_vals[:,n] = xrs_utilities.momtrans_inva(self.E0+self.eloss/1e3,self.E0,self.tth[n])
        else:
            for n in range( self.signals.shape[1] ):
                q_vals[:,n] = xrs_utilities.momtrans_au(self.E0+self.eloss/1e3,self.E0,self.tth[n])
        self.q_values = q_vals

        # return all q-values if a specific energy loss is given
        if energy_loss:
            ind = np.abs(self.eloss - energy_loss).argmin()
            return self.q_values[ind,:]

    def copy_edf_files( self, scan_numbers, dest_dir ):
        """ **copy_edf_files**

        Copies all EDF-files from given scan_numbers into given directory.

        Args:
          * scan_numbers (int or list) = Integer or list of integers defining the scan
            numbers of scans to be copied.
          * dest_dir             (str) = String with absolute path for the destination.

        """

        import shutil

        # make scan_numbers iterable
        numbers = []
        if not isinstance(scan_numbers,list):
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers
        fname = self.path + self.SPECfname

        # go through the scans, find the EDF-files and copy them
        for nscan in numbers:
            if use_PyMca:
                data, motors, counters = xrs_fileIO.PyMcaSpecRead(fname,nscan)
            else:
                data, motors, counters = xrs_fileIO.SpecRead(fname,nscan)

            for m in range(len(counters['ccdno'])):
                ccdnumber = counters['ccdno'][m]
                edfname   = self.path + self.EDFprefix + self.EDFname + "%04d" % ccdnumber + self.EDFpostfix
                shutil.copy2(edfname, dest_dir)

    def dump_spectrum_ascii( self, file_name, header='' ):
        """ **dump_spectrum_ascii**

        Stores the energy loss and signals in a txt-file.

        Args:
            filename (str): Path and filename to the file to be written.

        """

        data = np.zeros((len(self.eloss),self.signals.shape[1]*2+1))
        data[:,0] = self.eloss
        col = 1
        for ii in range(self.signals.shape[1]):
            data[:,col] = self.signals[:,ii]
            col += 1
            data[:,col] = self.errors[:,ii]
            col += 1

        np.savetxt( file_name, data, header=header )

    def dump_spectrum_hdf5( self, file_name, group_name, comment='' ):
        """ **dump_spectrum_hdf5**

        Writes the summed spectrum into an HDF5 file.

        Args:
            file_name  (str): Path and file name for the HDF5-file to be created.
            group_name (str): Group name under which to store status in the HDF5-file.
            comment    (str): Optional comment (no comment is default).

        """

        h5 = h5py.File(file_name,"a")
        h5.require_group(group_name)
        h5group =  h5[group_name]

        keys = ['energy', 'eloss', 'signals', 'errors']
        for key in keys:
            h5group[key] = getattr( self, key )

        h5group["comment"]  = comment
        h5.flush()
        h5.close()

    def dump_scans_ascii( self, scan_numbers, pre_fix, f_name, post_fix='.dat', header='' ):
        """ **dump_scans_ascii**

        Produce ASCII-type files with columns of energy, signal, and Poisson error.

        Args:
            scan_numbers (int or list): SPEC scan numbers of scans to be safed in ASCII format.
            pre_fix  (str): Path to directory where files should be written into.
            f_name   (str): Base name for the files to be written.
            post_fix (str): Extention for the files to be written (default is '.dat').

        """

        # make sure scan_numbers are iterable
        numbers = []
        if not isinstance(scan_numbers,list):
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers

        for number in numbers:
            scan_name = 'Scan%03d' % number
            if not scan_name in self.scans.keys():
                print ('Scan No. %03d is currently not loaded.'% number)
                return
            x = self.scans[scan_name].energy
            y = self.scans[scan_name].signals
            z = self.scans[scan_name].errors
            data = np.zeros((x.shape[0], y.shape[1]+z.shape[1]+1 ))
            data[:,0] = x
            data[:,1:1+y.shape[1]] = y
            data[:,1+y.shape[1]:1+y.shape[1]+z.shape[1]] = z
            file_name = pre_fix + f_name + '_' + scan_name + post_fix
            np.savetxt(file_name, data, header=header)

    def print_scan_length( self,scan_numbers ):
        """ **print_scan_length**

        Print out the numper of points in given scans.

        Args:
            scan_numbers (int or list): Scan number or list of scan numbers.

        """

        # make scan_numbers iterable
        numbers = []
        if not isinstance(scan_numbers,list):
            numbers.append(scan_numbers)
        else:
            numbers = scan_numbers

        # print out the number of points for all scan_numbers
        for ii in numbers:
            name = 'Scan%03d' % ii
            print( 'Length of scan %03d ' %ii + ' is ' + str(len(self.scans[name].energy)) + '.')

    def get_tths( self, rvd=None, rvu=None, rvb=None, rhl=None, rhr=None, rhb=None, order=[0,1,2,3,4,5] ):
        """ **get_tths**

        Calculates the scattering angles for all analyzer crystals based on
        the mean angle of the analyzer modules.

        Args:
              * rhl  (float): Mean scattering angle of the HL module (default is 0.0).
              * rhr  (float): Mean scattering angle of the HR module (default is 0.0).
              * rhb  (float): Mean scattering angle of the HB module (default is 0.0).
              * rvd  (float): Mean scattering angle of the VD module (default is 0.0).
              * rvu  (float): Mean scattering angle of the VU module (default is 0.0).
              * rvb  (float): Mean scattering angle of the VB module (default is 0.0).
              * order (list): List of integers (0-5) that describe the order of modules in which the
                ROIs were defined (default is VD, VU, VB, HR, HL, HB; i.e. [0,1,2,3,4,5]).

        """
        # reset all tth values
        self.tth   = []

        # try to grab from positions from SPEC-file (first scan in list)
        scan_name = list(self.scans.keys())[0]
        if not rvd:
            rvd = self.scans[scan_name].motors['RVD']
        if not rvu:
            rvu = self.scans[scan_name].motors['RVU']
        if not rvb:
            rvb = self.scans[scan_name].motors['RVB']
        if not rhr:
            rhr = self.scans[scan_name].motors['RHR']
        if not rhl:
            rhl = self.scans[scan_name].motors['RHL']
        if not rhb:
            rhb = self.scans[scan_name].motors['RHB']

        # horizontal modules
        # HL (motor name rhl)
        v_angles = self.TTH_OFFSETS1
        h_angles = self.TTH_OFFSETS2 + rhl
        HLtth    = []
        for n in range(len(h_angles)):
            HLtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)

        # HR (motor name rhr)
        v_angles = self.TTH_OFFSETS1
        h_angles = self.TTH_OFFSETS2 + rhr
        HRtth    = []
        for n in range(len(h_angles)):
            HRtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)

        # HB (motor name rhb)
        v_angles = self.TTH_OFFSETS1
        h_angles = self.TTH_OFFSETS2 + rhb
        HBtth = []
        for n in range(len(h_angles)):
            HBtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)

        # vertical modules
        # VD
        v_angles = self.TTH_OFFSETS2 + rvd
        h_angles = self.TTH_OFFSETS1
        VDtth    = []
        for n in range(len(h_angles)):
            VDtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)

        # VU
        v_angles = self.TTH_OFFSETS2 + rvu
        h_angles = self.TTH_OFFSETS1
        VUtth    = []
        for n in range(len(h_angles)):
            VUtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)

        # VB
        v_angles = self.TTH_OFFSETS2 + rvb
        h_angles = self.TTH_OFFSETS1
        VBtth    = []
        for n in range(len(h_angles)):
            VBtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)

        # list of TTH values
        tth = [VDtth, VUtth, VBtth, HRtth, HLtth, HBtth]

        # list all TTH values in one long list ordered by the 'order'-keyword
        for n in order:
            self.tth.extend(tth[n])

    def there_is_a_valid_roi_at( self,n ):
        """ **there_is_a_valid_roi_at**

        Checks if n is a valid ROI index.

        Args:
            n (int): Index to be checked.

        Returns:
            True, if n is a valid ROI index.

        """

        return n<len(self.roi_obj.indices) and len(self.roi_obj.indices[n])

    def get_pw_matrices( self, scan_numbers, method='pixel' ):
        """	**get_pw_matrices**

        Sums scans from pixelwise ROI integration for use in the pixel-wise ROI
        refinement.

        Args:
              * scan_numbers (int, list): Integer or list of scan numbers to be added.
              * method             (str): Keyword describing how to return the data
                (possible values are: 'pixel' for pixel-wise integration, 'column'
                for column-wise summation, and 'row' for row-wise summation.

        Returns:
            * Raw data in pixel-wise format.

        """
        # make scan_numbers iterable
        if isinstance(scan_numbers,list):
            scannums = scan_numbers
        elif isinstance(scan_numbers,int):
            scannums = [scan_numbers]
        else:
            print( 'Please provide keyword \'scan_numbers\' as integer or list of integers.' )
            return

        # make sure method is known
        if not method in ['pixel', 'column', 'row']:
            print('Unknown integration method. Use either \'pixel\', \'column\', or \'row\'')
            return

        # make sure scan exists
        for scannum in scannums:
            if not scannum in self.scan_numbers:
                self.load_scan(scannum, direct=True, method=method)

        # make sure raw_signals is existing
        if not self.scans['Scan%03d' % scannums[0]].raw_signals:
            self.get_raw_data( method=method )

        # one scan only
        if len(scannums)==1:
            scanname = 'Scan%03d' % scannums[0]
            raw_signals = self.scans[scanname].raw_signals # dict with raw_signals
            monitor     = self.scans[scanname].monitor

            # normalize data
            pw_matrices_norm = []
            for key in sorted(raw_signals):
                if method == 'pixel':
                    unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]*raw_signals[key].shape[2]))
                elif method in [ 'column' , 'row']:
                    unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]))
                else:
                    print('Method \''+method+'\' not supported, use either \'pixel\', \'column\', or \'row\'.')
                    return
                for ii in range(raw_signals[key].shape[0]):
                    if method == 'pixel':
                        unrav_mat[ii,:] = raw_signals[key][ii,:,:].ravel() / monitor[ii]
                    else:
                        unrav_mat[ii,:] = raw_signals[key][ii,:] / monitor[ii]
                pw_matrices_norm.append(unrav_mat)

        # multiple scans to be added
        else:
            # start with first scan
            scanname    = 'Scan%03d' % scannums[0]
            raw_signals = self.scans[scanname].raw_signals # dict with raw_signals
            monitor     = self.scans[scanname].monitor

            # normalize data (first scan)
            pw_matrices_norm = []
            for key in sorted(raw_signals):
                if method == 'pixel':
                    unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]*raw_signals[key].shape[2]))
                elif method in [ 'column' , 'row' ]:
                    unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]))
                else:
                    print('Method \''+method+'\' not supported, use either \'pixel\', \'column\', or \'row\'.')
                    return
                for ii in range(raw_signals[key].shape[0]):
                    if method == 'pixel':
                        unrav_mat[ii,:] = raw_signals[key][ii,:,:].ravel() / monitor[ii]
                    else:
                        unrav_mat[ii,:] = raw_signals[key][ii,:] / monitor[ii]
                pw_matrices_norm.append(unrav_mat)

            # add all other scans
            for ii in scannums[1:]:
                scanname    = 'Scan%03d' % ii
                raw_signals = self.scans[scanname].raw_signals # dict with raw_signals
                monitor     = self.scans[scanname].monitor

                for key,jj in zip(sorted(raw_signals), range(len(raw_signals))):
                    if method == 'pixel':
                        unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]*raw_signals[key].shape[2]))
                    elif method in [ 'column' , 'row']:
                        unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]))
                    for ii in range(raw_signals[key].shape[0]):
                        if method == 'pixel':
                            unrav_mat[ii,:] = raw_signals[key][ii,:].ravel() / monitor[ii]
                        else:
                            unrav_mat[ii,:] = raw_signals[key][ii,:] / monitor[ii]
                    pw_matrices_norm[jj] += unrav_mat
        return pw_matrices_norm

    def get_compensation_factor( self, scan_number, method='sum', roi_number=None ):
        """ **get_compensation_factor**

        Calculates the compensation factor from a given elastic line scan:
        - a pixel-wise center of mass for 'pixel' compensation.
        - a slope (eV/mm) for row-by-row compensation.
        - the center of mass for each ROI for no compensation.

        Args:
            scan_number (int): Scan number of elastic line scan to be used
                for finding the compensation factors.
            method      (str): Keyword describing what kind of compensation
                to be used. Can be \'sum\', \'row\', or \'pixel\'.
            roi_number  (int): ROI number (first ROI is Nr. 0) for which to
                calculate the line-by-line compensation factor.

        """

        # make sure a valid method is selected
        if not method in ['sum', 'row', 'pixel']:
            print( 'get_compensation_factor' )
            print ( 'Method \''+method+'\' not supported, use either \'sum\', \'row\', or \'pixel\'.')
            return

        # make sure there is a ROI defined
        if not self.roi_obj:
            print( 'Please set a ROI object first.' )
            return

        # simple sum: find center of mass for each ROI
        if method == 'sum':
            # reset values
            self.cenom_dict      = {}
            self.resolution_dict = {}

            # get EDF-files and pw_data
            self.load_scan(scan_number, direct=False, scan_type='elastic' )
            scan_name = 'Scan%03d' % scan_number
            self.scans[scan_name].get_raw_signals( self.roi_obj, method='sum' )

            # find CENOM of each ROI
            for key in sorted(self.scans[scan_name].raw_signals):
                cofm = xrs_utilities.find_center_of_mass(self.scans[scan_name].energy,self.scans[scan_name].raw_signals[key])
                self.cenom_dict[key] = cofm

            # set compensation factor
            self.comp_factor = 'sum'

        # pixel-by-pixel: find center of mass for each pixel in each ROI
        if method == 'pixel':
            # reset values
            self.cenom_dict      = {}
            self.resolution_dict = {}

            # get EDF-files and pw_data
            self.load_scan(scan_number, direct=False, scan_type='elastic' )
            scan_name = 'Scan%03d' % scan_number
            self.scans[scan_name].get_raw_signals( self.roi_obj, method='pixel' )

            # find CENOM for each pixel of each ROI
            for key in sorted(self.scans[scan_name].raw_signals):
                self.cenom_dict[key] = np.zeros_like(self.roi_obj.red_rois[key][1],dtype='float')
                for ii in range(self.cenom_dict[key].shape[0]):
                    for jj in range(self.cenom_dict[key].shape[1]):
                        cofm = xrs_utilities.find_center_of_mass(self.scans[scan_name].energy, self.scans[scan_name].raw_signals[key][:,ii,jj])
                        self.cenom_dict[key][ii,jj] = cofm

            # set compensation factor
            self.comp_factor = 'pixel'

        # line-by-line: find a RIXS-type compensation factor
        if method == 'row':
            # reset values
            self.cenom_dict      = {}
            self.resolution_dict = {}

            # which ROI should the compensation factor be calculated for:
            if not roi_number:
                roi_number = int(raw_input("Which ROI to calculate the factor (Python counting)? "))

            # get EDF-files and pw_data
            self.load_scan(scan_number, direct=False, scan_type='elastic' )
            scan_name = 'Scan%03d' % scan_number
            self.scans[scan_name].get_raw_signals( self.roi_obj, method='row' )

            # fit the response for each energy
            plt.ion()
            raw_signals  = self.scans[scan_name].raw_signals['ROI%02d' % roi_number]
            el_positions = np.zeros(raw_signals.shape[0] )
            for ii in range(el_positions.shape[0]):
                x = np.arange(raw_signals.shape[1])*self.PIXEL_SIZE+self.PIXEL_SIZE
                y = raw_signals[ii,:]
                try:
                    popt = optimize.curve_fit(math_functions.gauss_forcurvefit, x, y)[0]
                    g    = math_functions.gauss_forcurvefit(x,popt[0],popt[1],popt[2])
                    el_positions[ii] = (popt[1])
                except:
                    g = np.zeros_like(x)
                    el_positions[ii] = 0.0
                plt.cla()
                plt.plot(x,y,x,g)
                plt.draw()

            # fit the dispersion (eV/mm)
            plt.cla()

            # activate the zoom function already
            thismanager = plt.get_current_fig_manager()
            thismanager.toolbar.zoom()

            # make sure the energy vector is ascending
            sort   =  np.argsort(np.array(self.scans[scan_name].energy)*1e3)
            energy =  self.scans[scan_name].energy*1e3
            center = el_positions[sort]
            energy = energy[sort]

            # make user zoom into fig to define fitting-range
            plt.ion()
            plt.plot(center,energy)
            plt.xlabel('x_0 [mm]')
            plt.ylabel('energy [eV]')
            plt.draw()
            raw_input('Zoom in and press enter to continue.')
            limits = plt.axis()
            inds1  = np.where(np.logical_and(center >= limits[0], center <= limits[1]))[0]
            inds2  = np.where(np.logical_and(energy >= limits[2], energy <= limits[3]))[0]

            # make sure to set the correct window (vertical and horizontal)
            inds   = []
            for ind in inds1:
                if ind in inds2:
                    inds.append(ind)
            inds = np.array(inds)

            # do the fit
            fact = np.polyfit( center[inds], energy[inds], 1)

            # plot the result
            plt.cla()
            plt.plot(center,energy,center,np.polyval(fact, center), center[inds], energy[inds])
            plt.legend(['dispersion', 'fit', 'data-points used'])

            # assign compensation factor
            self.comp_factor = comp_factor = fact[0]
            print (' >>>> The compensation factor is: %0.6f [eV/mm].' %self.comp_factor )
            plt.ioff()


class Hydra_imaging(Hydra):
    """ **Hydra_imaging**
    """

    def __init__( self, path, SPECfname='hydra', EDFprefix='/edf/', EDFname='hydra_', \
                        EDFpostfix='.edf', en_column='sty', moni_column='izero' ):
        Hydra.__init__(self, path, SPECfname='hydra', EDFprefix='/edf/', EDFname='hydra_', \
                        EDFpostfix='.edf', en_column='sty', moni_column='izero')

    def load_scan( self, scan_numbers, scan_type='imaging', direct=True, scaling=None, method='column', scan_motor='STZ' ):
        """ **load_scan**

        Load a single or multiple scans.

        Note:
            When composing a spectrum later, scans of the same 'scan_type' will
            be averaged over. Scans with scan type 'elastic' or 'long' in their
            names are recognized and will be treated specially.

        Args:
              * scan_numbers (int or list): Integer or iterable of scan numbers to be loaded.
              * scan_type            (str): String describing the scan to be loaded (e.g. 'edge1' or 'K-edge').
              * direct           (boolean): Flag, 'True' if the EDF-files should be deleted after loading/integrating the scan.

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

        # keep track of the y- and z-axes
        y_scale = np.array([])
        z_scale = []

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


            # keep track of scales
            y_scale = scan.energy
            z_scale.append(scan.motors[scan_motor])

            # add it to the scans dict
            self.scans[scan_name] = scan

            # add the number to the list of scan numbers
            if not number in self.scan_numbers:
                self.scan_numbers.extend([number])

        self.y_scale = y_scale
        self.z_scale = np.sort(np.array(list(set(z_scale))))

    def load_reference_scan(self, scan_number, scan_type='image_ref', direct=True, scaling=None,):
        scan.load( self.path, self.SPECfname, self.EDFprefix, self.EDFname, self.EDFpostfix, number, \
                   direct=direct, roi_obj=self.roi_obj, scaling=scaling, scan_type=scan_type, \
                   en_column='stx', moni_column=self.moni_column, method=method, cenom_dict=self.cenom_dict )


    def get_compensation_factor( self, el_scan_numbers, scan_motor='sty',  plotting=False):
        """ **get_compensation_factor**

        Calculates the compensation factor for the case of imaging:

        Args:
            scan_number (int): Scan number of elastic line scan to be used
                for finding the compensation factors.
            method      (str): Keyword describing what kind of compensation
                to be used. Can be \'sum\', \'row\', or \'pixel\'.
            roi_number  (int): ROI number (first ROI is Nr. 0) for which to
                calculate the line-by-line compensation factor.

        """

        # make sure there is a ROI defined
        if not self.roi_obj:
            print( 'Please set a ROI object first.' )
            return

        # load the scans around the elastic line
        self.load_scan( el_scan_numbers, method='column', direct=True, scan_type='elastic')

        ## parse energy scale and signals for each ROI
        #energy  = {} #np.zeros(len(range(440,589)))
        #signals = {} #np.zeros(len(range(440,589)))
        #for roi_key in sorted(self.roi_obj.red_rois):
        #	energy[roi_key]  = np.zeros(len(el_scan_numbers))
        #	signals[roi_key] = np.zeros(len(el_scan_numbers))
        #	for scan_number, step in zip(el_scan_numbers, range(len(el_scan_numbers))):
        #		scan_name = 'Scan%03d'%scan_number
        #		energy[roi_key][step]  = self.scans[scan_name].motors['energy']
        #		signals[roi_key][step] = np.sum(self.scans[scan_name].raw_signals[roi_key])

        # parse energy scale
        energy = np.array([])
        for scan_key in sorted(self.scans):
            energy = np.append(energy, self.scans[scan_key].motors['energy'])
        energy = np.array(list(set(energy)))

        # parse signals
        signals = {}
        for roi_key, (pos,M) in sorted(self.roi_obj.red_rois.items()):
            signals[roi_key] = np.zeros( (len(energy), ) )
            for en,estep in zip(energy, range(len(energy))):
                for scan_key in sorted(self.scans):
                    if self.scans[scan_key].motors['energy'] == en:
                        signals[roi_key][estep] += np.sum(self.scans[scan_key].raw_signals[roi_key])


        # get centers of mass for each ROI and save in self.cenom_dict
        for key in signals:
            inds = np.argsort(energy)
            cofm = xrs_utilities.find_center_of_mass(energy[inds],signals[key][inds])
            self.cenom_dict[key] = cofm
            if plotting:
                plt.cla()
                plt.plot(energy[inds],signals[key][inds],'b-')
                plt.plot([cofm, cofm], [np.amin(signals[key]), np.amax(signals[key])], 'k-')
                plt.xlabel('energy [keV]')
                plt.ylabel('intensity [arb. units]')
                plt.title('CENOM for %s'%(key))
                plt.waitforbuttonpress()

        self.E0 = np.mean( [self.cenom_dict[key] for key in self.cenom_dict if self.cenom_dict[key]>0.0] )

    def interpolate_scans(self, step_motor='STZ'):
        """ **compensate_scans**

        Interpolate signals onto common energy-loss grid.
        """
        # define master energy-loss scale
        energy = np.array([])
        for scan_key in sorted(self.scans):
            energy = np.append(energy, self.scans[scan_key].motors['energy'])

        energy = np.sort(np.array(list(set(energy))))

        master_eloss = ( energy - self.E0 )*1.0e3
        self.eloss   = master_eloss

        # define master sty scale
        scan_scale_1 = self.scans[scan_key].energy

        # define master stz scale
        dmy = []
        for scan_key in self.scans:
            dmy.append(self.scans[scan_key].motors[step_motor])
        scan_scale_2 = np.array(list(set(dmy)))

        # concatenate scans to volume
        for roi_key, (pos,M) in sorted(self.roi_obj.red_rois.items()):
            self.raw_signals[roi_key] = np.zeros( (len(energy), len(scan_scale_1), M.shape[1], len(scan_scale_2)) )
            for en,estep in zip(energy, range(len(energy))):
                for scan_key in sorted(self.scans):
                    if self.scans[scan_key].motors['energy'] == en:
                        for zstep in range(len(scan_scale_2)):
                            if self.scans[scan_key].motors[step_motor] == scan_scale_2[zstep]:
                                self.raw_signals[roi_key][estep,:,:,zstep] = self.scans[scan_key].raw_signals[roi_key]

        # interpolate everything onto master energy-loss scale
        self.raw_signals_int = {}
        for roi_key in sorted(self.raw_signals):
            x = (energy - self.cenom_dict[roi_key])*1.0e3
            y = self.raw_signals[roi_key]
            f = interp1d(x, y, kind='linear', axis=0, bounds_error=False, fill_value=0.0)
            self.raw_signals_int[roi_key] = f(master_eloss)

class Fourc:
	"""Main class for handling RIXS data from ID20's high-resolution spectrometer 'Fourc'.

	This class is intended to read SPEC- and according EDF-files and perform dispersion
	compensations.

	Note:
		'Fourc' is the name of the high-energy-resolution spectrometer at ESRF's
		ID20 beamline. This class has been adopted specifically for this spectrometer.

		If you are using this program, please cite the following work:

		Sahle, Ch J., A. Mirone, J. Niskanen, J. Inkinen, M. Krisch, and S. Huotari.
		"Planning, performing and analyzing X-ray Raman scattering experiments."
		Journal of Synchrotron Radiation 22, No. 2 (2015): 400-409.

	Args:
		  * path        (str): Absolute path to directory holding the data.
		  * SPECfname   (str): Name of the SPEC-file ('rixs' is the default).
		  * EDFprefix   (str): Prefix for the EDF-files ('/edf/' is the default).
		  * EDFname     (str): Filename of the EDF-files ('rixs_' is the default).
		  * EDFpostfix  (str): Postfix for the EDF-files ('.edf' is the default).
		  * en_column   (str): Counter mnemonic for the energy motor ('energy' is the default).
		  * moni_column (str): Mnemonic for the monitor counter ('izero' is the default).
		  * EinCoor    (list): Coordinates, where to find the incident energy value in the
		    SPEC-file (default is [9,0])

	Attributes:
		  * path        (str): Absolute path to directory holding the data.
		  * SPECfname   (str): Name of the SPEC-file ('hydra' is the default).
		  * EDFprefix   (str): Prefix for the EDF-files ('/edf/' is the default).
		  * EDFname     (str): Filename of the EDF-files ('hydra_' is the default).
		  * EDFpostfix  (str): Postfix for the EDF-files ('.edf' is the default).
		  * en1_column  (str): Counter mnemonic for the energy motor ('anal energy' is the default).
		  * en2_column  (str): Counter mnemonic for the energy motor ('energy' is the default).
		  * moni_column (str): Mnemonic for the monitor counter ('izero' is the default).
		  * EinCoor    (list): Coordinates, where to find the incident energy value in the
		    SPEC-file (default is [9,0])
		  * scans        (dict): Dictionary holding all loaded scans.
		  * scan_numbers (list): List with scan number of all loaded scans.
		  * energy   (np.array): Array with the common energy scale.
		  * energy2  (np.array): Array with the common energy2 scale.
		  * signals  (np.array): Array with the signals for all analyzers (one column per anayzer).
		  * errors   (np.array): Array with the poisson errors for all analyzers (one column per anayzer).
		  * groups       (dict): Dictionary of groups of scans (instances of the 'scangroup' class,
		    such as 2 'elastic', or 5 'edge1', etc.).
		  * tth          (list): List of all scattering angles (one value for each ROI).
		  * resolution   (list): List of FWHM of the elastic lines (one for each analyzer).
		  * roi_obj  (instance): Instance of the roi_object class from the x0.055 # pixel size in mmrs_rois module defining all
		    ROIs for the current dataset (default is 'None').
		  * comp_factor (float): Compensation factor used for the dispersion correction.
		  * PIXEL_SIZE  (float): Pixel size of the used Maxipix detector (in mm).

	"""

	def __init__( self, path, SPECfname='rixs', EDFprefix='/edf/', EDFname='rixs_', \
						EDFpostfix='.edf', moni_column='izero', EinCoor='energy' ):

		self.path        = path
		self.SPECfname   = SPECfname

		if not os.path.isfile( path + SPECfname ):
			raise Exception( 'IOError! No such file or directory.' )

		self.EDFprefix   = EDFprefix
		self.EDFname     = EDFname
		self.EDFpostfix  = EDFpostfix
		#self.en1_column  = en1_column.lower()
		#self.en2_column  = en2_column.lower()
		self.en_column   = None
		self.moni_column = moni_column.lower()
		self.EinCoor     = EinCoor

		self.scans         = {}
		self.scan_numbers  = []
		self.energy        = []
		self.energy2       = []
		self.signals       = []
		self.errors        = []

		self.groups        = {}
		self.tth           = []
		self.resolution    = []
		self.cenom_dict    = {}
		self.raw_signals   = {}
		self.raw_errors    = {}

		self.roi_obj = []

		self.comp_factor   = 0.0
		self.comp_type     = None
		self.PIXEL_SIZE    = 0.055 # pixel size in mm

		print_citation_message( )

	def SumDirect( self,scan_numbers , clean_edf_stack=False ):
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

		# sum up all EDF-images from the given scans
		im_sum    = None
		en_column = None # uses first column in SPEC file as scanned motor
		for number in numbers:
			scan = xrs_scans.Scan()
			scan.load(self.path, self.SPECfname, self.EDFprefix, self.EDFname, self.EDFpostfix, number, \
					direct=False, roi_obj=None, scaling=None, scan_type='generic', \
					en_column=en_column, moni_column=self.moni_column, clean_edf_stack=clean_edf_stack)

			if im_sum is None:
				im_sum = np.zeros(scan.edfmats[0].shape ,"f")
			im_sum[:] += scan.edfmats.sum(axis=0)
		return im_sum

	def set_roiObj( self, roiobj ):
		""" **set_roiObj**

		Assign an instance of the 'roi_object' class to the current data set.

		Args:
			roiobj (instance): Instance of the 'roi_object' class holding all
				information about the definition of the ROIs.

		"""
		self.roi_obj = roiobj

	def load_scan( self, scan_numbers, direct=True, comp_factor=None, scan_type='generic', scaling=None, method='sum', rot_angles=None, clean_edf_stack=False ):
		""" **load_scan**

		Loads given scans and applies the dispersion compensation.

		Args:
			scan_numbers (int or list): Scan number(s) of scans to be loaded.
			direct           (boolean): Flag, if set to 'True', EDF-files are
				deleted after loading the scan (this is the default).
			comp_factor        (float): Compensation factor to be used. If 'None',
				the global compensation factor will be used. If provided, the global
				compensation factor will be overwritten.
			scan_type            (str): String describing the scan to be loaded.

		Note:
			If a compensation factor is passed to this function, the classes 'globel'
			compensation factor is overwritten.

        """

		# make sure scan_numbers are iterable
		if not isinstance(scan_numbers,list):
			numbers = []
			numbers.append(scan_numbers)
		else:
			numbers = scan_numbers

		# load scan/scans
		for number in numbers:

			# create a name for each scan
			scan_name = 'Scan%03d' % number

			# create an instance of the Scan class
			scan      = xrs_scans.Scan()

			# load scan, first column in SPEC file will be scanned motor
			scan.load( self.path, self.SPECfname, self.EDFprefix, self.EDFname, \
				self.EDFpostfix, number, direct=direct, roi_obj=self.roi_obj, \
				scaling=scaling, scan_type=scan_type, en_column=self.en_column, \
				moni_column=self.moni_column, method=method, cenom_dict=self.cenom_dict,\
				comp_factor=comp_factor,rot_angles=rot_angles, clean_edf_stack=clean_edf_stack )

			# assign one dictionary entry to each scan
			self.scans[scan_name] = scan
			if not number in self.scan_numbers:
				self.scan_numbers.extend([number])

			# save the incident energy
			try:
				self.scans[scan_name].Ein = scan.motors[self.EinCoor]
			except:
				self.scans[scan_name].Ein = None

	def get_compensation_factor( self, scan_number, method='sum', roi_number=None, rot_angles=None ):
		""" **get_compensation_factor**

		Calculates the compensation factor from a given elastic line scan:
		- a pixel-wise center of mass for 'pixel' compensation.
		- a slope (eV/mm) for row-by-row compensation.
		- the center of mass for each ROI for no compensation.

		Args:
			scan_number (int): Scan number of elastic line scan to be used
				for finding the compensation factors.
			method      (str): Keyword describing what kind of compensation
				to be used. Can be \'sum\', \'row\', or \'pixel\'.
			roi_number  (int): ROI number (first ROI is Nr. 0) for which to
				calculate the line-by-line compensation factor.

		"""

		# make sure there is a ROI defined
		if not self.roi_obj:
			print( 'Please set a ROI object first.' )
			return

		# simple sum: find center of mass for each ROI
		if method == 'sum':
			# reset values
			self.cenom_dict      = {}
			self.resolution_dict = {}

			# get EDF-files and pw_data
			elastic_scan = xrs_scans.Scan()
			elastic_scan.load( self.path, self.SPECfname, self.EDFprefix, \
								self.EDFname, self.EDFpostfix, scan_number, \
								direct=False, scan_type='elastic', \
								moni_column=self.moni_column, method='sum' )

			elastic_scan.get_raw_signals( self.roi_obj, method='sum' )

			# find CENOM of each ROI
			for key in elastic_scan.raw_signals:
				cofm = xrs_utilities.find_center_of_mass( elastic_scan.energy,elastic_scan.raw_signals[key] )
				self.cenom_dict[key] = cofm

			# set compensation factor
			self.comp_type = 'sum'

		# pixel-by-pixel: find center of mass for each pixel in each ROI
		elif method == 'pixel':
			# reset values
			self.cenom_dict      = {}
			self.resolution_dict = {}

			# get EDF-files and pw_data
			elastic_scan = xrs_scans.Scan()
			elastic_scan.load( self.path, self.SPECfname, self.EDFprefix, \
								self.EDFname, self.EDFpostfix, scan_number, \
								direct=False, scan_type='elastic', \
								moni_column=self.moni_column, method='pixel' )

			elastic_scan.get_raw_signals( self.roi_obj, method='pixel' )

			# find CENOM for each pixel of each ROI
			for key in sorted(elastic_scan.raw_signals):
				self.cenom_dict[key] = np.zeros_like(self.roi_obj.red_rois[key][1])
				for ii in range(self.cenom_dict[key].shape[0]):
					for jj in range(self.cenom_dict[key].shape[1]):
						cofm = xrs_utilities.find_center_of_mass(elastic_scan.energy, elastic_scan.raw_signals[key][:,ii,jj])
						self.cenom_dict[key][ii,jj] = cofm

			# set compensation factor
			self.comp_type = 'pixel'

		# line-by-line: find a RIXS-type compensation factor
		elif method == 'row':
			# reset values
			self.cenom_dict      = {}
			self.resolution_dict = {}
			self.comp_type = 'row'

			# which ROI should the compensation factor be calculated for:
			if not roi_number:
				roi_number = int(raw_input("Which ROI to calculate the factor (Python counting)? "))

			self.comp_factor = 0.0

			# get EDF-files and pw_data
			elastic_scan = xrs_scans.Scan()
			elastic_scan.load( self.path, self.SPECfname, self.EDFprefix, \
								self.EDFname, self.EDFpostfix, scan_number, \
								direct=False, scan_type='elastic', \
								moni_column=self.moni_column, method='row', rot_angles=rot_angles )

			elastic_scan.get_raw_signals( self.roi_obj, method='row' )

			# fit the response for each energy
			plt.ion()
			raw_signals  = elastic_scan.raw_signals['ROI%02d' % roi_number]
			el_positions = np.zeros(raw_signals.shape[0] )
			for ii in range(el_positions.shape[0]):
				x = np.arange(raw_signals.shape[1])*self.PIXEL_SIZE+self.PIXEL_SIZE
				y = raw_signals[ii,:]
				try:
					p0 = [ np.amax(y), xrs_utilities.find_center_of_mass(x,y), xrs_utilities.fwhm(x,y)[1]]
					popt = optimize.curve_fit(math_functions.gauss_forcurvefit, x, y, p0=p0)[0]
					g    = math_functions.gauss_forcurvefit(x,popt[0],popt[1],popt[2])
					el_positions[ii] = (popt[1])
				except:
					g = np.zeros_like(x)
					el_positions[ii] = 0.0
				plt.cla()
				plt.plot(x,y,x,g)
				plt.draw()

			# fit the dispersion (eV/mm)
			plt.cla()

			# activate the zoom function already
			thismanager = plt.get_current_fig_manager()
			thismanager.toolbar.zoom()

			# make sure the energy vector is ascending
			sort   =  np.argsort(np.array(elastic_scan.energy)*1e3)
			energy =  elastic_scan.energy*1e3
			center = el_positions[sort]
			energy = energy[sort]

			# make user zoom into fig to define fitting-range
			plt.ion()
			plt.plot(center,energy)
			plt.xlabel('x_0 [mm]')
			plt.ylabel('energy [eV]')
			plt.draw()
			raw_input('Zoom in and press enter to continue.')
			limits = plt.axis()
			inds1  = np.where(np.logical_and(center >= limits[0], center <= limits[1]))[0]
			inds2  = np.where(np.logical_and(energy >= limits[2], energy <= limits[3]))[0]

			# make sure to set the correct window (vertical and horizontal)
			inds   = []
			for ind in inds1:
				if ind in inds2:
					inds.append(ind)
			inds = np.array(inds)

			# do the fit
			fact = np.polyfit( center[inds], energy[inds], 1)

			# plot the result
			plt.cla()
			plt.plot(center,energy,center,np.polyval(fact, center), center[inds], energy[inds])
			plt.legend(['dispersion', 'fit', 'data-points used'])

			# assign compensation factor
			self.comp_factor = fact[0]
			print (' >>>> The compensation factor is: %0.6f [eV/mm].' %self.comp_factor )
			plt.ioff()

		# if method is unknown
		else:
			print ( 'Method \''+method+'\' not supported, use either \'sum\', \'row\', or \'pixel\'.')
			return

	def get_raw_data( self, method='sum', scaling=None ):
		""" **get_raw_data**

		Applies the ROIs to extract the raw signals from the
		EDF-files.

		"""

		# make sure there are some ROIs set
		if not self.roi_obj:
			print ( 'Did not find a ROI object, please set one first.' )
			return

		# call the get_raw_signals function for each scan
		for scan in self.scans:
			if len(self.scans[scan].edfmats):
				print ( "Integrating " + scan + " using method \'" + method + '\'.')
				self.scans[scan].get_raw_signals( self.roi_obj , method=method, scaling=scaling )

	def get_data( self, method='sum', scaling=None ):
		""" **get_data**

		Applies the ROIs to the EDF-files.

		This extracts the raw-data from the EDF-files subject to three
		different methods: 'sum' will simply sum up all pixels inside a
		ROI, 'row' will sum up over the dispersive direction, and 'pixel'
		will not sum at all.

		Args:
			method   (str): Keyword specifying the selected choice of data treatment:
				can be 'sum', 'row', or 'pixel'. Default is 'sum'.
			scaling (list): Optional scaling factors (one per ROI) to scale
				the data with.

		"""
		if not self.cenom_dict:
			print ( 'No compensation factor/center-of-mass info found, please provide first.' )
			return

		for scan in self.scans:
			# make sure there is raw data available
			if not self.scans[scan].raw_signals:
				self.get_raw_data( self, method=method, scaling=scaling )

			# apply compensation / get signals
			self.scans[scan].get_signals( method=method, cenom_dict=self.cenom_dict, comp_factor=self.comp_factor, scaling=scaling )

	def get_XES_spectrum( self, method='sum', interpolation=False ):
		""" **get_XES_spectrum**

		Constructs an emission spectrum based on the loaded single spectra.

		NOTE: this needs support for all methods ('sum', 'pixel', 'row')

		"""

		# cenom_dict for XES is zeros everywhere
		if method == 'sum':
			self.cenom_dict = {}
			for key in self.roi_obj.red_rois:
				self.cenom_dict[key] = 0.0

		# make sure there is data available
		for scan in self.scans:
			if not np.any(self.scans[scan].signals):
				self.get_data( method=method )

		# find the groups
		all_groups = xrs_scans.findgroups( self.scans )

		for group in all_groups:
			self.groups[group[0].get_type()] = xrs_scans.sum_scans_to_group( group, method=method, interp=interpolation )
			if method == 'sum':
				self.groups[group[0].get_type()].get_signals(method='sum', cenom_dict=self.cenom_dict )

		self.energy, self.signals, self.errors = xrs_scans.get_XES_spectrum( self.groups )

	def get_Ein_RIXS_map( self, scan_numbers, roi_number, logscaling=False, file_name=None):
		""" **get_Ein_RIXS_map**

		Returns a RIXS map that plots energy loss vs inciden energy for the
		specified scans.

		Args:
		        *  scan_numbers (int or list): SPEC scan numbers to be deleted.
                        *  logscaling       (boolean): If true numbers are returned on logarithmic scale
                           (False by default)
		        *  file_name            (str): Absolute path, if map should also be written into
			   an ascii-file.
		        *  roi_number           (int): ROI to use for creating the RIXS map (Python counting).

		"""

		# make sure scan_numbers are iterable
		numbers = []
		if not isinstance(scan_numbers,list):
			numbers.append(scan_numbers)
		else:
			numbers = scan_numbers

		# create the matrix
		scanname = 'Scan%03d' % numbers[0]
		pre_map  = np.zeros( ( len(self.scans[scanname].energy), len(numbers) ) )
		rixs_map = np.zeros( ( len(self.scans[scanname].energy), len(numbers) ) )
		e_transfer = np.flipud(self.scans[scanname].energy - self.scans[scanname].Ein)*1e3

		# save the incident energies
		e_incident = []

		# fill the matrix scan-by-scan
		for number,ii in zip(numbers, range(len(numbers))):
			scanname = 'Scan%03d' % number
			e_incident.append(self.scans[scanname].Ein)
			pre_map[:,ii] = np.interp( e_transfer, np.flipud(self.scans[scanname].energy - self.scans[scanname].Ein)*1e3, np.flipud(self.scans[scanname].signals[:,roi_number]))

		sort_inds = np.argsort(e_incident)
		e_incident_sort = []
		for ii in range(pre_map.shape[1]):
			rixs_map[:,ii] = pre_map[:,sort_inds[ii]]
			e_incident_sort.append(e_incident[sort_inds[ii]])

		return np.array(rixs_map), np.array(e_incident_sort), -np.flipud(np.array(e_transfer))


	def delete_scan( self, scan_numbers ):
		""" **delete_scan**

		Deletes scans of given scan numbers.

		Args:
			scan_numbers (int or list): SPEC scan numbers to be used for the RIXS map.

		"""

		# make sure scan_numbers are iterable
		numbers = []
		if not isinstance(scan_numbers,list):
			numbers.append(scan_numbers)
		else:
			numbers = scan_numbers

		# delete scans
		for number in numbers:
			scanname = 'Scan%03d' % number
			del(self.scans[scanname])
			self.scan_numbers.remove(number)

	def dump_scans_ascii( self, scan_numbers, pre_fix, f_name, post_fix='.dat', header='' ):
		""" **dum_scans_ascii**

		Produce ASCII-type files with columns of energy, signal, and Poisson error.

		Args:
		scan_numbers (int or list): SPEC scan numbers of scans to be safed in ASCII format.
		pre_fix  (str): Path to directory where files should be written into.
		f_name   (str): Base name for the files to be written.
		post_fix (str): Extention for the files to be written (default is '.dat').

		"""
		# make sure scan_numbers are iterable
		numbers = []
		if not isinstance(scan_numbers,list):
			numbers.append(scan_numbers)
		else:
			numbers = scan_numbers

		for number in numbers:
			scan_name = 'Scan%03d' % number
			if not scan_name in self.scans.keys():
				print ('Scan No. %03d is currently not loaded.'% number)
				return
			x = self.scans[scan_name].energy
			y = self.scans[scan_name].signals
			z = self.scans[scan_name].errors
			data = np.zeros((x.shape[0], y.shape[1]+z.shape[1]+1 ))
			data[:,0] = x
			data[:,1:1+y.shape[1]] = y
			data[:,1+y.shape[1]:1+y.shape[1]+z.shape[1]] = z
			file_name = pre_fix + f_name + '_' + scan_name + post_fix
			if len(header) == 0:
				header = ' Ein = ' + str(self.scans[scan_name].Ein) + ' keV'
			np.savetxt(file_name, data, header=header)

	def copy_edf_files( self, scan_numbers, dest_dir ):
		""" **copy_edf_files**

		Copies all EDF-files from given scan_numbers into given directory.

		Args:
			*       scan_numbers (int or list) = Integer or list of integers defining the scan
				numbers of scans to be copied.
			*       dest_dir             (str) = String with absolute path for the destination.

		"""

		import shutil

		# make sure scan_numbers is iterable
		numbers = []
		if not isinstance(scan_numbers,list):
			numbers.append(scan_numbers)
		else:
			numbers = scan_numbers
		fname = self.path + self.SPECfname

		# find EDF-file names and copy them
		for nscan in numbers:
			if use_PyMca:
				data, motors, counters = xrs_fileIO.PyMcaSpecRead(fname,nscan)
			else:
				data, motors, counters = xrs_fileIO.SpecRead(fname,nscan)

			for m in range(len(counters['ccdno'])):
				ccdnumber = counters['ccdno'][m]
				edfname   = self.path + self.EDFprefix + self.EDFname + "%04d" % ccdnumber + self.EDFpostfix
				shutil.copy2( edfname, dest_dir )

	def get_pw_matrices( self, scan_numbers, method='pixel' ):
		"""	**get_pw_matrices**

		Sums scans from pixelwise ROI integration for use in the pixel-wise ROI
		refinement.

		Args:
		        * scan_numbers (int, list): Integer or list of scan numbers to be added.
			* method             (str): Keyword describing how to return the data
			  (possible values are: 'pixel' for pixel-wise integration, 'column'
		          for column-wise summation, and 'row' for row-wise summation.

		Returns:
		       * Raw data in pixel-wise format.

		"""
		# make scan_numbers iterable
		if isinstance(scan_numbers,list):
			scannums = scan_numbers
		elif isinstance(scan_numbers,int):
			scannums = [scan_numbers]
		else:
			print( 'Please provide keyword \'scan_numbers\' as integer or list of integers.' )
			return

		# make sure method is known
		if not method in ['pixel', 'column', 'row']:
			print('Unknown integration method. Use either \'pixel\', \'column\', or \'row\'')
			return

		# make sure scan exists
		for scannum in scannums:
			if not scannum in self.scan_numbers:
				self.load_scan(scannum, direct=False, method=method)

		# make sure raw_signals is existing
		if not self.scans['Scan%03d' % scannums[0]].raw_signals:
			self.get_raw_data( method=method )

		# one scan only
		if len(scannums)==1:
			scanname = 'Scan%03d' % scannums[0]
			raw_signals = self.scans[scanname].raw_signals # dict with raw_signals
			monitor     = self.scans[scanname].monitor

			# normalize data
			pw_matrices_norm = []
			for key in sorted(raw_signals):
				if method == 'pixel':
					unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]*raw_signals[key].shape[2]))
				elif method in [  'column' , 'row']:
					unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]))
				else:
					print('Method \''+method+'\' not supported, use either \'pixel\', \'column\', or \'row\'.')
					return
				for ii in range(raw_signals[key].shape[0]):
					if method == 'pixel':
						unrav_mat[ii,:] = raw_signals[key][ii,:,:].ravel() / monitor[ii]
					else:
						unrav_mat[ii,:] = raw_signals[key][ii,:] / monitor[ii]
				pw_matrices_norm.append(unrav_mat)

		# multiple scans to be added
		else:
			# start with first scan
			scanname    = 'Scan%03d' % scannums[0]
			raw_signals = self.scans[scanname].raw_signals # dict with raw_signals
			monitor     = self.scans[scanname].monitor

			# normalize data (first scan)
			pw_matrices_norm = []
			for key in sorted(raw_signals):
				if method == 'pixel':
					unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]*raw_signals[key].shape[2]))
				elif method in [ 'column' , 'row']:
					unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]))
				else:
					print('Method \''+method+'\' not supported, use either \'pixel\', \'column\', or \'row\'.')
					return
				for ii in range(raw_signals[key].shape[0]):
					if method == 'pixel':
						unrav_mat[ii,:] = raw_signals[key][ii,:,:].ravel() / monitor[ii]
					else:
						unrav_mat[ii,:] = raw_signals[key][ii,:] / monitor[ii]
				pw_matrices_norm.append(unrav_mat)

			# add all other scans
			for ii in scannums[1:]:
				scanname    = 'Scan%03d' % ii
				raw_signals = self.scans[scanname].raw_signals # dict with raw_signals
				monitor     = self.scans[scanname].monitor

				for key,jj in zip(sorted(raw_signals), range(len(raw_signals))):
					if method == 'pixel':
						unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]*raw_signals[key].shape[2]))
					elif method in [  'column' ,'row'] :
						unrav_mat = np.zeros((raw_signals[key].shape[0], raw_signals[key].shape[1]))
					for ii in range(raw_signals[key].shape[0]):
						if method == 'pixel':
							unrav_mat[ii,:] = raw_signals[key][ii,:].ravel() / monitor[ii]
						else:
							unrav_mat[ii,:] = raw_signals[key][ii,:] / monitor[ii]
					pw_matrices_norm[jj] += unrav_mat
		return pw_matrices_norm


class read_id20:
    """
    Main class for handling raw data from XRS experiments on ESRF's ID20. This class
    is used to read scans from SPEC files and the according EDF-files, it provides access
    to all tools from the xrs_rois module for defining ROIs, it can be used to integrate
    scans, sum them up, stitch them together, and define the energy loss scale.
    INPUT:
      * absfilename   = path and filename of the SPEC-file
      * energycolumn  = name (string) of the counter for the energy as defined in the SPEC session (counter mnemonic)
      * monitorcolumn = name (string) of the counter for the monitor signals as defined in the SPEC session (counter mnemonic)
      * edfName       = name/prefix (string) of the EDF-files (default is the same as the SPEC-file name)
      * single_image  = boolean switch, 'True' (default) if all 6 detectors are merged in a single image,
        'False' if two detector images per point exist.
    """
    def __init__(self,absfilename,energycolumn='energy',monitorcolumn='kap4dio',edfName=None,single_image=True):
        self.scans = {} # was a dictionary before
        if absfilename is not None:
            if not os.path.isfile(absfilename):
                raise Exception('IOError! No such file, please check filename.')
            self.path          = os.path.split(absfilename)[0] + '/'
            self.filename = os.path.split(absfilename)[1]
            if not edfName:
                self.edfName = os.path.split(absfilename)[1]
            else:
                self.edfName = edfName
        self.single_image  = single_image
        self.scannumbers   = []
        self.EDF_PREFIXh   = 'edf/h_'
        self.EDF_PREFIXv   = 'edf/v_'
        self.EDF_PREFIX    = 'edf/'
        self.EDF_POSTFIX   = '.edf'
        self.DET_PIXEL_NUMx = 768
        self.DET_PIXEL_NUMy = 512
        self.DET_PIXEL_NUM  = 256
        # which column in the SPEC file to be used for the energy and monitor
        self.encolumn      = energycolumn.lower()
        self.monicolumn    = monitorcolumn.lower()
        # here are the attributes of the old rawdata class
        self.eloss    = []	# common eloss scale for all analyzers
        self.energy   = []	# common energy scale for all analyzers
        self.signals  = []	# signals for all analyzers
        self.errors   = []	# poisson errors
        self.qvalues  = []	# for all analyzers
        self.groups   = {}  # dictionary of groups (such as 2 'elastic', or 5 'edge1', etc.)
        self.cenom    = []  # list of center of masses of the elastic lines
        self.E0       = []  # energy value, mean value of self.cenom
        self.tth      = []  # list of scattering angles (one for each ROI)
        self.VDtth    = []
        self.VUtth    = []
        self.VBtth    = []
        self.HRtth    = []
        self.HLtth    = []
        self.HBtth    = []
        self.resolution   = []  # list of FWHM of the elastic lines for each analyzer
        self.signals_orig = []  # signals for all analyzers before interpolation
        self.errors_orig  = []  # errors for all analyzers before interpolation
        # TTH offsets from center of V and H modules
        # tth offsets in one direction (horizontal for V-boxes, vertical for H-boxes)
        self.TTH_OFFSETS1   = np.array([5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0])
        # tth offsets in the other direction (horizontal for H-boxes, vertical for V-boxes)
        self.TTH_OFFSETS2   = np.array([-9.71, -9.75, -9.71, -3.24, -3.25, -3.24, 3.24, 3.25, 3.24, 9.71, 9.75, 9.71])
        # input
        self.roi_obj = [] # an instance of the roi_object class from the xrs_rois module (new)



    def save_state_hdf5(self, filename, groupname, comment=""):
        import h5py
        h5 = h5py.File(filename,"a")

        h5.require_group(groupname)
        h5group =  h5[groupname]
        if(   "eloss" in h5group.keys() ):
            raise Exception(" Read data already present in  " + filename+":"+groupname)
        for key in [
    "eloss",
    "energy",
    "signals",
    "errors",
    "qvalues",
    ########### "groups",
    "cenom",
    "E0",
    "tth",
    "VDtth",
    "VUtth",
    "VBtth",
    "HRtth",
    "HLtth",
    "HBtth",
    "resolution",
    "signals_orig",
    "errors_orig"
            ]:
            h5group[key]  = getattr(self,key)
        h5group["comment"]  = comment
        h5.flush()
        h5.close()

    def load_state_hdf5(self, filename, groupname):
        import h5py

        print ( "filename  " , filename )
        print ( "groupname  " , groupname )

        h5 = h5py.File(filename,"r")

        h5group =  h5[groupname]

        chiavi = {"eloss":array,
                  "energy":array,
                  "signals":array,
                  "errors":array,
                  "qvalues":array,
                  ########### "groups":array,
                  "cenom":array,
                  "E0":float,
                  "tth":array,
                  "VDtth":array,
                  "VUtth":array,
                  "VBtth":array,
                  "HRtth":array,
                  "HLtth":array,
                  "HBtth":array,
                  "resolution":array,
                  "signals_orig":array,
                  "errors_orig":array
                  }

        for key in chiavi:
            # print key
            setattr(self,key,chiavi[key]( array(h5group[key])))

        h5.flush()
        h5.close()


    def save_state(self):
        d={}
        for key in [
    "eloss",
    "energy",
    "signals",
    "errors",
    "qvalues",
    "groups",
    "cenom",
    "E0",
    "tth",
    "VDtth",
    "VUtth",
    "VBtth",
    "HRtth",
    "HLtth",
    "HBtth",
    "resolution",
    "signals_orig",
    "errors_orig"
            ]:
            d[key]  = getattr(self,key)
        f=open("datas.pick","w")
        pickle.dump(d,f)
        f.close()

    def there_is_a_valid_roi_at(self,n):
         return n<len(self.roi_obj.indices) and len(self.roi_obj.indices[n])

    def read_just_first_scanimage(self,scannumber):
        fn = self.path + self.filename
        if use_PyMca == True:
            data, motors, counters, lables = xrs_fileIO.PyMcaSpecRead_my(fn,scannumber)
        else:
            data, motors, counters = xrs_fileIO.SpecRead(fn,scannumber)

        edfmat  =  xrs_fileIO.ReadEdf_justFirstImage(counters['ccdno'],
                                                     self.path,
                                                     self.EDF_PREFIX,
                                                     self.edfName, self.EDF_POSTFIX)
        return edfmat


    def readscan_new(self,scannumber,fromtofile=False):
        """
        Returns the data, motors, counter-names, and edf-files from the SPEC file defined when
        the xrs_read object was initiated.
        There should be an alternative that uses the PyMca module if installed.
        INPUT:
          * scannumber = number of the scan to be loaded
          * fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        # first see if scan can be loaded from npz-file
        if fromtofile:
            scanname = 'Scan%03d' % scannumber
            sub_path = os.path.split(self.path[:-1])[0]
            fname    = sub_path + '/scans/' + scanname + '.npz'
            return ReadScanFromFile(fname)
        else:
            print( 'Parsing EDF- and SPEC-files of scan No. %s' % scannumber )

        # load SPEC-file
        fn = self.path + self.filename
        if use_PyMca == True:
            data, motors, counters, lables = xrs_fileIO.PyMcaSpecRead_my(fn,scannumber)
        else:
            data, motors, counters = xrs_fileIO.SpecRead(fn,scannumber)

        # load EDF-files
        if not self.single_image:
            edfmats = xrs_fileIO.ReadEdfImages_TwoImages(counters['ccdno'], self.DET_PIXEL_NUMx, self.DET_PIXEL_NUMy, self.path, self.EDF_PREFIXh, self.EDF_PREFIXv, self.edfName, self.EDF_POSTFIX)
        else:
            edfmats = xrs_fileIO.ReadEdfImages(counters['ccdno'], self.DET_PIXEL_NUMx, self.DET_PIXEL_NUMy, self.path, self.EDF_PREFIX, self.edfName, self.EDF_POSTFIX)

        # add the scannumber to self.scannumbers, if not already present
        if not scannumber in self.scannumbers:
            self.scannumbers.extend([scannumber])

        # store scan in numpy zip-archive, if desired
        if fromtofile:
            scanname = 'Scan%03d' % scannumber
            sub_path = os.path.split(self.path[:-1])[0]
            fname = sub_path + '/scans/' + scanname
            xrs_fileIO.WriteScanToFile(fname,data,motors,counters,edfmats)

        return data, motors, counters, edfmats #<type 'list'> <type 'list'> <type 'dict'> <type 'numpy.ndarray'

    def readscan(self,scannumber,fromtofile=False):
        """
        Returns the data, motors, counter-names, and edf-files from the SPEC file defined when
        the xrs_read object was initiated.
        There should be an alternative that uses the PyMca module if installed.
        INPUT:
          * scannumber = number of the scan to be loaded
          * fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        # first see if scan can be loaded from npz-file
        if fromtofile:
            try:
                print( 'Trying to load scan from file.' )
                scanname = 'Scan%03d' % scannumber
                sub_path = os.path.split(self.path[:-1])[0]
                fname    = sub_path + '/scans/' + scanname + '.npz'
                scan     = np.load(fname)
                data     = list(scan['data'])
                motors   = list(scan['motors'])
                counters = scan['counters'].item()
                edfmats  = scan['edfmats']
                return data, motors, counters, edfmats
            except:
                print( 'Failed loading scan from file, will read EDF- and SPEC-file.' )
                pass

        # proceed with parsing edf- and SPEC-files, if scan could not be loaded from zip-archive
        if not fromtofile:
            print( 'Parsing EDF- and SPEC-files of scan No. %s' % scannumber )

        # load SPEC-file
        fn = self.path + self.filename
        # print( " READING ", fn, scannumber)
        data, motors, counters = xrs_utilities.specread(fn,scannumber)

        if not self.single_image:
            # initiate arrays for the edf-files
            edfmatsh = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy/2,self.DET_PIXEL_NUMx)))
            edfmatsv = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy/2,self.DET_PIXEL_NUMx)))
            edfmats  = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy,self.DET_PIXEL_NUMx)))
            # load edf-files
            for m in range(len(counters['ccdno'])):
                ccdnumber = counters['ccdno'][m]
                edfnameh   = self.path + self.EDF_PREFIXh + self.edfName + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
                edfnamev   = self.path + self.EDF_PREFIXv + self.edfName + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
                edfmatsh[m,:,:] = (xrs_utilities.edfread_test(edfnameh))
                edfmatsv[m,:,:] = (xrs_utilities.edfread_test(edfnamev))

                edfmats[m,0:self.DET_PIXEL_NUMy/2,:] = edfmatsv[m,:,:]
                edfmats[m,self.DET_PIXEL_NUMy/2:,:]  = edfmatsh[m,:,:]
        else:
            # initiate arrays for the edf-files
            #edfmats  = np.array(np.zeros((len(counters['ccdno']),self.DET_PIXEL_NUMy,self.DET_PIXEL_NUMx)))
            #edfmats = np.array([])
            # load edf-files
                        # summa=0
            for m in range(len(counters['ccdno'])):
                ccdnumber = counters['ccdno'][m]
                edfname   = self.path + self.EDF_PREFIX + self.edfName + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
                                # print " LEGGO ", edfname
                if m == 0:
                    edfshape = xrs_fileIO.EdfRead(edfname).shape
                    edfmats = np.zeros((  len(counters['ccdno']), edfshape[0] ,  edfshape[1] ))
                    edfmats[m,:,:] = (xrs_fileIO.EdfRead(edfname))
                else:
                    edfmats[m,:,:] = (xrs_fileIO.EdfRead(edfname))
                                # summa=summa+xrs_fileIO.EdfRead(edfname)
                        #print " SHAPE ", edfmats.shape
                        #spia = edfmats.sum(axis=0)
                        #from pylab import *
                        #spia=summa
                        #imshow(log(spia))
                        #show()

        # add the scannumber to self.scannumbers, if not already present
        if not scannumber in self.scannumbers:
            self.scannumbers.extend([scannumber])

        # store scan in numpy zip-archive, if desired
        if fromtofile:
            scanname = 'Scan%03d' % scannumber
            sub_path = os.path.split(self.path[:-1])[0]
            fname = sub_path + '/scans/' + scanname
            if os.path.exists(sub_path):
                print( 'trying to save file in numpy-archive.' )
                np.savez(fname, data=data, motors=motors, counters=counters, edfmats=edfmats)
            else:
                print( 'please create a folder ' + fname + ' first!' )
                pass

        return data, motors, counters, edfmats #<type 'list'> <type 'list'> <type 'dict'> <type 'numpy.ndarray'>

    def loadscan(self,scannumbers,scantype='generic',fromtofile=False):
        """
        Loads the files belonging to scan No. "scannumber" and puts it into an instance
        of the xrs_scan-class 'scan'. The default scantype is 'generic', later the scans
        will be grouped (and added) based on the scantype.
        INPUT:
          * scannumbers = integer or list of scannumbers that should be loaded
          * scantype    = string describing the scan to be loaded (e.g. 'edge1' or 'K-edge')
          * fromtofile  = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        # make sure scannumbers are iterable (list)
        if not isinstance(scannumbers,list):
            scannums = []
            scannums.append(scannumbers)
        else:
            scannums = scannumbers
        for number in scannums:
            scanname = 'Scan%03d' % number
            data, motors, counters, edfmats = self.readscan(number,fromtofile)
            # can assign some things here already (even if maybe redundant)
            monitor   = counters[self.monicolumn]
            monoangle = 1 # have to check for this later ... !!! counters['pmonoa']
            energy    = counters[self.encolumn]
            # create an instance of "scan" class for every scan
            onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
            # assign one dictionary entry to each scan
            self.scans[scanname] = onescan

    def loadloop(self,begnums,numofregions,fromtofile=False):
        """
        Loads a whole loop of scans based on their starting scannumbers and the number of single
        scans in the loop.
        INPUT:
          * begnums      = list of scannumbers of the first scans of each loop (is a list)
          * numofregions = number of scans in each loop (integer)

        """
        typenames = []
        for n in range(numofregions):
            typenames.append('edge'+str(n+1))
        numbers = []
        for n in range(len(begnums)):
            for m in range(numofregions):
                numbers.append((begnums[n]+m))
        typenames = typenames*len(begnums)
        for n in range(len(typenames)):
            thenumber = []
            thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
            self.loadscan(thenumber,typenames[n],fromtofile)

    def loadelastic(self,scann,fromtofile=False):
        """
        Loads a scan using the loadscan function and sets the scantype attribute to 'elastic'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='elastic')'.
        INPUT:
          * scann      = integer or list of integers
          * fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        self.loadscan(scann,'elastic',fromtofile)

    def loadlong(self,scann,fromtofile=False):
        """
        Loads a scan using the loadscan function and sets the scantype attribute to 'long'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='long')'.
        INPUT:
          * scann      = integer or list of integers
          * fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        self.loadscan(scann,'long',fromtofile)


    def set_roiObj(self,roiobj):
        self.roi_obj = roiobj



    def orderrois(self,arrangement='vertical',missing=None):
        """
        order the rois in an order provided such that e.g. autorois have the correct order
        """
        if not self.roi_obj:
            print( 'Please select some ROIs first!' )
            return

        # get some detector infos
        all_det_names = ['VD','VU','VB','HR','HL','HB']
        all_dets      = []
        for name in all_det_names:
            all_dets.append(xrs_utilities.maxipix_det(name,arrangement))

        # find the ROIs for each detector
        # go through all detectors and ROIs defined and see which one has centers in which detector
        for det in all_dets:
            det_x_min = det.get_pixel_range()[0]
            det_x_max = det.get_pixel_range()[1]
            det_y_min = det.get_pixel_range()[2]
            det_y_max = det.get_pixel_range()[3]
            for ii in range(len(self.roi_obj.indices)):
                x_mean = np.mean(self.roi_obj.x_indices[ii])
                y_mean = np.mean(self.roi_obj.y_indices[ii])
                if x_mean >= det_x_min and x_mean <= det_x_max and y_mean >= det_y_min and y_mean <= det_y_max:
                    det.roi_indices.append(self.roi_obj.indices[ii])
                    det.roi_x_indices.append(self.roi_obj.x_indices[ii])
                    det.roi_y_indices.append(self.roi_obj.y_indices[ii])
                    det.roi_x_means.append(x_mean)
                    det.roi_y_means.append(y_mean)

        # count, check with missing, break if something is wrong
        for det in all_dets:
            if not len(det.roi_indices) == 12:
                print( 'WARNING! Module ' + det.name + ' only has ' + '%d' % len(det.roi_indices) + ' ROIs defined, your numbering will be messed up!' )
        # rearrange according to 'arrangement' keyword
        if arrangement == 'vertical':
            verticalIndex = [0,3,6,9,1,4,7,10,2,5,8,11]
            for det in all_dets:
                # order from low to high y-mean
                roi_coords   = np.array(det.roi_indices)
                roi_y_means  = np.array(det.roi_x_means)
                sorting_inds = roi_y_means.argsort()
                roi_coords_increasing = [roi_coords[i] for i in sorting_inds]
                if det.get_det_name() in ['VD','VU','VB']:
                    # sort from high to low y-center
                    det.roi_indices = [roi_coords_increasing[i] for i in verticalIndex[::-1]]
                elif det.get_det_name() in ['HR','HL','HB']:
                    # sort from low to high y-center
                    det.roi_indices = [roi_coords_increasing[i] for i in verticalIndex]
                else:
                    print( 'Sorry, no such module.' )

        # reassign all roi_obj variables...
        allrois = []
        for det in all_dets:
            allrois.extend(det.roi_indices)

        self.roi_obj.indices = allrois
        self.roi_obj.roi_matrix     = xrs_rois.convert_inds_to_matrix(self.roi_obj.indices,(512,768))
        self.roi_obj.red_rois       = xrs_rois.convert_matrix_to_redmatrix(self.roi_obj.roi_matrix)
        self.roi_obj.x_indices      = xrs_rois.convert_inds_to_xinds(self.roi_obj.indices)
        self.roi_obj.y_indices      = xrs_rois.convert_inds_to_yinds(self.roi_obj.indices)
        self.roi_obj.masks          = xrs_rois.convert_roi_matrix_to_masks(self.roi_obj.roi_matrix)
        self.roi_obj.number_of_rois = np.amax(self.roi_obj.roi_matrix)

    def getrawdata(self):
        """
        Goes through all instances of the scan class and calls it's applyrois method
        to sum up over all rois.
        """
        if not np.any(self.roi_obj.indices):
            print( 'Please define some ROIs first.' )
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
            print( 'Please define some ROIs first.' )
            return
        for scan in self.scans:
            if len(self.scans[scan].edfmats):
                print ("integrating pixelwise "+scan)
                self.scans[scan].applyrois_pw(self.roi_obj.indices)


    def SumDirect(self,scannumbers):
        sum=None
        for number in scannumbers:
            data, motors, counters, edfmats = self.readscan(number)
            if sum is None:
                sum = np.zeros(edfmats[0].shape ,"f")
            sum[:] += edfmats.sum(axis=0)
        return sum

    def loadscandirect(self,scannumbers,scantype='generic',fromtofile=False,scaling=None):
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
            print( 'Please define some ROIs first' )
            return
        for number in scannums:
            scanname = 'Scan%03d' % number
            data, motors, counters, edfmats = self.readscan(number,fromtofile)
            # can assign some things here already (even if maybe redundant)
            monitor   = counters[self.monicolumn]
            monoangle = 1 # counters['pmonoa'] # this still needs checking
            energy    = counters[self.encolumn]
            # create an instance of "scan" class for every scan

            onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)

            onescan.applyrois(self.roi_obj.indices,scaling=scaling)

            print( 'Deleting -- EDF-files of Scan No. %03d' % number )
            onescan.edfmats = [] # delete the edfmats
            self.scans[scanname] = onescan

    def loadloopdirect(self,begnums,numofregions,fromtofile=False,scaling=None):
        """
        Loads a whole loop of scans based on their starting scannumbers and the number of single
        INPUT:
          * begnums      = list of scannumbers of the first scans of each loop (is a list)
          * numofregions = number of scans in each loop (integer)
          * fromtofile   = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)
          * scaling      = list of scaling factors to be applied, one for each ROI defined

        """
        typenames = []
        for n in range(numofregions):
            typenames.append('edge'+str(n+1))
        numbers = []
        for n in range(len(begnums)):
            for m in range(numofregions):
                numbers.append((begnums[n]+m))
        typenames = typenames*len(begnums)
        for n in range(len(typenames)):
            thenumber = []
            thenumber.append(numbers[n]) # make the scannumber into an interable list, sorry it's messy
            self.loadscandirect(thenumber,typenames[n],fromtofile,scaling=scaling)

    def loadelasticdirect(self,scann,fromtofile=False):
        """
        Loads a scan using the loadscan function and sets the scantype attribute to 'elastic'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='elastic')'.
        INPUT:
          * scann      = integer or list of integers
          * fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        self.loadscandirect(scann,'elastic',fromtofile)

    def loadlongdirect(self,scann,fromtofile=False,scaling=None):
        """
        Loads a scan using the loadscan function and sets the scantype attribute to 'long'.
        I.e. shorthand for 'obj.loadscan(scannumber,type='long')'.
        INPUT:
          * scann      = integer or list of integers
          * fromtofile = boolean flag, 'True' if the scan should be saved in a pickle-file (this is developmental)

        """
        self.loadscandirect(scann,'long',fromtofile,scaling=scaling)

    def deletescan(self,scannumbers):
        """
        Deletes scans from the class.
        INPUT:
          * scannumbers = integer or list of integers (SPEC scan numbers) to delete

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

    def getspectrum(self, include_elastic=False, absCounts=False):
        """
        Groups the instances of the scan class by their scantype attribute,
        adds equal scans (each group of equal scans) and appends them.
        INPUT:
          * include_elastic = boolean flag, skips the elastic line if set to 'False' (default)

        """
        # find the groups
        allgroups = xrs_scans.findgroups(self.scans)
        for group in allgroups:
            # self.makegroup(group)
            onegroup = xrs_scans.makegroup_nointerp(group, absCounts=absCounts)
            self.groups[onegroup.get_type()] = onegroup
        self.energy,self.signals,self.errors = xrs_scans.appendScans(self.groups,include_elastic)
        self.signals_orig = self.signals
        self.errors_orig  = self.errors

    def geteloss(self):
        """
        Defines the energy loss scale for all ROIs by finding the center of mass for each ROI's elastic line.
        Interpolates the signals and errors onto a commom energy loss scale. Finds the resolution (FWHM) of the
        'elastic' groups.
        """
        if not 'elastic' in self.groups:
            print( 'Please load/integrate at least one elastic scan first!' )
            return
        else:
            # reset values, in case this is run several times
            self.cenom = []
            self.resolution = []
            Origin=None
            valid_cenoms = []

            for n in range(len(self.roi_obj.indices)):
                # find the center of mass for each ROI

                cofm = xrs_utilities.find_center_of_mass(self.groups['elastic'].energy,self.groups['elastic'].signals_orig[:,n])
                self.cenom.append(cofm)
                if self.there_is_a_valid_roi_at(n)  :
                    valid_cenoms.append(cofm)
                    if Origin is None :
                        Origin = cofm

                # try finden the FWHM/resolution for each ROI
                FWHM,x0 = xrs_utilities.fwhm((self.groups['elastic'].energy - self.cenom[n])*1e3,self.groups['elastic'].signals_orig[:,n])
                # try:
                # 	FWHM,x0 = xrs_utilities.fwhm((self.groups['elastic'].energy - self.cenom[n])*1e3,self.groups['elastic'].signals_orig[:,n])
                # 	self.resolution.append(FWHM)
                # # append a zero if the FWHM routine fails
                # except:
                #     exc_type, exc_value, exc_traceback = sys.exc_info()
                #     print "*** print_tb:"
                #     traceback.print_tb(exc_traceback, limit=None, file=sys.stdout)
                #     print " need a more sofisticated way of finding the FWHM "
                #     self.resolution.append(0.0) # need a more sofisticated way of finding the FWHM
            self.E0 = np.mean(valid_cenoms)
            # define the first eloss scale as the 'master' scale for all ROIs
            self.eloss = (self.energy - cofm)*1e3 # energy loss in eV
            # define eloss-scale for each ROI and interpolate onto the 'master' eloss-scale
            for n in range(len(self.roi_obj.indices)):
                # inserting zeros at beginning and end of the vectors to avoid interpolation errors
                x = (self.energy-self.cenom[n])*1e3
                x = np.insert(x,0,-1e10)
                x = np.insert(x,-1,1e10)
                y = self.signals_orig[:,n]
                y = np.insert(y,0,0)
                y = np.insert(y,-1,0)
                f = interp1d(x,y, bounds_error=False,fill_value=0.0)
                self.signals[:,n] = f(self.eloss)
            # do the same for the errors
            for n in range(len(self.roi_obj.indices)):
                # inserting zeros at beginning and end of the vectors to avoid interpolation errors
                x = (self.energy-self.cenom[n])*1e3
                x = np.insert(x,0,-1e10)
                x = np.insert(x,-1,1e10)
                y = self.errors_orig[:,n]
                y = np.insert(y,0,0)
                y = np.insert(y,-1,0)
                f = interp1d(x,y, bounds_error=False,fill_value=0.0)
                self.errors[:,n] = f(self.eloss)

    def gettths(self,rvd=0.0,rvu=0.0,rvb=0.0,rhl=0.0,rhr=0.0,rhb=0.0,order=[0,1,2,3,4,5]):
        """
        Uses the defined TT_OFFSETS of the read_id20 class to set all scattering angles tth
        from the mean angle avtth of the analyzer modules.
        INPUT:
           * rhl = mean tth angle of HL module (default is 0.0)
           * rhr = mean tth angle of HR module (default is 0.0)
           * rhb = mean tth angle of HB module (default is 0.0)
           * rvd = mean tth angle of VD module (default is 0.0)
           * rvu = mean tth angle of VU module (default is 0.0)
           * rvb = mean tth angle of VB module (default is 0.0)
           * order = list of integers (0-5) which describes the order of modules in which the
             ROIs were defined (default is VD, VU, VB, HR, HL, HB; i.e. [0,1,2,3,4,5])


        """
        # reset all values, just in case mean angles are redefined (values are otherwise appended to existing values)
        self.VDtth = []
        self.VUtth = []
        self.VBtth = []
        self.HRtth = []
        self.HLtth = []
        self.HBtth = []
        self.tth   = []
        # horizontal modules
        # HL (motor name rhl)
        v_angles = self.TTH_OFFSETS1
        h_angles = self.TTH_OFFSETS2 + rhl
        for n in range(len(h_angles)):
            self.HLtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)
        # HR (motor name rhr)
        v_angles = self.TTH_OFFSETS1
        h_angles = self.TTH_OFFSETS2 + rhr
        for n in range(len(h_angles)):
            self.HRtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)
        # HB (motor name rhb)
        v_angles = self.TTH_OFFSETS1
        h_angles = self.TTH_OFFSETS2 + rhb
        for n in range(len(h_angles)):
            self.HBtth.append( np.arccos(np.cos(np.radians(h_angles[n]))*np.cos(np.radians(v_angles[n])))*180.0/np.pi)

        # vertical modules
        # VD
        v_angles = self.TTH_OFFSETS2 + rvd
        h_angles = self.TTH_OFFSETS1
        for n in range(len(h_angles)):
            self.VDtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)
        # VU
        v_angles = self.TTH_OFFSETS2 + rvu
        h_angles = self.TTH_OFFSETS1
        for n in range(len(h_angles)):
            self.VUtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)
        # VB
        v_angles = self.TTH_OFFSETS2 + rvb
        h_angles = self.TTH_OFFSETS1
        for n in range(len(h_angles)):
            self.VBtth.append( np.arccos(np.cos(np.radians(v_angles[n]))*np.cos(np.radians(h_angles[n])))*180.0/np.pi)

        # list of TTH values
        tth = [self.VDtth, self.VUtth, self.VBtth, self.HRtth, self.HLtth, self.HBtth]
        # list all TTH values in one long list ordered by the 'order'-keyword
        for n in order:
            self.tth.extend(tth[n])

    def getqvals(self,invangstr=False):
        """
        Calculates q-values from E0 and tth values in either atomic units (defalt) or
        inverse angstroms.
        """
        theqvals = np.zeros_like(self.signals)
        if invangstr:
            for n in range(len(self.signals[0,:])):
                theqvals[:,n] = xrs_utilities.momtrans_inva(self.E0+self.eloss/1e3,self.E0,self.tth[n])
        else:
            for n in range(len(self.signals[0,:])):
                theqvals[:,n] = xrs_utilities.momtrans_au(self.E0+self.eloss/1e3,self.E0,self.tth[n])
        self.qvalues = theqvals

    def getqvals_energy(self,energy):
        """
        Returns all q-values at a certain energy loss.
        INPUT:
          * energy = energy loss value for which all q-values are stored

        """
        ind = np.abs(self.eloss - energy).argmin()
        return self.qvalues[ind,:]

    def copy_edf_files(self,scannumbers,destdir):
        """
        Copies all edf-files from scan with scannumber or scannumbers into directory 'destdir'
        INPUT:
          * scannumbers = integer or list of integers defining the scannumbers from the SPEC file
          * destdir     = string with absolute path for the destination

        """
        import shutil
        numbers = []
        if not isinstance(scannumbers,list):
            numbers.append(scannumbers)
        else:
            numbers = scannumbers
        fn = self.path + self.filename
        if not self.single_image:
            for n in range(len(numbers)):
                data, motors, counters = xrs_utilities.specread(fn,numbers[n])
                for m in range(len(counters['ccdno'])):
                    ccdnumber = counters['ccdno'][m]
                    edfnameh   = self.path + self.EDF_PREFIXh + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
                    edfnamev   = self.path + self.EDF_PREFIXv + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
                    shutil.copy2(edfnameh, destdir)
                    shutil.copy2(edfnamev, destdir)
        if self.single_image:
            for n in range(len(numbers)):
                data, motors, counters = xrs_utilities.specread(fn,numbers[n])
                for m in range(len(counters['ccdno'])):
                    ccdnumber = counters['ccdno'][m]
                    edfname   = self.path + self.EDF_PREFIX + self.filename + '_' + "%04d" % ccdnumber + self.EDF_POSTFIX
                    shutil.copy2(edfname, destdir)


    def printlength(self,scannumbers):
        """
        Prints the number of energy points in a scan or a number of scans.
        INPUT:
          * scannumbers = integer or list of integers

        """
        numbers = []
        if not isinstance(scannumbers,list):
            numbers.append(scannumbers)
        else:
            numbers = scannumbers

        for i in numbers:
            name = 'Scan%03d' % i
            print( 'length of scan %03d ' %i + ' is ' + str(len(self.scans[name].energy)))

    def removeBackgroundRoi(self,backroinum,estart=None,estop=None):
        if not estart:
            estart = self.eloss[0]
        if not estop:
            estop = self.elsos[-1]
        if not self.tth:
            print( 'Please define the scattering angles first using the gettth method!' )
            return

        for ii in range(len(self.tth)):
            if ii != backroinum:
                inds       = np.where(np.logical_and(self.eloss>=estart, mpc96.eloss<=estop))
                expnorm    = np.trapz(self.signals[inds,ii],self.eloss[inds])
                backnorm   = np.trapz(self.signals[inds,backroinum],self.eloss[inds])
                background = self.signals[:,backroinum]/backnorm*expnorm
                self.signals[:,ii] -= background # subtract background roi

    def save_raw_data(self,filename):
        data = np.zeros((len(self.eloss),len(self.signals[0,:])))
        data[:,0]   = self.eloss
        data[:,1::] = self.signals
        np.savetxt(filename,data)

def animation(id20read_object,scannumber,logscaling=True,timeout=-1,colormap='jet'):
    """
    Shows the edf-files of a scan as a 'movie'.
    INPUT:
      * scannumber = integer/scannumber
      * logscaling = set to 'True' (default) if edf-images are to be shown on logarithmic-scale
      * timeout    = time in seconds defining pause between two images, if negative (default)
        images are renewed by mouse clicks
      * colormap   = matplotlib color scheme used in the display

    """
    if isinstance(scannumber,list):
        if len(scannumber)>1:
            print( 'this only works for a single scan, sorry' )
            return
        else:
            scannumber = scannumber[0]

    scanname = 'Scan%03d' % scannumber
    edfmats  = id20read_object.scans[scanname].edfmats
    scanlen  = np.shape(edfmats)[0]
    plt.ion()
    plt.clf()
    for n in range(scanlen):
        plt.clf()
        if logscaling:
            theimage = plt.imshow(np.log(edfmats[n,:,:]))
        else:
            theimage = plt.imshow(edfmats[n,:,:])
        plt.xlabel('detector x-axis [pixel]')
        plt.ylabel('detector y-axis [pixel]')
        if timeout<0:
            titlestring = 'Frame No. %d' % (n+1) + ' of %d' % scanlen + ', press key or mouse botton to continue'
            plt.title(titlestring)
        else:
            titlestring = 'Frame No. %d' % (n+1) + ' of %d' % scanlen + ', updating every %2.2f ' % timeout + ' seconds'
            plt.title(titlestring)
        theimage.set_cmap(colormap)
        plt.draw()
        plt.waitforbuttonpress(timeout=timeout)

def alignment_image(id20read_object,scannumber,motorname,filename=None):
    """
    Loads a scan from a sample position scan (x-scan, y-scan, z-scan), lets you choose a zoomroi and constructs a 2D image from this
    INPUT:
      * scannumber = number of the scan
      * motorname  = string that contains the motor name (must be the same as in the SPEC file)
      * filename   = optional parameter with filename to store the image


    """
    # load the scan
    data, motors, counters, edfmats = id20read_object.readscan(scannumber)

    # the scan motor
    position = counters[motorname.lower()]

    # define a zoom ROI
    image = xrs_utilities.sumx(edfmats)
    roi_finder_obj = roifinder_and_gui.roi_finder()
    roi_finder_obj.get_zoom_rois(image)

    # construct the image
    roixinds = roi_finder_obj.roi_obj.x_indices[0]
    roiyinds = roi_finder_obj.roi_obj.y_indices[0]

    # go through all edf files of the scan, sum over the height of the roi and stack the resulting lines into a matrix
    axesrange = [0,roiyinds[-1],position[-1],position[0]]
    theimage  = (np.sum(edfmats[:,np.amin(roixinds):np.amax(roixinds)+1,np.amin(roiyinds):np.amax(roiyinds)+1],axis=1))
    plt.close()
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.imshow(np.log(theimage),extent=axesrange)
    ax.set_aspect('auto')
    plt.xlabel('pixel along the beam')
    ylabelstr = motorname.lower() + ' position [mm]'
    plt.ylabel(ylabelstr)
    plt.show()

    # save the image, if a filename is provided
    if filename:
        from .xrs_imaging import LRimage
        f = open(filename, 'wb')
        yrange = np.arange(np.amin(roixinds),np.amax(roixinds)+1)
        theobject = LRimage(theimage, position, yrange)
        pickle.dump(theobject, f, protocol=-1)
        f.close()

def alignment_image_old(so,scan_number,motorname):
    """
    Loads a scan from a sample position scan (x-scan, y-scan, z-scan), lets you choose a zoomroi and constructs a 2D image from this
    INPUT:
      * scannumber = number of the scan
      * motorname  = string that contains the motor name (must be the same as in the SPEC file)
      * filename   = optional parameter with filename to store the image
    """
    # load the scan
    scan = xrs_scans.Scan()
    scan.load(so.path, so.SPECfname, so.EDFprefix, so.EDFname, so.EDFpostfix, scan_number)
    motors   = scan.motors
    counters = scan.counters
    edfmats  = scan.edfmats

    # the scan motor
    position = counters[motorname.lower()]

    # define a zoom ROI
    image = xrs_utilities.sumx(edfmats)
    roi_finder_obj = roifinder_and_gui.roi_finder()
    roi_finder_obj.get_zoom_rois(image)

    # construct the image
    roixinds = roi_finder_obj.roi_obj.x_indices[0]
    roiyinds = roi_finder_obj.roi_obj.y_indices[0]

    # go through all edf files of the scan, sum over the height of the roi and stack the resulting lines into a matrix
    axesrange = [0,roiyinds[-1],position[-1],position[0]]
    theimage  = (np.sum(edfmats[:,np.amin(roixinds):np.amax(roixinds)+1,np.amin(roiyinds):np.amax(roiyinds)+1],axis=1))
    plt.close()
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.imshow(np.log(theimage),extent=axesrange)
    ax.set_aspect('auto')
    plt.xlabel('pixel along the beam')
    ylabelstr = motorname.lower() + ' position [mm]'
    plt.ylabel(ylabelstr)
    plt.show()

def alignment_image_new( so, scan_number, log_scaling=True, cmap='Blues',
                             interpolation='nearest'):
    """
    Loads a scan from a sample position scan (x-scan, y-scan, z-scan), lets
    you choose a zoomroi and constructs a 2D image from this

    INPUT:
    so         = Hydra object
    scannumber = number of the scan

    """
    # load the scan
    scan = xrs_scans.Scan()
    scan.load(so.path, so.SPECfname, so.EDFprefix, so.EDFname, so.EDFpostfix, scan_number)
    motors   = scan.motors
    counters = scan.counters
    edfmats  = scan.edfmats

    # the scan motor
    position = scan.energy

    # define a zoom ROI
    image = xrs_utilities.sumx(edfmats)
    roi_finder = roifinder_and_gui.roi_finder()
    roi_finder.get_zoom_rois(image)

    # get the image
    scan.get_raw_signals( roi_finder.roi_obj, method='column' )
    if log_scaling:
        im       = np.log(scan.raw_signals['ROI00'])
        im_orig  = scan.raw_signals['ROI00']
    else:
        im       = scan.raw_signals['ROI00']

    im_shape = im.shape

    # plot the image
    axes_range = [0,im.shape[1], position[-1], position[0]]
    plt.close()
    fig = plt.figure( figsize=(7, 20) )

    # main plot
    ax0 = plt.subplot2grid((16, 16), (0, 0), colspan=12, rowspan=12)
    ax0.imshow( im, extent=axes_range, cmap=cmap, interpolation=interpolation )
    ax0.set_ylabel( '%s position [mm]'%scan.scan_motor.lower() )
    ax0.set_aspect( 'auto' )
    ax0.tick_params( left=True,
                    labelleft=True,
                    right=True,
                    bottom=True,
                    top=True,
                    labelbottom=False,
                    direction='in')

    # projection onto y
    ax1 = plt.subplot2grid((16, 16), (0, 12), colspan=4, rowspan=12)
    ax1.plot( np.sum(im_orig,axis=1), np.flipud(position), '-k' )
    ax1.set_ylim(ax1.get_ylim()[::-1])
    ax1.yaxis.set_label_position("right")
    ax1.set_ylabel( '%s position [mm]'%scan.scan_motor.lower() )
    ax1.tick_params( left=True,
                    labelleft=False,
                    right=True,
                    labelright=True,
                    bottom=True,
                    top=True,
                    labelbottom=False,
                    direction='in')

    # projection onto x
    ax2 = plt.subplot2grid((16, 16), (12, 0), colspan=12, rowspan=4)
    ax2.plot( np.sum(im_orig,axis=0), '-k' )
    ax2.set_xlim([ 0, im.shape[1]])
    ax2.tick_params( left=True,
                    labelleft=False,
                    right=True,
                    labelright=False,
                    bottom=True,
                    top=True,
                    labelbottom=True,
                    direction='in')
    ax2.set_xlabel( 'direction along beam [pixel]' )

    plt.show()



def get_scans_pw(id20read_object,scannumbers):
    """	**get_scans_pw**
    Sums scans from pixelwise ROI integration for use in the PW roi refinement.
    """
    if isinstance(scannumbers,list):
        scannums = scannumbers
    elif isinstance(scannumbers,int):
        scannums = [scannumbers]
    else:
        print('Please provide keyword \'scannumbers\' as integer or list of integers.')
        return
    if len(scannums)==1:
        scanname = 'Scan%03d' % scannums[0]
        pw_matrices = id20read_object.scans[scanname].signals_pw
        # normalize data
        pw_matrices_norm = []
        for matrix in pw_matrices:
            for ii in range(matrix.shape[1]):
                matrix[:,ii] /= id20read_object.scans[scanname].monitor
            pw_matrices_norm.append(matrix)
    else:
        scanname = 'Scan%03d' % scannums[0]
        pw_matrices = id20read_object.scans[scanname].signals_pw
        # normalize data
        pw_matrices_norm = []
        for matrix in pw_matrices:
            for ii in range(matrix.shape[1]):
                matrix[:,ii] /= id20read_object.scans[scanname].monitor
            pw_matrices_norm.append(matrix)

        for ii in scannums[1:]:
            scanname = 'Scan%03d' % ii
            for jj in range(len(pw_matrices)):
                pw_matrix = id20read_object.scans[scanname].signals_pw[jj]
                pw_matrices[jj]      += pw_matrix
                for kk in range(pw_matrix.shape[1]):
                    pw_matrix[:,kk] /= id20read_object.scans[scanname].monitor
                pw_matrices_norm[jj] += pw_matrix

    return pw_matrices_norm

class read_lerix:
    def __init__(self,exp_dir,elastic_name='elastic',nixs_name='nixs',wide_name='wide',energycolumn='energy',monitorcolumn='__enc'):
        self.scans         = {} # was a dictionary before
        self.path          = os.path.split(exp_dir)[0] + '/'
        self.monicolumn    = monitorcolumn
        self.encolumn      = energycolumn
        self.nixs_name, self.nixs_scans       = str.lower(nixs_name), []
        self.wide_name, self.wide_scans       = str.lower(wide_name), []
        self.elastic_name, self.elastic_scans = str.lower(elastic_name), []
        #split scans into NIXS and elastic and begin instance of XRStools scan class for each scan
        if not self.isValidDir(exp_dir):
            raise Exception('IOError! No such file, please check filename.')
        for file in self.sort_dir(self.path):
                scan_info = self.scan_info(file)
                if scan_info[2]=='elastic':
                    self.elastic_scans.append(file)
                if scan_info[2]=='nixs':
                    self.nixs_scans.append(file)
                if scan_info[2]=='wide':
                    self.wide_scans.append(file)
                else:
                    continue
        self.header_attrs  = {} #dictionary of useful information from scan files, inc. e0, comments, scan_time+date
        self.key           = {'Analyzer01':0, 'Analyzer02':1, 'Analyzer03':2,'Analyzer04':3,'Analyzer05':4,
                            'Analyzer06':5,'Analyzer07':6,'Analyzer08':7,'Analyzer09':8,'Analyzer10':9,'Analyzer11':10,'Analyzer12':11
                            ,'Analyzer13':12,'Analyzer14':13,'Analyzer15':14,'Analyzer16':15,'Analyzer17':16,'Analyzer18':17,'Analyzer19':18}
        self.scannumbers   = []
        self.scan_name     = []
        self.sample_name   = []
        self.eloss_avg     = [] #elastic eloss average
        self.signals_avg   = [] #elastic signals average used to plot analyzer resolutions at the end
        self.energy        = []
        self.signals       = []
        self.errors        = []
        self.is_checked    = [] #inserted later to save a list of the chosen analyzers after using .plot_data() save function
        self.tth           = []
        self.resolution    = {}
        self.E0            = []
        self.cenom         = []
        self.cenom_dict    = {}
        self.data          = {}
        #make cenom_dict which holds e0 and fwhm for all analyser
        for key in self.key.keys():
            self.cenom_dict[key] = {}
            self.data[key]       = {}
    ################################################################################
    # Get Ascii Info - parse a file and return the key details
    ################################################################################
    def getfloats(self, txt, allow_times=True):
        """
        function goes through a line and returns the line as a list of strings
        """
        words = [w.strip() for w in txt.replace(',', ' ').split()]
        # mktime = time.mktime
        for i, w in enumerate(words):
            val = None
            try:
                val = float(w)
            except ValueError:
                pass
        words[i] = val
        return(words)

    def colname(self, txt):
        """Function to replace bad characters with '_''s making a line of strings
        easier to handle."""
        return self.fixName(txt.strip().lower()).replace('.', '_')

    def isValidName(self, filename):
        """Function checks that a filename isn't in the list of reserved pythonic
        words. Returns corrected name or False"""
        if filename in RESERVED_WORDS:
            return False
        tnam = filename[:].lower()
        return NAME_MATCH(tnam) is not None

    def fixName(self, filename, allow_dot=True):
        if self.isValidName(filename):
            return filename
        if self.isValidName('_%s' % filename):
            return '_%s' % filename
        chars = []
        valid_chars = VALID_SNAME_CHARS
        if allow_dot:
            valid_chars = VALID_NAME_CHARS
        for s in filename:
            if s not in valid_chars:
                s = '_'
            chars.append(s)
        filename = ''.join(chars)
        # last check (name may begin with a number or .)
        if not self.isValidName(filename):
            filename = '_%s' % filename
        return filename

    def strip_headers(self, headers):
        #reorganise the headers and remove superfluous lines and commentchars
        header = []
        for hline in headers:
            hline = hline.strip().replace('\t', ' ')
            if len(hline) < 1:
                continue
            if hline[0] in COMMENTCHARS:
                hline = hline[1:].lstrip() #assumes reading l2r
            if len(hline) <1:
                continue
            header.append(hline)
        return(header)

    def separate_infile(self, text):
        """Function parses an 20ID ASCII file in reverse and separates it into
        headers, footers and data"""
        _labelline = None
        ncol = None
        dat, footers, headers = [], [], []
        try:
            text.reverse()
        except:
            text[::-1]
        section = 'FOOTER'
        for line in text:
            line = line.strip()
            if len(line) < 1: #remove any blank lines
                continue
            if section == 'FOOTER' and not None in self.getfloats(line):
                section = 'DATA'
            elif section == 'DATA' and None in self.getfloats(line):
                section = 'HEADER'
                _labelline = line
                if _labelline[0] in COMMENTCHARS:
                    _labelline = _labelline[1:].strip()
            if section == 'FOOTER': #reading footers but not using them currently
                footers.append(line)
            elif section == 'HEADER':
                headers.append(line)
            elif section == 'DATA':
                rowdat  = self.getfloats(line)
                if ncol is None:
                    ncol = len(rowdat)
                if ncol == len(rowdat):
                    dat.append(rowdat)
        return(headers, dat, footers)

    def pull_id20attrs(self, header):
        """Function takes headers of 20ID ASCII and parses it for key information
        for header_attrs - N.B. could be shortened by looping through a list of
        important key words rather than doing one by one."""
        bounds, steps, int_times = [], [], []
        header_attrs = {}
        line = -2
        #iterate through the header and pull out useful information and send it to header_attrs Dictionary
        for hhline in map(str.lower,header):
            line = line + 1 #counting to return the user comments which are on the next line
            try:
                if str(header[comment_line].strip()) == 'Scan config:':
                    header_attrs['User Comments'] = ""
                    pass
                else:
                    header_attrs['User Comments'] = str(header[comment_line].strip())
            except:
                pass
            if hhline.startswith('beamline'):
                words = hhline.split('beamline',1)
                header_attrs['beamline'] = str(words[1].strip())
            elif hhline.startswith('e0'):
                if ':' in hhline:
                    words = hhline.split(':',1)
                    header_attrs[words[0]] = float(words[1].strip(' ').split(' ',1)[0])
                elif '=' in hhline:
                    words = hhline.split('=',1)
                    header_attrs[words[0]] = float(words[1].strip(' ').split(' ',1)[0])
            elif hhline.startswith('user comment'):
                comment_line = line
            elif "scan time" in hhline:
                #search for scan date and time see: https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
                try:
                    words = hhline.split('scan time',1)
                    header_attrs['scan_time'] = datetime.strptime(words[1].strip(), '%H hrs %M min %S sec.').time()
                    header_attrs['scan_date'] = datetime.strptime(words[0].split('panel',1)[1].strip().strip(';'), '%m/%d/%Y  %I:%M:%S %p').date()
                except:
                    continue
            elif "scan bounds" in hhline:
                words = hhline.split('scan bounds',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        bounds.append(float(i))
                    except:
                        pass
                header_attrs['scan_bounds'] = bounds
            elif "scan step(s)" in hhline:
                words = hhline.split('scan step(s)',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        steps.append(float(i))
                    except:
                        pass
                header_attrs['scan_steps'] = steps
            elif "integration times" in hhline:
                words = hhline.split('integration times',1)
                for i in words[1].strip(':').split(' '):
                    try:
                        int_times.append(float(i))
                    except:
                        pass
                header_attrs['int_times'] = int_times
        return(header_attrs)

    def get_col_headers(self, header):
        col_headers = []
        for i in self.colname(header[0]).split('___'): #need three _ to work
            if not i:
                continue
            col_headers.append(i.strip('_'))
        return(col_headers)

    def scan_info(self, f):
        """get the scan number, name, type and file extention from the title of
        the scan assuming typical format e.g. elastic.0001, nixs.0001
        returns:
        [0] -> scan number (e.g. 0001)
        [1] -> scan name (e.g. nixs0001)
        [2] -> scan_type (e.g. nixs)
        [3] -> file name (e.g. nixs.0001)"""
        f = os.path.basename(f) #this allows both directories and files to be passed to get scan_info
        fn,fext = os.path.splitext(f)
        if str.lower(fn)==str.lower(self.nixs_name):
            scan_type = 'nixs'
        elif str.lower(fn)==str.lower(self.elastic_name):
            scan_type = 'elastic'
        elif str.lower(fn)==str.lower(self.wide_name):
            scan_type = 'wide'
        else:
            print(""">> LERIX >> WARNING  \n You have probably called the scan_info
            function without specifying a correct \n <class>.nixs/wide/elastic_name if you
            are calling scan_info manually - you can change this by setting:\n\n
            <class>.nixs_name = '<nixs_name>'""")
            sys.exit()
        scan_number = fext.lstrip('.')
        scan_number = int(scan_number)
        scan_name = scan_type + '%04d' %scan_number
        return(scan_number, scan_name, scan_type, f)

    def sort_dir(self, dir):
        """Returns a list of directory contents after filtering out scans without
        the correct format or size e.g. 'elastic.0001, nixs.0001 '"""
        dir_scans = []
        for file in os.listdir(dir):
            file_lc = str.lower(file)
            fn,fext = os.path.splitext(file_lc)
            if not file_lc.startswith('.'):
                    if fext.lstrip('.').isdigit():
                        if not os.stat(dir + '/' + file).st_size > 8000:
                            print("{} {}".format(">> >> Warning!! skipped empty scan (<8KB): ", file))
                            continue
                        elif not os.stat(dir + '/' + file).st_size < MAX_FILESIZE:
                            print("{} {}".format(">> >> Warning!! skipped huge scan (>100MB): ", file))
                            continue
                        else:
                            if fn==self.nixs_name:
                                dir_scans.append(file)
                            elif fn==self.elastic_name:
                                dir_scans.append(file)
                            elif fn==self.wide_name:
                                dir_scans.append(file)
        sorted_dir = sorted(dir_scans, key=lambda x: os.path.splitext(x)[1])
        return sorted_dir

    def isValidDir(self,dir):
        """Show that the scan directory is valid, that the directory holds a scan
        with the correct elastic name, nixs name and then let the user know if it
        has not found a wide scan. Returns True if valid directory."""
        if not os.path.isdir(dir):
            print('Check the directory you have supplied')
            return False
        elif not os.path.isfile(dir+'/'+self.elastic_name+'.0001'):
            print("The directory you supplied does not have a elastic.0001 file!!! \n If your elastic scan has a different name, please specify as: 'elastic_name'")
            return False
        elif not os.path.isfile(dir+'/'+self.nixs_name+'.0001'):
            print("The directory you supplied does not have a NIXS.0001 file!!! \n If your raman scan has a different name, please specify as: 'NIXS_name'")
            return False
        elif not os.path.isfile(dir+'/'+self.wide_name+'.0001'):
            print("No wide scans found. Continuing...")
            return True
        else:
            return True

    ################################################################################
    # Read Scan
    ################################################################################
    def update_cenom(self, analyzers="all"):
        """Internal Function to get the centre of mass of the elastic peak and
        the E0 for each elastic scan using XRStools"""
        self.cenom = []
        if not self.cenom_dict: #check that the cenom-dict is populated first.
            print('Cenom Dictionary is empty. Please load elastics first!')
            return(0)
        if analyzers == "all":
            print("Running 'Update_Cenom' Script for All analysers")
            analyzers = sorted(self.key.keys())
        elif type(analyzers) is list:
            print("Running 'Update_Cenom' Script for analysers:", analyzers)
            tmp = []
            for i in analyzers:
                i = i - 1
                tmp.append(sorted(self.key.keys())[i])
            analyzers = tmp
        else:
            print("list of analysers for E0 calculation must be a list type.")
            return(False)
        for analyzer in analyzers:
            avg_fwhm, avg_cenom = [],[]
            for scan in self.elastic_scans:
                scan = self.scan_info(scan)[1] # change from elastic.0001 to elastic0001
                avg_cenom.append(self.cenom_dict[analyzer][scan]['e0'])
                avg_fwhm.append(self.cenom_dict[analyzer][scan]['fwhm'])
            # annoying bit of code because nanmean() of a list of just nan's returns a warning, this step avoids ugly (yet harmless) warnings.
            if np.all(np.isnan(avg_cenom)):
                avg_cenom = np.nan
            else:
                avg_cenom = np.nanmean(np.array(avg_cenom)) # was just nanmean
            if np.all(np.isnan(avg_fwhm)):
                avg_fwhm = np.nan
            else:
                avg_fwhm = np.nanmean(np.array(avg_fwhm)) # was just nanmean
            self.cenom_dict[analyzer].update({'average': {'fwhm': avg_fwhm, 'e0': avg_cenom}})
            self.cenom.append(avg_cenom / 1e3) #divide by thousand to go from eV to keV
        self.E0 = np.nanmean(np.array(self.cenom)) # was just nanmean
        print("{} {}".format("E0 was found to be (keV): ", self.E0))
        print("{} {}".format("Average FWHM for the elastics is (eV): ", avg_fwhm))
        for i in range(len(analyzers)):
            if np.isnan(self.cenom[i]):
                print(analyzers[i], 'Elastic peak is less than 100 counts, setting to average e0')
                self.cenom[i] = self.E0
            else:
                continue
        #self.resolution['Resolution'] = round(np.mean(resolution),3)

    def readscan_20ID(self, file):
        """Read an ID20-type ASCII file and return header attributes and data as
        a dictionary. Takes a file path.
        header_attrs -> int_times, scan_steps, scan_bounds, e0, comments, beamline,
                        scan_time, scan_date
        data         -> dictionary of np.array (float64) with callable column names

        New XRStools uses scan class:
        # create an instance of "scan" class for every scan
        onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
        # assign one dictionary entry to each scan
        self.scans[scanname] = onescan
        edfmats: 2D images from pixel counting detector (not required here)
        number: scan number
        energy: energy axis
        monitor: i0 axis
        counters: interpreted as the detector channels
        motors: motor positions (not required)
        data: full data matrix
        scantype: elastic, nixs or long
        """
        import pandas as pd
        scan_info = self.scan_info(file)
        qixs_list = []
        f = open(file, "r") #read starts here
        text = f.read()
        text = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
        headers, dat, footers = self.separate_infile(text)
        try:
            dat = [map(list,zip(*dat))[i][::-1] for i in range(len(dat[1]))] # this function does the inverse ([::-1]) transposition of the dat object, doesn't seem to work in windows
        except:
            dat = [list(map(list,zip(*dat)))[i][::-1] for i in range(len(dat[1]))]
        names = self.get_col_headers(self.strip_headers(headers)) #returns a list of names in the order found in the data file.
        data = pd.DataFrame(np.array(dat).T,columns = np.array(names).T, dtype='float64') #returns a pandas array with the data arranged into labelled columns
        for column in sorted(data.columns): #sort by name so that analyzers are in correct (numerical) order
            if not column.rfind('i0') == -1:
                tmp_monitor = np.array(data[column].values)
            if not column.rfind('__enc') == -1:
                tmp_energy = np.array(data[column].values)
            if not column.rfind('qixs') == -1:
                qixs_list.append(column)
        tmp_signals = np.array(data[qixs_list].values)
        tmp_errors = np.sqrt(np.absolute(tmp_signals))
        scan_attrs = self.pull_id20attrs(self.strip_headers(headers)) #get scan_attrs
        if scan_info[2]=='elastic':
            for analyzer in self.key.keys(): #The analyzer channels in the scan ASCII
                # check counts are high-enough, using XIA filters avoids broadening FWHM
                if 200 <= np.max(tmp_signals[:,self.key[analyzer]]) <= 4000:
                    try:
                        fwhm, cenom = xrs_utilities.fwhm(tmp_energy,tmp_signals[:,self.key[analyzer]])
                    except:
                        print('Elastic scan is wrong shape, skipping')
                        fwhm, cenom = np.nan, np.nan
                else:
                    fwhm, cenom = np.nan, np.nan
                self.cenom_dict[analyzer].update({scan_info[1]: {'fwhm': fwhm, 'e0': cenom}})
        elif scan_info[2]=='nixs' or scan_info[2]=='wide':
            #create empty array with shape energy.v.signals
            eloss = np.zeros(tmp_signals.shape)
            self.tth = list(range(9,180,9)) #assign tth to self
            try:
                e_zero = self.E0 * 1e3 # convert e0 back to eV to perform subtraction
                tmp_eloss = np.subtract(tmp_energy,e_zero)
            except:
                print('>> No elastic <class>.E0 found! Make sure elastic scans have been loaded first')
                e_zero = scan_attrs['e0'] * 1e3 # convert e0 back to eV to perform subtraction
                tmp_eloss = np.subtract(tmp_energy,e_zero)
        # create an instance of "scan" class for every scan
        edfmats, motors = [],scan_attrs #no 2D pixel detector at LERIX
        number          = scan_info[0]
        energy          = tmp_energy
        monitor         = tmp_monitor
        counters        = tmp_signals
        data            = data
        scantype        = scan_info[2]
        onescan = xrs_scans.scan(edfmats,number,energy,monitor,counters,motors,data,scantype)
        self.scans[scan_info[1]] = onescan
        if scan_info[2]=='nixs' or scan_info[2]=='wide':
            self.scans[scan_info[1]].eloss = tmp_eloss
            self.scans[scan_info[1]].signals = np.divide(tmp_signals.T,monitor).T #transpose seems to be necessary, but don't know why?
            self.scans[scan_info[1]].errors = tmp_errors
        f.close()

    ################################################################################
    # Begin the reading
    ################################################################################
    def load_elastics(self,exp_dir=None,scans='all',analyzers='all'):
        """Function to load scan data from a typical APS 20ID Non-Resonant inelastic
        X-ray scattering experiment. With data in the form of elastic.0001, allign.0001
        and NIXS.0001. Function reteurns the averaged energy loss, signals, errors, E0
        and 2theta angles for the scans in the chosen directory."""
        scann = []
        if exp_dir is None:
            exp_dir = self.path
        if scans is 'all':
            chosen_scans = self.elastic_scans
        elif isinstance(scans,list):
            scann[:] = [x - 1 for x in scans] #scan 1 will be the 0th item in the list
            chosen_scans = [self.nixs_scans[i] for i in scann]
        else:
            print("scans must be list of scan numbers (e.g. [1,2,3]) or all")
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading elastic scan: ", file))
            self.readscan_20ID(exp_dir + '/' + file)
        self.update_cenom(analyzers)

    def load_nixs(self,exp_dir=None,scans='all',analyzers='all'):
        """Blah Blah"""
        scann = []
        if exp_dir is None:
            exp_dir = self.path
        if scans is 'all':
            chosen_scans = self.nixs_scans
        elif isinstance(scans,list):
            scann[:] = [x - 1 for x in scans] #scan 1 will be the 0th item in the list
            chosen_scans = [self.nixs_scans[i] for i in scann]
        else:
            print("scans must be list of scan numbers (e.g. [1,2,3]) or all")
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading NIXS scan: ", file))
            self.readscan_20ID(exp_dir + '/' + file)
        #average the data over the chosen scans
        self.energy   = np.array([self.scans[self.scan_info(i)[1]].energy  for i in chosen_scans]).mean(axis=0)
        self.signals  = np.array([self.scans[self.scan_info(i)[1]].signals for i in chosen_scans]).mean(axis=0)
        self.eloss    = np.array([self.scans[self.scan_info(i)[1]].eloss   for i in chosen_scans]).mean(axis=0)
        self.errors   = np.array([self.scans[self.scan_info(i)[1]].errors  for i in chosen_scans]).mean(axis=0)

    def load_wides(self,exp_dir=None,scans='all',analyzers='all'):
        """Blah Blah"""
        scann = []
        if exp_dir is None:
            exp_dir = self.path
        if scans is 'all':
            chosen_scans = self.wide_scans
        elif isinstance(scans,list):
            scann[:] = [x - 1 for x in scans] #scan 1 will be the 0th item in the list
            chosen_scans = [self.nixs_scans[i] for i in scann]
        else:
            print("scans must be list of scan numbers (e.g. [1,2,3]) or all")
        for file in chosen_scans:
            scan_info = self.scan_info(file)
            print("{} {}".format("Reading Wide scan: ", file))
            self.readscan_20ID(exp_dir + '/' + file)

    def plot_data(self,analyzer=False):
        """<classObj>.plot_data() Function that can be called to plot the eloss
        data for each channel and build an average by clicking a button.
        Does not require matplotlib >2.1"""
        import matplotlib.pyplot as plt
        from matplotlib.widgets import CheckButtons, Button, Cursor
        channels = []
        for analyzer in self.resolution:
            if analyzer.startswith('Analyzer'):
                if self.resolution[analyzer] < 1.0:
                    channels.append(int(analyzer.lstrip('Analyzer'))-1)
        data = np.average(self.signals[:,channels],axis=1)
        fig, ax = plt.subplots()
        ax.plot(self.eloss, data, lw=2)
        ax.set_xlabel('Energy Loss (eV)')
        ax.set_ylabel('S(q,w) [1/eV]')
        ax.set_title('Plotting Raman Analysers')
        plt.subplots_adjust(left=0.3)
        checkbuttonaxis = plt.axes([0.02, 0.15, 0.2, 0.8])
        anlabels, anvals = list(self.key), (False,)*len(list(self.key))
        anstates = dict(zip(anlabels,anvals))
        analyzers = CheckButtons(checkbuttonaxis, anlabels, anvals)
        buttonaxis = plt.axes([0.01, 0.01, 0.3, 0.09])
        bStatus  = Button(buttonaxis,'Save Averaged Analyzers')

        def onclick(label):
            """Tell the user what they have clicked - also good for de-bugging
            """
            anstates[label] = not anstates[label]
            print('un'*(not anstates[label]) + 'checked %s' %label)
            func()

        def savebutton(val):
            import pandas as pd
            import sys
            from PyQt4.QtGui import QApplication, QWidget, QFileDialog
            if not self.is_checked:
                print('please select your chosen analysers first!')
            else:
                print('selected analysers (python counting):  ', self.is_checked)
                save_signals = np.average(self.signals[:,self.is_checked],axis=1)
                save_errors = np.average(self.errors[:,self.is_checked],axis=1)
                df = pd.DataFrame(list(zip(self.eloss,save_signals,save_errors)), columns=['eloss','signals','errors'])
                print(df)
                try:
                    w = QWidget() #bug 'start a Qapp before a QPaintDevice py36 mac'
                    filename = str(QFileDialog.getSaveFileName(w, 'Save Analyzer Average','Result.csv'))
                    df.to_csv(filename,sep=',',na_rep='nan')
                    print('Saved as: ',filename)
                except:
                    print("{} {}".format(">> Warning >>", "file save was unsuccessful"))

        def func():
            ax.clear()
            self.is_checked = []
            for ii in anlabels:
                if anstates[ii]:
                    self.is_checked.append(self.key[ii])
            ax.plot(self.eloss, np.average(self.signals[:,self.is_checked],axis=1),lw=2)
            cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
            ax.autoscale(True)
            ax.set_xlabel('Energy Loss (eV)')
            ax.set_title('Plotting Raman Analysers')
            plt.draw()

        bStatus.on_clicked(savebutton)
        analyzers.on_clicked(onclick)
        cursor = Cursor(ax, useblit=False, color='red', linewidth=2)
        plt.show()

    def save_H5(self,H5name='20ID_APS_data.H5'):
        #if the user asks, call function to write all info to H5 file
        H5path = self.path + H5name
        if not os.path.isdir(os.path.dirname(H5path)):
            print('H5 path directory does not exist!')
        if os.path.isfile(H5path):
            H5file = h5py.File(saveloc, "a")
        else:
            H5file = h5py.File(saveloc, "w")
        g = H5file.create_group(H5name) #H5 subgroup with the name of the sample
        H5_ela = g.create_group('elastic') #H5 subgroup for elastics
        H5_xrs = g.create_group('XRS')     #H5 subgroup for NIXS
        all_scans = self.elastic_scans + self.nixs_scans + self.wide_scans
        for file in all_scans:
            scan_info = self.scan_info(file)
            if scan_info[2] == 'elastic':
                h5group = H5_ela.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                h5group.create_dataset("cenoms",data=self.scans[scan_info[1]].cenom)
            elif scan_info[2]=='nixs':
                h5group = H5_xrs.create_group(scan_info[1])
                h5group.create_dataset("energy",data=self.scans[scan_info[1]].energy)
                h5group.create_dataset("signals",data=self.scans[scan_info[1]].signals)
                h5group.create_dataset("eloss",data=self.scans[scan_info[1]].eloss)
                h5group.create_dataset("errors",data=self.scans[scan_info[1]].errors)
                #h5group.create_dataset("tth",data=self.scans[scan_info[1]].tth)
        g.create_dataset("energy",data=self.energy)
        g.create_dataset("signals",data=self.signals)
        g.create_dataset("eloss",data=self.eloss)
        g.create_dataset("errors",data=self.errors)
        g.create_dataset("tth",data=self.tth)
        g.create_dataset("Mean Resolutions", data=np.array(self.resolution.items()))
        #Never forget to close an open H5 file!!!
        H5file.close()
