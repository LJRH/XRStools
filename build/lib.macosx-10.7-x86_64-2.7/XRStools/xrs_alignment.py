from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import six
from six.moves import range
from six.moves import zip
#!/usr/bin/python
# Filename: xrs_alignment.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group and contains practical functions, 
# most of which are translated from Matlab functions from the University of
# Helsinki Electronic Structure Laboratory.
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

import numpy as np
from . import xrs_scans, xrs_read, roifinder_and_gui, math_functions
from scipy import interpolate, optimize
from matplotlib import pylab as plt

scan72_motornames = ['vdtx1','vdtx2','vdtx3','vdtx4','vdtx5','vdtx6','vdtx7','vdtx8','vdtx9','vdtx10','vdtx11','vdtx12', \
                     'vutx1','vutx2','vutx3','vutx4','vutx5','vutx6','vutx7','vutx8','vutx9','vutx10','vutx11','vutx12', \
                     'vbtx1','vbtx2','vbtx3','vbtx4','vbtx5','vbtx6','vbtx7','vbtx8','vbtx9','vbtx10','vbtx11','vbtx12', \
                     'hrtx1','hrtx2','hrtx3','hrtx4','hrtx5','hrtx6','hrtx7','hrtx8','hrtx9','hrtx10','hrtx11','hrtx12', \
                     'hltx1','hltx2','hltx3','hltx4','hltx5','hltx6','hltx7','hltx8','hltx9','hltx10','hltx11','hltx12', \
                     'hbtx1','hbtx2','hbtx3','hbtx4','hbtx5','hbtx6','hbtx7','hbtx8','hbtx9','hbtx10','hbtx11','hbtx12']

def optimize_analyzer_focus(path, SPECfname, EDFprefix, EDFname, EDFpostfix, roi_obj, scan_number):
    """Returns position for all 72 TX motors that optimize the analyzer foci.

    Args:
        path        (str): Absolute path to location of the SPEC-file.
        SPECfname   (str): SPEC-file name.
        EDFprefix   (str): Prefix to where EDF-files are stored.
        EDFname     (str): Base name of the EDF-files.
        EDFpostfix  (str): Post-fix used for the EDF-files.
        roi_obj (roi_obj): ROI object of the xrs_rois class.
        Scan_number (int): Scan number of the 72-motor scan.

    Returns:
        Dictionary with motorname - position pairs.

    """
    TX_positions = {}

    scan72 = xrs_scans.Scan()
    scan72.load(path, SPECfname, EDFprefix, EDFname, EDFpostfix, scan_number)

    edfmats     = scan72.edfmats
    for ii in range(len(scan72_motornames)):
        motor_scale = scan72.counters[scan72_motornames[ii]]
        try:
            TX_pos      = findAnalyzerFocus(edfmats, motor_scale, roi_obj, ii)
            TX_positions[scan72_motornames[ii]] = TX_pos
        except:
            print ('Fit failed for ROI No. %d.'%ii)

    return TX_positions

def fit_foci_2d(edfmats, roi_obj):
    """ **fit_foci_2d**

    Finds FWHMs of 2D Gaussians of the content of each ROI for a given stack of
    EDF-images.

    Args:

    edfmats (np.array): Stack of EDF-matrices from a d72scan of all translation
                        motors of the analyzer crystals.
    roi_obj (obj): ROI object of the xrs_rois class.

    Returns:

    sigma_1 (dict): Dictionary containing FWHMs in the vertical direction.
    sigma_2 (dict):    Dictionary containing FWHMs in the horizontal direction.
    """
    sigma_1 = {} # FWHM along dimension 1 (vertical)
    sigma_2 = {} # FWHM along dimension 2 (horizontal)
    for key, (pos, M) in six.iteritems(roi_obj.red_rois):
        S     = M.shape
        inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))
        ind   = 0
        sigma_1[key] = np.zeros(edfmats.shape[0])
        sigma_2[key] = np.zeros(edfmats.shape[0])
        for ii in range(len(edfmats)):
            sub_mat = edfmats[ii, inset[0], inset[1]] * (M/M.max())
            x = np.array(list(range(sub_mat.shape[0])))
            y = np.array(list(range(sub_mat.shape[1])))
            xx, yy = np.meshgrid(y, x)
            initial_guess = (np.amax(sub_mat),x.mean(),y.mean(),0.5,0.5,0)
            try:
                popt, pcov = optimize.curve_fit(math_functions.flat2DGaussian, (xx, yy), sub_mat.ravel(), p0=initial_guess)
                if np.abs(popt[3]) < 10.0:
                    sigma_1[key][ii] = popt[3]
                else:
                    sigma_1[key][ii] = 0.0
                if np.abs(popt[4]) < 10.0:
                    sigma_2[key][ii] = popt[4]
                else:
                    sigma_2[key][ii] = 0.0
            except:
                sigma_1[key][ii] = 0.0
                sigma_2[key][ii] = 0.0    
            ind += 1
    return sigma_1, sigma_2

def fit_foci_1d(edfmats, roi_obj):
    """ **fit_foci_1d**

    Finds FWHMs of 2D Gaussians of the content of each ROI for a given stack of
    EDF-images.

    Args:

    edfmats (np.array): Stack of EDF-matrices from a d72scan of all translation
                        motors of the analyzer crystals.
    roi_obj (obj): ROI object of the xrs_rois class.

    Returns:

    sigma_1 (dict): Dictionary containing FWHMs in the vertical direction.
    sigma_2 (dict):    Dictionary containing FWHMs in the horizontal direction.
    """
    sigma_1 = {} # FWHM along dimension 1 (vertical)
    sigma_2 = {} # FWHM along dimension 2 (horizontal)
    for key, (pos, M) in six.iteritems(roi_obj.red_rois):
        S     = M.shape
        inset = (slice( pos[0], pos[0]+(S[0]) ), slice( pos[1], pos[1]+(S[1]) ))
        ind   = 0
        sigma_1[key] = np.zeros(edfmats.shape[0])
        sigma_2[key] = np.zeros(edfmats.shape[0])
        for ii in range(len(edfmats)):
            sub_mat = edfmats[ii, inset[0], inset[1]] * (M/M.max())
            y0 = np.sum(sub_mat,axis=0)
            y1 = np.sum(sub_mat,axis=1)
            initial_guess0 = (np.where(y0 == y0.max())[0][0], 1.0, 1.0, y0.max(), 0.0)
            initial_guess1 = (np.where(y1 == y1.max())[0][0], 1.0, 1.0, y1.max(), 0.0)
            try:
                popt0, pcov0 = optimize.curve_fit(math_functions.pearson7_forcurvefit, np.arange(len(y0)), y0, p0=initial_guess0)
                popt1, pcov1 = optimize.curve_fit(math_functions.pearson7_forcurvefit, np.arange(len(y1)), y1, p0=initial_guess1)
                if np.abs(popt0[1]) < 10.0:
                    sigma_1[key][ii] = popt0[1]
                else:
                    sigma_1[key][ii] = 0.0
                if np.abs(popt1[1]) < 10.0:
                    sigma_2[key][ii] = popt1[1]
                else:
                    sigma_2[key][ii] = 0.0
            except:
                print (key)
                sigma_1[key][ii] = 0.0
                sigma_2[key][ii] = 0.0    
            ind += 1
    return sigma_1, sigma_2





#for name,key in zip(scan72_motornames,sorted(roifinder.roi_obj.red_rois)):
#    motor_scale[key] = scan72.counters[name]

def findBestFocus( sigma_1, sigma_2, counters, margin=3.0, verbose=False ):
    for key, ii in zip(sorted(sigma_1), list(range(len(sigma_1)))):
        motor_scale = counters[scan72_motornames[ii]]
        # find the best compromise/minimum for sigma_1 and sigma_2
        x = np.array(motor_scale)
        y1 = np.array(sigma_1[key])
        y2 = np.array(sigma_2[key])
        if show_fits:
            plt.cla()
            plt.title(key)
            plt.plot(x,y1,'-o')
            plt.plot(x,y2,'-o')
            plt.plot(x, np.polyval(np.polyfit(x,y1,2),x))
            plt.plot(x, np.polyval(np.polyfit(x,y2,2),x))
            plt.xlabel('motor position [mm]')
            plt.ylabel('FWHM [pixels]')
            plt.legend(['dim1', 'dim2', 'fit1', 'fit2'])
            plt.waitforbuttonpress()
        actual_pos = motor_scale[len(motor_scale)//2]
        min1_ind = np.where(np.polyval(np.polyfit(x,y1,2),x) == np.amin(np.polyval(np.polyfit(x,y1,2),x)))[0]
        min2_ind = np.where(np.polyval(np.polyfit(x,y2,2),x) == np.amin(np.polyval(np.polyfit(x,y2,2),x)))[0]
        min1 = motor_scale[min1_ind]
        min2 = motor_scale[min2_ind]
        if np.abs(min1 - min2) <= margin:
            min_pos = np.mean([min1, min2])
        else:
            min_pos = None
        if min_pos and np.abs(min_pos - actual_pos) <= margin:
            if verbose:
                print(key + ' :')
                print('current motor position is: %0.4f'%actual_pos )
                print('optimum focus for estimated to:')
                print('umv ' + scan72_motornames[ii] + ' %0.4f # current position is: %0.4f '%(min_pos,actual_pos))
                print('\n')
            else:
                print('umv ' + scan72_motornames[ii] + ' %0.4f # current position is: %0.4f '%(min_pos,actual_pos))



#np.polyval(np.polyfit(x,y2,2),x)


def findAnalyzerFocus(edfmats,motor_scale,roi_obj,roiNumber):
    """ **findAnalyzerFocus**

    Returns motor position that optimizes the analyzer focus subject to a 2D Gaussian fit.

    Args: 
    -----
    edfmats (np.array): 3D Numpy array containing the EDF-matrices of a scan.
    motorScale (np.array): Motor positions along the scan (analyzer tx-scan).
    roi_obj (xrs_rois.roi_object): ROI object, defined from the scan, should have some margins around the spots.
    roiNumber (int): Number indicating which ROI to be optimized.

    Returns:
    --------
    optPos (float): Motor position that optimizes the focus of the analyzer in question.

    """
    xmin = min(roi_obj.x_indices[roiNumber])
    xmax = max(roi_obj.x_indices[roiNumber])
    ymin = min(roi_obj.y_indices[roiNumber])
    ymax = max(roi_obj.y_indices[roiNumber])
    x = np.arange(xmin,xmax)
    y = np.arange(ymin,ymax)
    xx, yy = np.meshgrid(y, x)

    sigma_1 = [] # FWHM along dimension 1
    sigma_2 = [] # FWHM along dimension 2

    # go through all images and fit the 2D Gaussian
    for ii in range(edfmats.shape[0]):
        initial_guess = (np.amax(edfmats[ii,xmin:xmax,ymin:ymax]),(ymin+ymax)/2.,(xmin+xmax)/2.,1.0,1.0,0)
        popt, pcov = optimize.curve_fit(math_functions.flat2DGaussian, (xx, yy), edfmats[ii,xmin:xmax,ymin:ymax].ravel(), p0=initial_guess)
        sigma_1.append(popt[3])
        sigma_2.append(popt[4])

    # find the best compromise/minimum for sigma_1 and sigma_2
    x = np.array(motor_scale)
    y1 = np.array(sigma_1)
    y2 = np.array(sigma_2)

    # popt will have (x,amp,x0,fwhm)
    popt_1, pcov_1 = optimize.curve_fit( gauss_forcurvefit, x, y1 )
    popt_2, pcov_2 = optimize.curve_fit( gauss_forcurvefit, x, y2 )

    return np.mean([popt_1[2], popt_2[2]])



def fouc_det_focus(path, scan_number, SPECfname='rixs', EDFprefix='/edf/', EDFname='rixs_', EDFpostfix='.edf'):
    """ **fouc_det_focus**
    Returns best focus for FOURC spectrometer.

    Args:
        path        (str): Absolute path to location of the SPEC-file.
        roi_obj (roi_obj): ROI object of the xrs_rois class.
        scan_number (int): Scan number of the 72-motor scan.
        SPECfname   (str): SPEC-file name.
        EDFprefix   (str): Prefix to where EDF-files are stored.
        EDFname     (str): Base name of the EDF-files.
        EDFpostfix  (str): Post-fix used for the EDF-files.

    Returns:
        Optimized dtx and dtz position.

    """

    # create xrs_read object
    fourc_obj = xrs_read.Fourc(path,SPECfname=SPECfname, EDFprefix=EDFprefix, EDFname=EDFname, EDFpostfix=EDFpostfix)
    im4rois   = fourc_obj.SumDirect(scan_number)

    # creat ROI object
    roifinder = roifinder_and_gui.roi_finder()
    roifinder.get_zoom_rois(im4rois)

    # create scan object, load scan, use first ROI defined
    scan = xrs_scans.Scan()
    scan.load(path, SPECfname, EDFprefix, EDFname, EDFpostfix, scan_number)
    scan.get_raw_signals(roifinder.roi_obj, method='row')

    # a2scan motors
    dtx = scan.counters['detector x']
    dtz = scan.counters['detector z']

    # fit width of ROI at each point of the a2scan
    roi_shape = scan.raw_signals['ROI00'].shape
    fwhms = []
    for ii in range(roi_shape[0]):
        y  = scan.raw_signals['ROI00'][ii,:]
        x  = np.arange(len(y))
        p0 = ( y.max(), x[np.where(y==y.max())[0][0]], roi_shape[1]/3 )
        try:
            popt, pcov = optimize.curve_fit(math_functions.gauss_forcurvefit, x, y, p0=p0)
            y_fit = math_functions.gauss_forcurvefit(x, popt[0], popt[1], popt[2])
            plt.cla()
            plt.plot(x, y, '-ok')
            plt.plot(x, y_fit,'-r')
            plt.xlabel('detector pixel')
            plt.ylabel('intensity')
            plt.legend(['data', 'Gaussian fit'])
            plt.hold(True)
            plt.draw()
            #plt.waitforbuttonpress()
            plt.pause(0.01)
        except:
            print('aaaaaaaaaaahhhhhhhhh')
            popt = np.zeros((3,))
        fwhms.append(popt[2])

    # fit quadartic function to all FWHMs
    try:
        fit1 = np.polyval( np.polyfit(dtx, fwhms, 2), dtx)
        fit2 = np.polyval( np.polyfit(dtz, fwhms, 2), dtz)
    except:
        print ('whaoooo')
        fit1 = np.zeros_like(dtx)
        fit2 = np.zeros_like(dtz)

    # plot results, should be in one figure with double x-axis
    plt.figure()
    plt.plot(dtx, fwhms, '-o')
    plt.plot(dtx, fit1, '-')
    #plt.hold(True)
    plt.draw()

    plt.figure()
    plt.plot(dtz, fwhms, '-o')
    plt.plot(dtz, fit2, '-')
    #plt.hold(True)
    plt.show()

    # return/print out optimal positions
    dtx_min = dtx[np.where(fit1==fit1.min())[0][0]]
    dtz_min = dtz[np.where(fit1==fit1.min())[0][0]]
    print('Minimum position for dtx = %6.4f'%dtx_min )
    print('Minimum position for dtz = %6.4f'%dtz_min )
    return dtx_min, dtz_min
    










