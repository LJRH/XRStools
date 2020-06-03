from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .XRStools import xrs_read, xrs_extraction, xrs_rois, roifinder_and_gui, xrs_utilities, math_functions
from pylab import *
from six.moves import range
ion()
import numpy as np

# set some input variables
# except for the 'path', all are default (see xrs_read.Hydra class) 
# and the user should not need to set them for a standard experiment
path        = '/data/id20/inhouse/data/run4_16/run1_ihr/'
path        = '/media/christoph/Seagate Expansion Drive/data/run4_16/run1_ihr/'
SPECfname   = 'hydra'
EDFprefix   = '/edf/'
EDFname     = 'hydra_'
EDFpostfix  = '.edf'
en_column   = 'energy'
moni_column = 'izero'

# in case the ROI should be saved:
#ROIpath = '/media/christoph/Seagate Expansion Drive/data/run4_16/run1_ihr/rois/hbn_1_72_rois_ordered_mini-example.h5'
ROIpath = None

# file name for ASCII output / HDF5 file
ASCIIpath = None
HDF4path  = None

################
# hexagonal BN #
################

hbn = xrs_read.Hydra(path, SPECfname=SPECfname, EDFprefix=EDFprefix, EDFname=EDFname, \
                        EDFpostfix=EDFpostfix, en_column=en_column, moni_column=en_column)

#
# here the ROI-GUI could/should be used
#

image4roi =  hbn.SumDirect( [475] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)
if ROIpath:
    roifinder.roi_obj.writeH5(ROIpath)

# set the ROI obj and load scans
hbn.set_roiObj(roifinder.roi_obj)
hbn.load_scan( 475, scan_type='elastic', direct=True, scaling=None )

# the long scan needs some scaling (unfortunately a different factor for
# each ROI (here ROIs in groups of three))
scaling = np.zeros(72)
scaling[list(range(3))] = 126.439/2.2
scaling[list(range(3,6))] = 126.439/1.8
scaling[list(range(6,9))] = 126.439/2.9
scaling[list(range(9,12))] = 126.439/8.8
scaling[list(range(12,15))] = 251.
scaling[list(range(15,18))] = 241.
scaling[list(range(18,21))] = 238.1
scaling[list(range(21,24))] = 201.0
scaling[list(range(24,27))] = 180.1
scaling[list(range(27,30))] = 149.7
scaling[list(range(30,33))] = 115.2
scaling[list(range(33,36))] = 109.31

hbn.load_scan(438, scan_type='long', direct=True, scaling=scaling )
hbn.load_scan([476,479], scan_type='B', direct=True, scaling=None )
hbn.load_scan([477,480], scan_type='N', direct=True, scaling=None )

# stitch the spectrum together
hbn.get_spectrum(include_elastic=False, abs_counts=False)

# check/visualize the result
ion()
plot(hbn.eloss, np.sum(hbn.signals[:,list(range(33,36))],axis=1))

# save the spectrum as ASCII-file
if ASCIIpath:
    hbn.dump_signals(ASCIIpath)

# or alternatively/additionally, the whole instance could be saved 
# in HDF5 (BUT I could not get it working yet)
if HDF4path:
    hbn.save_state_hdf5( filename, groupname, comment="" )










