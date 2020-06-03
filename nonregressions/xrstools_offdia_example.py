from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .XRStools import xrs_read, theory, xrs_extraction, xrs_rois, roifinder_and_gui, xrs_utilities, math_functions
from pylab import *
from six.moves import range

ion()
import numpy as np
from scipy import interpolate, signal, integrate, constants, optimize, ndimage
from .XRStools import ixs_offDiagonal,xrs_utilities

##############################################################################
# off-focus resolution tests
##############################################################################
repertorio = "/home/christoph/data/ihr_sep15/"

############
# elastic line
############
offdia = ixs_offDiagonal.offDiagonal(repertorio + 'rixs',scanMotor='srz',monitorName='kaprixs',edfName=None,armLength=1.0)

image4roi =  offdia.SumDirect( [442] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)

offdia.set_roiObj(roifinder.roi_obj)
offdia.loadRockingCurve(list(range(477,757,3)),direct=True)
offdia.loadRockingCurve(list(range(478,758,3)),direct=True)

offdia.stitchRockingCurves()

offdia.offDiaDataSets[1].normalizeSignals()

offdia.offDiaDataSets[1].alignRCmonitor()

ion()
imshow(offdia.offDiaDataSets[1].alignedRCmonitor)

offdia.offDiaDataSets[1].removeElastic()
offdia.offDiaDataSets[1].windowSignalMatrix(5.0,50)









