from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
reload(xes_read)
from .XRStools import xes_read, xrs_rois, roifinder_and_gui
import numpy as np
from pylab import *
ion()

xest = xes_read.read_id20('/home/christoph/data/hc2270_xes/rixs')

# ROI
image4roi =  xest.SumDirect( [71] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)

# set ROI object, load scans
xest.set_roiObj(roifinder.roi_obj)
xest.loadscandirect(70,'part1')
xest.loadscandirect(71,'part2')

xest.getXESspectrum()


xest.getrawdata()









