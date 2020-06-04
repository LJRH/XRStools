from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .XRStools import xrs_read,xrs_imaging, xrs_rois, roifinder_and_gui, rixs_read
from pylab import *
ion()
import numpy as np

data_home = '/home/christoph/data/rixs_example/'
#data_home = '/data/id20/inhouse/data/run3_15/run6_ihr/'
rixs = rixs_read.read_id20(data_home+'rixs1',energyColumn='Anal Energy',monitorColumn='kap4dio')

# ROI
image4roi =  rixs.SumDirect( [237] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)

#
rixs.set_roiObj(roifinder.roi_obj)
rixs.loadscan(237,'elastic')









