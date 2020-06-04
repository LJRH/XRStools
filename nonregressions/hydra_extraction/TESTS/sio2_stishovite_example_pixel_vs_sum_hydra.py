from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .XRStools import xrs_read, roifinder_and_gui
from pylab import *
import numpy as np

path = '/media/christoph/Seagate Expansion Drive/data_nonregressions/hydra_extraction/'

###################################
# sio2 stishovite - no compensation
###################################
sio2_sum = xrs_read.Hydra(path)

image4roi =  sio2_sum.SumDirect( [30] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)
sio2_sum.set_roiObj(roifinder.roi_obj)

sio2_sum.load_scan(31, direct=True, scan_type='long', method='sum')
sio2_sum.load_scan(32, direct=True, scan_type='edge1', method='sum')
sio2_sum.load_scan(33, direct=True, scan_type='edge2', method='sum')
sio2_sum.load_scan(30, direct=True, scan_type='elastic', method='sum')

sio2_sum.get_spectrum_new(method='sum',include_elastic=True)

#ion()
#plot(sio2_sum.eloss, sio2_sum.signals)

###############################################
# sio2 stishovite - pixel-by-pixel compensation
###############################################
sio2_pixel = xrs_read.Hydra(path)

sio2_pixel.set_roiObj(roifinder.roi_obj)

sio2_pixel.load_scan(31, direct=True, scan_type='long', method='pixel')
sio2_pixel.load_scan(32, direct=True, scan_type='edge1', method='pixel')
sio2_pixel.load_scan(33, direct=True, scan_type='edge2', method='pixel')
sio2_pixel.load_scan(30, direct=True, scan_type='elastic', method='pixel')

sio2_pixel.get_spectrum_new( method='pixel', include_elastic=True )


plot(sio2_sum.eloss, sio2_sum.signals[:,0])
plot(sio2_pixel.eloss, sio2_pixel.signals[:,0])
legend(['sum, ROI00', 'pixel, ROI00'])
xlabel('energy loss [eV]')
ylabel('intensity [arb. units]')
show()











