from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


#!/usr/bin/python
# Filename: oneD_imaging_example.py

from .XRStools import xrs_read,xrs_imaging, xrs_rois, roifinder_and_gui
from pylab import *
from six.moves import range
from six.moves import zip
ion()
import numpy as np

#data_home = '/home/christoph/data/ihr_jul15/'
data_home = '/data/id20/inhouse/data/run3_15/run6_ihr/'

##### STY scan
fey = xrs_imaging.oneD_imaging(data_home+'hydra',monitorcolumn='kapraman',energycolumn='sty')

image4roi =  fey.SumDirect( [270] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)
fey.set_roiObj(roifinder.roi_obj)

fey.loadscan_2Dimages(270,scantype='sty')

# try making a 3D image from ROI 0
s = np.zeros((np.shape(fey.twoDimages['Scan270'][0].matrix)[0],np.shape(fey.twoDimages['Scan270'][0].matrix)[1],len(fey.twoDimages)))
for ii,number in zip(list(range(s.shape[-1])),[270]):
    scanname = 'Scan%03d' % number
    s[:,:,ii] = fey.twoDimages[scanname][0].matrix

from mayavi import mlab
# s is a 3D np.ndarray
src = mlab.pipeline.scalar_field(s)
mlab.pipeline.volume(src,vmin=10.0, vmax=20.0)
mlab.show()


##### STX scan
fex = xrs_imaging.oneD_imaging(data_home+'hydra',monitorcolumn='kapraman',energycolumn='sty')

image4roi =  fex.SumDirect( [273] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)
fex.set_roiObj(roifinder.roi_obj)

fex.loadscan_2Dimages(273,scantype='stx')

# try making a 3D image from ROI 0
s = np.zeros((np.shape(fex.twoDimages['Scan273'][0].matrix)[0],np.shape(fex.twoDimages['Scan273'][0].matrix)[1],len(fex.twoDimages)))
for ii,number in zip(list(range(s.shape[-1])),[273]):
    scanname = 'Scan%03d' % number
    s[:,:,ii] = fex.twoDimages[scanname][0].matrix

from mayavi import mlab
# s is a 3D np.ndarray
#src = mlab.pipeline.scalar_field(s)
#mlab.pipeline.volume(src,vmin=1.0, vmax=20.0)
#mlab.show()


src = mlab.pipeline.scalar_field(s)
mlab.pipeline.iso_surface(src, contours=[s.min()+0.05*s.ptp(), ], opacity=0.6)
#mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
mlab.show()









