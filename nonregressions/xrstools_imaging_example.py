from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import range
from six.moves import zip
#!/usr/bin/python
# Filename: oneD_imaging_example.py
__doc__="""
ici la doc  ::

  python xrstools_imaging_example.py


"""
from XRStools import xrs_read,xrs_imaging, xrs_rois, roifinder_and_gui
import sys
import pylab
# from pylab import *

def main():
    pylab.ion()
    import numpy as np


    # data_home = '/home/christoph/data/ihr_jul15/'
    data_home = '/data/id20/inhouse/data/run3_15/run6_ihr/'
    gasket = xrs_imaging.oneD_imaging(data_home+'hydra',monitorcolumn='kapraman',energycolumn='sty')

    image4roi =  gasket.SumDirect( [372] )
    roifinder = roifinder_and_gui.roi_finder()
    roifinder.get_zoom_rois(image4roi)
    gasket.set_roiObj(roifinder.roi_obj)

    # gasket.loadscan_2Dimages(range(372,423),scantype='sty')
    gasket.loadscan_2Dimages(list(range(372,374)),scantype='sty')

    # try making a 3D image from ROI 0
    s = np.zeros((np.shape(gasket.twoDimages['Scan372'][0].matrix)[0],np.shape(gasket.twoDimages['Scan372'][0].matrix)[1],len(gasket.twoDimages)))
    for ii,number in zip(list(range(s.shape[-1])),list(range(372,423))):
        scanname = 'Scan%03d' % number
        s[:,:,ii] = gasket.twoDimages[scanname][0].matrix


    from mayavi import mlab

    #x, y, z = np.ogrid[-10:10:20j, -10:10:20j, -10:10:20j]

    #s = np.sin(x*y*z)/(x*y*z)

    # s is a 3D np.ndarray
    src = mlab.pipeline.scalar_field(s)
    mlab.pipeline.iso_surface(src, contours=[s.min()+0.05*s.ptp(), ], opacity=0.6)
    #mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
    mlab.show()

    # s is a 3D np.ndarray
    src = mlab.pipeline.scalar_field(s)
    mlab.pipeline.volume(src,vmin=1000.0, vmax=2000.0)
    mlab.show()

if(sys.argv[0][-12:]!="sphinx-build"):
    main()













