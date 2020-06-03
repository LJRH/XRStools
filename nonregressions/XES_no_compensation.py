from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
###################################################
# generic example for extraction of a XES spectrum
# using the Fourc spectrometer at ID20 with bent
# crystal analyzers
###################################################

from .XRStools import xrs_read, roifinder_and_gui
import numpy as np
import pylab as pl

path_to_data = '/data/id20/inhouse/XRStools_nonregression_data/xes_no_compensation/'
path_to_data = '/media/christoph/Seagate Expansion Drive/data/hc2270_xes/'

xes = xrs_read.Fourc( path_to_data,SPECfname='rixs', moni_column='kaprixs', EinCoor='energy' )

# ROI
image4roi =  xes.SumDirect( [71] )
roifinder = roifinder_and_gui.roi_finder()
roifinder.get_zoom_rois(image4roi)

# pass ROI object to main class
xes.set_roiObj(roifinder.roi_obj)

# load scans
xes.load_scan(70, direct=True, scan_type='part1', method='sum')
xes.load_scan(71, direct=True, scan_type='part2', method='sum')

# construct emission spectrum
xes.get_XES_spectrum()

# plot results
pl.plot(xes.energy, xes.signals)
pl.xlabel('energy [keV]')
pl.ylabel('intensity [arb. units]')
pl.show()







