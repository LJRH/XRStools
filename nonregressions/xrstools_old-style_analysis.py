from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .XRStools import xrs_read, theory, extraction, xrs_rois, roifinder_and_gui
import numpy as np
from pylab import *
ion()

lif = xrs_read.read_id20('/home/christoph/data/ihr_may15/rixs',monitorcolumn='kaprixs')
lif.loadelastic(59)
lif.set_roiObj(roifinder_and_gui.get_zoom_rois(lif.scans,59))

lif.loadelasticdirect([59])
lif.loadloopdirect([60,63,66,69,72,75,78],2)
lif.loadlongdirect(58)

lif.







