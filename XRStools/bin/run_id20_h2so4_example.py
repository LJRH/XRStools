from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .xrstools import xrs_read, theory, extraction
from pylab import *
import numpy as np
ion()
h2so4 = xrs_read.read_id20('/home/christoph/data/ch3914/orig/raman',energycolumn='energy',monitorcolumn='kap4dio')
#h2so4 = xrs_read.read_id20('xrstools/things/h2so4_test_files/raman',energycolumn='energy',monitorcolumn='kap4dio')
h2so4.loadelastic([908,914])
h2so4.getzoomrois(908,numrois=12)

h2so4.loadloopdirect([909,915],5)

h2so4.getrawdata()
h2so4.getspectrum()
h2so4.geteloss()

plot(h2so4.eloss,np.sum(h2so4.signals,axis=1))







