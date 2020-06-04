from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: superresolution.py

import numpy as np
from pylab import *
from six.moves import range
ion()

from .xrstools import superresolution

reload(superresolution)

sr = superresolution.imageset()

sr.loadhe3070('xrstools/things/he3070-zscan.mat')
#sr.loadkimberlite('xrstools/things/kimberlite1.mat')
#sr.loadkimberlite('xrstools/things/he3070-zscan_likesimo.mat')
sr.estimate_xshifts()

#sr.interpolate_xshift_images(2)
sr.interpolate_xshift_images(4,whichimages=[0,1,2,4,5,6,7,8])

sr.plotSR()

for ii in range(9):
    figure()
    sr.plotLR(ii)









