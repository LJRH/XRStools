from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .XRStools import xrs_ComptonProfiles as CP
import numpy as np

from pylab import *
ion()

filename = '/home/christoph/sources/XRStools/data/ComptonProfiles.dat'

cp = CP.AtomProfile('Si',filename)

E0 = 9.7
twotheta = 35.0

cp.get_elossProfiles(E0, twotheta)







