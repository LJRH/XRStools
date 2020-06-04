from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from pylab import *
ion()
from .xrstools import id20_imaging
reload(id20_imaging)
test = id20_imaging.imaging('xrstools/things/licl_test_files/raman')
test.loadscan(90)
test.getlinrois(90,numrois=1)
test.image_from_line(90,'energy_cc')


a = test.image_from_line(90,'energy_cc')
from pylab import *
ion()
imshow(np.squeeze(a))









