from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .xrstools import xrs_read

reload(xrs_read)

foo = xrs_read.read_p01('/home/christoph/data/p01jun13/')
foo.loadelastic(268)
foo.getautorois(268)
#foo.getzoomrois(298)

foo.loadedge([276])#137,138,143,145
#foo.getzoomrois(192,numrois=3)

foo.getrawdata()
foo.getspectrum()

import pylab
pylab.plot(foo.energy,foo.signals[:,0],'-')
pylab.show()








