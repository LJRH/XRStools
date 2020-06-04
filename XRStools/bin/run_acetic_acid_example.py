from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .xrstools import xrs_read
from pylab import *
from six.moves import range
ion()

# C K-edge 180 C, 3 bars 
aaC180 = xrs_read.read_id20('/home/csahle/data/ch3898/raman',energycolumn='energy',monitorcolumn='kap4dio')
aaC180.loadelastic([150])
aaC180.getlinrois([150],numrois = 3)

aaC180.loadloop([135,143],4)

aaC180.getrawdata()
aaC180.getspectrum()
aaC180.geteloss()

for n in range(len(aaC180.signals[0,:])):
    plot(aaC180.eloss,aaC180.signals[:,n],'-')








