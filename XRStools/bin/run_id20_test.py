from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .xrstools import xrs_read, theory, extraction
from pylab import *
from six.moves import range

t = xrs_read.read_id20('/home/christoph/data/ch3898/raman',energycolumn='energy',monitorcolumn='kap4dio')
t.loadelastic(259)
t.getlinrois(259,numrois=3)

# t.loadloop([210,221,232,243],3)
t.loadscan([260],'edge1')
t.loadscan([263],'edge2')
t.loadscan([264],'edge3')
t.loadscan([265],'edge4')
t.loadscan([266],'edge5')
t.loadscan([267],'edge6')
t.loadscan([268],'edge7')
t.loadlong([258])

t.getrawdata()
t.getspectrum()
t.geteloss()

t.gettths(rhl=45.0,rhr=45.0,rhb=143.38,rvd=-45.0,rvu=87.66,rvb=121.88)

# plot the experimental data
for n in range(len(t.signals[0,:])):
    plot(t.eloss,t.signals[:,n],'-')

show()

#for n in range(len(t.signals[0,:])):
#    plot(t.energy,t.signals[:,n],[t.cenom[n],t.cenom[n]],[0,1]))

from .xrstools import xrs_read, theory, extraction
from pylab import *
ion()
t = xrs_read.read_id20('/home/csahle/data/ch3898/raman',energycolumn='energy',monitorcolumn='kap4dio')
t.loadelastic(259)
t.loadloop([265],4) 
t.loadlong(258)
t.getlinrois(259,numrois=12)
t.getrawdata()
t.getspectrum()
t.geteloss()








