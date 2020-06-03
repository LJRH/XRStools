from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from scipy import interpolate
from .xrstools import xrs_read, theory, extraction
from pylab import *
from six.moves import range

# try loading old Si data (this is more complicated since this is not id16 data)
counters=[3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19]
data = []
for counter in counters:
    try:
        nr = '%02d' % counter
        fn = 'xrstools/things/Si/fig_raw_si_' + nr + '.dat'
        data.append(np.loadtxt(fn))
    except:
        print( 'Det ' + nr + ' not found!')

tth = np.array([26, 32, 48, 63, 69, 79, 87, 96, 106, 115, 123, 135, 145, 157, 175])+8

# initiate data
eloss   = data[0][:,0]
signals = np.zeros((len(eloss),len(counters)))
errors  = np.zeros(np.shape(signals))
E0      = 9.893 # np.zeros((len(counters),1))+9.893

# spline onto first detector
for n in range(len(counters)):
    f = interpolate.interp1d(data[n][:,0],data[n][:,1], bounds_error=False, fill_value=0.0)
    signals[:,n] = f(eloss)
    f = interpolate.interp1d(data[n][:,0],data[n][:,2], bounds_error=False, fill_value=0.0)
    errors[:,n]  = f(eloss)

# initiate id20 instance (from a random Spec-file, since the init checks if the specified file actually exists)
si = xrs_read.read_id16('xrstools/things/licl_test_files/raman')

si.eloss = eloss
si.signals = signals
si.errors  = errors
si.E0 = E0
si.tth = tth

#pylab.plot(si.eloss,si.signals,'.')
#pylab.show()

# initiate HFspectrum instance from the theory module for low-q extraction without core asymmetry correction:
hf = theory.HFspectrum(si,['Si'],[1.0],correctasym=[[1.5]],correctasym_pertth=[0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5])

# initiate instance of extraction class for background removal 
# (using the instance of the read_id20 and the HFspectrum class):
extr = extraction.extraction(si,hf)

# apply energy dependent corrections to the experimental data (absorption, self-absorption, relativistic scattering cross section)
extr.energycorrect(list(range(15)),10,2.3,0.2)

ion()
plot(extr.eloss,extr.signals[:,14],extr.eloss,extr.J[:,14])


# example of removing a pearson function from low-q data (here, you have to play around with
# the fitting regions quite a bit, and difficult already for whichq=2)
extr.removecoreppearson2(0,[60,98],[105,350])
extr.removecoreppearson2(1,[80,98],[105,250])


# high q valence extraction
extr.extractval(14,linrange1=[40,95],linrange2=[850,1200],mirror=True)
extr.getallvalprof(14,smoothgval=5.0,stoploop=False)
extr.remvalenceprof(14)

# example of using the quick valence removal
extr.remquickval(13,[700.0,1500.0],[98.0,115.0],10.0)

plot(extr.eloss,extr.signals[:,0],extr.eloss,extr.J[:,0],extr.eloss,extr.C[:,0])
show()










