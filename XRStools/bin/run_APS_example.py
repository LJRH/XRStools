from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from scipy import interpolate
from .xrstools import xrs_read, theory, extraction
from pylab import *
from six.moves import range

# try loading old Si data (this is more complicated since this is not id16 data)
counters=[1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19]
data = np.loadtxt('xrstools/things/Ba24Si100_rawdata_APS11.dat')

tth = np.array(list(range(9,171,9)))

# initiate data
eloss   = data[:,0]
signals = data[:,1:]
errors  = np.sqrt(np.absolute(data[:,1:]))
E0      = 7.912

# initiate id20 instance
basi = xrs_read.read_id16('xrstools/things/licl_test_files/raman')

basi.eloss = eloss
basi.signals = np.absolute(signals)
basi.errors  = errors
basi.E0 = E0
basi.tth = tth

# initiate HFspectrum instance from the theory module
hf = theory.HFspectrum(basi,['Ba24Si100'],[1.0],correctasym=[[0.0,0.0]])

# initiate instance of extraction class for background removal
# (using the instance of the read_id20 and the HFspectrum class):
extr = extraction.extraction(basi,hf)

# apply energy dependent corrections to the experimental data (absorption, self-absorption, relativistic scattering cross section)
extr.energycorrect(list(range(17)),10,2.3,0.1)


# for a quick estimate of the core profile and if the edge is not too far into the compton profile:
#extr.removeelastic(range(13,17),[0,50],[800,850],overwrite=True,stoploop=False)
#extr.remquickval(range(13,17),[400,1000],[90,120],20)

# for extraction of the valence profile:
extr.removeelastic(list(range(13,17)),[0,50],[800,850],overwrite=True,stoploop=False)










