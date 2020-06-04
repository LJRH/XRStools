from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from .xrstools import xrs_read, theory, extraction
import pylab
from six.moves import range

# create an instance of the read_id16 class from the xrs_read module
licl = xrs_read.read_id16('xrstools/things/licl_test_files/raman')
# load a scan (Spec- and edf-file) into your class (here: an elastic scan)
licl.loadelastic(90)
# use the scan you loaded (No. 90) for finding ROIs automatically
licl.getautorois(90)
# other options to find ROIs are 
# getzoomedgerois(90,energyval)
# getlinedgerois(90,energyval)
# getzoomrois(90,numrois=9)
# getlinrois(90)

# load more scans, scan sequence of 2 loops (starting at scan No. 165 and 170) of 5 scans each
licl.loadloop([165,170],5)

#licl.loadscan([165,170],'edge1')
# load  a wide overview scan
#licl.loadscan([166,171],'edge2')
#licl.loadscan([167,172],'edge3')
#licl.loadscan([168,173],'edge4')
#licl.loadscan([169,174],'edge5')


licl.loadlong(157)

# apply the ROIs to all scans
licl.getrawdata()
licl.getspectrum()

# find the center of mass of the elastic line (define energy loss scale)
licl.geteloss()
# set the scattering angle
licl.gettths(35)
# calculate the momentum transfers
licl.getqvals()

# plot the experimental data
for n in range(len(licl.signals[0,:])):
    pylab.plot(licl.eloss,licl.signals[:,n],'.')

pylab.axis([-50.0,620.0,0.0,0.005])
pylab.show()

# create an instance of the HFspectrum class from the theory module
hf = theory.HFspectrum(licl,['H2O'],[1.0],correctasym=[[0.0,0.0]])
# create an instance of the extraction class from the extraction 
# module (using the instances of the read_id20 and HFspectrum classes)
extr = extraction.extraction(licl,hf)

# remove a constant from the raw data and scale the data to the Hartree-Fock edge jump
extr.removeconstpcore(list(range(9)),[527.0,534.5],[545.0,588.0],weights=[1,1],stoploop=False)
# average over all 9 spectra
extr.averageqs(list(range(9)))

# plot the result
pylab.plot(extr.eloss,extr.sqwav)
pylab.show()










