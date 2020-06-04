from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from scipy import interpolate
from .xrstools import xrs_read, theory, extraction
from .xrstools.helpers import *
from pylab import *
from six.moves import range
ion()

# try loading old Si data
counters=[1, 2, 3, 4, 5, 6, 7, 8, 9]
data = np.loadtxt('xrstools/things/diamond_data_E016keV_tth123.dat')
inds = np.where(np.logical_and(data[:,0]>=6.0,data[:,0]<=7.5))[0]
data[inds,1::] = 0

# remove some glitches in the data
glitchindex1 = np.where(data[:,0]==1147.0999999999999)[0][0]
data = np.delete(data, glitchindex1, 0)
glitchindex2 = np.where(data[:,0]==397.10000000000002)[0][0]
data = np.delete(data, glitchindex2, 0)

# scale some data points
inds = np.where(data[:,0]>1136.08)
data[inds,5::] = data[inds,5::] + 0.0004

clf();plot(data[:,0],data[:,1::])

meantth = 133
# tth = np.array([meantth-13.0, meantth-6.5, meantth-6.5, meantth+0.0, meantth+0.0, meantth+7.5, meantth+7.5, meantth+14.0, meantth+14.0])
tth = np.array([115, 121, 121, 129, 129, 137, 137, 143, 143])

# initiate data
eloss   = data[:,0]
signals = data[:,1:]
errors  = np.sqrt(np.absolute(data[:,1:]))
E0      = 15.8

# reset the eloss scale
inds = np.where(np.logical_and(eloss>=-10.0,eloss<=10.0))[0]
for ii in range(len(data[0,1:])):
    cenom = np.trapz(signals[inds,ii]*eloss[inds],eloss[inds])
    maxi  = signals[inds,ii] == np.amax(signals[inds,ii])
    #print cenom
    escale = eloss-eloss[maxi]
    signals[:,ii] = np.interp(eloss, escale, signals[:,ii])

#plot(eloss,signals)
pz = e2pz(eloss/1e3+E0,E0,tth[0])[0]


# initiate id20 instance (from a random Spec-file, since the init checks if the specified file actually exists)
dia = xrs_read.read_id16('xrstools/things/licl_test_files/raman')

dia.eloss   = eloss
dia.signals = np.absolute(signals)
dia.errors  = errors
dia.E0      = E0
dia.tth     = tth
dia.cenom   = np.zeros(len(tth))+15.8

# initiate cprofiles instance
hf = theory.HFspectrum(dia,['C'],[1.0],correctasym=[[0.0]])

#plot(hf.eloss,hf.J[:,0],dia.eloss,dia.signals[:,0])

# initiate instance of extraction class for background removal:
#reload(extraction)
extr = extraction.extraction(dia,hf,prenormrange=[20,np.inf])
plot(extr.eloss,extr.signals[:,0],extr.eloss,extr.J[:,0])

extr.extractval_test(8,linrange1=[340,470],linrange2=[1500,2500])

extr.getallvalprof(8,smoothgval=50.0,stoploop=False)
extr.remvalenceprof_test(list(range(9)),eoffset=20.0)

ion()
clf()
plot(extr.eloss,extr.sqw[:,0],extr.eloss,extr.C[:,0])

# save some things for gnuplot
thedata = np.zeros((len(extr.eloss),5))
thedata[:,0] = extr.eloss
thedata[:,1] = extr.signals[:,8]-0.0002
thedata[:,2] = extr.C[:,8]
thedata[:,3] = np.interp(extr.eloss,extr.eloss,extr.valence[:,8])
thedata[:,4] = np.interp(extr.eloss,extr.eloss,extr.valasymmetry[:,8]*8)
np.savetxt('/home/csahle/Dropbox/tool_paper/figures/analysis/val_extraction_SCVA_det9.dat',thedata)

thesqw = np.zeros((len(extr.eloss),10))
thesqw[:,0] = extr.eloss
theval = np.zeros((len(extr.eloss),10))
theval[:,0] = extr.eloss
thedata = np.zeros((len(extr.eloss),10))
thedata[:,0] = extr.eloss
for ii in range(1,11):
    thesqw[:,ii] = extr.sqw[:,ii]
    theval[:,ii] = extr.valence[:,ii]
    thedata[:,ii] = extr.signals[:,ii]-extr.background[:,ii]

np.savetxt('/home/csahle/Dropbox/tool_paper/figures/analysis/val_extraction_sqw_alldet.dat',thesqw)
np.savetxt('/home/csahle/Dropbox/tool_paper/figures/analysis/val_extraction_val_alldet.dat',theval)
np.savetxt('/home/csahle/Dropbox/tool_paper/figures/analysis/val_extraction_data_alldet.dat',thedata)









