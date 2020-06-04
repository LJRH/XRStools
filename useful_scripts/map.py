from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import h5py
import numpy as np
from six.moves import filter

"SPECTRA.h5:/myexp_84_106"


f=h5py.File("SPECTRA.h5","r")

datagroup = f["/myexp_84_106/2/"]

def reorder(sl):
    sl=  sorted(  sl     , key = lambda x:  int(list(filter(str.isdigit, str(x) ))) )
    return sl

e0keys = reorder([k for k in  datagroup.keys() if k[:3]=="E0_"] )
ekeys  = reorder([k for k in  datagroup.keys() if k[:3]=="ene"] )
skeys  = reorder([k for k in  datagroup.keys() if k[:4]=="spec"])

print( skeys)

print( ekeys)

print( e0keys)

first_spectra = datagroup[skeys[0] ][:]
elen =  len(first_spectra)


map = np.zeros([ len(e0keys)  ,  elen ])

for i,k in enumerate(skeys):
    spectra = datagroup[k ][:]
    map[i,:] = spectra

print( datagroup[ekeys[0] ][:]- datagroup[ e0keys[0] ].value*1000)

import pylab
pylab.imshow(map)
pylab.show()

    











