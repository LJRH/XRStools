from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
import glob

fs = glob.glob("../XRStools/scripts/*")
print( fs)

for f in  fs:
    n = os.path.basename(f)
    os.system("copy tipo  %s.bat"%n)







