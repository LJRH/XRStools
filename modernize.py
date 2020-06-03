from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
fl=open("del.txt","r")
for l in fl:
    os.system("python-modernize -w %s"% l.strip())
