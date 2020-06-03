from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
fl=open("del.txt","r")
for l in fl:
    f=open(l,"r").read()
    print f[:100]
