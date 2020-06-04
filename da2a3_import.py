from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

fl=open("del.txt","r")
for l in fl:
    fn = l.strip()
    s=open(fn,"r").read()
    print( " FILE : ", fn)

    sl = s.split("\n")
    for nl,ll in enumerate(sl):
        llorig = ll
        # ll=ll.strip()
        if len(ll.strip()) and ll.strip()[0]=="#":
            continue
        # ll_sp = ll.split()
        status = 0
        if "import" in ll:
            print( "%10d"% nl,"  " ,      ll)

