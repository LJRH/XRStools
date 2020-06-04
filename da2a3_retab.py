from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



fl=open("del.txt","r")
for l in fl:
    
    fn = l.strip()
    print( fn)
    s=open(fn,"r").read()

    sl = s.split("\n")
    for nl,ll in enumerate(sl):
        if len(ll) and ll[0]!="\"" and "\t" in ll:
            print( fn , ":", ll)
            # ll = ll.replace("\t","    ")
            # sl[nl] = ll

    if False : # or "Wizard.py" in fn:
        f=open(fn,"w")
        for l in sl:
            f.write(l+"\n")
