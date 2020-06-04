from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

tochange={}
f=open("div2.txt","r").read()
fl = f.split("\n")
for l in fl:
    if "FILE :" in l:
        fn= l[7:].strip()
        print( fn)
    else:
        if "//" in l:
            ln = int(l[:10])
            l = l[14:]
            print( l,)
            tochange[(fn,ln)]=l


fl=open("del.txt","r")
for l in fl:
    fn = l.strip()
    s=open(fn,"r").read()

    sl = s.split("\n")
    for nl,ll in enumerate(sl):
        if (fn,nl) in tochange:
            sl[nl] = tochange[(fn,nl)]

    if True : # or "Wizard.py" in fn:
        f=open(fn,"w")
        for l in sl:
            f.write(l+"\n")
