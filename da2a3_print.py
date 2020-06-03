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
        ll=ll.strip()
        if len(ll) and ll[0]=="#":
            continue
        ll_sp = ll.split()
        if "print" in ll_sp:
            i = ll_sp.index("print")
            if i<len(ll_sp) and ll_sp[i+1][0] =="(":
                continue
            if len(ll_sp) and ll_sp[0]=="will":
                continue
            i = llorig.find("print")
            llnew = llorig[:i+5]+"("+llorig[i+5:]+")"
            print( "                 ", nl,"            " ,      llnew)
            sl[nl]=llnew
    if True or "qfile.py" in fn:
        f=open(fn,"w")
        for l in sl:
            f.write(l+"\n")

