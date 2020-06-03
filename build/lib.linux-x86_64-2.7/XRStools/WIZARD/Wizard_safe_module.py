from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import h5py
import tokenize    
import math, pylab, numpy
import numpy as np

class funct_Readline:
    def __init__(self,s):
        self.sl=s.split("\n")
    def __call__(self):
        if len(self.sl)==0:
            raise StopIteration
        res = self.sl[0]
        self.sl = self.sl[1:]
        return res

def python_code_is_safe(s):
    rl=funct_Readline(s)
    tk = tokenize.generate_tokens(rl)

    for t in tk:
        t=t[1]

        if t in ["import","exec"]:
            return 0
        if t[:1]=="_":
            return 0
        if t[:4]=="exec":
            return 0

    return 1

def make_obj( h5group , h5, groupname     ):
    runit=None
    if isinstance(h5group,h5py._hl.dataset.Dataset):
        pos = groupname.rfind("python_")
        if pos>=0 and "/" not in groupname[pos:]:
            func_name = str(groupname[pos:][len("python_"):])
            groupname= str(groupname[:pos])
            h5group = h5[groupname]
            runit=func_name
    if isinstance(h5group,h5py._hl.dataset.Dataset):
        oggetto = h5group.value
    else:
        dictio={}
        for key in h5group:
            if isinstance(h5group[key],h5py._hl.dataset.Dataset):
                dictio[key]= h5group[key].value
                if key[:len("python_")]=="python_":
                    func_name = key[len("python_"):]
                    sfunc = h5group[key].value
                    print( " SFUNC " , sfunc)
                    s="def %s(self):\n"%func_name
                    sl = sfunc.split("\n")
                    for l in sl:
                        s+="  "+l+"\n"
                    s+="foo=%s\n"%func_name
                    if python_code_is_safe(s):
                        print( " ESEGUO " , s)
                        resdic = {"math":math,"numpy":np,"np":np,"pylab":pylab}
                        
                        exec(s , resdic, resdic)
                        foo = resdic["foo"]
                        dictio[func_name] = foo
                        if func_name == runit:
                            runit = foo
                    else:
                        print( "cowardly refusing to execute suspicious python code ", s)
        oggetto = type('MyObject',(object,),dictio)()

    return runit, oggetto








