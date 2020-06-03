from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

fl=open("del.txt","r")
for l in fl:
    fn = l.strip()
    s=open(fn,"r").read()
    s = """from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
""" +s
    open(fn,"w").write(s )
