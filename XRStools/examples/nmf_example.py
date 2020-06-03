from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/env python

import numpy as np

w1 = np.array([[1,2,3],[4,5,6]])
h1 = np.array([[1,2],[3,4],[5,6]])
w2 = np.array([[1,1,3],[4,5,6]])
h2 = np.array([[1,1],[3,4],[5,6]])

v = np.dot(w1,h1)

(wo,ho) = nmf(v, w2, h2, 0.001, 10, 10)
print( wo)
print( ho)







