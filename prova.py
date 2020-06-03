from  PyQt4 import Qt, QtCore
import PyQt4

def dum():
    print (" gestate qlineedit")
    return {}

class Pippo:
    def __init__(self):
        self.ql=Qt.QLineEdit()
        self.ql.__getstate__= dum
        pass
    # def __getstate__(self):
    #     print(" in getstate ")
    #     return {}
    
app=Qt.QApplication([])

import copy
pippo=Pippo()
a=copy.deepcopy(pippo)
print (a)
print (a.ql)
