from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from  PyQt4 import Qt, QtCore, QtGui
# from .. import ui


from XRStools.installation_dir import installation_dir

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]


# @ui.UILoadable
class scansTable(QtGui.QWidget) :
    def __init__(self, parent=None):
        
        super( scansTable , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "scansTable.ui" ), self)
        
        qb = QtGui.QPushButton("+")
        qb.setFixedWidth(20)
        sb = QtGui.QSpinBox( )
        sb.setMinimum(1)
        sb.setMaximum(1000)
        self.gridLayout.addWidget(sb ,0,0)
        self.gridLayout.addWidget(qb,0,1)
        self.gridLayout.addWidget(QtGui.QPushButton("-"),0,2)
        

if __name__ =="__main__":
    app=Qt.QApplication([])
    w = scansTable()
    w.show()
    app.exec_()
    
        
