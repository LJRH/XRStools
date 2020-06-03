from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
from PyQt4 import QtGui,QtCore

class   MyTree( QtGui.QTreeView):
    def __init__(self, *args):
        QtGui.QTreeView.__init__(self, *args)
        self.clicked[QtCore.QModelIndex].connect(self.itemClicked)

    def itemClicked(self, modelIndex):
        event ="itemDoubleClicked"
        print( event, modelIndex)
        
        index = self.selectedIndexes()[0]
        pippo = index.model().data(index)

        print( pippo)
        print( str(pippo))
        print( dir(pippo))
        print( pippo.toString())
    
        pippo = index.model().filePath(index)
        print( str(pippo))

class Myview(QtGui.QMainWindow):
    def __init__(self,parent=None):
        QtGui.QMainWindow.__init__(self)
        model = QtGui.QFileSystemModel()
        model.setNameFilters(["winfo_*.py"])
        self.model=model
        model.setRootPath('C:\Myfolder')
        model.setNameFilterDisables(False)
        # view = QtGui.QTreeView()
        view = MyTree()
        view.setModel(model)
        view.setRootIndex(model.index("C:\\Users")) 
        self.setCentralWidget(view)


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    myview = Myview()
    myview.show()

    sys.exit(app.exec_())







