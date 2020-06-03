from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui

# from .. import ui
import collections

from XRStools.installation_dir import installation_dir

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]

# @ui.UILoadable
class acquisition(QtGui.QWidget) :
    def __init__(self, parent=None):

        
        super( acquisition , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path ,  "acquisition.ui" ), self)

        
        self.checkBox_includeElastic.setChecked(True)
        extraActions = collections.OrderedDict([("Browse",self.browse)] )
        self.lineEdit_outputPrefix.contextMenuEvent  = Functor_contextMenuEvent(self.lineEdit_outputPrefix , extraActions = extraActions )
	
    def set_selection( self, selection ):
        if selection[3] is None:
            return
        self.comboBox_method.setCurrentIndex(self.comboBox_method.findText( selection[0]) )
        self.spinBox_scan.setValue(   selection[1] ) 
        self.checkBox_includeElastic.setChecked( selection[2]   )
        self.lineEdit_outputPrefix.setText(selection[3])
    def get_selection( self ):
        res = [     str(self.comboBox_method.currentText())       ,   self.spinBox_scan.value() , self.checkBox_includeElastic.isChecked() , str(self.lineEdit_outputPrefix.text())  ]
        return res

    def browse(self):
        filename = self.lineEdit_outputPrefix.text()
        if len(filename):                 
            filename =  QtGui.QFileDialog.getSaveFileName(None, "select", filename)
        else:
            filename = QtGui.QFileDialog.getSaveFileName(None, "select", )

        if isinstance(filename, tuple):
            filename = filename[0]

        filename=str(filename)
        if len(filename):
            self.lineEdit_outputPrefix.setText(filename)



    
if __name__ =="__main__":
    app=Qt.QApplication([])
    w = acquisition()
    w.set_selection(  ["pixel" ,611  ]   )
    w.show()
    app.exec_()

    print(" SELECTION =")
    print(w.get_selection())
            
class Functor_contextMenuEvent:
    def __init__(self,  ql, extraActions):
        self.ql=ql
        self.extraActions = extraActions
        

        
    def __call__(self, event):

        if not hasattr(self.ql,"createStandardContextMenu"):
            self.mem=[]
        else:
            menu = self.ql.createStandardContextMenu()
            self.mem=[]
            menu.insertSeparator(menu.actions()[0])
            menu.insertSeparator(menu.actions()[0])
            for na, aa in list(self.extraActions.items())[::-1]:
                testAction = Qt.QAction( na   ,None);
                self.mem.append(testAction)
                menu.insertAction(  menu.actions()[0]   , testAction)
                testAction.triggered.connect( aa )
            menu.exec_(event.globalPos());
            del menu
