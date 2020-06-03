from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from  PyQt4 import Qt
from silx.gui import qt as Qt
import os

from XRStools.installation_dir import installation_dir
my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]



import PyMca5.PyMcaIO.specfile as specfile




# this will use the ui module to load
# the structure specified in localfilesdialog.ui

class localfilesdialog(Qt.QDialog):
    def __init__(self, user_input, parent=None):
           
        super(localfilesdialog , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "localfilesdialog.ui" ), self)

        


        
        
        self.BrowseSpec_pushButton.clicked.connect(self.__onBrowseSpec)
        self.BrowseImage_pushButton.clicked.connect(self.__onBrowseFile)
        self.SpecFileName_lineEdit.textChanged.connect(self.__onChangeSpec)
        self.ScanNumber_spinBox.setMaximum(-1)
        
        if "sf" in user_input:
            self.SpecFileName_lineEdit.setText(  user_input["sf"]  )   
        if "fn" in user_input:
            self.FileName_lineEdit.setText(  user_input["fn"]  )   
        
    def __onBrowseSpec(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename is not None:
            filename=str(filename)


        
        self.SpecFileName_lineEdit.setText(filename)

    def __onChangeSpec(self):
        filename =  str(self.SpecFileName_lineEdit.text())
        print( filename)

        try:
            s=specfile.Specfile(filename)
        except:
            s=None
        if s is not None:
            ns = len(s)
            self.ScanNumber_spinBox.setMinimum(0)
            self.ScanNumber_spinBox.setMaximum(ns)
        else:
            self.ScanNumber_spinBox.setMaximum(-1)

# s[450].alllabels()
# s[450]["ccdno"]
# s[424].datacol("ccdno")
# s[425].datacol("ccdno")
# ls
# history




   
    def __onBrowseFile(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        if isinstance(filename, tuple):
            filename = filename[0]
        
        self.FileName_lineEdit.setText(filename)
   



if __name__=="__main__":
    app=Qt.QApplication([])
    w = localfilesdialog()
    w.show()
    app.exec_()







