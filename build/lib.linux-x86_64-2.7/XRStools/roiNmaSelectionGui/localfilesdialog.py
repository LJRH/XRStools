from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from  PyQt4 import Qt
from silx.gui import qt as Qt

import os
import os.path




from XRStools.installation_dir import installation_dir
my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]



import PyMca5.PyMcaIO.specfilewrapper  as specfile
import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget


class hdf5dialog(Qt.QDialog):
    def __init__(self, parent=None):
           
        super( hdf5dialog, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "hdf5dialog.ui" ), self)

 

        
       



# this will use the ui module to load
# the structure specified in localfilesdialog.ui

class localfilesdialog(Qt.QDialog):
    def __init__(self, user_input, parent=None):


        super( localfilesdialog, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "localfilesdialog.ui" ), self)

 

        
        
        self.BrowseSpec_pushButton.clicked.connect(self.__onBrowseSpec)
        self.BrowseImage_pushButton.clicked.connect(self.__onBrowseFile)
        self.SpecFileName_lineEdit.textChanged.connect(self.__onChangeSpec)

        if "sf" in user_input:
            self.SpecFileName_lineEdit.setText(  os.path.dirname(user_input["sf"])  )   
        if "roi" in user_input:
            self.FileName_lineEdit.setText(  user_input["roi"]  )   
        
        # self.ScanNumber_spinBox.setMaximum(-1)

    def __onBrowseSpec(self):
        filename =  Qt.QFileDialog.getExistingDirectory()
        if isinstance(filename, tuple):
            filename = filename[0]
        
        # getOpenFileName()
        self.SpecFileName_lineEdit.setText(filename)

    def __onChangeSpec(self):
        filename =  str(self.SpecFileName_lineEdit.text())
        print( filename)

        try:
            s=specfile.Specfile(filename)
        except:
            s=None
        # if s is not None:
        #     ns = len(s)
        #     self.ScanNumber_spinBox.setMinimum(0)
        #     self.ScanNumber_spinBox.setMaximum(ns)
        # else:
        #     self.ScanNumber_spinBox.setMaximum(-1)

# s[450].alllabels()
# s[450]["ccdno"]
# s[424].datacol("ccdno")
# s[425].datacol("ccdno")
# ls
# history




   
    def __onBrowseFile(self):


        filename =  Qt.QFileDialog.getOpenFileName(None,'Open hdf5 file with rois',filter="hdf5 (*h5)\nall files ( * )"  )
        
        if isinstance(filename, tuple):
            filename = filename[0]

        filename=str(filename)
        print( filename)
        if len(filename):
            storage=[None]
            def mySlot(ddict):
                name = ddict["name"]
                storage[0]=name
            # browse
            self.__hdf5Dialog = hdf5dialog()
            self.__hdf5Dialog.setWindowTitle('Select a Group  containing roi_definitions by a double click')
            self.__hdf5Dialog.mainLayout =  self.__hdf5Dialog.verticalLayout_2
            fileModel = HDF5Widget.FileModel()
            fileView = HDF5Widget.HDF5Widget(fileModel)
            hdf5File = fileModel.openFile(filename)
            shiftsDataset = None
            fileView.sigHDF5WidgetSignal.connect(mySlot)
            self.__hdf5Dialog.mainLayout.addWidget(fileView)
            self.__hdf5Dialog.resize(400, 200)

            ret = self.__hdf5Dialog.exec_()
            print( ret)
            hdf5File.close()
            
            if ret:
                print( " Obtained " )
                name =  storage[0]        
                self.FileName_lineEdit.setText(filename+":"+name)
   



if __name__=="__main__":
    app=Qt.QApplication([])
    w = localfilesdialog()
    w.show()
    app.exec_()







