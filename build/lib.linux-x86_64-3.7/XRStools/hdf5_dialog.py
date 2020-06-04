from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from XRStools.installation_dir import installation_dir
my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]



class hdf5dialog(Qt.QDialog):
    def __init__(self, parent=None):
        
        super(hdf5dialog , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "hdf5dialog.ui" ), self)


def hdf5dialog():
    filename =  Qt.QFileDialog.getOpenFileName(None,'Open hdf5 file with rois',filter="hdf5 (*h5)\nall files ( * )"  )
        
    if isinstance(filename, tuple):
        filename = filename[0]    
    
    filename=str(filename)
    if len(filename):
        import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget
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
            return filename +":" + name
        else:
            return filename +":" + name
    else:
        return filename
            







