from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui


import collections
from . import hdf5dialog
import os


from XRStools.installation_dir import installation_dir

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]



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


# @ui.UILoadable
class source(QtGui.QWidget) :
    def __init__(self, parent=None):

            
        super( source, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "source.ui" ), self)


        extraActions = collections.OrderedDict([("Browse",self.browse),("pymcaview" ,self.pymcaview)] )
        self.lineEdit_specfileName.contextMenuEvent  = Functor_contextMenuEvent(self.lineEdit_specfileName , extraActions = extraActions )

        extraActions = collections.OrderedDict([("Browse",self.browse_hdf5),("pymcaview" ,self.pymcaview)] )
        self.lineEdit_roiAddress.contextMenuEvent  = Functor_contextMenuEvent(self.lineEdit_roiAddress , extraActions = extraActions )
        self.load_user_input  = {}

    def update_user_input(self, uinput):
        self.load_user_input.update(uinput)
        if "sf" in self.load_user_input :
            self.lineEdit_specfileName.setText(self.load_user_input["sf"])
        if "roi_spectral" in self.load_user_input :
            self.lineEdit_roiAddress.setText(self.load_user_input["roi_spectral"])
    def browse(self):
        filename = self.lineEdit_specfileName.text()
        if len(filename):                 
            filename =  QtGui.QFileDialog.getOpenFileName(None, "select", filename)
        else:
            filename = QtGui.QFileDialog.getOpenFileName(None, "select", )
            
        if isinstance(filename, tuple):
            filename = filename[0]
            
        filename=str(filename)
        if len(filename):
            self.lineEdit_specfileName.setText(filename)
            
    def pymcaview(self):
        filename = str(self.lineEdit_specfileName.text())
        os.system("pymca -f '%s' &"%filename)


    def browse_hdf5(self):
        filename = str(self.lineEdit_roiAddress.text())
        fn, dg = hdf5_filedialog(filename, 1)
        if fn is not None:
            filename=str(fn)
            if dg is not None:
                filename = filename +":"+str(dg)
            self.lineEdit_roiAddress.setText(filename)


    def set_selection(self, selection):

        self.lineEdit_specfileName.setText(str(selection[0]))
        self.lineEdit_roiAddress.setText(str(selection[1]))
            
    def get_selection(self):
        selection = []
        selection.append(str(self.lineEdit_specfileName.text()))
        selection.append(str(self.lineEdit_roiAddress.text()))
        return selection

def split_hdf5_address(dataadress):
    pos = dataadress.rfind(":")
    if ( pos==-1):
        return None
    filename, groupname = dataadress[:pos], dataadress[pos+1:]
    return filename, groupname 

def hdf5_filedialog(hint, fileme=1, myaction = None, modal=1):
    sphint= split_hdf5_address(hint)
    groupname="/"
    if sphint is not None and len(sphint):
        hint , groupname = sphint
    if fileme in [1]:
        if len(hint):
            # QtGui.QApplication.instance().processEvents()
            # dialog = QWS.QFileDialog(parent)
            # dialog.setFileMode(QWS.QFileDialog.ExistingFile)
            # dialog.setNameFilter("hdf5 (*h5)\nall files ( * )"        )
            # dialog.selectFile(hint)
            # dialog.show()
            # if (dialog.exec()):
            #     filename = dialog.selectedFiles()
            
            filename =  QtGui.QFileDialog.getOpenFileName( None,'Open hdf5 file and then choose groupname', hint,filter="hdf5 (*h5)\nall files ( * )"  )
        else:
            filename =  QtGui.QFileDialog.getOpenFileName(None,'Open hdf5 file and then choose groupname',filter="hdf5 (*h5)\nall files ( * )"  )
    else:
        if len(hint):
            filename =  QtGui.QFileDialog.getSaveFileName(None,'Open hdf5 file and then choose groupname', hint,filter="hdf5 (*h5)\nall files ( * )"  )
        else:
            filename =  QtGui.QFileDialog.getSaveFileName(None,'Open hdf5 file and then choose groupname',filter="hdf5 (*h5)\nall files ( * )"  )

    if isinstance(filename, tuple):
        filename = filename[0]
            
    filename=str(filename)
    if len(filename) :
        ret =None
        if os.path.exists(filename):
            storage=[None]
            window = hdf5dialog.hdf5dialog(filenames = [filename], storage = storage, groupname = groupname)
            window.setModal(True)
            ret = window.exec_()
        if ret:
            name =  storage[0]
            return filename, name
        else:
            # return filename , None
            return None , None
    else:
        
        return None,None

    

        
            

if __name__ =="__main__":
    app=Qt.QApplication([])
    w = source()
    w.show()
    w.set_selection([   "/data/id20/inhouse/data/run5_17/run7_ihr/hydra",   "/mntdirect/_scisoft/users/mirone/WORKS/Christoph/XRStoolsSuperResolution/sandbox/spectral_myroi.h5"   ])
    app.exec_()
    print("===== SELECTION =============== " )
    print(w.get_selection())
    
        
