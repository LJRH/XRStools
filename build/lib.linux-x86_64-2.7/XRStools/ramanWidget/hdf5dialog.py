#!/usr/bin/env python
# coding: utf-8
# /*##########################################################################
#
# Copyright (c) 2016-2017 European Synchrotron Radiation Facility
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
"""Qt Hdf5 widget examples
"""

import logging
import sys
import tempfile
import numpy
import h5py
# from .. import ui


from XRStools.installation_dir import installation_dir
import os

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]

logging.basicConfig()
_logger = logging.getLogger("customHdf5TreeModel")
"""Module logger"""

from silx.gui import qt
import silx.gui.hdf5
from silx.gui.data.DataViewerFrame import DataViewerFrame
from silx.gui.widgets.ThreadPoolPushButton import ThreadPoolPushButton
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel 



# class Hdf5TreeModel(originalH5model):
#     f2close=[]
#     def appendFile(self, filename):
#         try:
#             h5file = silx_io.open(filename)
#             f2close = f2close+[h5file]
#             self.insertH5pyObject(h5file, row=-1)
#         except IOError:
#             print ("File '%s' can't be read.", filename)
#             raise
#     def chiudi(self):
#         for obj in self.f2close:
#             obj.close()

def readme_to_dict(README):
    res={}
    sl=README.split("\n")
    for l in sl:
        ll=l.split(":")
        if len(ll)>1:
            nm = ll[0]
            doc=""
            for t in ll[1:]:
                doc=doc+t+":"
            res[nm.strip()]=doc
    return res


    
class CustomTooltips(qt.QIdentityProxyModel):
    """Custom the tooltip of the model by composition.
    It is a very stable way to custom it cause it uses the Qt API. Then it will
    not change according to the version of Silx.
    But it is not well integrated if you only want to add custom fields to the
    default tooltips.
    """

    def data(self, index, role=qt.Qt.DisplayRole):
        if role == qt.Qt.ToolTipRole:

            # Reach information from the node
            sourceIndex = self.mapToSource(index)
            sourceModel = self.sourceModel()
            originalTooltip = sourceModel.data(sourceIndex, qt.Qt.ToolTipRole)
            originalH5pyObject = sourceModel.data(sourceIndex, Hdf5TreeModel.H5PY_OBJECT_ROLE)

            # We can filter according to the column
            if sourceIndex.column() == Hdf5TreeModel.TYPE_COLUMN:
                return super(CustomTooltips, self).data(index, role)

            # Let's create our own tooltips
            template = u"""<html>
            <dl>
            <dt><b>Original</b></dt><dd>{original}</dd>
            <dt><b>Parent name</b></dt><dd>{parent_name}</dd>
            <dt><b>Name</b></dt><dd>{name}</dd>
            <dt><b>Power of 2</b></dt><dd>{pow_of_2}</dd>
            <dt><b>Help </b></dt><dd>{README}</dd>

            </dl>
            </html>
            """
           
            try:
                README = originalH5pyObject.parent["README"].value
            except:
                README = "NA"
            
            try:
                data = originalH5pyObject[()]
                if data.size <= 10:
                    result = data ** 2
                else:
                    result = "..."
            except Exception:
                result = "NA"
                
            name=originalH5pyObject.name
            ipos = name.rfind("/")
            if ipos!=-1:
                name=name[ipos+1:]

            
            print(" README ", README )
            print (" name ", name)
            
            if name =="README":
                pass
            else:
                dico = readme_to_dict(README)
                print(" dico " , dico)
                if name in dico:
                    README = dico[name]
                else:
                    README = ""
                
            info = dict(
                original=originalTooltip,
                parent_name=originalH5pyObject.parent.name,
                name=originalH5pyObject.name,
                pow_of_2=str(result),
                README = README
            )


            
            return template.format(**info)

        elif ( not self.RO and role == qt.Qt.BackgroundRole):
            
            
            s=self.mled.toPlainText()
            sl =s.split("\n")

            sourceIndex = self.mapToSource(index)
            sourceModel = self.sourceModel()

            myn = sourceModel.data(sourceIndex, Hdf5TreeModel.H5PY_OBJECT_ROLE )
            sl = self.selected_h5_objects
            if(myn in sl):
                return qt.QColor(255, 0, 0)



        
        return super(CustomTooltips, self).data(index, role)


_file_cache = {}


def get_hdf5_with_all_types():
    global _file_cache
    ID = "alltypes"
    if ID in _file_cache:
        return _file_cache[ID].name

    tmp = tempfile.NamedTemporaryFile(prefix=ID + "_", suffix=".h5", delete=True)
    tmp.file.close()
    h5 = h5py.File(tmp.name, "w")

    g = h5.create_group("arrays")
    g.create_dataset("scalar", data=10)
    g.create_dataset("list", data=numpy.arange(10))
    base_image = numpy.arange(10**2).reshape(10, 10)
    images = [base_image,
              base_image.T,
              base_image.size - 1 - base_image,
              base_image.size - 1 - base_image.T]
    dtype = images[0].dtype
    data = numpy.empty((10 * 10, 10, 10), dtype=dtype)
    for i in range(10 * 10):
        data[i] = images[i % 4]
    data.shape = 10, 10, 10, 10
    g.create_dataset("image", data=data[0, 0])
    g.create_dataset("cube", data=data[0])
    g.create_dataset("hypercube", data=data)
    g = h5.create_group("dtypes")
    g.create_dataset("int32", data=numpy.int32(10))
    g.create_dataset("int64", data=numpy.int64(10))
    g.create_dataset("float32", data=numpy.float32(10))
    g.create_dataset("float64", data=numpy.float64(10))
    g.create_dataset("string_", data=numpy.string_("Hi!"))
    # g.create_dataset("string0",data=numpy.string0("Hi!\x00"))

    g.create_dataset("bool", data=True)
    g.create_dataset("bool2", data=False)
    h5.close()

    _file_cache[ID] = tmp
    return tmp.name


         
# @ui.UILoadable
class hdf5dialog(qt.QDialog):
    """
    This window show an example of use of a Hdf5TreeView.
    The tree is initialized with a list of filenames. A panel allow to play
    with internal property configuration of the widget, and a text screen
    allow to display events.
    """
    f2close = None
    def __del__(self):
        
        if self.f2close is not None:
            
            self.f2close.close()
        
    def __init__(self,filenames=None, storage=None,  parent=None, groupname = "",RO=True):

        super(hdf5dialog  , self).__init__(parent)
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path ,  "hdf5dialog.ui" ), self)
   
        
        self.RO=RO
    # def closeEvent(self,event):
    #     if self.hdf5File is not None:
    #         self.hdf5File.close()
    #     event.accept()

        
        # qt.QMainWindow.__init__(self)

        self.storage=storage
        
        
        self.setWindowTitle("Silx HDF5 widget example")

        self.__asyncload = False
        self.__treeview = silx.gui.hdf5.Hdf5TreeView(self)
        """Silx HDF5 TreeView"""

        self.__sourceModel = self.__treeview.model()
        """Store the source model"""

        self.__text = qt.QTextEdit(self)
        """Widget displaying information"""

        self.__dataViewer = DataViewerFrame(self)
        
        vSpliter = qt.QSplitter(qt.Qt.Vertical)
        vSpliter.addWidget(self.__dataViewer)
        vSpliter.addWidget(self.__text)
        vSpliter.setSizes([10, 0])

        spliter = qt.QSplitter(self)
        spliter.addWidget(self.__treeview)
        spliter.addWidget(vSpliter)
        spliter.setStretchFactor(1, 1)

        main_panel = qt.QWidget(self)
        layout = qt.QVBoxLayout()
        layout.addWidget(spliter)
        layout.addWidget(self.createTreeViewConfigurationPanel(self, self.__treeview))
        layout.setStretchFactor(spliter, 1)
        main_panel.setLayout(layout)

        # self.setLayout( qt.QVBoxLayout())

        self.mainLayout =  self.verticalLayout_2
        self.mainLayout.addWidget(main_panel)

        hl = qt.QWidget()
        hl.setLayout(qt.QHBoxLayout())
        self.varname = qt.QLineEdit("Set_here_variable_name")
        hl.layout().addWidget(self.varname)
        self.verticalLayout.addWidget(hl)
        if not RO:
        
            button = qt.QPushButton("Add Selected to the list")
            button.clicked.connect(self.__AddToList)
            self.verticalLayout.layout().addWidget(button)

            button = qt.QPushButton("Remove Selected from the list")
            button.clicked.connect(self.__RemoveFromList)
            self.verticalLayout.layout().addWidget(button)



            hl = qt.QWidget()
            hl.setLayout(qt.QHBoxLayout())
            button = qt.QPushButton("Exit Cancelling Selecteds")
            button.clicked.connect(self.__WriteCancel)
            hl.layout().addWidget(button)

            # button = qt.QPushButton("Exit Keeping Selecteds")
            # button.clicked.connect(self.__WriteKeep)
            # hl.layout().addWidget(button)

            button = qt.QPushButton("Discard")
            button.clicked.connect(self.reject)
            hl.layout().addWidget(button)

            self.verticalLayout.addWidget(hl)


            self.buttonBox.clear()

            self.mled = qt.QTextEdit()
            self.selected_h5_objects = []
            self.mled.setReadOnly(True)
            self.verticalLayout.layout().addWidget(self.mled)

        
        # self.layout().addWidget(main_panel)
        # append all files to the tree
        self.filenames = filenames
        for file_name in filenames:
            f = h5py.File(file_name,"r")
            self.f2close = f
            # self.__treeview.findHdf5TreeModel().appendFile(file_name)
            self.__treeview.findHdf5TreeModel().insertH5pyObject(f)

            
        self.__treeview.activated.connect(self.displayData)
        self.__treeview.doubleClicked.connect(self.toConsole)
        self.__treeview.clicked.connect(self.toConsole1)
        self.__useCustomLabel()

        
        ## ce qui est ci dessous ca va marcher pour le premier fichier
        # si non il faudra chercher le fichier


         
        
       # h5 = self.__treeview.model().data(indice, Hdf5TreeModel.H5PY_OBJECT_ROLE)
        if groupname !="":

            if 0 and hasattr(   self.__treeview, "setSelectedH5Node" ):
                print(" DDDDDDDDDDDDDDDDDDDDISPONIBILE PUOI CANCELLARE CODICE PER SELECTION\n"*1, f[groupname], self.__treeview.__class__.__module__)
                self.__treeview.setSelectedH5Node(f[groupname] )
            else:
                indice = self.__treeview.model().index(0,0,qt.QModelIndex())  # parce que on va au premier fichier
                self.__treeview.expand(indice) 
                groupname.replace("//","/")
                gn_l = groupname.split("/")
                i=None
                for tn in gn_l:
                    if tn=="":
                        continue
                    nl = [      ]
                    for k in range( self.__treeview.model().rowCount(indice)    ):
                        ind_tmp =   self.__treeview.model().index(k,0,indice)
                        nl.append(   str(self.__treeview.model().data(ind_tmp)    ) )
                    if tn not in  nl:
                        break
                    i = nl.index( tn  )
                    indice = self.__treeview.model().index(i,0,indice)
                    print (self.__treeview.model().rowCount(indice))
                    self.__treeview.expand(indice)
                    # h5=h5[tn]

                if i is not None:
                    self.__treeview.setCurrentIndex(indice)
    def __AddToList(self):
        s=self.mled.toPlainText()
        sl =s.split("\n")

        selecteds = list(self.__treeview.selectedH5Nodes())
        ns = ""
        for sel in selecteds:
            if str(sel.name) not in sl:
                ns = ns+ sel.name +"\n"
                self.selected_h5_objects.append(sel.h5py_object)
        s=str(s)+ ns
        self.mled.setPlainText(s)
        
    def __RemoveFromList(self):
        s=self.mled.toPlainText()
        sl =s.split("\n")

        selecteds = list(self.__treeview.selectedH5Nodes())
        for sel in selecteds:
            if str(sel.name)  in sl:
                i = sl.index( str(sel.name)  )
                del sl[i]
                self.selected_h5_objects.remove(sel.h5py_object)
        s=""
        for t in sl:
            s=str(s)+ t+"\n"
        self.mled.setPlainText(s)

        
    def __WriteCancel(self):
        msgBox = qt.QMessageBox ()
        msgBox.setText("You are going change the hdf5 file");
        msgBox.setInformativeText(" Do you want to proceed?");
        msgBox.setStandardButtons(qt.QMessageBox.Ok |  qt.QMessageBox.Cancel);
        msgBox.setDefaultButton(qt.QMessageBox.Cancel);
        ret = msgBox.exec_();
        res = (ret==qt.QMessageBox.Ok)
        
        if res:
            s  = self.mled.toPlainText()
            sl = s.split("\n")

            self.storage[0]=sl
            
            # target = h5py.File(self.filenames[0],"a")

            # for l in sl:
                
            #     if not  len(l) : continue
            #     if not subcontained(l,sl ):
            #         print( " CANCELLO :", l)
            #         del target[l]
            # target.flush()
            # target.close()
            # target = None
            self.accept()
        else:
            pass
    
    def __WriteKeep(self):
  
            self.reject()


    


        
    def toConsole1(self, index):
        originalH5pyObject = self.__treeview.model().data(index, Hdf5TreeModel.H5PY_OBJECT_ROLE)
        self.storage[0] = originalH5pyObject.name
        
    def toConsole(self, index):
        

        # sourceIndex = self.mapToSource(index)
        # sourceModel = self.sourceModel()
        
        originalH5pyObject = self.__treeview.model().data(index, Hdf5TreeModel.H5PY_OBJECT_ROLE)
        self.storage[0] = originalH5pyObject.name

        # originalH5pyObject = sourceModel.data(sourceIndex, Hdf5TreeModel.H5PY_OBJECT_ROLE)
        
        vn = str(self.varname.text())
        if len(vn)==0:
            vn="grabbed"
        self.consoleAction (originalH5pyObject,  vn)

    def consoleAction(self,h5group, vn ):

        from . import ShadokWidget
        from . import  Wizard_safe_module

        groupname = h5group.name
        h5 = h5group.file
        runit, oggetto = Wizard_safe_module.make_obj( h5group ,h5,  groupname     )


        if runit is not None:
            runit(oggetto)

        update = {vn:oggetto}
        if runit is None:
            ShadokWidget.ShadokWidget.LaunchConsole()
            ShadokWidget.ShadokWidget.my_console[0].updateNamespace(update)
            ShadokWidget.ShadokWidget.my_console[0].showMessage("# %s is now in namespace "% list(update.keys())[0])


    def displayData(self):
        """Called to update the dataviewer with the selected data.
        """
        selected = list(self.__treeview.selectedH5Nodes())
        if len(selected) == 1:
            # Update the viewer for a single selection
            data = selected[0]
            # data is a hdf5.H5Node object
            # data.h5py_object is a Group/Dataset object (from h5py, spech5, fabioh5)
            # The dataviewer can display both
            self.__dataViewer.setData(data)

    def __fileCreated(self, filename):
        if self.__asyncload:
            self.__treeview.findHdf5TreeModel().insertFileAsync(filename)
        else:
            self.__treeview.findHdf5TreeModel().insertFile(filename)

    def __hdf5ComboChanged(self, index):
        function = self.__hdf5Combo.itemData(index)
        self.__createHdf5Button.setCallable(function)

    def __edfComboChanged(self, index):
        function = self.__edfCombo.itemData(index)
        self.__createEdfButton.setCallable(function)

    def __useCustomLabel(self):
        customModel = CustomTooltips(self.__treeview)
        customModel.RO = self.RO
        customModel.setSourceModel(self.__sourceModel)
        
        if not self.RO :
            customModel.mled = self.mled            
            customModel.selected_h5_objects = self.selected_h5_objects
            
        self.__treeview.setModel(customModel)

    def __useOriginalModel(self):
        self.__treeview.setModel(self.__sourceModel)

    def createTreeViewConfigurationPanel(self, parent, treeview):
        """Create a configuration panel to allow to play with widget states"""
        panel = qt.QWidget(parent)
        panel.setLayout(qt.QHBoxLayout())

        # content = qt.QGroupBox("Create HDF5", panel)
        #  content.setLayout(qt.QVBoxLayout())
        # panel.layout().addWidget(content)

        # combo = qt.QComboBox()
        # combo.addItem("Containing all types", get_hdf5_with_all_types)
        # combo.activated.connect(self.__hdf5ComboChanged)
        # content.layout().addWidget(combo)

        # button = ThreadPoolPushButton(content, text="Create")
        # button.setCallable(combo.itemData(combo.currentIndex()))
        # button.succeeded.connect(self.__fileCreated)
        # content.layout().addWidget(button)

        # self.__hdf5Combo = combo
        # self.__createHdf5Button = button

        # content.layout().addStretch(1)

        self.varname.hide()
        if 0:
            self.varname = qt.QLineEdit("Set_here_variable_name")
            panel.layout().addWidget(self.varname)        
            button = qt.QPushButton("Delete")
            button.clicked.connect(self.__Cancel)
            panel.layout().addWidget(button)

        # option = qt.QGroupBox("Custom model", panel)
        # option.setLayout(qt.QVBoxLayout())
        # panel.layout().addWidget(option)

        # button = qt.QPushButton("Original model")
        # button.clicked.connect(self.__useOriginalModel)
        # option.layout().addWidget(button)

        # button = qt.QPushButton("Custom tooltips by composition")
        # button.clicked.connect(self.__useCustomLabel)
        # option.layout().addWidget(button)

        # option.layout().addStretch(1)

        panel.layout().addStretch(1)
        return panel


def main(filenames):
    """
    :param filenames: list of file paths
    """
    app = qt.QApplication([])
    sys.excepthook = qt.exceptionHandler
    window = Hdf5TreeViewExample(filenames)
    window.show()
    result = app.exec_()
    # remove ending warnings relative to QTimer
    app.deleteLater()
    sys.exit(result)


if __name__ == "__main__":
    main(sys.argv[1:])
