from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from PyQt4 import Qt, QtGui, QtCore
import sys
from six.moves import range
if sys.platform =="win32":
    rootdir = "c:\\users\\aless\\src\\xrstoolssuperresolution\\xrstools\\WIZARD\\methods"
else:
    rootdir = "/scisoft/users/mirone/WORKS/Christoph/XRStoolsSuperResolution/XRStools/WIZARD/methods"
    
   
class MyLoop(QtCore.QEventLoop):
    def __init__(self):
        QtCore.QEventLoop.__init__(self)

    def quit(self,path):
        print( " loaded " , path)
        QtCore.QEventLoop.quit(self)

My_loop = MyLoop()
        
class   MyTree( QtGui.QTreeView):
    def __init__(self, *args):
        QtGui.QTreeView.__init__(self, *args)
        self.clicked[QtCore.QModelIndex].connect(self.itemClicked)

    def itemClicked(self, modelIndex):

            event ="itemDoubleClicked"
            print( event, modelIndex)
            
            index = self.selectedIndexes()[0]


            # if not( sys.platform =="win32"):
            #     index = self.model().mapToSource(index )
            pippo = index.model().data(index)
            
            print( pippo)
            
            pippo = index.model().filePath(index)
            print( str(pippo))
            

class MethodsProxyModel(Qt.QAbstractProxyModel):
    def __init__(self, parent=None, model=None, rootIndexes=[]):
        self.rootIndexes = rootIndexes
        self.mapping = {}
        self.model=model
        Qt.QSortFilterProxyModel.__init__(self, model)
        # self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self.fixModel()
        # self.emit(QtCore.SIGNAL("layoutChanged()"))
        self.setSourceModel(model )
    

    def fixModel(self):
        print( " fix model ")
        self.mapping.clear()
       
        for ir, root in enumerate(self.rootIndexes):
            print( root)
            
            sourceIndex = self.model.index(root)
 

            print( " sourceIndex.internalPointer()  " , sourceIndex.internalId())
 
            proxy = self.createIndex(ir, 0, sourceIndex.internalId())


            tmp = Qt.QPersistentModelIndex(sourceIndex)
            self.mapping[tmp] = proxy            



            print( " STACK AGGIUNGO ", sourceIndex.data().toString(), "  con rowcount ", self.model.rowCount(sourceIndex  ) )
            todo = list(range(  self.model.rowCount(sourceIndex  )))
            stack = [ ( iter(todo) , sourceIndex    )   ]
            
            print( stack)
            
            while len(stack):
                try:
                    it, si = stack[-1]
                    ir = next(it)
                    for ic in range(4):
                        print( ic)
                        sourceIndex = self.model.index( ir,ic, si)

                        if sourceIndex and sourceIndex != Qt.QModelIndex():

                            print( " AGGIUNGO per ir ic ", ir, ic,  " ITEM sourceIndex   ", sourceIndex.data().toString())

                            proxy = self.createIndex(ir, ic, sourceIndex.internalId())
                            tmp = Qt.QPersistentModelIndex(sourceIndex)
                            self.mapping[tmp] = proxy

                            if ic==0:
                                # if sourceIndex.model().isDir(sourceIndex):
                                if self.model.isDir(sourceIndex):
                                    percorso = self.model.filePath( sourceIndex  )
                                    self.model.setRootPath(percorso)
                                    My_loop.exec_()
                                
                                print( " STACK AGGIUNGO ", sourceIndex.data().toString(), "  con rowcount ", self.model.rowCount(sourceIndex  ) )
                                todo = list(range(  self.model.rowCount(sourceIndex  )))
                                stack.append( ( iter(todo) , sourceIndex    )   )
                except StopIteration:
                    stack = stack[:-1]

    def columnCount( self, index ):
        print( " column count " )
        return 4
    
    def rowCount( self, parent ):
        print( " row count A")
        sourceParent = Qt.QModelIndex()
        if parent.isValid():  
            sourceParent = self.mapToSource(parent)


        if  sourceParent ==Qt.QModelIndex():
            print( " row count B", 1)

            return 1
            
        count = self.model.rowCount(sourceParent)
        
        print( " row count C", count)

        return count


     
    def mapToSource(self, proxyIndex):
        print( "DB1 map to source ")
        index = proxyIndex
        coords = []
        while index != QtCore.QModelIndex():
            coords.append( [ index.row(), index.column()  ]  )
            index = index.parent()


        coords.reverse()
        print( "DB1 coords ", coords)


        if len(coords)==0:
            print( "DB1 RITORNO 0")
            return QtCore.QModelIndex()
        r,c = coords[0]
        if c>0 or r>= len(self.rootIndexes):
            print( "DB1 RITORNO 0b None")
            return QtCore.QModelIndex()
            # return None

        # index = self.rootIndexes[c]
        index = self.model.index(self.rootIndexes[0])
        print( "DB1  RADICE ", self.model.filePath(index))

        for r,c in coords[1:]:
            index = index.child(r,c)
            print( "DB1  in profondita ", self.model.filePath(index))

        print( "DB1  RITORNO ", index)
        print( "DB1  RITORNO ", self.model.filePath(index))
        return index
        

    def mapFromSource(self, sourceIndex):
        print( "DB1 map from source ")
        if not sourceIndex.isValid():
            print( "DB1 map from source  non valido  ")
            return Qt.QModelIndex()

        print( "DB1 map from sourc filePath  ", self.model.filePath(sourceIndex))
        pi = Qt.QPersistentModelIndex(sourceIndex)
        if pi in self.mapping:
            print( "DB1 mapfrom source ben in mapping")
            return self.mapping[pi]
        else:
            print( " DB1 non c'e'")
            return Qt.QModelIndex()
        


#     def mapFromSource(self, sourceIndex):
#         print " in map to proxy "
#         index = sourceIndex
#         coords = []
#         while index not in [ None, QtCore.QModelIndex()] and index not in self.rootIndexes  :
#             coords.append( [ index.row(), index.column()  ]  )
#             index = index.parent()
#         if index in [None, QtCore.QModelIndex()]:
#             return QtCore.QModelIndex()

#         r,c =  self.rootIndexes.index(index)    , 0
#         coords.append([r,c])
#         coords.reverse()

#         index = QtCore.QModelIndex()
#         for r,c in coords:
#             index = index.child(r,c)
#         return index


    def index(self, row, column, parent):
        print( " IndEX , parent", parent.data().toString())
        
        sourceParent= Qt.QModelIndex()
        if parent.isValid() :
            
            sourceParent = self.mapToSource(parent)
            print( "DEB2 good ", sourceParent.model().filePath(sourceParent), row, column)
        else:
            print( "DEB2 fault")
            if (row,column)==(0,0):
                sourceParent = self.model.index(self.rootIndexes[row])
                return self.mapFromSource(sourceParent)
            else:
                sourceParent = None
                return sourceParent
            #sourceParent = self.model.invisibleRootItem()


        print( " IndEX , parent", sourceParent.model().filePath(sourceParent))

            
        sind = self.model.index(row, column, sourceParent     )
        return self.mapFromSource(sind)

if(1):
        
        app = Qt.QApplication(sys.argv)

        widget = Qt.QWidget()
        
        button = Qt.QPushButton("Make only A1 + 'A' children visible", widget);
        tree_view = MyTree(widget)
        
        lay = Qt.QVBoxLayout(widget)
        lay.addWidget(button)
        
        
        # button.connect(button, QtCore.SIGNAL("clicked()"), proxyModel.hideEverythingButA1AndChildren)
        lay.addWidget(tree_view)


        methods_model = QtGui.QFileSystemModel()
        # methods_model.setNameFilters(["winfo_*.py"])
        # methods_model.setNameFilterDisables(False)


        
        methods_model.setRootPath(rootdir)



        methods_model.connect(methods_model,QtCore.SIGNAL("directoryLoaded(QString)"),  My_loop.quit )
        My_loop.exec_()



        folder = rootdir

        extrawpaths=[]
        if  ( 1) :
            # self.methods_proxy.setSourceModel(self.methods_model)
            rindexes=[]
            wpaths =extrawpaths+ [folder]
            for p in wpaths:
                rindexes.append(methods_model.index(p.lower()))

            # methods_proxy = MethodsProxyModel( model=methods_model, rootIndexes = rindexes )
            print( " CREO PROXY")

            methods_proxy = MethodsProxyModel( model=methods_model, rootIndexes = wpaths )
            # methods_proxy.fixModel()



            print( " CREO OK ")
            
            # tree_view.setRootIndex(Qt.QModelIndex()) 
 

            # self.methods_model.setRootPath("\\")
            tree_view.setModel(methods_proxy)
            # self.tree_view.expandAll()
        else: 
            methods_model.setRootPath(folder)
            tree_view.setModel(methods_model)
            tree_view.setRootIndex(methods_model.index(folder)) 
           
        tree_view.header().setStretchLastSection(False);
        tree_view.header().setResizeMode(0, Qt.QHeaderView.Stretch);  

        widget.show()
        app.processEvents(QtCore.QEventLoop.AllEvents) 


        app.exec_()

 

if (0):
    
    def displayMethodCreator(self, ):
        metodo = str(self.method_selection_list.currentText())
        interfaceDesc = copy.deepcopy( self.InterfaceDescription[metodo]  )        
        view = mainwindow2( interfaceDesc, self.viewsTab  )
        self.viewsTab.addTab(view, "Go create : %s"% metodo)


        
    def displayMethodCreatorFromTree(self, count=[0] ):
        count[0]+=1

        
        index = self.tree_view.selectedIndexes()[0]
        # if not (sys.platform == "win32"):
        #     index = self.tree_view.model().mapToSource(index )
        pippo = index.model().data(index)

        print( pippo)
    
        pippo = index.model().filePath(index)
        pippo = str(pippo)
        foo = imp.load_source('pippo'+str(count[0]),pippo)
        interfaceDesc  = copy.deepcopy( foo.getMethod(self.options) ) 
        interfaceDesc ["tobeimported"] = [   'pippo'+str(count[0]), pippo   ]
        view = mainwindow2( interfaceDesc, self.viewsTab  )
        self.viewsTab.addTab(view, "Go create : %s"% os.path.basename(  pippo  ))









