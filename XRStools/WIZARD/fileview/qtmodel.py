from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
    #include <QtGui>
    
import os
import sys
from PyQt4 import Qt, QtCore, QtGui
from six.moves import range

list = []
# QList<QStandardItem *> list;
     
#    class SortProxy : public QAbstractProxyModel
#    {
#    Q_OBJECT




print( " SORTYPROXY " )

class SortProxy(Qt.QAbstractProxyModel):

    def __init__(self, parent=0):
        Qt.QAbstractProxyModel.__init__(self, parent)
        self.hideThem = False
        # self.hideThem = True

        self.mapping = {}
        self.mapping2key = {}
        self.proxySourceParent  = {}
            
        self.fixModel()
        
    # SortProxy(QObject *parent = 0) : QAbstractProxyModel(parent), hideThem(false)
    # {
    # fixModel();
    # }
    
    def rowCount(self, parent):
        sourceParent = Qt.QModelIndex()
        if parent.isValid():  
            sourceParent = self.mapToSource(parent)
        count = 0
        for key,value in self.proxySourceParent.items():
            if  value == sourceParent:
                count+=1
        return count
    
    # int rowCount(const QModelIndex &parent) const
    # {
    #     QModelIndex sourceParent;
    #     if (parent.isValid())        sourceParent = mapToSource(parent);
    #     int count = 0;
    #     QMapIterator<QPersistentModelIndex, QPersistentModelIndex> it(proxySourceParent);
    #     while (it.hasNext()) {
    #             it.next();
    #             if (it.value() == sourceParent)
    #             count++;
    #     }
    #     return count;
    # }

    def columnCount(self, index):
        return 2
    
    
    # int columnCount(const QModelIndex &) const
    # {
    # return 1;
    # }

    def index(self, row, column, parent):
        sourceParent= Qt.QModelIndex()
        if parent.isValid() : sourceParent = self.mapToSource(parent)
        for key, value in self.proxySourceParent.items():
            if value == sourceParent and key.row() == row and key.column() == column:
                return key
        return Qt.QModelIndex()
     
    # QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const
    # {
    #     QModelIndex sourceParent;
    #     if (parent.isValid())
    #     sourceParent = mapToSource(parent);
    #     QMapIterator<QPersistentModelIndex, QPersistentModelIndex> it(proxySourceParent);
    #     while (it.hasNext()) {
    #             it.next();
    #             if (it.value() == sourceParent && it.key().row() == row &&
    #                 it.key().column() == column)
    #             return it.key();
    #     }
    #     return QModelIndex();
    # }


    def parent(self, child):
        mi = self.proxySourceParent[ child  ]
        if mi.isValid():
            return self.mapFromSource(mi)
        return Qt.QModelIndex()
    
    #     QModelIndex parent(const QModelIndex &child) const
    #     {
    #     QModelIndex mi = proxySourceParent.value(child);
    #     if (mi.isValid())
    #     return mapFromSource(mi);
    #     return QModelIndex();
    #     }
    
    def mapToSource( self,proxyIndex):

        if not proxyIndex.isValid():
            # print " RITORNO INVALIDO "
            return Qt.QModelIndex()
        for a,b in  self.mapping.items():
            # print a,b, ( b == proxyIndex)    , proxyIndex
            if b == proxyIndex:
                # print " RITORNO ", a
                return Qt.QModelIndex(a)

        print( " NON TROVATO ")
        return Qt.QModelIndex()


    #     QModelIndex mapToSource(const QModelIndex &proxyIndex) const
    #     {
    #         if (!proxyIndex.isValid())   return QModelIndex();
    #         return mapping.key(proxyIndex);
    #         }
    

    def mapFromSource(self, sourceIndex):
        if not sourceIndex.isValid():
            return Qt.QModelIndex()
        if sourceIndex in self.mapping:
            return self.mapping[sourceIndex]
        else:
            return Qt.QModelIndex()
        

    #     QModelIndex mapFromSource(const QModelIndex &sourceIndex) const{
    #         if (!sourceIndex.isValid())
    #         return QModelIndex();
    #         return mapping.value(sourceIndex);
    #         }
    

    def  hideEverythingButA1AndChildren(self):
        self.hideThem = not self.hideThem
        self.emit(QtCore.SIGNAL("layoutAboutToBeChanged()"))
        self.fixModel()
        self.emit(QtCore.SIGNAL("layoutChanged()"))



#     public slots:
#     void hideEverythingButA1AndChildren()
#     {
#         hideThem = !hideThem;
#         // Now we set up the proxy <-> source mappings
#         emit layoutAboutToBeChanged();
#         fixModel();
#         emit layoutChanged();
#        }
        
    

    
    def fixModel(self):
        # if self.amodel is not None:
        #     print " OK "
        #     raise
        self.mapping.clear()
        self.mapping2key.clear()
        self.proxySourceParent.clear()
        oldrow=0
        for si in list:
            
            if (self.hideThem):

                # if ( (not si.text().startsWith(Qt.QString("A"))) or  (not si.parent())):
                #    continue
                if ( (not si.text().startsWith(Qt.QString("A")))):
                    continue
                print( si.text(), si.row(), si.column())


                row = si.row()
                if row-oldrow>1:
                    row = oldrow+1

                proxy = self.createIndex(row, si.column(), si.index().internalPointer())
     
                oldrow = si.row()

                tmp = Qt.QPersistentModelIndex(si.index())
                self.mapping[tmp] = proxy
                self.mapping2key [ proxy ] = tmp
                
                sourceParent = Qt.QModelIndex()
                if (si.parent()) :
                    # if (si.parent().parent()) :
                    sourceParent = Qt.QPersistentModelIndex(si.parent().index())

                self.proxySourceParent[proxy] = sourceParent




            else :
                # if ( (not si.text().startsWith(Qt.QString("A"))) ):
                #     break
                proxy = self.createIndex(si.row(), si.column(), si.index().internalPointer())

                tmp = Qt.QPersistentModelIndex(si.index())
                self.mapping[tmp] = proxy
                self.mapping2key [ proxy ] = tmp
                
                sourceParent = Qt.QModelIndex()
                if (si.parent()) :
                    sourceParent = Qt.QPersistentModelIndex(si.parent().index())
                self.proxySourceParent[proxy] = sourceParent


print( " Tree " )



global proxyModel 

class Tree(Qt.QTreeView):

    def __init__(self, parent=None):
        global proxyModel 

        Qt.QTreeView.__init__(self, parent)



        self.amodel = None
        if 1:
            methods_model = QtGui.QFileSystemModel()
            rootdir = "/scisoft/users/mirone/WORKS/Christoph/XRStoolsSuperResolution/XRStools/WIZARD/methods"
            methods_model.setRootPath(rootdir)
            
            class MyLoop(QtCore.QEventLoop):
                def __init__(self):
                    QtCore.QEventLoop.__init__(self)

                def quit(self,path):
                    print( " loaded " , path)
                    QtCore.QEventLoop.quit(self)
            loop = MyLoop()
            methods_model.connect(methods_model,QtCore.SIGNAL("directoryLoaded(QString)"),  loop.quit )
            loop.exec_()


            self.amodel = methods_model




        
        sourceModel = Qt.QStandardItemModel(self)

        parentA = sourceModel.invisibleRootItem()

        for i in range(2):
            self.itemA = Qt.QStandardItem(Qt.QString("AAA %d"%i));
            self.itemAb = Qt.QStandardItem(Qt.QString("Abis %d"%i));
            parentA.appendRow([  self.itemA, self.itemAb])
            parentA = self.itemA
            list.append(self.itemA)
            list.append(self.itemAb)


        self.itemA = Qt.QStandardItem(Qt.QString("A 2"))
        parentA.appendRow(self.itemA)
        list.append(self.itemA)

        self.itemA3 = Qt.QStandardItem(Qt.QString("A 3"))
        list.append(self.itemA3)
        parentA.appendRow(self.itemA3)

        self.itemNonA = Qt.QStandardItem(Qt.QString("Non A"))
        list.append(self.itemNonA)
        parentA.appendRow(self.itemNonA)


        self.itemA4 = Qt.QStandardItem(Qt.QString("A 4"))
        list.append(self.itemA4)
        parentA.appendRow(self.itemA4)

        # parentB = sourceModel.invisibleRootItem()
        # self.itemB = Qt.QStandardItem(Qt.QString("B %d"%100 ))
        # parentB.appendRow(self.itemB)
        # list.append(self.itemB);


        parentB = sourceModel.invisibleRootItem()

        self.itemB = Qt.QStandardItem(Qt.QString("B %d"%0 ))
        parentB.appendRow(self.itemB)
        parentB = self.itemB;
        list.append(self.itemB)

        self.itemB = Qt.QStandardItem(Qt.QString("B %d"%1 ))
        parentB.appendRow(self.itemB)
        list.append(self.itemB)

        self.itemB = Qt.QStandardItem(Qt.QString("B %d"%2 ))
        parentB.appendRow(self.itemB)
        parentB = self.itemB;
        list.append(self.itemB)


        
        if(1):





            parentB = sourceModel.invisibleRootItem()

   

            for i in range(3): 
                self.itemB = Qt.QStandardItem(Qt.QString("B %d"%i ))
                parentB.appendRow(self.itemB)
                parentB = self.itemB;
                list.append(self.itemB);


            parentC = sourceModel.invisibleRootItem()

            for i in range(3): 
                self.itemC = Qt.QStandardItem(Qt.QString("C %d"%i ))
                parentC.appendRow(self.itemC)
                parentC = self.itemC
                list.append(self.itemC)

     
        proxyModel =  SortProxy(self)
        proxyModel.setSourceModel(sourceModel)
        self.setModel(proxyModel)
        self.expandAll()
     
     
    #include "main.moc"
def main()     :
    print( " qui " )
    app = Qt.QApplication(sys.argv)

    widget = Qt.QWidget()

    button = Qt.QPushButton("Make only A1 + 'A' children visible", widget);
    tree = Tree(widget)

    lay = Qt.QVBoxLayout(widget)
    lay.addWidget(button)


    button.connect(button, QtCore.SIGNAL("clicked()"), proxyModel.hideEverythingButA1AndChildren)
    lay.addWidget(tree)
    widget.show()
    app.exec_()
    print( " cucu ")
    
main()







