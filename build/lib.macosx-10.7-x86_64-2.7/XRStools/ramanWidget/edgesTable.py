from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui
import os


# from .. import ui
from .. import xrs_utilities
import collections



from XRStools.installation_dir import installation_dir
my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]


deletion_mode = 1
    

# @ui.UILoadable
class edgesTable(QtGui.QWidget) :
    def __init__(self, parent=None):
            
        super( edgesTable, self).__init__(parent)
        print(my_relativ_path) 
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "edgesTable.ui" ), self)

        

        self.gridLayout.setHorizontalSpacing(0)

        self.abstract=[]
        self.abstract_edges=[]

        self.render_abstract()
        
    def clear_layout(self):
        
        for l in self.bottoni_info.keys():
                self.gridLayout_edges.removeWidget(l)
                if deletion_mode:
                    l.hide()
                    l.setParent(None)
                else:
                    l.deleteLater()
                del l
                
        for l in self.abstract_edges:
            print( " RIMUOVO =====================", l[0].text())
            for t in  l:
                self.gridLayout_edges.removeWidget(t)
                if deletion_mode:
                    t.hide()
                    t.setParent(None)
                else:
                    t.deleteLater()
                del t
                
        for l in self.abstract:
            for t in  l:
                
                self.gridLayout.removeWidget(t)
                if deletion_mode:                
                    t.hide()
                    t.setParent(None)
                else:
                    t.deleteLater()
                del t


    def get_abstractN(self):
        res=[]
        for line in self.abstract:
            if len(line):
                aline = [    (str(line[0].text()))   ,      str(line[1].text())    ]
            res.append(aline)
            
        res_edges = {}
        for line in self.abstract_edges:
            if len(line):
                dl = collections.OrderedDict()
                el = str(line[0].text())[:-1] ## removing ":"
                res_edges[el]=dl
                for tok in line[1:]:
                    dl[ tok.text  ] = tok.isChecked()
            
        return res, res_edges
            
    def render_abstract(self):
        self.bottoni_info={}

        for irow,line in enumerate(self.abstract):
            for icol,tok in enumerate(line):
                self.gridLayout.addWidget(tok ,irow,icol)
                tok.show()
 
        qb = QtGui.QPushButton("+line",None)
        self.bottoni_info[ qb  ]="addrow"
        qb.clicked.connect(self.modifica)
        self.gridLayout.addWidget(qb,len(self.abstract),  0)
        self.plusLine = qb

        qb = QtGui.QPushButton("-line",None)
        self.bottoni_info[ qb  ]="removerow"
        qb.clicked.connect(self.modifica)
        self.gridLayout.addWidget(qb,len(self.abstract),  1)
        
        for irow,line in enumerate(self.abstract_edges):
            for icol,tok in enumerate(line):
                self.gridLayout_edges.addWidget(tok ,irow,icol)
                tok.show()


        
    def vLine(self):
        line = QtGui.QFrame()
        line.setGeometry(Qt.QRect())
        line.setFrameShape(QtGui.QFrame.VLine);
        line.setFrameShadow(QtGui.QFrame.Sunken);
        return line

    def set_abstractN(self,abst_stoichio, elements_edges ):
        self.clear_layout()
        self.abstract=[]
        self.abstract_edges=[]

        for line in abst_stoichio:
            if len(line):
                qlf  = QtGui.QLineEdit()
                qlf.setToolTip("Weight")
                qlf.setText(str( line[1]))
                qle =  QtGui.QLineEdit()
                qle.setToolTip("Formula")
                qle.setText( str(line[0]))
                qle.editingFinished.connect(self.onFormulaChange)
                al=[qle,qlf]
                # for ipos,val in enumerate(line[2:]):
                #     cb = MyCheckbox(str(ipos  ) )
                #     cb.setChecked(val>0)
                #     al.append(cb)
                self.abstract.append(al)
                
        for el, edglist in elements_edges.items():
            print("========================= " ,  el,)
            print("========================= " ,   edglist)
            ql = QtGui.QLabel(str(el)+":")
            f = Qt.QFont ( "Arial", 20, Qt.QFont.Bold);
            ql.setFont( f);
            al=[ql]
            
            for edg in edglist:
                cb = MyCheckbox(edg)
                cb.setChecked( edglist[edg]>0  )
                al.append( cb  )
            self.abstract_edges.append(al)

        self.render_abstract()
        
    def modifica(self):
        print( " RICEVUTO DA ", self.bottoni_info[self.sender()])
        if self.bottoni_info[self.sender()] == "addrow":
            abstr_stoichio, edges = self.get_abstractN()
            abstr_stoichio.append([   "formula here",1.0]    )
            self.set_abstractN(  abstr_stoichio, edges     )
            return
        if self.bottoni_info[self.sender()] == "removerow":
            abstr_stoichio, edges = self.get_abstractN()
            abstr_stoichio = abstr_stoichio[:-1]
            self.set_abstractN(abstr_stoichio, edges )
            return

    def onFormulaChange(self):
        # abstN, edges = self.get_abstractN()
        abstN, edges = self.get_selection()
        print("  abstN, edges  " ,  abstN, edges )
        self.set_selection( abstN, edges    )
        pass

        
    def set_selection(self, stoichio_arg, edges_arg):
        abst_stoichio=[]
        elements_edges = {}
        for formu,f in stoichio_arg:
            print( formu, f)
            abst_stoichio.append( [ formu,f      ]  )
            elements,weights = xrs_utilities.parseformula( formu)
            for el in elements:
                if xrs_utilities.element(el) is None:
                    break
            else:
                for el in elements:
                    if not (el in elements_edges) :
                        elements_edges[el] = collections.OrderedDict()
                    z=xrs_utilities.element(el)
                    if z<= 35:
                        edges =  ['pz', 'total', 'K',       'L1',      'L23',     'M1',      'M23',     'M45',     'N1',      'N23'    ]
                    else:
                        edges = ['pz', 'total', 'K', 'L1', 'L2', 'L3', 'M1', 'M2', 'M3', 'M4', 'M5', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'O1', 'O2', 'O3', 'O4', 'O5', 'P1', 'P2', 'P3']
                    for ed in edges:
                        elements_edges[el][ed]=0
                    
        for el in elements_edges:
            if el in edges_arg:
                elist = edges_arg[el]
                for ed in elist:
                    if ed in elements_edges[el]:
                        elements_edges[el][ed] = 1
        print(" CHIAMO SET ABSTRACT ", abst_stoichio, elements_edges  )
        self.set_abstractN(abst_stoichio, elements_edges     )
        
    def get_selection(self):
        abstN, edges = self.get_abstractN()
        stoichio_arg=abstN
        edges_arg = {}
        print ( " EDGES ", edges)
        for atom, dictio in edges.items():
            l=[]
            for ed,val in dictio.items():
                if val:
                    l.append(ed)
            if len(l):
                edges_arg[atom]=l
                
        return stoichio_arg, edges_arg


class MyCheckbox(QtGui.QWidget):
    def __init__(self, text):
        Qt.QWidget.__init__(self, None)
        self.setLayout(QtGui. QVBoxLayout()  )
        self.cb = QtGui.QCheckBox(self)
        self.layout().addWidget(self.cb)
        self.layout().setContentsMargins(2,0,2,0)
        self.layout().addWidget(QtGui.QLabel(text))
        self.text = text
        


        sp = Qt.QSizePolicy()
        sp.setVerticalPolicy(Qt.QSizePolicy.Fixed) #, Qt.QSizePolicy.Expanding
        sp.setVerticalStretch(0)
        self.setSizePolicy(sp)
        
    def setChecked(self,val):
        self.cb.setChecked(val)
    def isChecked(self):
        return self.cb.isChecked()
            

if __name__ =="__main__":
    app=Qt.QApplication([])
    w = edgesTable()
    # w.set_selection(  [[1.0, "H2O"] ] ,  {'O':['K'] }  )
    w.set_selection(  [[ "H2O",1.0],["Formula Here",1.0]] ,  {'O':['K'] }  )
    w.show()
    app.exec_()    
    print(" SELECTION ")
    print( w.get_selection())
     
