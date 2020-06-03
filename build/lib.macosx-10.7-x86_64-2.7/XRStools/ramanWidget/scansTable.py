from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui

from XRStools.installation_dir import installation_dir

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]



deletion_mode = 1
    
# @ui.UILoadable
class scansTable(QtGui.QWidget) :
    def __init__(self, parent=None):
            
        super( scansTable, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "scansTable.ui" ), self)

        


        self.abstract=[]

        self.render_abstract()
        
    def clear_layout(self):
        for t in list(self.bottoni_info.keys()):
            print( " levo ", t) 
            self.gridLayout.removeWidget(t)
            if deletion_mode:
                t.hide()
                t.setParent(None)
            else:
                t.deleteLater()
            del t
        # self.    bottoni_info={}
        for l in self.abstract:
            for t in  l:
                self.gridLayout.removeWidget(t)
                if deletion_mode:
                    t.hide()
                    t.setParent(None)
                else:
                    t.deleteLater()
                del t

        self.abstract=[]

        
    def get_abstractN(self):
        res=[]
        done = False
        
        if len(self.abstract):
            minS = 100000
            maxS=-1
            for line in self.abstract:
                if len(line):
                    aline = [ str(line[1].text())    ]
                    for t in line[3:]:
                        done= True
                        aline.append(  t.value() )
                        if t.value()> maxS:
                            maxS = t.value()
                        if t.value()< minS:
                            minS = t.value()
                    res.append(aline)
        if not done :
            minS = 1
            maxS= 100000
        
        return res, minS, maxS
            
    def render_abstract(self):
        self.bottoni_info={}

        maxlen=0
        for l in self.abstract:
            if  len(l)> maxlen:
                maxlen = len(l)
        for irow,line in enumerate(self.abstract):
            for icol,tok in enumerate(line):
                self.gridLayout.addWidget(tok ,irow,icol)
                tok.show()
                
            qb = QtGui.QPushButton("+",None)
            qb.setFixedWidth(20)
            self.bottoni_info[ qb  ]=+1, irow
            qb.clicked.connect(self.modifica)
            self.gridLayout.addWidget(qb,irow,  maxlen)
            qb = QtGui.QPushButton("-",None)
            qb.setFixedWidth(20)
            self.bottoni_info[ qb  ]=-1, irow
            qb.clicked.connect(self.modifica)
            self.gridLayout.addWidget(qb,irow,  maxlen+1)
        
        qb = QtGui.QPushButton("+line",None)
        self.plusLine = qb
        # qb.setFixedWidth(20)
        self.bottoni_info[ qb  ]="addrow"
        qb.clicked.connect(self.modifica)
        self.gridLayout.addWidget(qb,len(self.abstract),  0)
        qb = QtGui.QPushButton("-line",None)
        # qb.setFixedWidth(20)
        self.bottoni_info[ qb  ]="removerow"
        qb.clicked.connect(self.modifica)
        self.gridLayout.addWidget(qb,len(self.abstract),  1)
        print ( " aggiunto ", qb)
        line = self.vLine()
        self.gridLayout.addWidget(line ,len(self.abstract),  2)
        self.bottoni_info[ line ]=0


    def vLine(self):
        line = QtGui.QFrame()
        line.setGeometry(Qt.QRect())
        line.setFrameShape(QtGui.QFrame.VLine);
        line.setFrameShadow(QtGui.QFrame.Sunken);
        return line
    def set_abstractN(self, aN):
        self.clear_layout()
        self.abstract=[]
        for line in aN:
            if len(line):
                ql  = QtGui.QLabel("")
                qle =  QtGui.QLineEdit()
                qle.setText( line[0])
                al=[ql,qle, self.vLine() ]
                for N in line[1:]:
                    sb = QtGui.QSpinBox( )
                    sb.setMinimum(1)
                    sb.setMaximum(100000)
                    sb.setValue(N)
                    al.append(sb)
                self.abstract.append(al)
        self.render_abstract()
        
    def modifica(self):
        print( " RICEVUTO DA ", self.bottoni_info[self.sender()])
        if self.bottoni_info[self.sender()] == "addrow":
            resN, minS, maxS = self.get_abstractN()
            self.set_abstractN(resN+[["name here"]])
            
            return
        if self.bottoni_info[self.sender()] == "removerow":
            resN, minS, maxS = self.get_abstractN()
            self.set_abstractN(resN[:-1])
            return

        todo, irow = self.bottoni_info[self.sender()]
        if todo>0 :
            res, minS, maxS = self.get_abstractN()
            res[irow].append(minS)
            self.set_abstractN(res)
        else:
            res, minS, maxS = self.get_abstractN()
            
            if(len(res[irow])>2):
                del res[irow][-1]
                self.set_abstractN(res)
            
    def set_selection(self, abst):
        return self.set_abstractN(abst)
    
    def get_selection(self):
        return self.get_abstractN()
            

if __name__ =="__main__":
    app=Qt.QApplication([])
    w = scansTable()
    w.set_selection(  [["elastic" , 611,615,619,623 ] ,   [  'ok1', 612,616,620,624 ]  ]   )
    w.show()
    app.exec_()
    print(" SELECTION =")
    print(w.get_selection())
   
        
