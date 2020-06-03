from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from  PyQt4 import Qt, QtCore, QtGui
import os
from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui



# from .. import ui


from XRStools.installation_dir import installation_dir

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]



# @ui.UILoadable
class subsetTable(QtGui.QWidget) :
    def __init__(self, parent=None):

            
        super( subsetTable, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "subsetTable.ui" ), self)


        
        self.gridLayout.setHorizontalSpacing(0)

        self.abstract=[]

        self.render_abstract()
        
    def clear_layout(self):
        for t in list(self.bottoni_info.keys()):
            print( " levo ", t) 
            self.gridLayout.removeWidget(t)
            t.hide()
            t.setParent(None)
            del t
        for l in self.abstract:
            for t in  l:
                self.gridLayout.removeWidget(t)
                t.hide()
                t.setParent(None)
                del t

    def get_abstractN(self):
        res=[]
        for line in self.abstract:
            if len(line):
                aline = [    (str(line[0].text()))   ,      str(line[1].text())    ]
                for t in line[3:]:
                    aline.append(  t.isChecked() )
            res.append(aline)
        return res
            
    def render_abstract(self):
        self.bottoni_info={}

        maxlen=0
        for irow,line in enumerate(self.abstract):
            for icol,tok in enumerate(line):
                self.gridLayout.addWidget(tok ,irow,icol)
                tok.show()
 
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
                qlf  = QtGui.QLineEdit()
                qlf.setToolTip("Scaling")
                qlf.setText(str( line[0]))
                qle =  QtGui.QLineEdit()
                qle.setToolTip("SubSet Name")
                qle.setText( str(line[1]))
                

                al=[qlf,qle, self.vLine() ]
                for ipos,val in enumerate(line[2:]):
                    cb = MyCheckbox(str(ipos//12 +1)+"-"+str(ipos%12+1) )
                    cb.setChecked(val>0)
                    al.append(cb)
                self.abstract.append(al)
        self.render_abstract()
        
    def modifica(self):
        print( " RICEVUTO DA ", self.bottoni_info[self.sender()])
        if self.bottoni_info[self.sender()] == "addrow":
            abstrN = self.get_abstractN()
            abstrN.append([ 1.0,  "name here"]+[False]*72    )
            self.set_abstractN(abstrN)
            return
        if self.bottoni_info[self.sender()] == "removerow":
            abstrN = self.get_abstractN()[:-1]
            print ( " ECCO abstr ", abstrN)
            self.set_abstractN(abstrN)
            return

    def set_selection(self, sele_arg):
        sele=[]
        for line in sele_arg:
            l=[line[0] , line[1] ]+[0]*72
            for t in line[2:]:
                l[2+t]=True
            sele.append(l)
        self.set_abstractN(sele)
        
    def get_selection(self):
        sele = self.get_abstractN()
        sele_arg=[]
        for line in sele:
            l=[float(line[0]),line[1] ]
            for pos, t in enumerate(line[2:]):
                if t:
                    l.append(pos)
            sele_arg.append(l)
        return sele_arg


class MyCheckbox(QtGui.QWidget):
    def __init__(self, text):
        Qt.QWidget.__init__(self, None)
        self.setLayout(QtGui. QVBoxLayout()  )
        self.cb = QtGui.QCheckBox(self)
        self.layout().addWidget(self.cb)
        self.layout().setContentsMargins(2,0,2,0)
        self.layout().addWidget(QtGui.QLabel(text))
    def setChecked(self,val):
        self.cb.setChecked(val)
    def isChecked(self):
        return self.cb.isChecked()
            

if __name__ =="__main__":
    app=Qt.QApplication([])
    w = subsetTable()
    w.set_selection(  [[2.0, "lowq" , 0,1,2,3 ] ,   [ 2.0,  'middleq', 12,13,14,15 ]  ]   )
    w.show()
    app.exec_()
    
    print(" SELECTION ")
    print( w.get_selection())
     
