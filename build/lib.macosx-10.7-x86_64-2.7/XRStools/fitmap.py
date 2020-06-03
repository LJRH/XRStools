from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/env python

# embedding_in_qt4.py --- Simple Qt4 application embedding matplotlib canvases
#
# Copyright (C) 2005 Florent Rougon
#               2006 Darren Dale
#
# This file is an example program for matplotlib. It may be used and
# modified with no restriction; raw copies as well as modified versions
# may be distributed without limitation.

from __future__ import unicode_literals
import sys
import os
import random
import PyQt4
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends import qt_compat
from matplotlib.colors import LogNorm
import matplotlib.patches
import math
import h5py 
from six.moves import range

use_pyside = qt_compat.QT_API == qt_compat.QT_API_PYSIDE

from PyQt4 import QtGui, QtCore

from numpy import arange, sin, pi
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


import numpy as np
import glob


progname = os.path.basename(sys.argv[0])
progversion = "0.1"


class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.fig=fig
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)


        #
        FigureCanvas.__init__(self, fig)

        self.setParent(parent)

        self.compute_initial_figure()

        
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)


        self.cidpress = fig.canvas.mpl_connect(
            'button_press_event', self.on_press)
        
        self.cidmotion = fig.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)


    def compute_initial_figure(self):
        pass
    def on_press(self, event):
        pass
    def on_motion(self, event):
        pass


class Plot2D(MyMplCanvas):
    """2D Plot with interactions."""

    def change_action(self, i):
        print( "change action ", i)
        print( self.combobox.currentText())
        self.Mask_annotate.changeAction( i)
        for t in self.Points_annotate:
            t.changeAction( i)
            
    def renorm(self):
        mask = self.Mask_annotate.mask
        vmax = (self.img[mask]).max()
        print( " DA RINORMALIZZARE ", vmax)
        self.im.set_norm(LogNorm(vmax = vmax)) 

        self.draw()

    def write_all(self):
        error = np.array(self.error)
        merr = error.max()
        self.Mask_annotate.update_mask()
        error[True  - self.Mask_annotate.mask ] = merr*1.0e10
        file = h5py.File("fitinit.h5","w")
        file["error"] = error
        file["data"] =  self.img
        curve_group = file.require_group("curves")

        
        for i,t in enumerate(self.Points_annotate):
            rects = t.confirmed_rects
            points=[p for r,p in rects.items()]
            if len(p):
                curve_group[str(i)] = np.array(points)
        rects = self.Mask_annotate.confirmed_rects
        masks = [ r for m,r in   rects.items()     ]
        file["masks"] = np.array(masks)
        file.close()
        
    def read_all(self):
        
        file = h5py.File("fitinit.h5","r")
        masks = file["masks"][:]

        for x0,x1,y0,y1 in  masks:
            rect = matplotlib.patches.Rectangle((x0,y0),(x1-x0), (y1-y0),  facecolor='None', edgecolor='green')
            self.Mask_annotate.ax.add_patch(rect)
            self.Mask_annotate.confirmed_rects[ rect   ] =  [ x0,x1,y0,y1]

        self.Mask_annotate.update_mask()

            
        curve_group = file["curves"]
        for i in range(5):
            points = curve_group[str(i)][:]
            anno = self.Points_annotate[i]
            for x0,y0 in points:
                print( " riaggiungo in " , x0,y0)
                rect = matplotlib.patches.Circle((  x0,y0  ), radius=1)
                anno.confirmed_rects[ rect   ] =  [ x0,y0]
                anno.ax.add_patch(rect )
            anno.fig.canvas.draw()
  



            
    def compute_initial_figure(self):
        file_list = glob.glob("scan_*.txt")
        file_list.sort()

        img=[]
        error = []
        for f in file_list:
            sc = np.loadtxt(  f )
            img.append( sc[:,-2]  )
            error.append( sc[:,-1]  )
            print( len(img[-1]))

        img   = np.array(img   ).T
        error = np.array(error ).T

    
        self.img = img
        self.error = error
        
        self.im = self.axes.imshow(img,  norm=LogNorm(),  #vmax=1000),
                              aspect='auto', origin='lower',
                              interpolation='nearest')
        # print " IMIMIMIM ", im.norm
        # print dir(im)
        self.Mask_annotate = Annotate(self.fig, self.axes,0)
        self.Points_annotate = []
        for i in range(5):
            self.Points_annotate.append( Annotate(self.fig, self.axes,i+1) ) 

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        pass
        # print " ##################################### "
        # print dir(event)
        # print " BUTTON ", event.button, event.xdata, event.ydata

    def on_motion(self, event):
        if event.xdata is not None:
            ix,iy = int(round(event.xdata)), int(round(event.ydata))
            self.label.setText( "X=%10d\nY=%10d\nZ=%10.5e " %(ix,iy ,  self.img[iy,ix])  )

        pass
        # print " ##################################### "
        # print " BUTTON move", event.button, event.xdata, event.ydata
        

class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""

    def compute_initial_figure(self):
        t = arange(0.0, 3.0, 0.01)
        s = sin(2*pi*t)
        self.axes.plot(t, s)

def  readjust(  x0,x1,y0,y1    ) :

    tx0 = math.ceil(min(x0,x1)+0.5 )-0.5
    tx1 = math.floor(max(x0,x1)+0.5)-0.5
    ty0 = math.ceil(min(y0,y1)+0.5)-0.5
    ty1 = math.floor(max(y0,y1)+0.5 )-0.5
    
    return  tx0,tx1,ty0,ty1
    
class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""

    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(1000)

    def compute_initial_figure(self):
        self.axes.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        l = [random.randint(0, 10) for i in range(4)]

        self.axes.plot([0, 1, 2, 3], l, 'r')
        self.draw()

class Annotate(object):
    def __init__(self, fig, axes, myaction):
        self.fig = fig
        self.ax = axes
        self.rect = matplotlib.patches.Rectangle((0,0), 1, 1, facecolor='None', edgecolor='green')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.ax.add_patch(self.rect)
        
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.pressed = False
        self.confirmed_rects = {}

        self.LX0,self.LX1 = self.ax.get_xlim()
        self.LY0,self.LY1 = self.ax.get_ylim()

        self.LX0 = math.ceil(self.LX0 +0.5 )-0.5
        self.LX1 = math.floor(self.LX1+0.5)-0.5
        self.LY0 = math.ceil(self.LY0+0.5)-0.5
        self.LY1 = math.floor(self.LY1  +0.5 )-0.5
            
        self.mask = np.zeros([ self.LY1-self.LY0  , self.LX1-self.LX0  ], dtype = bool)

        self.currentAction = 0
        self.myaction=myaction

    def changeAction(self,i):
        
        self.currentAction = i
        isactive=False
        if self.currentAction==self.myaction:
            isactive=True

        for rect in self.confirmed_rects:
            if not isactive:
                rect.set_facecolor('red')
                rect.set_edgecolor('blue')
            else:
                rect.set_facecolor('blue')
                rect.set_edgecolor('red')
                
        self.fig.canvas.draw()
        
        
    def update_mask(self):
        if self.myaction==0:
            self.mask[:]= True
            for rect, (x0,x1,y0,y1) in self.confirmed_rects.items():
                self.mask[   int(math.ceil(y0)):int(math.ceil(y1)) , int(math.ceil(x0)):int(math.ceil(x1)) ]   = False
        
    def on_press(self, event):
        if self.currentAction!=self.myaction  : return
        # print 'press'
        if event.xdata  is None:
            return
        print( " BUTTON ", event.button)
        if event.button==3:
            
            for rect in list(self.confirmed_rects):
                contains, attrd = rect.contains( event )
                print( contains, attrd)
                if contains:
                    print( " RIMUOVO " , self.confirmed_rects[rect])
                    del self.confirmed_rects[rect]
                    rect.remove()
            self.fig.canvas.draw()
            if self.myaction==0:
                self.update_mask()
            return 
        
        self.x0 = event.xdata
        self.y0 = event.ydata    
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.oldx=None
        self.oldy=None


        if self.myaction==0:
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('dashed')
            self.fig.canvas.draw()
            self.pressed = True
        else:
            if event.button==1:
                print( " VORREI AGGIUNGER UN PUNTO ")
                self.rect = matplotlib.patches.Circle((  self.x0,self.y0  ), radius=1)
                self.confirmed_rects[ self.rect   ] =  [ self.x0,self.y0]
                self.ax.add_patch(self.rect )
                self.fig.canvas.draw()
            elif event.button==2:
                print( "oppure iniziare a spostarlo se bottone di mezzo ")
                self.rect=None
                self.pressed=False
                for rect in list(self.confirmed_rects):
                    contains, attrd = rect.contains( event )
                    if contains:
                        self.rect = rect
                        self.pressed=True
                
        
    def on_motion(self,event):
        if self.currentAction!=self.myaction  :
            return
        
        
        if not self.pressed :
            return
        
        if  event.xdata is None:
            return

        self.oldx=self.x1
        self.oldy=self.y1

        self.x1 = event.xdata
        self.y1 = event.ydata
        
        print( self.x1)
        print( self.y1)
        # pc = self.ax.transData.transform(np.vstack([[self.x1],[self.y1]]).T)
        # print pc.T


        if self.myaction ==0:
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            self.rect.set_linestyle('dashed')
            self.fig.canvas.draw()
        else:
            if event.button==2:
                print( " provo  aspostare ", self.rect, "  su ",  (self.x1, self.y1))
                self.rect.center = (self.x1, self.y1)
                self.confirmed_rects[self.rect] = (self.x1, self.y1)
                self.fig.canvas.draw()
            

        
    def on_release(self, event):
        if self.myaction !=0:
            return
        if self.currentAction!=self.myaction  : return
        self.pressed = False
        print( 'release')
        
        if event.xdata is not None:
            self.x1 = event.xdata
            self.y1 = event.ydata
        elif self.oldx is not None:
            X0,X1 = self.ax.get_xlim()
            Y0,Y1 = self.ax.get_ylim()
            if abs(self.oldx-self.x1 )>  abs(X1-self.x1 ):
                self.x1 = X1
            elif abs(self.oldx-self.x1 )>  abs(X0-self.x1 ):
                self.x1 = X0
            if abs(self.oldy-self.y1 )>  abs(Y1-self.y1 ):
                self.y1 = Y1
            elif abs(self.oldy-self.y1 )>  abs(Y0-self.y1 ):
                self.y1 = Y0
            
        if self.x0 is None:
            return

        self.x0,self.x1,self.y0,self.y1 = readjust(  self.x0,self.x1,self.y0,self.y1    ) 
        print( " READJUSTED ", self.x0,self.x1,self.y0,self.y1)

        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.rect.set_linestyle('solid')

        if self.x0==self.x1 or self.y0==self.y1:
            self.rect.remove()
            ret_val = [ None]*4
        else:
            ret_val = [ self.x0,self.x1,self.y0,self.y1]

            self.confirmed_rects[ self.rect   ] =  [ self.x0,self.x1,self.y0,self.y1]


        self.fig.canvas.draw()

        
        self.rect = matplotlib.patches.Rectangle((0,0), 1, 1, facecolor='None', edgecolor='green')
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.ax.add_patch(self.rect)

        self.update_mask()

        
        return ret_val


class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.file_menu = QtGui.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtGui.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.help_menu.addAction('&About', self.about)

        self.main_widget = QtGui.QWidget(self)

        hlayout =  QtGui.QHBoxLayout(self.main_widget)
        self.panel_widget = QtGui.QWidget(self)
        self.plot_widget  = QtGui.QWidget(self)

        hlayout.addWidget(self.panel_widget)
        hlayout.addWidget(self.plot_widget )

        self.vl_cmds =  QtGui.QVBoxLayout(self.panel_widget)
        self.label = QtGui.QLabel("pippo")
        
        self.vl_cmds.addWidget(self.label)

        
        if 0:
            l = QtGui.QVBoxLayout(self.main_widget)
            sc = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100)
            dc = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=100)
            l.addWidget(sc)
            l.addWidget(dc)

        if 1:
            l = QtGui.QVBoxLayout(self.plot_widget)
            # sc = MyStaticMplCanvas(self.plot_widget, width=5, height=4, dpi=100)
            # dc = MyDynamicMplCanvas(self.plot_widget, width=5, height=4, dpi=100)
            sc = Plot2D(self.plot_widget, width=5, height=4, dpi=100)
            sc.label = self.label
            l.addWidget(sc)
            # l.addWidget(dc)

        but = QtGui.QPushButton('Renorm', self)
        but.clicked.connect(sc.renorm)
        self.vl_cmds.addWidget(but)


        combobox = QtGui.QComboBox( self)
        
        combobox.insertItems( 0,  ["Work On Masks","Work On N 0","Work On N 1","Work On N 2","Work On N 3","Work On N 4"] )
        
        # comboBox.connect(comboBox, QtCore.SIGNAL("activated(int)"), self.textChangedFilter)
        combobox.activated.connect(sc.change_action )

        sc.combobox = combobox
        
        self.vl_cmds.addWidget(combobox)


        write_button = QtGui.QPushButton('Write', self)
        write_button.clicked.connect(sc.write_all)
        self.vl_cmds.addWidget(write_button)
        
        read_button = QtGui.QPushButton('Read', self)
        read_button.clicked.connect(sc.read_all)
        self.vl_cmds.addWidget(read_button)
        
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("All hail matplotlib!", 2000)

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QtGui.QMessageBox.about(self, "About",
                                """embedding_in_qt4.py example
Copyright 2005 Florent Rougon, 2006 Darren Dale

This program is a simple example of a Qt4 application embedding matplotlib
canvases.

It may be used and modified with no restriction; raw copies as well as
modified versions may be distributed without limitation."""
                                )


qApp = QtGui.QApplication(sys.argv)

aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()







