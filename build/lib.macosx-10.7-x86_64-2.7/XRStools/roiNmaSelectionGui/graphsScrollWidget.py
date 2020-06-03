from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from  PyMca5.PyMcaGui  import MaskImageWidget  as sole_MaskImageWidget
# from  PyQt4 import Qt, QtCore, QtGui
from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui
import os

import numpy as np
import string

from PyMca5.PyMcaGraph.Plot import Plot
from six.moves import range
from silx.gui.plot.PlotWindow import Plot1D  , Plot2D
from silx.gui.plot import PlotWidget, actions, items
from silx.gui.plot.MaskToolsWidget import MaskToolsDockWidget
from silx.gui.plot.Profile import ProfileToolBar
from silx.gui.plot.PlotTools import PositionInfo
from silx.gui import qt
from sklearn.decomposition import NMF

from . import nnma
import time



from XRStools.installation_dir import installation_dir

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]



#

# from PyMca5.PyMcaGraph.backends.OpenGLBackend import OpenGLBackend
# Plot.defaultBackend = OpenGLBackend 
# Plot.defaultBackend = OpenGLBackend 


#
# DA FARE
#
#   Lanciare con selezione default componente 0   OK 
#
#   Inizializzare la maschera a zero         OK
#   renderla persistente                     OK
#
#   Usare i radio button per verificare nella gestione degli eventi OK
#
#   Vedere Che le componenti siano solo sulla roi iniziale   OK
#
#   Mettere a zero i dati iniziali con mask ioniziale      OK
#     
#   Bloccare la selezione dove la maschera iniziale e' 1   OK
#
#   Automatizzare tutto con scelte prestabilite :  usare prima comp, trheshold a una frazione  OK
#
#   Slider Trasparenza   OK 
#   Slider dimensione  OK 
#
#   Automatizzazione globale con scelte prestabilite ( bisognera conservare creati in liste) OK 
#
#  Salvare ROI
#  Salvare separatamente stato  con selezione N comp, numero di comp maschera threshold
#  Rilettura


class MyPlot1D(Plot1D):
    def __init__(self, parent=None):
        super(MyPlot1D, self).__init__(parent)  # , backend = "gl")
        self.shint = 400
        
        self.setSizePolicy( Qt.QSizePolicy.Fixed, Qt.QSizePolicy.Fixed  ) # ou maximum        
    def sizeHint(self ) :
        return Qt.QSize( self.shint, self.shint)
    
    def setSizeHint(self, val ) :
        self.shint = val
        self.updateGeometry()
        
    


class MyPlot(PlotWidget):

    def __init__(self, parent=None):
        super(MyPlot, self).__init__(parent)  # , backend = "gl")
        
        self.resetZoomAction = actions.control.ResetZoomAction(self)
        self.addAction(self.resetZoomAction)
        
        self.colormapAction = actions.control.ColormapAction(self)
        self.addAction(self.colormapAction)

        toolbar = qt.QToolBar('Plot', self)
        toolbar.addAction(self.resetZoomAction)
        toolbar.addAction(self.colormapAction)
        self.addToolBar(toolbar)

        self.profile = ProfileToolBar(plot=self)
        self.addToolBar(self.profile)

        #self.maskToolsDockWidget = MaskToolsDockWidget(
        #    plot=self, name='Mask')
        #self.maskToolsDockWidget.hide()
        #self.addDockWidget(qt.Qt.BottomDockWidgetArea, self.maskToolsDockWidget)
        
        posInfo = [
            ('X', lambda x, y: x),
            ('Y', lambda x, y: y),
            ('Data', self._getImageValue)]
        self.positionWidget = PositionInfo(plot=self, converters=posInfo)
        self.statusBar().addWidget(self.positionWidget)

        self.getDefaultColormap().setName('temperature')  # 'temperature'....viridis

        self.transp = 100
        self.shint = 400
        
        self.setSizePolicy( Qt.QSizePolicy.Fixed, Qt.QSizePolicy.Fixed  ) # ou maximum

        
    #def getSelectionMask(self):
    #    return self.maskToolsDockWidget.getSelectionMask()
    def sizeHint(self ) :
        return Qt.QSize( self.shint, self.shint)
    
    def setSizeHint(self, val ) :
        self.shint = val
        self.updateGeometry()
    

    def reset_SelectionMask(self, transp=100):
        self.transp=transp
        self.setSelectionMask(  self.getSelectionMask() )
        
    def setSelectionMask(self, mask):
        #return bool(self.maskToolsDockWidget.setSelectionMask(mask))
        image = np.zeros(mask.shape + (4,), dtype=np.uint8)
        image[:, :, :3] = (255, 0, 0)   # rgb color
        image[:, :, 3][mask == 1] = self.transp  # transparency
        self.addImage(image, legend='mask', z=2, replace=False, info=mask)
        
    def getSelectionMask(self,):
        #return bool(self.maskToolsDockWidget.setSelectionMask(mask))
        return self.getImage(legend='mask').getInfo()

    def _getImageValue(self, x, y):
        """Get status bar value of top most image at position (x, y)
        :param float x: X position in plot coordinates
        :param float y: Y position in plot coordinates
        :return: The value at that point or '-'
        """
        value = '-'
        valueZ = -float('inf')
        mask = 0
        maskZ = -float('inf')

        mask=" "
        for image in self.getAllImages():
            data = image.getData(copy=False)
            isMask = image.getLegend() == 'mask'
            
            if isMask:
                data = image.getInfo()
            #    zIndex = maskZ
            #else:
            zIndex = valueZ
            if image.getZValue() >= zIndex:
                # This image is over the previous one
                ox, oy = image.getOrigin()
                sx, sy = image.getScale()
                row, col = (y - oy) / sy, (x - ox) / sx
                if row >= 0 and col >= 0:
                    # Test positive before cast otherwise issue with int(-0.5) = 0
                    row, col = int(row), int(col)
                    if (row < data.shape[0] and col < data.shape[1]):
                        v, z = data[row, col], image.getZValue()
                        if not isMask:
                            value = v
                            valueZ = z
                        else:
                            mask = v
                            maskZ = z
        if maskZ > valueZ and mask > 0:
            return value, "Masked"
        return value


    
# @ui.UILoadable
class pickup_choice(QtGui.QWidget):
    def __init__(self, parent):
        
            
        super( pickup_choice, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "pickup_choice.ui" ), self)


        
        self.choosebypointvalue.setChecked(True)
        self.visualizefromroi.setChecked(True)
# @ui.UILoadable
class NNMF_control(QtGui.QWidget):
    def __init__(self, parent,data_tobeanalyzed):
      
        super(NNMF_control , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "NNMF_control.ui" ), self)

  
        self.data_tobeanalyzed = data_tobeanalyzed
        self.pushButton_byNofComps.clicked.connect(self.onNNMF_button_clik)
        self.pushButton_percentage.clicked.connect(self.onNNMF_prec_button_clik)
        self. lineEdit_percentage.hide()
        key , data, image , (corner, roi_mask) = data_tobeanalyzed
        self.data=data
        self.image=image
        self.roi_mask = roi_mask

    def threshold(self):# , val=None):
        print("THRESHOLD")
        # if val is None:
        val =  float(str( self.pickup_choice_obj.threshold_value.text()))
        self.threshold_byval(val)


    def threshold_byval(self, val):
        mask = np.less(val*self.scp.max() ,self.scp )
        self.plot_comp_spatial.setSelectionMask(  mask   )
        self.pickup_choice_obj.threshold_value.setText( str(val) )
        
    def reset(self, val):
        mask = np.less(self.scp.max()+1 ,self.scp )
        self.plot_comp_spatial.setSelectionMask(  mask   )

        
    def onNNMF_button_clik(self):
        ncomps = int(str(self. lineEdit_byNofComps.text()   ))
        self.NNMF_by_comps(ncomps)

        
    def onNNMF_prec_button_clik(self):
        ncomps = int(str(self. lineEdit_byNofComps.text()   ))
        self.NNMF_by_percentage(ncomps)
    ##  I segnali corrono sulla dimensione lenta
    def NNMF_by_comps(self, ncomps):

        self.cslider.setMinimum(0)
        self.cslider.setMaximum(ncomps-1)
        self.cslider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.cslider.setTickInterval(1)
        
        self.lineEdit_byNofComps.setText(str(ncomps) )
        key , data, image , (corner, roi_mask) = self.data_tobeanalyzed

        data = data * roi_mask.flatten()
        
        param={}
        print (" NNMF per ", data.shape, str(data))
        start=time.time()
        #X,Y,obj,count,converged = nnma.FastHALS(data, ncomps, eps=1e-5, verbose=0,
        #                                        maxcount=1000, **param)
        # X,Y,obj,count,converged = nnma.ALS(np.array(data), ncomps, eps=1e-5, verbose=0,
        #                                         maxcount=1000, **param)
        # print (" OK " )
        # assert(converged)
        # print("obj = %E  count=%5d  converged=%d  TIME=%.2f secs" % \
        #              (obj,count, converged, time.time()-start))

        model = NMF(n_components=ncomps, init='random', random_state=0)
        X = model.fit_transform((data*self.roi_mask.flatten()).T)
        Y = model.components_
        
        contributed_area  =  X.sum(axis=0) * Y.sum(axis=-1)
        order  = np.argsort(   contributed_area      )[::-1]
        X=X[:,order]
        Y=Y[order,:]
        contributed_area = contributed_area[order]
        print ("contributed_area ", contributed_area)
        self.weightplot.addCurve(x=np.arange( contributed_area.shape[0])+1, y=contributed_area, legend="Component weigth", replace=True)
        # self.weightplot.addCurve(x=np.arange( contributed_area.shape[0])+1, y=contributed_area[::-1], legend="Component weigth", replace=False, axis="left")

        self.X = X
        self.Y = Y
        for icomp in range(  Y.shape[0]  ):
            c = Y[icomp,:]
            self.compsplot.addCurve(x=np.arange( c.shape[0])+1, y=c, legend="Component #%d"%icomp, replace=(icomp==0))

        
        # return
        if(0):
            distances=[]
            for i in range(1,ncomps+1):
                # Xp,Yp,obj,count,converged = nnma.FastHALS(data, i, eps=1e-5, verbose=0,
                #                                          maxcount=1000, **param)
                # Xp,Yp,obj,count,converged = nnma.ALS(np.array(data), i, eps=1e-5, verbose=0,
                #                                           maxcount=1000, **param)
                model = NMF(n_components=i, init='random', random_state=0)
                Xp = model.fit_transform(data.T)
                Yp = model.components_
                dist = np.abs( data.T - np.dot(  Xp,Yp   )       ).sum()
                distances.append(dist)
            self.plot_comp_spatial.addCurve(x=np.arange( contributed_area.shape[0])+1, y=distances, legend="distance", replace=False, yaxis="right")
            
        for where in ["left", "right"]:
            ax=self.weightplot.getYAxis(axis=where)
            vmin,vmax = ax.getLimits()
            ax.setLimits(0,vmax)

        self.cslider.setValue(0)

    def onComponentSelected(self, previous, legend):
        if legend==previous:
            return
        key , data, image , (corner, roi_mask) = self.data_tobeanalyzed
        print(" SELEZIONATO ", previous, legend)
        icomp = int(''.join(filter(str.isdigit, str(legend) )))

        self.cslider.setValue(icomp)


        
    def NNMF_by_percentage(self, ncomps):
        key , data, image , (corner, roi_mask) = self.data_tobeanalyzed
        data = data * roi_mask.flatten()

        model = NMF(n_components=ncomps, init='random', random_state=0)
        
        distances=[]
        for i in range(1,ncomps+1):
            model = NMF(n_components=i, init='random', random_state=0)
            Xp = model.fit_transform(data.T)
            Yp = model.components_
            
            dist = np.abs( data.T - np.dot(  Xp,Yp   )       ).sum()
            distances.append(dist)
        distances=np.array(distances)
        self.weightplot.addCurve(x=np.arange( distances.shape[0])+1, y=distances, legend="distance", replace=False, axis="right")
        for where in ["left", "right"]:
            ax=self.weightplot.getYAxis(axis=where)
            vmin,vmax = ax.getLimits()
            ax.setLimits(0,vmax)

    def cslider_valuechange(self, icomp):
        key , data, image , (corner, roi_mask) = self.data_tobeanalyzed

        scp = self.X[:,icomp]
        scp.shape = roi_mask.shape

        oldmask = self.plot_comp_spatial .getImage(legend='mask')
        if oldmask is not None:
            oldmask = oldmask.getInfo()

        
        self.plot_comp_spatial.addImage( scp   )
        self.scp = scp
        
        self.plot_comp_label.setText( key + ": Spatial component #%d"%icomp  )

        self.compsplot.setActiveCurve("Component #%d"%icomp)
  
        if oldmask is not None:
            self.plot_comp_spatial.setSelectionMask(oldmask)
            
    def plotSignalHandler(self, signal):
        key , data, image , (corner, roi_mask) = self.data_tobeanalyzed

        if self.pickup_choice_obj.visualizefromroi.isChecked():
            data = data*self.roi_mask.flatten()
        
        if signal["event"]=="mouseMoved":
            ix=int(signal["x"])
            iy=int(signal["y"])

            if(ix<0 or ix>= roi_mask.shape[1] ):
                return
            if(iy<0 or iy>= roi_mask.shape[0] ):
                return
            
            if self.pickup_choice_obj.choose_byline.isChecked():
                ipos = iy*self.roi_mask.shape[1]+np.arange(  self.roi_mask.shape[1])
                spectra = (self.data[:,ipos]).sum(axis=-1)
            elif self.pickup_choice_obj.choose_bycolumn.isChecked():
                ipos =  np.arange(  self.roi_mask.shape[0])*self.roi_mask.shape[1]+ix
                spectra = (self.data[:,ipos]).sum(axis=-1)
            else:
                ipos = iy*self.roi_mask.shape[1]+ix
                spectra = self.data[:,ipos]
            
            self.plot_spectra.addCurve(  x = np.arange( spectra.shape[0])  , y= spectra , legend = "pixel spectra", replace=True)
        elif signal["event"]=="mouseClicked" and signal["button"]=="left":
            
            if self.pickup_choice_obj.choosebypointvalue.isChecked():
                return
            
            ix=int(signal["x"])
            iy=int(signal["y"])
            
            if(ix<0 or ix>= roi_mask.shape[1] ):
                return
            if(iy<0 or iy>= roi_mask.shape[0] ):
                return
            
            if self.pickup_choice_obj.choose_bypoint.isChecked():
                mask = self.plot_comp_spatial.getSelectionMask()
                mask[iy,ix] = (1-mask[iy,ix])*self.roi_mask[iy,ix]
                self.plot_comp_spatial.setSelectionMask(mask)
                
            elif self.pickup_choice_obj.choose_byline.isChecked():
                mask = self.plot_comp_spatial.getSelectionMask()
                if mask[iy,:].sum():
                    mask[iy,:] = 0
                else:
                    mask[iy,:] = self.roi_mask[iy,:]
                self.plot_comp_spatial.setSelectionMask(mask)
                
            elif self.pickup_choice_obj.choose_bycolumn.isChecked():
                mask = self.plot_comp_spatial.getSelectionMask()
                if mask[:,ix].sum():
                    mask[:,ix] = 0
                else:
                    mask[:,ix] = self.roi_mask[:,ix]
                self.plot_comp_spatial.setSelectionMask(mask)
  
    
class graphsScrollWidget(QtGui.QScrollArea ):

    def set_size(self, val):
        print ( " in set size ", val)

        # rcount = self.scrolledGrid.rowCount()
        # ccount = self.scrolledGrid.colCount()

        # for i in range(rcount):
        #     self.scrolledGrid.        
        
        count = self.scrolledGrid.count()
        for i in range(count):
            item = self.scrolledGrid.itemAt(i)
            if not ( isinstance(item.widget(), MyPlot)  or     isinstance(item.widget(), MyPlot1D)   ) :
                continue
            item.widget().setSizeHint( val)
            # item.setGeometry( Qt.QRect(0,0,val,val))

    def setComps(self):
        
        nofcomps = int(str(self.gc.nofcomps.text()))
        ccomp    = int(str(self.gc.choosedcomp.text()))
        
        for nnmf in self.nnmf_list       :
            nnmf.NNMF_by_comps( nofcomps ) 
            nnmf.cslider_valuechange( ccomp ) 


    def setThreshold(self):
        threshold = float(str(self.gc.threshold.text()))
        for nnmf in self.nnmf_list       :
            nnmf.threshold_byval( threshold ) 
            
    def set_transp(self, val):
        for p in self.rois_w_list+self.refined_rois_w_list:
            p.reset_SelectionMask(val)
        
        count = self.scrolledGrid.count()
        for i in range(count):
            item = self.scrolledGrid.itemAt(i)
            if not isinstance(item.widget(), MyPlot):
                continue
            print(i, item.widget())
            print(i, item.widget().height() )

    def setSelectionMaskS(self, masksDict, metadata):
        for w, tba, nnmf in zip(self.refined_rois_w_list, self.tobeanalyzed_items  , self.nnmf_list ):
            key , data, image , (corner, roi_mask) = tba

            if key in masksDict:
                print(" PER \n"*10)
                print(" KEY ", key)
                print( "mask ",masksDict[key])
                corners, mask = masksDict[key]
                w.setSelectionMask(mask)
                
            if key in metadata:
                mdata = metadata[key]
                if "globalThreshold" in mdata :
                    self.gc.threshold.setText(str( mdata[ "globalThreshold"]  ) )
                if "globalNofComps" in mdata :
                    self.gc.nofcomps.setText(str( mdata["globalNofComps"] ))
                if "globalChoosedComp" in mdata :
                    self.gc.choosedcomp.setText(str( mdata["globalChoosedComp"]) )

                if "threshold" in mdata:
                    nnmf.pickup_choice_obj.threshold_value.setText(  str( mdata["threshold"])  )
                if "nofcomps" in mdata:
                    nnmf.lineEdit_byNofComps.setText(  str( mdata["nofcomps"])  )

                    if 1:
                        nofcomps = int( str( mdata["nofcomps"]))
                        if nofcomps>0 and nofcomps<100:
                            nnmf.NNMF_by_comps(nofcomps)

                        if "icomp" in mdata:
                            icomp = int( str( mdata["icomp"]) ) 
                            nnmf.cslider_valuechange(icomp)
                    # except:
                    #     pass


                    # @@@@@ QUA SI PERDE INFORMAZIONE GEOMETRIA INIZIALZE
    def getSelectionMaskS(self, mask=None):
        metadata = {}
        metageo  = {}
        for w, tba, nnmf in zip(self.refined_rois_w_list, self.tobeanalyzed_items  , self.nnmf_list ):
            key , data, image , (corner, roi_mask) = tba
            if mask is None:
                mask = np.zeros(image.shape )
            wmask = w.getSelectionMask()
            meta  = {}
            inum =    int(''.join(filter(str.isdigit, str(key) )))+1
            meta["globalThreshold"]=     str(self.gc.threshold.text())
            meta["globalNofComps"]=     str(self.gc.nofcomps.text())
            meta["globalChoosedComp"]=     str(self.gc.choosedcomp.text())
            meta["threshold"]           =  str(nnmf.pickup_choice_obj.threshold_value.text())
            meta["nofcomps"]           = str(nnmf.lineEdit_byNofComps.text() )
            meta["icomp"]              = nnmf.cslider.value()
            metadata[key]= meta
            mask[ corner[0]:corner[0]+roi_mask.shape[0]               , corner[1]:corner[1]+roi_mask.shape[1]   ] = inum*wmask
            metageo[inum] = corner, roi_mask.shape
        return mask, metadata, metageo
            
    def __init__(self, parent, tobeanalyzed_items, initial_ncomps=2):        
        QtGui.QScrollArea.__init__(self, parent )
        scrolledWidget = QtGui.QWidget()
        scrolledGrid   = QtGui.QGridLayout()

        self.rois_w_list = []
        self.refined_rois_w_list = []
        self.nnmf_list = []
        self.tobeanalyzed_items  = tobeanalyzed_items
        for irow, (key , data, image , (corner, roi_mask) ) in enumerate(tobeanalyzed_items):
            icol=0
            
            scrolledGrid.addWidget( QtGui.QLabel( key + ": Commands for components" ), irow*2, icol)
            tobeanalyzed = (key , data, image , (corner, roi_mask) )
            
            nnmf_ctrl = NNMF_control(None, tobeanalyzed)
            self.nnmf_list.append(nnmf_ctrl)

            scrolledGrid.addWidget(nnmf_ctrl , irow*2+1, icol)
            
            icol+=1
            
            scrolledGrid.addWidget( QtGui.QLabel( key + ": Region with original selection " ), irow*2, icol);
            roiplot = MyPlot()
            roiplot.getYAxis().setInverted(True)
            self.rois_w_list.append(roiplot)            
            
            roiplot.addImage(  image[ corner[0]:corner[0]+ roi_mask.shape[0],   corner[1]:corner[1]+ roi_mask.shape[1]   ]   )
            roiplot.setSelectionMask(  roi_mask   )
            scrolledGrid.addWidget(roiplot ,irow*2+1, icol)


            icol+=1
            
            scrolledGrid.addWidget( QtGui.QLabel( key + ": Tot. w. per c. & res. dist." ), irow*2, icol);
            
            weightplot=  MyPlot1D()#control=True)
            weightplot.getFitAction().setVisible(True)
            scrolledGrid.addWidget(weightplot ,  irow*2+1, icol);
            nnmf_ctrl.weightplot = weightplot


            icol+=1
            
            scrolledGrid.addWidget( QtGui.QLabel( key + ": Components" ), irow*2, icol);
            
            compsplot=  MyPlot1D()#control=True)
            compsplot.getFitAction().setVisible(True)
            scrolledGrid.addWidget(compsplot ,  irow*2+1, icol);
            nnmf_ctrl.compsplot = compsplot
            compsplot.sigActiveCurveChanged.connect(   nnmf_ctrl.onComponentSelected    )




            icol +=1

            cslider  = QtGui.QSlider(Qt.Qt.Vertical)
            scrolledGrid.addWidget(cslider ,  irow*2+1, icol);
            
            cslider.valueChanged.connect(nnmf_ctrl.cslider_valuechange)
            
            icol +=1

            plot_comp_spatial =  MyPlot()
            self.refined_rois_w_list.append(plot_comp_spatial) 
            
            plot_comp_spatial.getYAxis(axis="left").setInverted(True)

            clab =  QtGui.QLabel( key + ": Spatial component")
            scrolledGrid.addWidget( clab, irow*2, icol )
            scrolledGrid.addWidget(plot_comp_spatial ,  irow*2+1, icol);

            icol +=1


            plot_spectra =  MyPlot()
            scrolledGrid.addWidget( plot_spectra,  irow*2+1, icol);

            icol+=1
            
            pickup_choice_obj = pickup_choice(None)
            scrolledGrid.addWidget( pickup_choice_obj,  irow*2+1, icol);
            
            plot_comp_spatial.sigPlotSignal.connect( nnmf_ctrl.plotSignalHandler    )

            
            
            nnmf_ctrl.plot_comp_spatial = plot_comp_spatial
            nnmf_ctrl.plot_comp_label = clab
            nnmf_ctrl.cslider = cslider
            nnmf_ctrl.plot_spectra = plot_spectra
            nnmf_ctrl.pickup_choice_obj = pickup_choice_obj
            
            pickup_choice_obj.threshold.clicked.connect( nnmf_ctrl.threshold )
            pickup_choice_obj.pushButton_reset.clicked.connect( nnmf_ctrl.reset )
            
            self.nnmf_ctrl = nnmf_ctrl
            nnmf_ctrl.NNMF_by_comps(initial_ncomps)
            nnmf_ctrl.cslider_valuechange(0)
            nnmf_ctrl.threshold_byval(  0.2 )

        if(0):
            scrolledGrid.addWidget( QtGui.QLabel("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), 0, 0);
            for i in range(80):
                scrolledGrid.addWidget( QtGui.QLabel(" B "), 1, i);
            scrolledGrid.addWidget( QtGui.QLabel("C"), 2, 0);
            scrolledGrid.addWidget( QtGui.QLabel("D"), 3, 0);
            scrolledGrid.addWidget( QtGui.QLabel("E"), 4, 0); 

            myplot=  MyPlot1D()#control=True)
            myplot.getFitAction().setVisible(True)
            x=np.arange(0,10,0.1)
            y=np.sin(x)
            myplot.addCurve(x=x, y=y, legend="pippo", replace=True)
            scrolledGrid.addWidget(myplot , 5, 0);


        
        self.scrolledGrid = scrolledGrid
        scrolledWidget.setLayout(scrolledGrid)
        self.setWidget(scrolledWidget)




