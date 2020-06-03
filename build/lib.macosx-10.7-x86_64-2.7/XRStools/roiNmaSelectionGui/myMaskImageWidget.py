from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from  PyMca5.PyMcaGui  import MaskImageWidget  as sole_MaskImageWidget
# from  PyQt4 import Qt, QtCore
from silx.gui import qt as Qt 
from silx.gui import qt as QtCore
import numpy 
import string

from PyMca5.PyMcaGraph.Plot import Plot
from six.moves import range

# from PyMca5.PyMcaGraph.backends.OpenGLBackend import OpenGLBackend
# Plot.defaultBackend = OpenGLBackend 
# Plot.defaultBackend = OpenGLBackend 

class MaskImageWidget(sole_MaskImageWidget.MaskImageWidget):
    changeTagOn=False
    def _graphSignal(self, ddict, ownsignal = None):
 
        if self.changeTagOn:
            if self.__selectionMask is not None and ddict['event']=="mouseClicked" and ddict['button']=="middle" :
                x,y = int(ddict["x"]), int(ddict["y"])

                y, x = sole_MaskImageWidget.convertToRowAndColumn(x , y,
                                             self.__imageData.shape,
                                             xScale=self._xScale,
                                             yScale=self._yScale,
                                             safe=True)

                id_target = self.__selectionMask[y,x]
                if id_target and id_target!= self._nRoi:
                    mask_target =  (self.__selectionMask == id_target  )
                    mask_swap   =  (self.__selectionMask== self._nRoi)
                    self.__selectionMask [mask_target] = self._nRoi
                    self.__selectionMask [mask_swap] = id_target
                    emitsignal = True
                    if emitsignal:
                        self.plotImage(update = False)
                        self._emitMaskChangedSignal()
                    return
        # super(MaskImageWidget, self)._graphSignal(ddict, ownsignal) 
        sole_MaskImageWidget.MaskImageWidget._graphSignal(self, ddict, ownsignal)


    def dragEnterEvent(self,event):
        print( dir(event.mimeData()))
        print( list(event.mimeData().formats()))
        print( event.mimeData().text())

        model = Qt.QStandardItemModel()
        model.dropMimeData(event.mimeData(), QtCore.Qt.CopyAction, 0,0, Qt.QModelIndex())
        print( model)

        if event.mimeData().hasFormat('application/x-qabstractitemmodeldatalist'):
            event.acceptProposedAction()

            bytearray = event.mimeData().data('application/x-qabstractitemmodeldatalist')
            data_items = decode_data(bytearray)
            print( data_items)

        if event.mimeData().hasFormat("text/plain"):
            print( " OK ")
            event.acceptProposedAction()
 
    def dropEvent(self, e):

   

        localpos = self.graph.getWidgetHandle().mapFromGlobal( Qt.QCursor().pos()   )
        x,y = localpos.x(), localpos.y()

        x,y = self.graph.pixelToData(x,y ) 
        print( "POSITION ",x,y)
    
        mask = self.getSelectionMask()
        
        Ct = mask[int(y), int(x)]
        print( " VALORE MASCHERA ", Ct)

        if Ct:
            print( str(e.mimeData().text()))
            Cc = int( str(e.mimeData().text()))
            # bytearray = e.mimeData().data('application/x-qabstractitemmodeldatalist')
            # data = decode_data(bytearray)
            # print data
            # Cc = (2-data[0])*4+data[1] + 1
            # print Cc
            zonet = (mask==Ct)
            zonec = (mask==Cc)
            mask[zonet]=Cc
            mask[zonec]=Ct
            self.setSelectionMask(mask)
            self.annotateSpots()


    def annotateSpots(self, a_ids = None, offset=None):

        self.graph.clearMarkers()
        mask = self.getSelectionMask().astype("i")

        nspots = mask.max()
        for i in range(1,nspots+1):
            m = (mask==i).astype("f")
            msum=m.sum()
            if msum:
                ny,nx = m.shape
                px= (m.sum(axis=0)*numpy.arange(nx)).sum()/msum
                py= (m.sum(axis=1)*numpy.arange(ny)).sum()/msum

                print( " ##################################    ", px, py)
                extra_info = ""
                if a_ids is not None:
                    if offset is None:
                        extra_info = "(A%d)"% a_ids[i-1]
                    else:
                        extra_info = "(A[%d],ROI%02d )"% (a_ids[i-1], offset+i-1)

                # res=self.graph.insertMarker( px, py, "legend"+str(i)+extra_info, "%d"%(i)+extra_info, color='black', selectable=False, draggable=False, searchFeature=True, 
                #                            xytext = (-20, 0),
                #                            textcoords = 'offset points',
                #                            ha = 'right', va = 'bottom',
                #                            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.4),
                #                            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
                
                res=self.graph.addMarker( px, py, "legend"+str(i)+extra_info, "%d"%(i)+extra_info, color='black')




def decode_data( bytearray):
        
    data = []
    item = {}
    
    ds = QtCore.QDataStream(bytearray)
    while not ds.atEnd():
            
        row = ds.readInt32()
        column = ds.readInt32()

        return row, column

        print( row, column)

        map_items = ds.readInt32()
        for i in range(map_items):
               
            key = ds.readInt32()
            
            value = QtCore.QVariant()
            ds >> value
            item[QtCore.Qt.ItemDataRole(key)] = value
               
        data.append(item)
           
    return data







