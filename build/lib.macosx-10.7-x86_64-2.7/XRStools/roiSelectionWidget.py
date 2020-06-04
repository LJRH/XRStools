from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from  PyQt4 import Qt, QtCore
from silx.gui import qt as Qt
from silx.gui import qt as QtCore
# from . import ui
import os
import PyMca5.PyMcaIO.EdfFile as edf
import PyMca5

import PyMca5.PyMcaIO.specfilewrapper as specfile
from . import localfilesdialog
from . import myMaskImageWidget

import numpy
import string
import six
from six.moves import range
from six.moves import zip

from . import xrs_rois
from . import xrs_utilities

from . import spotdetection
from . import match
import pickle
import h5py
from  XRStools.installation_dir import installation_dir

SKIP_WARNING = False

FASTDEBUG=False
try:
    import PyTango
except:
    print( " Could not load PyTango")


def h5_assign_force(h5group, name, item):
    if name in h5group:
        del h5group[name]
    h5group[name] = item


# @ui.UILoadable
class spotdetectioncontrol(Qt.QDockWidget):
    def __init__(self, parent,flag,  detectionCallBack=None,  fatteningCallBack=None,
                 thresholdCallBack=None,
                 annotateMaskCallBack=None,    relabeliseMaskCallBack=None,
                 resetMask           = None):


              
        super(spotdetectioncontrol , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" ,  "spotdetectioncontrol.ui" ), self)
        
        # Qt.QDialog.__init__(self, parent, flag)
        # self.loadUi()  # load the ui file

        self.detectionButton .setToolTip("Run a Canny edge detection. If threshold entry contains a float >0 and <1,\n then all pixel <threshold*maxval are set to zero beforehand." )
        self.HoughDetection  .setToolTip("Run a specialised  Canny edge detection for LINES. If threshold entry contains a float >0 and <1 ,\nthe  all pixel <threshold*maxval are set to zero beforehand." )
        self.thresholdButton .setToolTip("Use this for existing ROIs : ROIs is regenerated starting from its maximum  and expanding till threshold*maxvalue"  )
        
        # self.detectionButton.clicked.connect(detectionCallBack)
        self.detectionCallBack = detectionCallBack
        self.detectionButton.clicked.connect(self.my_detectionCallBack)


        
        self.HoughDetection.clicked.connect(self.my_HoughCallBack)


        
        # self.fatteningSpinBox.valueChanged[int].connect(fatteningCallBack)

        self.fatteningCallBack = fatteningCallBack
        self.inflateButton.clicked.connect(self.callcallback)


        self.thresholdCallBack = thresholdCallBack
        self.thresholdButton.clicked.connect(self.callThrecallback)


        self.annotateButton.clicked.connect(annotateMaskCallBack)
        self.relabeliseButton.clicked.connect(relabeliseMaskCallBack)
        self.resetButton.clicked.connect(resetMask)

        self.geo_informations = None
        
    def my_detectionCallBack(self):
        thr_s =  str(self.thresholdEdit.text())
        self.detectionCallBack(thr_s = thr_s)
        
    def my_HoughCallBack(self):
        thr_s =  str(self.thresholdEdit.text())
        self.detectionCallBack(thr_s = thr_s, Hough=True)
        

    def callThrecallback(self):
        print( self.thresholdEdit.text())
        print( str(self.thresholdEdit.text()))
        # value= string.atof(str(self.thresholdEdit.text()))
        value= float(str(self.thresholdEdit.text()))
        self.thresholdCallBack(value)

    def callcallback(self):
        value= self.fatteningSpinBox.value()
        self.fatteningCallBack(value)
   
# @ui.UILoadable
class imageview(Qt.QWidget):
    all_layouts=[]
    def __init__(self, parent=None, isglobal=False, layoutsNames =None):

              
        super(imageview, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" ,  "imageview.ui" ), self)

        
        # Qt.QWidget.__init__(self, parent)
        # self.loadUi()  # load the ui file
        
        print(  layoutsNames  ) 
        if not isglobal:
            self.registeringLayoutComboBox.addItems( layoutsNames     )
            self.all_layouts.append(self.registeringLayoutComboBox)
        else:
            self.registeringLayoutComboBox.addItems( layoutsNames )
            self.registeringLayoutComboBox.currentIndexChanged.connect(self.changeAll)
    
    def changeAll(self,index):
        for t in self.all_layouts:
            t.setCurrentIndex(index)

            
   
# @ui.UILoadable
class spotregistrationcontrol(Qt.QDockWidget):
    def __init__(self, parent,flag,  globregistrationCallBack=None,
                 registrationCallBack=None):


              
        super( spotregistrationcontrol, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" ,  "spotregistrationcontrol.ui" ), self)
   
        # Qt.QDialog.__init__(self, parent, flag)
        # self.loadUi()  # load the ui file
        
        self.globregistrationButton.clicked.connect(globregistrationCallBack)
        self.registrationButton.clicked.connect(registrationCallBack)

        self.tableWidget.setRowCount(3)
        self.tableWidget.setColumnCount(4)
        count=1
        for n in range(1,12+1):
            self.tableWidget.setItem( 2-(n-1)//4,(n-1)%4,  Qt.QTableWidgetItem("%d"%n))

        self.tableWidget.resizeColumnsToContents()

        def _mousePressEvent( event):
            """ in case of problems, generate the code for spotregistrationcontrol from ui
                and reimplement it as a class method """

            print( " evento at ", event.pos())
            print( "  widget ",   self.tableWidget.itemAt( event.pos()))
            print(    self.tableWidget.itemAt( event.pos()).text())


            drag      =  Qt.QDrag(   self.tableWidget)
            mimeData  =  Qt.QMimeData()

            mimeData.setText( self.tableWidget.itemAt( event.pos()).text()   );
            drag.setMimeData(mimeData)
            ## iconPixmap = Qt.QPixmap(_cross_data)
            ## drag.setPixmap(iconPixmap)

            dropAction = drag.exec_();


        self.tableWidget.mousePressEvent = _mousePressEvent
        

class MyTableModel(Qt.QAbstractTableModel):
    def __init__(self, parent, mylist, header, *args):
        Qt.QAbstractTableModel.__init__(self, parent, *args)
        self.mylist = mylist
        self.header = header
    def rowCount(self, parent):
        return len(self.mylist)
    def columnCount(self, parent):
        return len(self.mylist[0])
    def data(self, index, role):
        if not index.isValid():
            return None
        # elif role != Qt.DisplayRole:
        #     return None
        return self.mylist[index.row()][index.column()]
    def headerData(self, col, orientation, role):
        # if orientation == Qt.Horizontal and role == Qt.DisplayRole:
        #     return self.header[col]
        return None
    # def sort(self, col, order):
    #     """sort table by given column number col"""
    #     self.emit(SIGNAL("layoutAboutToBeChanged()"))
    #     self.mylist = sorted(self.mylist,
    #         key=operator.itemgetter(col))
    #     if order == Qt.DescendingOrder:
    #         self.mylist.reverse()
    #     self.emit(SIGNAL("layoutChanged()"))


    
# @ui.UILoadable
class mainwindow(Qt.QMainWindow):
    user_input_signal = QtCore.pyqtSignal(object)
    def __init__(self, parent=None, labelformat="ROI%02d",layout="2X3-12"):
      
        super(mainwindow, self).__init__(parent)

        Qt.loadUi(  os.path.join(  installation_dir,"resources" ,  "mainwindow.ui" ), self)

        self.geo_informations=None
        self.isOK=0
        self.labelformat = labelformat

        self.layout=layout
        self.image=None
        


        # self.actionLoad.triggered.connect(self.LoadRemote )
        self.actionSelectScanFiles.triggered.connect(self.LoadLocal)
        self.actionSpot_detection.triggered.connect(self.CreateSpotDetectionDockWidget)
        self.actionGlobalSpotDetection.triggered.connect(self.CreateGlobalSpotDetectionDockWidget)
        self.showIsData = True
        # self.actionShowDatas.triggered.connect(self.showDatas)
        # self.actionShowMasks.triggered.connect(self.showMasks)
        self.actionShowDatas.triggered.connect(self.showToggle)
        # self.actionShowMasks.triggered.connect(self.showToggle)

        self.actionRegistration.triggered.connect(self.CreateRegistrationWidget)


        self.actionWrite_mask_on_file.triggered.connect(  self.write_mask_on_file  )
        self.actionLoad_mask_from_file.triggered.connect(  self.read_mask_from_file  )
        ## self.actionRemote_load.triggered.connect(  self.remoteMaskload  )
        self.actionRemote_load.triggered.connect(  self.remoteMaskload  )
        self.actionPush_mask_remotely.triggered.connect(  self.PushMask  )
        self.actionWrite_masksDict_on_file.triggered.connect(  self.write_masksDict_on_file  )
        self.actionLoad_masksDict_from_file.triggered.connect(  self.load_masksDict_from_file  )
        
        
        self.actionExit.triggered.connect(self.Exit)
        self.actionConfirm_And_Exit.triggered.connect(self.confirm_and_Exit)

        self.actionLoad_image_from_hdf5.triggered.connect(self.load_image_from_hdf5)

        self.actionLoad_maskDict_from_hdf5.triggered.connect(self.load_maskDict_from_hdf5)
        self.actionWrite_maskDict_to_hdf5.triggered.connect(self.write_maskDict_to_hdf5)

        self.load_user_input = {}

    def Exit(self):
        self.isOK=0
        self.close()
    def confirm_and_Exit(self):
        self.isOK=1
        self.close()

    def remoteMaskload(self):
        ## 
        pathtolima = "id20/xrs/mpx-ram"
        dev = PyTango.DeviceProxy(pathtolima ) 

        string =dev.getRoifrompickle()
        masks = pickle.loads(string)
        masksDict = self.masks2MasksDict(masks)

        self.load_masksDict(masksDict)
            

    def PushMask(self):
        ## 
 
        masks = self.getMasks()

        correspondance = self.getLabelCorrespondance()

        new_masks = []
        for m in masks:
            new_masks.append(  ( correspondance[m[0]-1], m[1],m[2])   )
        masks = new_masks
    
        
        print( " de mainWindow je vais pusher : " , masks)
        print( " ========================================" )
        masks_string = pickle.dumps(masks)
        
        pathtolima = "id20/xrs/mpx-ram"
        dev = PyTango.DeviceProxy(pathtolima ) 
        # open("mask","w").write(masks_string)
        dev.setRoifromPickle( masks_string  )

    def getMasks(self) :
        totnofrois, gmask = self.recomposeGlobalMask()
        globalMask = self.mws[0].getSelectionMask().astype("i")
        masks = [] 
        nrois = globalMask.max()
        nummaxroi = (globalMask).max()
        masks = []
        for n in range(1,totnofrois +1):
            rawindices = numpy.nonzero( globalMask == n )
            if(len(rawindices)):
                if(len(rawindices[0])):
                    Y1 = rawindices[0].min()
                    Y2 = rawindices[0].max()+1
                    X1 = rawindices[1].min()
                    X2 = rawindices[1].max()+1
                    submask = (globalMask[Y1:Y2, X1:X2 ] == n).astype("i")
                    # submask = (globalMask[Y1:Y2, X1:X2 ] / n).astype("i"
                    masks.append(  ( n, [Y1,X1], submask  )  )
            else:
                masks.append(  ( n, [-1,-1], numpy.array([[1]])  )  )
        return masks

    # def getMasks(self) :
    #     totnofrois = self.recomposeGlobalMask()
    #     globalMask = self.mws[0].getSelectionMask().astype("i")
    #     masks = [] 
    #     nrois = globalMask.max()
    #     nummaxroi = (globalMask).max()
    #     masks = []
    #     for n in range(1, nummaxroi+1):
    #         rawindices = numpy.nonzero( globalMask == n )
    #         if(len(rawindices)):
    #             if(len(rawindices[0])):
    #                 Y1 = rawindices[0].min()
    #                 Y2 = rawindices[0].max()+1
    #                 X1 = rawindices[1].min()
    #                 X2 = rawindices[1].max()+1
    #                 submask = (globalMask[Y1:Y2, X1:X2 ] == n).astype("i")
    #                 masks.append(  ( n, [Y1,X1], submask  )  )
    #     return masks

    def getRoiObj(self):
        
        roiob=xrs_rois.roi_object()
        if self.image is not None:
            roiob.load_rois_fromMasksDict( self.getMasksDict(), newshape = self.image.shape )
            roiob.input_image = self.image
        return roiob


    def write_maskDict_to_hdf5(self, fn=None):
        if(fn is None) :
            filename =  (Qt.QFileDialog.getSaveFileName())
            if isinstance(filename, tuple):
                filename = str(filename[0])
        else:
            filename = fn
            
        if filename is None:
            return
        
        self.load_user_input ["roi"]= str(filename)+":/roi_from_selector/"

        self.user_input_signal.emit( self.load_user_input  )
        
        self.saveMaskDictOnH5( str(filename)+":/roi_from_selector/"  )

    ## w4r.load_maskDict_from_givenhdf5andgroup(roiaddress[0], roiaddress[1] )
    def  saveMaskDictOnH5( self, roiaddress, masktype ="ROI" , filterMask = None ) :
        
        if not ( isinstance(roiaddress , list) or isinstance(roiaddress, tuple) ):
            roiaddress =   xrs_utilities.split_hdf5_address( roiaddress )
                 
        h5=h5py.File(roiaddress[0],'a')
        image4roi = self.image
        if  masktype=="ROI":
            h5.require_group(roiaddress[1]+"/rois_definition")
            h5group =  h5[roiaddress[1]+"/rois_definition"]
        else:
            h5.require_group(roiaddress[1])
            h5group =  h5[roiaddress[1]]
        README="""image : the image upon which the rois have been selected
        rois_dict : a directory contain for each selected roi a subdirectory with the rois mask and its bottomleft corener position.
        """
        h5_assign_force(h5group , "image" , image4roi   )
        h5_assign_force(h5group , "README" , README   )
        if masktype=="ROI":
            masks = self.getMasksDict()
            h5group.require_group("rois_dict")
            h5group=h5group["rois_dict"]
            # remove all previous ROIS from h5group
            for key in h5group:
                del h5group[key]
                
            ### filterMask = None
            xrs_rois.write_rois_toh5(h5group, self.getMasksDict(), filterMask=filterMask )
        else:
            totnrois, filtermask  =  self.recomposeGlobalMask()
            filtermask = np.equal(0,filtermask).astype("i")
            h5_assign_force( h5group,"filter",  filtermask)
        h5.flush()
        h5.close()





    
    def getMasksDict(self) :
        masks = self.getMasks()
        return self.masks2MasksDict(masks)

    def masks2MasksDict(self,masks):
        masksDict={}
        for m in masks:
            if(m[1][0]>=0 and m[1][1]>=0  ):
                masksDict[ self.labelformat% (m[0]-1)  ]=[m[1],m[2]]
        return masksDict
        


    def recomposeGlobalMask(self):
      offset=0
      globalMask = self.mws[0].getSelectionMask().astype("i")
      for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
          localmask = numpy.less(0,mw.getSelectionMask() ) 
          globalMask[geo] =  offset*localmask +  mw.getSelectionMask()
          offset += nofrois
      self.mws[0].setSelectionMask(globalMask  )
      return offset, globalMask

    
    def decomposeGlobalMask(self):
      offset=0
      globalMask = self.mws[0].getSelectionMask().astype("i")
      for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
          localmask = numpy.less(0,globalMask[geo]  ) 
          Mask =    globalMask[geo] - offset *localmask
          offset += nofrois
          mw.setSelectionMask(Mask  )
          

    def write_mask_on_file(self):
        filename =  Qt.QFileDialog.getSaveFileName()
        if isinstance(filename, tuple):
            filename = filename[0]
        print( filename)
        if filename is not None:
            filename=str(filename)
            self.recomposeGlobalMask()
            globalMask = self.mws[0].getSelectionMask().astype("i")
            ef = edf.EdfFile( filename, "w+")
            ef.WriteImage( {}, globalMask )

    def write_masksDict_on_file(self):
        filename =  Qt.QFileDialog.getSaveFileName()
        if isinstance(filename, tuple):
            filename = filename[0]

        
        print( filename)
        if filename is not None:
            filename=str(filename)
            masksDict = self.getMasksDict()
            filename=str(filename)
            f = open(filename, 'wb')
            pickle.dump(masksDict ,  f)
            f.close()
            
    def load_masksDict_from_file(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        if isinstance(filename, tuple):
            filename = filename[0]

        print( filename)
        if filename is not None:
            filename=str(filename)
            f = open(filename, 'rb')
            masksDict = pickle.load( f)
            f.close()
            self.load_masksDict( masksDict)


    def load_masksDict(self, masksDict):
            self.recomposeGlobalMask()
            mask  = self.mws[0].getSelectionMask().astype("i")
            mask[:]=0
            mask =   convert_redmatrix_to_matrix(masksDict,mask, offsetX=0, offsetY=0)
            self.mws[0].setSelectionMask(mask)
            self.decomposeGlobalMask()
            self.annotateAllMasksCallBack()


    def read_mask_from_file(self):
        filename =  Qt.QFileDialog.getOpenFileName()
        if isinstance(filename, tuple):
            filename = filename[0]

        if filename is not None:
            filename=str(filename)
            ef = edf.EdfFile( filename, "r")
            mask = ef.GetData(0)
            self.recomposeGlobalMask()
            self.mws[0].setSelectionMask(mask)
            self.decomposeGlobalMask()



 
    def detectionCallBack(self, thr_s="", Hough=False):
        print( " in  detectionCallBack, Hough " , Hough)
        itab =  self.viewsTab.currentIndex()
        if itab==0:
            return 
        
        globalMask = self.mws[0].getSelectionMask().astype("i")
        roiroiMask = self.roiroiw.getSelectionMask().astype("i")
        
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:itab-1], self.mws[1:itab]):
            offset += nofrois

        (name,geo,nofrois), mw   =     self.names_geos_nofrois[itab-1], self.mws[itab]
        self.detectSpotsSubMask( name,geo,nofrois, mw  , globalMask , roiroiMask, offset,thr_s=thr_s , Hough=Hough  )
        self.mws[0].setSelectionMask( globalMask )


    def get_geo(self):
        if self.geo_informations is None:
            subset_infos = xrs_rois.get_geo_informations( (self.image.shape+(self.layout,) ) )
        else:
            subset_infos = self.geo_informations


        dl = subset_infos["analyser_nIDs"]
        dl1k = list(dl[list(dl.keys())[0]].keys())
        if len(dl1k)==1:
            for t in imageview.all_layouts:
                t.setCurrentIndex(1)
                
        return subset_infos


    def getLabelCorrespondance(self):
        res = []
        itab=1
        # print " CORRESPONDANCE "
        # print self.mws[1:]
        # print self.names_geos_nofrois[:]
        
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
            subset_infos = self.get_geo()
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]
 
            # res.extend(list(numpy.array(a_ids) +len(res) ))
            res.extend(list(numpy.array(range(1,1+12) ) +len(res) ))
            # print " res "  ,res
            itab+=1
        return res
    
    def annotateOneMaskCallBack(self):
        itab =  self.viewsTab.currentIndex()     
        if itab>0:
            subset_infos = self.get_geo()

            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

            self.mws[itab].annotateSpots( a_ids, self.get_offset(  itab ) )

    def annotateAllMasksCallBack(self):
        itab=1
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):

            subset_infos = self.get_geo()
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

            mw.annotateSpots( a_ids, self.get_offset(itab) )
            itab+=1

    def get_offset(self,itab):
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:itab-1], self.mws[1:itab]):
            offset += nofrois   
        return offset

    def relabeliseAllMasksCallBack(self):
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
            self.relabeliseSpots( mw, nofrois, name , geo, offset)

    def relabeliseOneMaskCallBack(self):
        itab =  self.viewsTab.currentIndex()     
        if itab==0: return
        offset=0
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:itab-1], self.mws[1:itab]):
            offset += nofrois
        (name,geo,nofrois), mw  =    self.names_geos_nofrois[itab-1], self.mws[itab]
        self.relabeliseSpots( mw , nofrois, name, geo, offset)

    def relabeliseSpots(self,mw, nofrois, name, geo, offset):

        mask = mw.getSelectionMask( ).astype("i")
        mask = (mask>0).astype("i")
        newmask = spotdetection.relabelise(mask,mask, nofrois)
        self.checkNspots(newmask.max(), nofrois , name) 
        mw.setSelectionMask( newmask )

        globalMask = self.mws[0].getSelectionMask().astype("i")
        globalMask[geo] = newmask
        self.mws[0].setSelectionMask( globalMask )

    def resetOneMask(self):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        mw = self.mws[itab]
        mw.graph.clearMarkers()
        mask = mw.getSelectionMask( ).astype("i")
        mask[:]=0
        mw.setSelectionMask( mask)

    def resetAllMasks(self):
        ret = self.warnForGloablChange()
        print( ret)
        if ret:
            for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois[:], self.mws[1:]):
                mw.graph.clearMarkers()
                mask = mw.getSelectionMask( ).astype("i")
                mask[:]=0
                mw.setSelectionMask( mask)



    def threshold(self,itab,value ):
        globalMask = self.mws[0].getSelectionMask().astype("i")
        name,geo,nofrois =   self.names_geos_nofrois[itab-1]

        mw = self.mws[itab]
        mask = mw.getSelectionMask( )
        print( mask.sum())

        data = self.image[geo]

        mask = spotdetection.threshold_mask(mask, data ,  value  ) 
        mw.setSelectionMask(mask  )

        globalMask[geo] = mask
        self.mws[0].setSelectionMask( globalMask )
        
        print( mask.sum())


    def localThresholdCallBack(self,value):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        self.threshold(itab,value)

    def globalThresholdCallBack(self,value):
        ret = self.warnForGloablChange()
        print( ret)
        if ret:
            for itab in range(1, len(self.mws )  ) :
                self.threshold(itab,value)


    def fatten(self, itab, value  ) :
        globalMask = self.mws[0].getSelectionMask().astype("i")
        name,geo,nofrois =   self.names_geos_nofrois[itab-1]

        mw = self.mws[itab]
        mask = mw.getSelectionMask( )
        if value>0:
            mask = spotdetection.grow_mask(mask, 1+2*value  )
        else:
            mask = spotdetection.shrink_mask(mask, 1-2*value  )
        mw.setSelectionMask(mask  )

        globalMask[geo] = mask
        self.mws[0].setSelectionMask( globalMask )


    def fatteningCallBack(self,value):
        itab =  self.viewsTab.currentIndex()
        if itab==0: return 
        self.fatten(itab,value)

    def GlobalfatteningCallBack(self,value):
        print( " in GlobalfatteningCallBack  ", value  )
        ret = self.warnForGloablChange()
        print( ret)
        if ret:
            for itab in range(1, len(self.mws )  ) :
                self.fatten(itab,value)


        
    def checkNspots(self,nspots, nofrois,name ) :
        if nspots != nofrois:
            if not SKIP_WARNING :
                msgBox = Qt.QMessageBox ()
                msgBox.setText("Warning: found %d spots instead of expected %d "%(nspots , nofrois ) );
                msgBox.setInformativeText("For detector %s " %  name);
                msgBox.setStandardButtons(Qt.QMessageBox.Ok );
                msgBox.setDefaultButton(Qt.QMessageBox.Cancel);
                ret = msgBox.exec_();

    def    detectSpotsSubMask( self, name,geo,nofrois, mw   , globalMask,roiroiMask,  offset ,thr_s="", Hough=False):
        itab =  self.viewsTab.currentIndex()
        if itab: 
            subset_infos = self.get_geo() 
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]

                
        else:
            a_ids = None

        submatrix = mw.getImageData()
        subroiroi = roiroiMask[geo]
        tval = -1
        try:
            tval = float(thr_s)
        except:
            tval=-1


        geos = self.get_geo()
        dl    = geos["analyser_nIDs"]
        dll   = dl[list(dl.keys())[0]]
        ll = dll[list(dll.keys())[0]]
        nrois = len(ll)
        
        mask  = spotdetection.get_spots_mask( submatrix,subroiroi,
                                              median_size = 5 , tval=tval, nofroi=nrois, Hough=Hough ) 
    
        self.checkNspots(mask.max(), nofrois , name) 

        mw.setSelectionMask( mask)
        globalMask[geo]=mask+offset*numpy.less(0,mask)



        mw.annotateSpots( a_ids , self.get_offset(itab)  )
        
    def GlobaldetectionCallBack(self, warn=True,  thr_s="", Hough=False):
        print( " GlobaldetectionCallBack ",  warn)
        if 1 or warn:
            print( " in Global  detectionCallBack" )
            ret = self.warnForGloablChange()
        else:
            ret = True
        if ret:
            offset=0
            globalMask = self.mws[0].getSelectionMask().astype("i")
            roiroiMask = self.roiroiw.getSelectionMask().astype("i")
         
            for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois, self.mws[1:]):
                self.detectSpotsSubMask( name,geo,nofrois, mw  , globalMask ,roiroiMask,  offset ,thr_s=thr_s, Hough=Hough)
                offset += nofrois
            self.mws[0].setSelectionMask( globalMask )


    def showToggle(self):
        if self.showIsData :
            self.showMasks()
            self.showIsData  = not self.showIsData 
        else:
            self.showDatas()
            self.showIsData  = not self.showIsData 
            
    def showMasks(self):
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois, self.mws[1:]):
            mask = mw.getSelectionMask().astype("i")
            mw.setImageData(mask  , xScale=(0.0, 1.0), yScale=(0., 1.))

    def showDatas(self):
        Data = self.mws[0].getImageData()
        for (name,geo,nofrois), mw in    zip(self.names_geos_nofrois, self.mws[1:]):
            d  = Data[geo]
            mask = mw.getSelectionMask().astype("i")
            if mask.sum():
                mm = (d*mask).max()
                d=numpy.minimum(mm,d )
            mw.setImageData(d  , xScale=(0.0, 1.0), yScale=(0., 1.))


    def warnForGloablChange(self):
        if not SKIP_WARNING:
            msgBox = Qt.QMessageBox ()
            msgBox.setText("You are going to recalculate the GLOBAL mask");
            msgBox.setInformativeText("This will reset all modifications to local masks. Do you want to proceed?");
            msgBox.setStandardButtons(Qt.QMessageBox.Ok |  Qt.QMessageBox.Cancel);
            msgBox.setDefaultButton(Qt.QMessageBox.Cancel);
            ret = msgBox.exec_();
            return ret==Qt.QMessageBox.Ok
        else:
            return True
        

    def CreateSpotDetectionDockWidget(self):
        w = spotdetectioncontrol(self, QtCore.Qt.Widget,detectionCallBack=self.detectionCallBack,  
                                 fatteningCallBack=self.fatteningCallBack,
                                 thresholdCallBack  =   self.localThresholdCallBack, 
                                 annotateMaskCallBack=self.annotateOneMaskCallBack,
                                 relabeliseMaskCallBack=self.relabeliseOneMaskCallBack,
                                 resetMask           = self.resetOneMask
                                 )
        self.addDockWidget ( QtCore.Qt.LeftDockWidgetArea, w )
        # w.setAllowedAreas (QtCore.Qt.AllDockWidgetAreas)
        w.show()

    def globregistration(self):
        print( " in globregistration ")
        for itab in range(1,len(self.mws)):
            self.registerTab(itab)

    def registration(self):
        itab =  self.viewsTab.currentIndex()     
        print( self.layouts[itab].currentText())
        self.registerTab(itab)

    def registerTab(self, itab):
        if itab>0:
            
            subset_infos = self.get_geo() 
            if "3x4" in str(self.layouts[itab].currentText() ):
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]]["3x4"]
            else:
                a_ids = subset_infos["analyser_nIDs"][ subset_infos["subnames"][itab-1]] ["Vertical"]
                
            name,geo,nofrois = self.names_geos_nofrois[itab-1]
            self.registerSpots( self.mws[itab], self.layouts[itab].currentText(), name, nofrois , a_ids)


            self.mws[itab].annotateSpots(a_ids , self.get_offset(itab))

    def registerSpots( self,  mw, layoutString , name, nofrois, a_ids) :

        mask = mw.getSelectionMask().astype("i")
        newmask = spotdetection.relabelise(mask,mask, nofrois)
        self.checkNspots(newmask.max(), nofrois , name) 
        mask = newmask
        nspots = mask.max()
        spots = []
        for i in range(1,nspots+1):
            zone = (mask==i)
            mask[zone]=i+100
            m = zone.astype("f")
            msum=m.sum()
            if msum:
                ny,nx = m.shape
                px= (m.sum(axis=0)*numpy.arange(nx)).sum()/msum
                py= (m.sum(axis=1)*numpy.arange(ny)).sum()/msum
                spots.append((py,px,i))

        self.checkNspots(len(spots), nofrois , name) 
            
        print( type(str(layoutString)))
        if "3x4" in str(layoutString):
            
            positions = [ [px,py]   for (py,px,i) in spots ]
            print( str( numpy.array(positions)))
            choices = match.register(numpy.array(positions)  )
            

            newspots = []
            for (cx,cy),(y,x,i) in zip(choices, spots ):
                print( x,y, cx, cy)
                newspots.append((y,x,i,  int((cy*4+cx)+1)  ) )
        else:
            spots.sort()
            newspots = []
            for k,(y,x,i) in enumerate( spots):
                newspots.append((y,x,i,k+1) ) 

        print( " NEWSPOTS ", newspots)

        for (y,x,i,k) in newspots:
            zone = (mask==(i+100))
            mask[zone] = a_ids[k-1]
                
        mw.setSelectionMask(mask)
            
    def CreateRegistrationWidget(self):

        w = spotregistrationcontrol(self, QtCore.Qt.Widget,globregistrationCallBack=self.globregistration,  
                                    registrationCallBack=self.registration
                                )

        self.addDockWidget ( QtCore.Qt.LeftDockWidgetArea, w )
        w.show()

    def CreateGlobalSpotDetectionDockWidget(self):
        w = spotdetectioncontrol(self, QtCore.Qt.Widget,detectionCallBack=self.GlobaldetectionCallBack, 
                                 fatteningCallBack=self.GlobalfatteningCallBack,
                                 thresholdCallBack  =   self.globalThresholdCallBack, 
                                 annotateMaskCallBack=self.annotateAllMasksCallBack,
                                 relabeliseMaskCallBack=self.relabeliseAllMasksCallBack,
                                 resetMask           = self.resetAllMasks
                                 )
        w.detectionButton.setText("GLOBAL detection")
        w.HoughDetection.setText("GLOBAL Hough")
        self.addDockWidget ( QtCore.Qt.LeftDockWidgetArea, w )
        # w.setAllowedAreas (QtCore.Qt.AllDockWidgetAreas)
        self.globalSpotDectionWidget = w
        w.show()

    def set_roiob(self, roiob):
        self.roiob=roiob

    def LoadRemote(self):

        if self.roiob is None:
            mb = qt.QMessageBox();
            mb.setText("No roi-object has been associated to the roi manager. Cannot Load.");
            mb.exec_();
            return
        self.showImage(self.roiob)


    def create_viewWidget(self, image = None, nofrois = None, changeTagOn=False, isglobal=False , layoutsNames = None) :
        view = imageview(self, isglobal, layoutsNames =layoutsNames)        
        maskW  =myMaskImageWidget.MaskImageWidget(self, aspect=True,
                                                  profileselection=True,
                                                  maxNRois=nofrois )
        maskW.setY1AxisInverted(1)

        maskW.setAcceptDrops(True)


        maskW.setDefaultColormap(2, logflag=True)

        maskW.changeTagOn = changeTagOn
        maskW.setImageData(image  , xScale=(0.0, 1.0), yScale=(0., 1.))
        maskW.setSelectionMask(image*0 , plot=True)
        view.roiContainerWidget.layout().addWidget(maskW)
        return view, maskW
            
    def showImage(self, image, geo_informations=None):
        self.image=image
        if geo_informations is None:
            geo_informations = self.get_geo()  
        self.geo_informations = geo_informations

        self.image=image

        subset_infos = geo_informations

        self.viewsTab.clear()       
        totNofrois = subset_infos["nofrois"]
        nofrois = totNofrois//len( subset_infos["subgeos"]  )
        
        # "analyser_nIDs": {"DETECTOR":{ "Vertical": V1234} }

        layoutsNames = list(list(subset_infos["analyser_nIDs"].items())[0][1].keys())
        layoutsNames.sort(reverse=True)
        
        self.mws=[]
        self.layouts=[]
        view, mw = self.create_viewWidget(image = image, nofrois = totNofrois , changeTagOn = False, isglobal=True  , layoutsNames = layoutsNames)
        
        self.viewsTab.addTab(view, "Global")
        self.mws.append(mw)
        self.layouts.append(view.registeringLayoutComboBox)

        self.names_geos_nofrois = list(zip(subset_infos["subnames"],  subset_infos["subgeos"], [nofrois]*len(subset_infos["subgeos"])  ))

        for name,geo,nofr in self.names_geos_nofrois :
            view, mw = self.create_viewWidget(image = image[     geo ], nofrois = nofrois, changeTagOn = True , layoutsNames = layoutsNames )
            
            self.viewsTab.addTab(view, name)
            self.mws.append(mw)
            self.layouts.append(view.registeringLayoutComboBox)
            

        view, roiroiw = self.create_viewWidget(image = image, nofrois = 1 , isglobal=True, changeTagOn = False  , layoutsNames = layoutsNames)
        
        self.viewsTab.addTab(view, "ROI of ROIS")
        self.roiroiw = roiroiw
        self.layouts.append(view.registeringLayoutComboBox)



            
    def LoadLocal(self, sf=None, fn=None, ns=None):
        if not FASTDEBUG and sf is None :
            w = localfilesdialog.localfilesdialog(self.load_user_input)
            result = w.exec_()
            if not result:
                return
            sf =  str(w.SpecFileName_lineEdit.text())
            fn =  str(w.FileName_lineEdit.text())
            ns =  w.ScanNumber_spinBox.value()
        elif sf is not None:
            pass
        else:
            sf = "/data/id20/inhouse/data/run2_18/run5_ihr/hydra"
            fn = "/data/id20/inhouse/data/run2_18/run5_ihr/edf/hydra_0000.edf"
            ns =  189


        self.load_user_input = {  "sf":sf, "fn":fn    }

        self.user_input_signal.emit( self.load_user_input  )

        
        template = getTemplateName( fn )

        s=specfile.Specfile(sf)
        Scan = s[ns-1]
        numbers = Scan.datacol("ccdno")

        roiob=None

        print( numbers)
        imagesum=numpy.array([0.0])
        for n in numbers:

            image = [ edf.EdfFile(tok%n,"r").GetData(0)   for tok in template ] 
            image = numpy.concatenate(image, axis=0)
            # if roiob is None:
            #     roiob = xrs_rois.roi_object() 
            #     shape = image.shape
            #     roiob.prepare_rois(  len(numbers) ,  "%dx%d"%shape ) 
            # roiob.process(image)
            imagesum=imagesum+image
            
        self.geo_informations=None
        self.showImage(imagesum)
        self.roiob = roiob

        if False and FASTDEBUG:
            self.GlobaldetectionCallBack(warn=False)


    def load_maskDict_from_hdf5(self):
        filename =  Qt.QFileDialog.getOpenFileName(None,'Open hdf5 file with rois',filter="hdf5 (*h5)\nall files ( * )"  )
        if isinstance(filename, tuple):
            filename = filename[0]
        if filename is None:
            return
        
        filename=str(filename)
        
        print( filename)
        if len(filename):
            import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget
            storage=[None]
            def mySlot(ddict):
                name = ddict["name"]
                storage[0]=name
            # browse
            self.__hdf5Dialog = hdf5dialog()
            self.__hdf5Dialog.setWindowTitle('Select a Group  containing roi_definitions by a double click')
            self.__hdf5Dialog.mainLayout =  self.__hdf5Dialog.verticalLayout_2
            fileModel = HDF5Widget.FileModel()
            fileView = HDF5Widget.HDF5Widget(fileModel)
            hdf5File = fileModel.openFile(filename)
            shiftsDataset = None
            fileView.sigHDF5WidgetSignal.connect(mySlot)
            self.__hdf5Dialog.mainLayout.addWidget(fileView)
            self.__hdf5Dialog.resize(400, 200)

            ret = self.__hdf5Dialog.exec_()
            print( ret)
            hdf5File.close()
            
            if ret:
                print( " Obtained " )
                name =  storage[0]
                print( name)
                file= h5py.File(filename,"r")
                datagroup = file[name]
                masks={}
                xrs_rois.load_rois_fromh5(datagroup,masks)
                file.close()
                self.load_masksDict(masks)

                
    def loadMaskDictFromH5(self,  filename ):
        self.load_maskDict_from_givenhdf5andgroup(filename)


    def load_maskDict_from_givenhdf5andgroup(self, filename, gname=None):

        if gname is None:
            filename, gname  =    xrs_utilities.split_hdf5_address(filename )

        
        file= h5py.File(filename,"r")
        datagroup = file[gname]
        masks={}
        xrs_rois.load_rois_fromh5(datagroup,masks)
        file.close()
        self.load_masksDict(masks)
        


                
    def load_image_from_hdf5(self):
        print( " load " )
        filename =  Qt.QFileDialog.getOpenFileName()
        if isinstance(filename, tuple):
            filename = filename[0]

        if filename is None:
            return
        
        print( " OK " )
        filename=str(filename)
        print( filename)
        if len(filename):
            import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget

            storage=[None]

            def mySlot(ddict):
                name = ddict["name"]
                storage[0]=name
                print( " MY SLOT " )
                print( name)

            # browse
            self.__hdf5Dialog = hdf5dialog()
            self.__hdf5Dialog.setWindowTitle('Select your data set by a double click')
            self.__hdf5Dialog.mainLayout =  self.__hdf5Dialog.verticalLayout_2
            fileModel = HDF5Widget.FileModel()
            fileView = HDF5Widget.HDF5Widget(fileModel)
            hdf5File = fileModel.openFile(filename)
            shiftsDataset = None
            fileView.sigHDF5WidgetSignal.connect(mySlot)
            self.__hdf5Dialog.mainLayout.addWidget(fileView)
            self.__hdf5Dialog.resize(400, 200)
            
            # self.__hdf5Dialog.setModal(True)
            # self.__hdf5Dialog.show()

            ret = self.__hdf5Dialog.exec_()

            hdf5File.close()
            
            if ret:
                print( " Obtained " )
                name =  storage[0]
                file= h5py.File(filename,"r")
                image4roi = file[name][:]
                file.close()

                self.showImage(image4roi)


# @ui.UILoadable
class hdf5dialog(Qt.QDialog):
    def __init__(self, parent=None):

        super(hdf5dialog  , self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" ,  "hdf5dialog.ui" ), self)
        
        # Qt.QDialog.__init__(self, parent)
        # self.loadUi()  # load the ui file
       

def convert_redmatrix_to_matrix( masksDict,mask, offsetX=0, offsetY=0):
    for key, (pos,M)  in six.iteritems(masksDict):
        num=int("".join([c for c in key if c.isdigit()]))
        S = M.shape
        inset =    (slice(offsetY+pos[0]  , offsetY+pos[0]+S[0]   ), slice(  offsetX+pos[1]  , offsetX+pos[1]+S[1] ) )
        M=numpy.less(0,M)
        mask[  inset   ][M>0] =  (num+1)*M[M>0]
    return mask


def getTemplateName(name):
    dirname=os.path.dirname(str(name))
    name=os.path.basename(str(name))
    ls = len(name)
    fine = None
    inizio=None

    for n in list(range(ls))[::-1]:
        if fine is None and name[n] in "1234567890":
            fine = n+1
        if name[n] in  "1234567890":
            inizio=n
        if fine is not None and name[n] not in "1234567890":
            break
        
    print( name)
    print( inizio)
    print( fine)
        
    name = name[:inizio]+"%" + ("0%d"%(fine-inizio)+"d"+name[fine:])

    print( name)

    if name[:2] in ["h_", "v_"]:
        return [dirname+"/"+"h_"+name[2:],dirname+"/"+"v_"+name[2:]]
    else:
        return [dirname+"/"+name]

_cross_data = "\
\x00\x00\x02\x71\
\x89\
\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d\x49\x48\x44\x52\x00\
\x00\x00\x10\x00\x00\x00\x10\x08\x06\x00\x00\x00\x1f\xf3\xff\x61\
\x00\x00\x00\x06\x62\x4b\x47\x44\x00\xff\x00\xff\x00\xff\xa0\xbd\
\xa7\x93\x00\x00\x00\x09\x70\x48\x59\x73\x00\x00\x0d\xd7\x00\x00\
\x0d\xd7\x01\x42\x28\x9b\x78\x00\x00\x00\x07\x74\x49\x4d\x45\x07\
\xde\x07\x10\x0d\x3b\x14\x28\x69\xf4\x66\x00\x00\x01\xfe\x49\x44\
\x41\x54\x38\xcb\x85\x93\xcf\x6a\x14\x41\x10\xc6\x7f\x5d\xdd\x59\
\x77\x3d\x86\x44\x14\xd4\xa0\x41\xfc\x17\xc5\x9b\x20\x1e\x66\x86\
\x3d\xa9\x2f\x10\x04\x4f\x7a\x51\x1f\xc3\x47\x10\x05\x51\x3c\x08\
\x82\x4f\x90\xcb\x30\x33\x37\x11\xd4\x8d\xa2\x17\x49\x02\x1a\x3c\
\x18\x04\x0f\x2e\xe8\x6c\xa6\xbb\x3c\xa4\x27\x4c\x16\xc5\x86\xa6\
\xa1\xaa\xbe\xaf\xbf\xae\xfe\xca\x0c\x87\x49\x8f\x9d\x25\x40\x88\
\xa7\xc6\xfd\xbf\x98\x97\x08\xb6\x31\xe9\x00\x35\x46\x4d\x9e\x57\
\x13\xdf\xc8\x29\xdf\xc8\x25\xef\x65\x2e\xcf\xab\x7a\xba\x0e\xb0\
\xd2\x61\xb4\x80\xcf\xf3\xaa\x6e\xb6\xed\xbd\x34\xc9\x14\x18\x01\
\x25\xca\x66\x9a\x64\x3f\x7c\x23\x57\x22\x91\x6f\x55\x99\xf8\x04\
\x07\x78\x2c\xea\x6b\x59\x03\x8e\xf0\xef\xf5\xb8\xac\x8a\x5b\xc3\
\x61\x32\x00\x42\xab\xc0\xe7\x79\x55\xfb\x5a\x56\x80\x59\xeb\xc2\
\x00\xb8\xb3\x07\x66\xb8\x0c\x9c\x00\x6e\xa6\x49\xb6\x1c\x55\x9b\
\xb6\x39\x92\x26\xd9\x19\x20\x03\xdc\xb1\xa5\xad\xba\xbf\x7f\xfb\
\x51\x4b\x62\x0c\x4b\x22\xe1\x2d\x30\x89\x74\x4f\xe3\x53\x42\x4b\
\x00\x70\x3d\x9e\xbd\xb5\xd5\x83\x1b\x46\x34\x64\xcb\x1f\x1e\x18\
\xd1\x79\x23\x61\x23\x78\x39\x04\x7c\x8e\x35\xfb\xd2\x24\x3b\x07\
\x88\x00\xe2\x7a\x7e\x02\x1c\xdf\x15\x0b\x0b\xbf\xc6\xbd\xd1\xfa\
\xea\x01\x23\xa2\x3f\x83\x97\x59\x60\x7d\xaa\x17\x8b\x80\x0a\x10\
\xfc\xb6\xed\x01\x5b\x9d\xa4\x07\xae\x6d\x7e\x9a\xeb\xe7\x79\x55\
\x97\x55\xf1\x15\x78\x38\x45\x30\x6e\xcd\x23\xaa\x04\xe0\x7d\x4c\
\xd4\xc6\xe8\x61\x6b\xc3\xf7\x10\xe4\x74\x9a\x64\xef\x00\xca\xaa\
\xb8\x0d\x3c\x6b\xd1\xbd\x41\x53\xec\x69\x62\x59\x15\x4f\xe2\xcd\
\xa3\xa2\x2c\xbf\x79\x2f\x17\x50\xde\x00\xe7\xd3\x24\x7b\x1d\x49\
\x6e\x44\xfc\x0b\x0d\xc6\x01\xa1\xf5\x81\x45\x51\x1f\xe4\x22\x4a\
\x05\x7c\x01\x8e\x4e\x49\xfe\x18\x2f\x58\xb4\x2e\xcc\x77\xbf\x71\
\xc7\x89\x06\x71\xd6\xbf\x04\xd2\xbf\x80\x01\xce\x02\x7d\x23\xba\
\x80\xee\xce\x84\xee\xb1\xb1\x62\xac\x75\xe1\x55\x59\x15\x06\xb8\
\x0b\x3c\x07\x72\xe0\x3e\x70\xb5\xac\x8a\x93\x22\x3a\xc6\x60\x76\
\xfd\x35\x1c\x26\xfd\xce\x94\xd1\x25\xb4\x2e\x04\x37\xe3\xc3\xe4\
\xb7\x13\x55\x33\xd3\x99\x81\xb6\xce\xb8\x18\xb4\x31\xd0\xfa\xa0\
\x01\xc4\x37\x22\xbe\x11\xe9\xc6\xa6\xea\xfc\x1f\xdb\x37\xde\x59\
\x25\x68\xce\x04\x00\x00\x00\x00\x49\x45\x4e\x44\xae\x42\x60\x82\
\
"

def launch4MatPlotLib(layout=None, im4roi=None, manageQApp = True):
    
    if manageQApp:
        dostop=0
        app = Qt.QApplication.instance()
        if app is None:
            app=Qt.QApplication([])
            dostop = 1 
    w4r = mainwindow(layout=layout)
    w4r.showImage( im4roi )
    w4r.show()
    
    if manageQApp:
        app.exec_()
        if dostop:
            app.quit()

            
    return w4r



if __name__=="__main__":
    app = Qt.QApplication.instance()
    print( " APP IS " , app)
    if app is None:
        app=Qt.QApplication([])
    w = mainwindow()
    w.show()
    app.exec_()








