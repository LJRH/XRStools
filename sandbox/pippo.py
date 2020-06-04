from XRStools.roiNmaSelectionGui import roiNmaSelectionWidget

# import  XRStools.roiSelectionWidget as roiNmaSelectionWidget
from  PyQt4 import Qt, QtCore
app=Qt.QApplication([])
w4r = roiNmaSelectionWidget.mainwindow()
w4r.show()

sf = '/data/id20/inhouse/data/run5_17/run7_ihr/'
fn = "myroi.h5:/datas/ROI"
ns_s =  [613,617,621,625]
w4r.load_rois(sf, fn , ns_s )
w4r.loadMaskDictFromH5("pippo.h5:salva_qui" )


app   .exec_()
w4r.saveMaskDictOnH5( "pippo.h5:salva_qui"   ) 
md, meta=w4r.getMasksDict(withmeta=True)
print( md)
print ( " META ")
print( meta)

roiob = w4r.getRoiObj()
print ( " roiobj " )
print ( roiob)


w4r = None
