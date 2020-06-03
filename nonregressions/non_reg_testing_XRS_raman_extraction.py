from __future__ import print_function
#####
# this is important with old version of matplotlib  to avoid
# that matplotlib import pyqt4 while silx imports pyqt5
#
try:
    import PyQt4.QtCore
    print("   QT4 ")
except:
    pass

import os
import numpy as np
# from  pylab import  *
# ion()

import  XRStools
import  XRStools.ramanWidget
import  XRStools.ramanWidget.MainWindow

import XRStools.roiSelectionWidget

XRStools.roiSelectionWidget.SKIP_WARNING = True


def check_results(fn_a, fn_b):
    a = np.loadtxt(fn_a)
    b = np.loadtxt(fn_b)
    amin = a[0, 0]
    amax = a[-1,0]
    bmin = b[0, 0]
    bmax = b[-1,0]
    c0 = max(amin,bmin)
    c1 = min(amax,bmax)
    try:
        assert(  abs( (  c1-c0    )    /(amax-amin)      )   >0.9 )
        assert(  abs( (  c1-c0    )    /(bmax-bmin)      )   >0.9 )
    except:
        print (  "  amin,amx, bmib, bmax, c0, c1 " ,amin,amx, bmib, bmax, c0, c1  ) 
        raise



    ref = np.interp( a[:,0],  b[:,0],  b[:,1]      )

    try:
        assert(    (abs( a[:,1] -  ref)[np.less( c0, a[:,0])* np.less(  a[:,0], c1)]   ).sum()  <1.0e-5    )
    except:
        raise





#####################################################################
# for programmatic testing
from silx.gui.utils.testutils import QTest, TestCaseQt
from silx.gui import qt

config_file = None
config_file = "conf_raman_nonreg.yaml"


def qWait( ms=None):
    if ms is None:
        ms = cls.DEFAULT_TIMEOUT_WAIT

    if qt.BINDING in ('PySide', 'PySide2'):
        # PySide has no qWait, provide a replacement
        timeout = int(ms)
        endTimeMS = int(time.time() * 1000) + timeout
        while timeout > 0:
            _qapp.processEvents(qt.QEventLoop.AllEvents,
                                maxtime=timeout)
            timeout = endTimeMS - int(time.time() * 1000)
    else:
        QTest.qWait(ms )

##
########################################################################

class Pippo(TestCaseQt):
    def runTest(self):
        pass

_qapp = qt.QApplication.instance() or qt.QApplication([])

print( " CHIAMO " )
wtest = XRStools.ramanWidget.MainWindow.main( manageQApp = False)
qWait( ms=10)

if config_file is  None:
    ntab = 0
    QTest.mouseClick(  wtest.tabWidget.tabBar(),
                       qt.Qt.LeftButton ,
                       pos=qt.QPoint(
                           wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                           wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
    )

    qWait( ms=40)

    wtest.rsw.LoadLocal( sf="/data/id20/inhouse/data/run5_17/run7_ihr/hydra", fn="/data/id20/inhouse/data/run5_17/run7_ihr/edf/hydra_12828.edf", ns=611)

    print( " PIPPO " )

    qWait( ms=40)

    for act in wtest.rsw.menuActions .actions():
        if  str(act.objectName()) == "actionGlobalSpotDetection":
            act.trigger()
    QTest.mouseClick(  wtest.rsw.globalSpotDectionWidget.detectionButton ,qt.Qt.LeftButton, pos= qt.QPoint(3,3))

    wtest.rsw.write_maskDict_to_hdf5("spatialroi.h5")

    qWait( ms=40)


    ntab = 1
    QTest.mouseClick(  wtest.tabWidget.tabBar(),
                       qt.Qt.LeftButton ,
                       pos=qt.QPoint(
                           wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                           wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
    )

    qWait( ms=40)

    uinp =  wtest.sprsw.load_user_input
    print (uinp)


    # uinp["sf"] = "/data/id20/inhouse/data/run5_17/run7_ihr/" c'e' gia
    # uinp["roi"] = "spatialroi.h5:/datas/ROI" c'e' gia
    uinp["ns_s"] = [613,617,621,625]

    wtest.sprsw.LoadLocalOption(uinp)

    qWait( ms=40)
    wtest.showMaximized()
    qWait( ms=40)

    for ntab in range( 1,7):
        QTest.mouseClick(  wtest.sprsw.viewsTab.tabBar(),
                           qt.Qt.LeftButton ,
                           pos=qt.QPoint(
                               wtest.sprsw.viewsTab.tabBar().tabRect(ntab).x()+wtest.sprsw.viewsTab.tabBar().tabRect(ntab).width()/2 ,
                               wtest.sprsw.viewsTab.tabBar().tabRect(ntab).y()+   wtest.sprsw.viewsTab.tabBar().tabRect(ntab).height()/2 )
        )

        qWait( ms=40)
        mw = wtest.sprsw.mws[ntab]
        mw.gc.nofcomps.setText("3")
        mw.gc.choosedcomp.setText("0")
        mw.gc.threshold.setText("0.2")
        qWait( ms=40)
        mw.gc.pushButton_calccomps.clicked.emit(True)
        qWait( ms=40)
        mw.gc.pushButton_calccomps.clicked.emit(True)
        mw.gc.pushButton_threshold.clicked.emit(True)
        qWait( ms=40)

    wtest.sprsw.write_maskDict_to_hdf5_option("sp_roi.h5")

    QTest.qWait(40)

    ntab = 2
    QTest.mouseClick(  wtest.tabWidget.tabBar(),
                       qt.Qt.LeftButton ,
                       pos=qt.QPoint(
                           wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                           wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
    )

    qWait( ms=40)
else:
    wtest.loadConfigurationOption(config_file)

ntab = 3
QTest.mouseClick(  wtest.tabWidget.tabBar(),
                   qt.Qt.LeftButton ,
                   pos=qt.QPoint(
                       wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                       wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
)

qWait( ms=40)
myw = wtest.tabWidget.currentWidget().widget()
myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[0][1].setText("lowq")
l = myw.abstract[0][3:3+24]+ myw.abstract[0][36+3:36+3+24]
for t in l:
    t.cb.setChecked(True)
myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[1][1].setText("mediumq")
l = myw.abstract[1][3+24:3+24+24]
for t in l:
    t.cb.setChecked(True)

myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[2][1].setText("highq")

l = myw.abstract[2][3+60:3+60+12]
for t in l:
    t.cb.setChecked(True)

qWait( ms=40)
ntab = 4

QTest.mouseClick(  wtest.tabWidget.tabBar(),
                   qt.Qt.LeftButton ,
                   pos=qt.QPoint(
                       wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                       wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
)

qWait( ms=40)
myw = wtest.tabWidget.currentWidget().widget()
print(myw)

myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[0][1].setText("elastic")
qWait( ms=40)
myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[1][1].setText("ok0")
qWait( ms=40)
myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[2][1].setText("ok1")

qWait( ms=40)
myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[3][1].setText("ok2")


qWait( ms=40)
myw.plusLine.clicked.emit(True)
qWait( ms=40)
myw.abstract[4][1].setText("ok3")
qWait( ms=40)


fois = [4,1,4,4,4]
for row in range(5):
    for i in range(fois[row]):
        for k in myw.bottoni_info.keys():
            if myw.bottoni_info[k] == (+1,row):
                k.clicked.emit(True)
                qWait( ms=40)
                break
    for f in range(fois[row]):
        if row ==1 :
            myw.abstract[row][3+f].setValue(628)
        else:
            myw.abstract[row][3+f].setValue(611 +f*4 +row-(row>1))      


qWait( ms=40)


ntab = 5
QTest.mouseClick(  wtest.tabWidget.tabBar(),
                   qt.Qt.LeftButton ,
                   pos=qt.QPoint(
                       wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                       wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
)

qWait( ms=40)
myw = wtest.tabWidget.currentWidget()
myw.plusLine.clicked.emit(True)
qWait( ms=40)

print(myw.abstract)
myw.abstract[0][0].setText("H2O")
myw.abstract[0][0].editingFinished.emit()
qWait( ms=40)


abstr_stoichio, edges = myw.get_abstractN()

print("   abstr_stoichio      " , abstr_stoichio)
print("   edges      " , edges)

edges["O"]["K"]=1

myw.set_abstractN(abstr_stoichio, edges)

abstr_stoichio, edges
# myw.abstract[0][1].setText("lowq")
# l = myw.abstract[0][3:3+24]+ myw.abstract[0][36+3:36+3+24]
# for t in l:
#     t.cb.setChecked(True)
# myw.plusLine.clicked.emit(True)
# qWait( ms=40)
# myw.abstract[1][1].setText("mediumq")
# l = myw.abstract[1][3+24:3+24+24]
# for t in l:
#     t.cb.setChecked(True)

# myw.plusLine.clicked.emit(True)
# qWait( ms=40)
# myw.abstract[2][1].setText("highq")

qWait( ms=40)




ntab = 6
QTest.mouseClick(  wtest.tabWidget.tabBar(),
                   qt.Qt.LeftButton ,
                   pos=qt.QPoint(
                       wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                       wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
)

myw = wtest.tabWidget.currentWidget()

myw.spinBox_scan.setValue(611)
qWait( ms=40)


myw.pushButton.clicked.emit(True)

qWait( ms=40)



for ntab in [7,8,9] :
    QTest.mouseClick(  wtest.tabWidget.tabBar(),
                       qt.Qt.LeftButton ,
                       pos=qt.QPoint(
                           wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                           wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
    )

    myw = wtest.tabWidget.currentWidget()
    qWait( ms=40)
    rois = myw.plot.getCurvesRoiDockWidget().getRois()
    rois["range1"].setFrom(100.0)
    rois["range1"].setTo  (500.0)
    myw.lineEdit_hfshift .setText("-5.0")
    if ntab==9 :
        qWait( ms=40)   
        myw.pushButton_guess.clicked.emit(True)
        myw.lineEdit_A0.setText("300.0")
        myw.lineEdit_A0.editingFinished.emit()

    qWait( ms=40)    
    myw.pushButton_guess.clicked.emit(True)
    qWait( ms=40)
    rois["range2"].setFrom(550.0)
    rois["range2"].setTo  (594.0)
    myw.pushButton_fit.clicked.emit(True)
    rois["Output"].setFrom(457.0)
    rois["Output"].setTo  (577.0)
    rois["Norm"].setFrom(516.0)
    rois["Norm"].setTo  (553.0)
    qWait( ms=40)

ntab = 6
QTest.mouseClick(  wtest.tabWidget.tabBar(),
                   qt.Qt.LeftButton ,
                   pos=qt.QPoint(
                       wtest.tabWidget.tabBar().tabRect(ntab).x()+wtest.tabWidget.tabBar().tabRect(ntab).width()/2 ,
                       wtest.tabWidget.tabBar().tabRect(ntab).y()+   wtest.tabWidget.tabBar().tabRect(ntab).height()/2 )
)
myw = wtest.tabWidget.currentWidget()
myw.lineEdit_outputPrefix.setText("non_reg_output_gui_raman")
myw.pushButton_saveAnalysis.clicked.emit(True)
qWait( ms=40)


check_results("non_reg_output_gui_raman_lowq.txt","non_reg_output_gui_raman_reference_lowq.txt")
check_results("non_reg_output_gui_raman_mediumq.txt","non_reg_output_gui_raman_reference_mediumq.txt")
check_results("non_reg_output_gui_raman_highq.txt","non_reg_output_gui_raman_reference_highq.txt")
print(" OK " )
qWait( ms=400000)
