from __future__ import print_function
#####
# this is important with old version of matplotlib  to avoid
# that matplotlib import pyqt4 while silx imports pyqt5
#
try:
    import PyQt4.QtCore
except:
    pass

import os
import numpy as np
# from  pylab import  *
# ion()
from XRStools import xrs_read, roifinder_and_gui, xrs_extraction, roiSelectionWidget



#####################################################################
# for programmatic testing
from silx.gui.utils.testutils import QTest, TestCaseQt
from silx.gui import qt


AUTOMATIC_TEST = 1


if AUTOMATIC_TEST :
    _qapp = qt.QApplication.instance() or qt.QApplication([])


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

lowq = list(range(24))
lowq.extend(list(range(36,60)))
medq = list(range(24,36))
highq = list(range(60,72))

########################################################################

data_path = '/data/id20/inhouse/data/run5_17/run7_ihr/'
save_path = './'
rois_path = './'

########################################################################
# liquid water 
########################################################################
lw = xrs_read.Hydra(data_path)

# ROI definition
im4roi = lw.SumDirect([611]) # elastic on the sample

##################################################################################################################
# manageQApp is set to false because we are programmatically generating events by using the QTest class from silx
# which manages qapp by itself.
# In case of a normal usage, dont pass the manageQApp  argument

if AUTOMATIC_TEST  :
    delay=1000
    roiSelectionWidget.SKIP_WARNING = True
    w4r = roiSelectionWidget.launch4MatPlotLib(im4roi=im4roi, layout = "2X3-12", manageQApp = False )

    for act in w4r.menuActions .actions():
        if  str(act.objectName()) == "actionGlobalSpotDetection":
            act.trigger()
    QTest.qWait(delay)
    QTest.mouseClick(  w4r.globalSpotDectionWidget.detectionButton ,qt.Qt.LeftButton, pos= qt.QPoint(3,3))
    QTest.qWait(delay)


    # # Select the second Tab of the MainWindow
    # QTest.mouseClick( w4r.viewsTab.widget(3)   , qt.Qt.LeftButton, qt.Qt.NoModifier )
#  w4r.viewsTab.tabBar().tabBarClicked.emit(3)


    QTest.mouseClick(  w4r.viewsTab.tabBar(),
                       qt.Qt.LeftButton ,
                       pos=qt.QPoint(
                           w4r.viewsTab.tabBar().tabRect(3).x()+w4r.viewsTab.tabBar().tabRect(3).width()/2 ,
                           w4r.viewsTab.tabBar().tabRect(3).y()+   w4r.viewsTab.tabBar().tabRect(3).height()/2 )
                   )
    QTest.qWait(delay)
    #  alternativa w4r.viewsTab.setCurrentIndex(3)
    QTest.qWait(delay)

    myw =  w4r.mws[3]
    myw.show()
    myw.graph.setFocus()

    pos1 = myw.graph.dataToPixel( 115,250  )
    pos2 = myw.graph.dataToPixel( 115+40,250-5  )

    QTest.mouseMove( myw.graph )
    QTest.qWait(delay)
    QTest.mouseMove( myw.graph ,  pos= qt.QPoint(  pos1[0]  ,pos1[1]),  delay=-1)
    QTest.qWait(delay)
    myw.graph.onMousePress( pos1[0], pos1[1], "left" )
    # QTest.mousePress( myw.graph ,  qt.Qt.LeftButton, pos= qt.QPoint(  pos1[0]  ,pos1[1]), delay=10)
    QTest.qWait(delay)
    QTest.mouseMove( myw.graph   ,  pos= qt.QPoint(  pos2[0] ,pos2[1]), delay=10)
    QTest.qWait(delay)
    # QTest.mouseRelease( myw.graph  ,  qt.Qt.LeftButton, pos= qt.QPoint(  pos2[0]  ,pos2[1]), delay=10)
    myw.graph.onMouseRelease( pos2[0], pos2[1], "left" )
    QTest.qWait(delay)
    QTest.mouseClick(  w4r.globalSpotDectionWidget.relabeliseButton ,qt.Qt.LeftButton, pos= qt.QPoint(3,3))
    QTest.qWait(delay)    
    QTest.mouseClick(  w4r.globalSpotDectionWidget.annotateButton ,qt.Qt.LeftButton, pos= qt.QPoint(3,3))
    QTest.qWait(delay)    
    for act in w4r.menuActions .actions():
        if  str(act.objectName()) == "actionRegistration":
            act.trigger()

    QTest.qWait(delay)    
    for act in w4r.menuFIle .actions():
        if  str(act.objectName()) == "actionConfirm_And_Exit":
            act.trigger()

    QTest.qWait(delay)    


else:
    w4r = roiSelectionWidget.launch4MatPlotLib(im4roi=im4roi, layout = "2X3-12", manageQApp = True )

roi = w4r.getRoiObj()
roi.writeH5(os.path.join(rois_path,'ROI_widget_roi.H5'))

##########################################################################
# load data, use sum-algorithm
##########################################################################
lw.set_roiObj(roi) 
lw.get_compensation_factor(611, method='sum')
lw.load_scan([611], method='sum', direct=True, scan_type='elastic')
lw.load_scan([612], method='sum', direct=True, scan_type='ok1')
lw.load_scan([613], method='sum', direct=True, scan_type='ok2')
lw.load_scan([614], method='sum', direct=True, scan_type='ok3')
lw.get_spectrum_new(method='sum', include_elastic=True)


lw.get_tths(rvd=28.0, rvu=28.0, rvb=65.0, rhr=30.0, rhl=30.0, rhb=143.0, order=[0, 1, 2, 3, 4, 5])
lw_ex = xrs_extraction.edge_extraction(lw,['H2O'],[1.0],{'O':['K']})

# O edge low-q
lw_ex.analyzerAverage(lowq, errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[250.0,534.0],[570.0,600.0],weights=[2,1],HFcore_shift=-5.0, guess= [-1.07743447e+03, 8.42895443e+02, 4.99035465e+01, 3193e+01, -3.80090286e-07, 2.73774370e-03, 5.11920401e+03],scaling=1.2)
lw_ex.save_average_Sqw(os.path.join(save_path,'h2o_sum_lq.dat'), emin=00.0, emax=610.0, normrange=[520.,600.])


def check_results(fn_a, fn_b):
    if os.path.exists(  fn_b ):
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
            assert(    (abs( a[:,1] -  ref)[np.less( c0, a[:,0])* np.less(  a[:,0], c1)*   np.less(  a[:,0], 100) ]   ).sum()  <2000.0     )
            assert(    (abs( a[:,1] -  ref)[np.less( c0, a[:,0])* np.less(  a[:,0], c1)*   np.less( 100,  a[:,0]) ]   ).sum()  <10.0     )
        except:
            print (   a[:,1] -  ref  )
            print(   (abs( a[:,1] -  ref)[np.less( c0, a[:,0])* np.less(  a[:,0], c1)  ]   ).sum()   )
            raise

check_results( os.path.join(save_path,'h2o_sum_lq.dat') , "h2o_sum_lq_ref.dat"   )


# O edge med-q
lw_ex.analyzerAverage(medq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[300.0,534.0],[570.0,600.0], weights=[2,1], HFcore_shift=-5.0, guess=[-1.39664220e+03 ,  1.03655696e+03 ,  7.67728511e+02,   7.30355600e+02,  7.93995221e-04,  -4.76580011e-01,  -1.37652621e+03], scaling=1.2)
lw_ex.save_average_Sqw(save_path+'/h2o_sum_mq.dat', emin=0.0, emax=610.0, normrange=[520.0,600.0])


check_results( os.path.join(save_path,'h2o_sum_mq.dat') , "h2o_sum_mq_ref.dat"   )


# O edge high-q
lw_ex.analyzerAverage(highq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[52.0,534.0],[570.0,600.0],weights=[2,1], guess=[ 3.40779687e+02, 2.57030454e+02, 1.27747244e+03, 4.55875194e-01, -8.59501907e-06, 1.39969288e-02, 2.60071705e+00], HFcore_shift=-5.0,scaling=3.55)
lw_ex.save_average_Sqw(save_path+'/h2o_sum_hq.dat', emin=0.0, emax=600.0, normrange=[520.0,600.0])



check_results( os.path.join(save_path,'h2o_sum_hq.dat') , "h2o_sum_hq_ref.dat"   )



########################################################################
# liquid water, use pixel-algorithm
########################################################################
lw = xrs_read.Hydra(data_path)

# # ROI definition
# im4roi = lw.SumDirect([611]) # elastic on the sample 
# w4r = roiSelectionWidget.launch4MatPlotLib(im4roi=im4roi, layout = "2X3-12")
# roi = w4r.getRoiObj()
# roi.writeH5('//data/id20/inhouse/data/run5_17/run7_ihr/rois_nr/ROI_widget_roi.H5')

# load data
lw.set_roiObj(roi) 
lw.get_compensation_factor(611, method='pixel')
lw.load_scan([611], method='pixel', direct=True, scan_type='elastic')
lw.load_scan([612], method='pixel', direct=True, scan_type='ok1')
lw.load_scan([613], method='pixel', direct=True, scan_type='ok2')
lw.load_scan([614], method='pixel', direct=True, scan_type='ok3')
lw.get_spectrum_new(method='pixel', include_elastic=True)


lw.get_tths(rvd=28.0, rvu=28.0, rvb=65.0, rhr=30.0, rhl=30.0, rhb=143.0, order=[0, 1, 2, 3, 4, 5])
lw_ex = xrs_extraction.edge_extraction(lw,['H2O'],[1.0],{'O':['K']})

# O edge low-q
lw_ex.analyzerAverage(lowq, errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[250.0,534.0],[570.0,600.0],weights=[2,1],HFcore_shift=-5.0, guess= [-1.07743447e+03, 8.42895443e+02, 4.99035465e+01, 3193e+01, -3.80090286e-07, 2.73774370e-03, 5.11920401e+03],scaling=1.2)
lw_ex.save_average_Sqw(save_path+'/h2o_pixel_lq.dat', emin=00.0, emax=610.0, normrange=[520.,600.])

check_results( os.path.join(save_path,'h2o_pixel_lq.dat') , "h2o_pixel_lq_ref.dat"   )



# O edge med-q
lw_ex.analyzerAverage(medq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[300.0,534.0],[570.0,600.0], weights=[2,1], HFcore_shift=-5.0, guess=[-1.39664220e+03 ,  1.03655696e+03 ,  7.67728511e+02,   7.30355600e+02,  7.93995221e-04,  -4.76580011e-01,  -1.37652621e+03], scaling=1.2)
lw_ex.save_average_Sqw(save_path+'/h2o_pixel_mq.dat', emin=0.0, emax=610.0, normrange=[520.0,600.0])

check_results( os.path.join(save_path,'h2o_pixel_mq.dat') , "h2o_pixel_mq_ref.dat"   )


# O edge high-q
lw_ex.analyzerAverage(highq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[52.0,534.0],[570.0,600.0],weights=[2,1], guess=[ 3.40779687e+02, 2.57030454e+02, 1.27747244e+03, 4.55875194e-01, -8.59501907e-06, 1.39969288e-02, 2.60071705e+00], HFcore_shift=-5.0,scaling=3.55)
lw_ex.save_average_Sqw(save_path+'/h2o_pixel_hq.dat', emin=0.0, emax=600.0, normrange=[520.0,600.0])

check_results( os.path.join(save_path,'h2o_pixel_hq.dat') , "h2o_pixel_hq_ref.dat"   )
