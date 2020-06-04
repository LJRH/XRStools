
import numpy as np
import subprocess
from  pylab import  *
### ion()
from XRStools import xrs_read, roifinder_and_gui, xrs_extraction, xrs_rois

from  PyQt4 import Qt, QtCore
app=Qt.QApplication([])

import matplotlib
matplotlib.use("Qt4Agg")

from XRStools  import roiSelectionWidget
from XRStools.roiNmaSelectionGui import roiNmaSelectionWidget
from XRStools import roiSelectionWidget


lowq = list(range(24))
lowq.extend(range(36,60))
medq = list(range(24,36))
highq = list(range(60,72))

########################################################################
# set the filename and path (this is the O K-edge of liquid water)
# path = '/data/id20/inhouse/data/run5_17/run7_ihr/'
# scans = [613,617,621,625]

path = '/data/id20/inhouse/data/run1_18/run4_es689/'
scans = list(range(63,79,2))


########################################################################

########################################################################
# define ROIs

USE_SWISSKNIFE = 1

if 1 :  # create the initial spatial mask

    if USE_SWISSKNIFE:
        s="""
        create_rois :
           expdata :  "{path}/hydra"  # this points to a spec file
           scans   : {scans}                   # a list containing one or more scans 
                                                 # for the elastic part.
                                                 # They will be summed to an image
                                                 # and rois will then be extracted from this image.
           roiaddress : "{target}"  # the target destination for rois
        """
        s=s.format(path=path, scans = scans, target =  "myroi.h5:/datas/ROI"  )
        open("input.yaml","w").write(s)
        comando = 'XRS_swissknife input.yaml' 
        p1 = subprocess.Popen(args=comando.split( " ") ,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)        
        ms,errors=  p1.communicate()
        if p1.returncode:
            print ( " ERRORS ")
            print (  errors.decode("utf-8")  )
                    
            raise Exception("The create_rois instruction was stopped with error message ", errors.decode("utf-8"))

    else:
        lw = xrs_read.Hydra(path)
        im4roi = lw.SumDirect(scans)
        w4r = roiSelectionWidget.mainwindow(layout="2X3-12")
        w4r.showImage(im4roi  )
        ####### TO RELOAD A PREVIOUS NNMA SELECTION inside the widget  ########################
        ## w4r.loadMaskDictFromH5( "spectral_myroi.h5:/datas/ROI" )
        ###################################################################
        ####### TO RELOAD A PREVIOUS NNMA SELECTION inside the widget  ########################
        # w4r.loadMaskDictFromH5( "myroi.h5:/datas/ROI" )
        ###################################################################
        w4r.show()
        app.exec_()
        assert(w4r.isOK)
        w4r.saveMaskDictOnH5( "myroi.h5:/datas/ROI"  )
        w4r = None

if 1:  # create the mask by spectral analisys


    if USE_SWISSKNIFE:

        target = "spectral_myroi.h5:/datas/ROI"
        s="""
        create_spectral_rois :
           expdata :  "{path}/"   # this points to a directory containing a spec file or to the specfile itself 
                                    # You can simply pass the directory. By default the file hydra will be selected
           scans   : {scans}               # a list containing one or more scans 
                                           # They will be summed up  to form an image to display 
                                           # and analysed in the depth direction for the spectral part to refine
                                           # the rois by NNMA decomposition
           spatial_roiaddress  : {spatial_roi}   # a previous ROIs spatial selection                 
           spectral_roiaddress : "{target}"      # the target destination for rois
        """
        s=s.format(path=path, scans = scans,  spatial_roi =  "myroi.h5:/datas/ROI" ,  target= "spectral_myroi.h5:/datas/ROI")
        
        open("input.yaml","w").write(s)
        comando = 'XRS_swissknife input.yaml' 
        p1 = subprocess.Popen(args=comando.split( " ") ,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)        
        ms,errors=  p1.communicate()
        if p1.returncode:
            print ( " ERRORS ")
            print (  errors.decode("utf-8")  )
            raise Exception("The create_spectral_rois instruction was stopped with error message ", errors.decode("utf-8"))
        
        myroi = xrs_rois.load_rois_fromh5_address("spectral_myroi.h5:/datas/ROI")

    else:
        spectral_w4r = roiNmaSelectionWidget.mainwindow()
        sf = path
        fn = "myroi.h5:/datas/ROI"
        ns_s =  scans
        spectral_w4r.load_rois(sf, fn , ns_s )
        spectral_w4r.show()
        ####### TO RELOAD A PREVIOUS NNMA SELECTION inside the widget  ########################
        # spectral_w4r.loadMaskDictFromH5( "spectral_myroi.h5:/datas/ROI" )
        ###################################################################
        app   .exec_()
        assert(spectral_w4r.isOK)
        spectral_w4r.saveMaskDictOnH5( "spectral_myroi.h5:/datas/ROI"  )
    
        myroi=spectral_w4r.getRoiObj()
        spectral_w4r = None
    
else:  # simply reuse an already written nnma mask
    
    myroi = xrs_rois.load_rois_fromh5_address("spectral_myroi.h5:/datas/ROI")



print  ( list(myroi.red_rois.keys()  ) ) 
roinums = [ int(''.join(filter(str.isdigit, str(key) )))           for key in     myroi.red_rois.keys()     ]

roinums.sort()
lowq = [ i for i in range(len(roinums)) if roinums[i] in lowq ]
medq = [ i for i in range(len(roinums)) if roinums[i] in medq ]
highq = [ i for i in range(len(roinums)) if roinums[i] in highq ]

lw = xrs_read.Hydra(path)
lw.set_roiObj(myroi)

########################################################################
# load the spectra (something is going wrong during the normalization,
# so I am adding this scaling parameter, I think it happens when loading the
# scans that the monitor variable is multiplied by a wrong factor)


lw.get_compensation_factor(62, method='sum') 
lw.load_scan([62], method='sum', direct=True, scan_type='elastic')
lw.load_scan(scans, method='sum', direct=True, scan_type='ok0')

lw.get_spectrum_new(method='sum', include_elastic=True)

########################################################################
# setting the scattering angles
lw.get_tths(rvd=28.0, rvu=28.0, rvb=65.0, rhr=30.0, rhl=30.0, rhb=143.0, order=[0, 1, 2, 3, 4, 5])

########################################################################
# starting the subtraction of the background
# 1. average over some crystals
# 2. subtract a Pearson function
# 3. write the spectra into a txt file

lw_ex = xrs_extraction.edge_extraction(lw,['H2O'],[1.0],{'O':['K']})

# O edge low-q
lw_ex.analyzerAverage(lowq, errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[250.0,532.0],[550.0,600.0],weights=[2,1],HFcore_shift=-5.0, guess= [-1.07743447e+03, 8.42895443e+02, 4.99035465e+01, 3193e+01, -3.80090286e-07, 2.73774370e-03, 5.11920401e+03],scaling=1.32)
lw_ex.save_average_Sqw('water_OK_nnmf_lq.dat', emin=00.0, emax=610.0, normrange=[520.,600.])

# O edge med-q
lw_ex.analyzerAverage(medq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[300.0,532.0],[550.0,600.0], weights=[2,1], HFcore_shift=-5.0, guess=[-1.39664220e+03 ,  1.03655696e+03 ,  7.67728511e+02,   7.30355600e+02,  7.93995221e-04,  -4.76580011e-01,  -1.37652621e+03], scaling=1.34)
lw_ex.save_average_Sqw('water_OK_nnmf_mq.dat', emin=0.0, emax=610.0, normrange=[520.0,600.0])

# O edge high-q
lw_ex.analyzerAverage(highq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[52.0,532.0],[550.0,600.0],weights=[2,1], guess=[ 3.40779687e+02, 2.57030454e+02, 1.27747244e+03, 4.55875194e-01, -8.59501907e-06, 1.39969288e-02, 2.60071705e+00], HFcore_shift=-5.0,scaling=3.55)
lw_ex.save_average_Sqw('water_OK_nnmf_hq.dat', emin=0.0, emax=600.0, normrange=[520.0,600.0])

