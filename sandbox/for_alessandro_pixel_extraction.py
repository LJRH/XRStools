import matplotlib
matplotlib.use("Qt4Agg")

import numpy as np
from  pylab import  *
ion()
from XRStools import xrs_read, roifinder_and_gui, xrs_extraction

lowq = range(24)
lowq.extend(range(36,60))
medq = range(24,36)
highq = range(60,72)

########################################################################
# set the filename and path (this is the O K-edge of liquid water)
path = '/data/id20/inhouse/data/run5_17/run7_ihr/'
lw = xrs_read.Hydra(path)
########################################################################

########################################################################
# define ROIs
roifinder = roifinder_and_gui.roi_finder()
#im4roi = lw.SumDirect([613,617,621,625])
#roifinder.get_zoom_rois(im4roi)
#roifinder.roi_obj.writeH5(path+'rois/lw_OK_raw_bigroi.h5')
roifinder.roi_obj.loadH5(path+'rois/lw_OK_raw_bigroi.h5')
#roifinder.roi_obj.loadH5(path+'rois/lw_OK_MFrefined.h5')
lw.set_roiObj(roifinder.roi_obj)


########################################################################
# load the spectra (something is going wrong during the normalization,
# so I am adding this scaling parameter, I think it happens when loading the
# scans that the monitor variable is multiplied by a wrong factor)
scaling = np.zeros(72)
scaling[lowq] = 4.3
scaling[medq] = 4.3
scaling[highq]= 4.4

method = 'pixel'
lw.get_compensation_factor(611, method=method)
#lw.load_scan([628], method=method, direct=True, scan_type='long', scaling=scaling)
lw.load_scan([611,615,619,623], method=method, direct=True, scan_type='elastic')
lw.load_scan([612,616,620,624], method=method, direct=True, scan_type='ok1')
lw.load_scan([613,617,621,625], method=method, direct=True, scan_type='ok2')
lw.load_scan([614,618,622,626], method=method, direct=True, scan_type='ok3')
lw.get_spectrum_new(method=method, include_elastic=True)

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
lw_ex.save_average_Sqw('water_OK_pixel_lq.dat', emin=00.0, emax=610.0, normrange=[520.,600.])

# O edge med-q
lw_ex.analyzerAverage(medq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[300.0,532.0],[550.0,600.0], weights=[2,1], HFcore_shift=-5.0, guess=[-1.39664220e+03 ,  1.03655696e+03 ,  7.67728511e+02,   7.30355600e+02,  7.93995221e-04,  -4.76580011e-01,  -1.37652621e+03], scaling=1.34)
lw_ex.save_average_Sqw('water_OK_pixel_mq.dat', emin=0.0, emax=610.0, normrange=[520.0,600.0])

# O edge high-q
lw_ex.analyzerAverage(highq,errorweighing=False)
lw_ex.removeCorePearsonAv('O','K',[52.0,532.0],[550.0,600.0],weights=[2,1], guess=[ 3.40779687e+02, 2.57030454e+02, 1.27747244e+03, 4.55875194e-01, -8.59501907e-06, 1.39969288e-02, 2.60071705e+00], HFcore_shift=-5.0,scaling=3.55)
lw_ex.save_average_Sqw('water_OK_pixel_hq.dat', emin=0.0, emax=600.0, normrange=[520.0,600.0])

