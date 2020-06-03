from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: id20_imaging.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import io
from itertools import groupby
from scipy.interpolate import Rbf, RectBivariateSpline
from scipy.optimize import leastsq, fmin
import scipy
import scipy.ndimage
from . import xrs_rois

from .helpers import *
#from xrs_read import rois

from . import xrs_read

import PyMca5.PyMcaIO.specfilewrapper as SpecIO
from six.moves import map
from six.moves import range
from six.moves import zip

class oneD_imaging(xrs_read.read_id20):
    """ **oneD_imaging**
    Class to construct images using the 1D piercing mode.
    """
    def __init__(self,absfilename,energycolumn='sty',monitorcolumn='kapraman', monitor_divider = 1.0,
                     edfName=None,single_image=True,
                     sumto1D = 1, recenterings=None ):
        self.sumto1D = sumto1D
        try:
            self.path     = os.path.split(absfilename)[0] + '/'
            self.filename = os.path.split(absfilename)[1]
        except IOError:
            print('IOError! No such SPEC file, please check SPEC-filename.')
        if not edfName:
            self.edfName = os.path.split(absfilename)[1]
        else:
            self.edfName = edfName
        self.single_image  = single_image
        self.scannumbers   = []
        self.EDF_PREFIXh   = 'edf/h_'
        self.EDF_PREFIXv   = 'edf/v_'
        self.EDF_PREFIX    = 'edf/'
        self.EDF_POSTFIX   = '.edf'
        self.DET_PIXEL_NUMx = 768
        self.DET_PIXEL_NUMy = 512
        self.DET_PIXEL_NUM  = 256
        self.encolumn       = energycolumn.lower()
        self.monicolumn     = monitorcolumn.lower()
        self.monitor_divider = monitor_divider
        self.eloss    = []    # common eloss scale for all analyzers
        self.energy   = []    # common energy scale for all analyzers
        self.signals  = []    # signals for all analyzers
        self.errors   = []    # poisson errors
        self.qvalues  = []    # for all analyzers
        self.groups   = {}  # dictionary of groups (such as 2 'elastic', or 5 'edge1', etc.)
        self.cenom    = []  # list of center of masses of the elastic lines
        self.E0       = []  # energy value, mean value of self.cenom
        self.tth      = []  # list of scattering angles (one for each ROI)
        self.VDtth    = []
        self.VUtth    = []
        self.VBtth    = []
        self.HRtth    = []
        self.HLtth    = []
        self.HBtth    = []
        self.resolution   = []  # list of FWHM of the elastic lines for each analyzer
        self.signals_orig = []  # signals for all analyzers before interpolation
        self.errors_orig  = []  # errors for all analyzers before interpolation
        # TTH offsets from center of V and H modules
        # tth offsets in one direction (horizontal for V-boxes, vertical for H-boxes)
        self.TTH_OFFSETS1   = np.array([5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0, 5.0, 0.0, -5.0])
        # tth offsets in the other direction (horizontal for H-boxes, vertical for V-boxes)
        self.TTH_OFFSETS2   = np.array([-9.71, -9.75, -9.71, -3.24, -3.25, -3.24, 3.24, 3.25, 3.24, 9.71, 9.75, 9.71]) 
        # input
        self.roi_obj = [] # an instance of the roi_object class from the xrs_rois module (new)
        # images
        self.twoDimages = {} # dictionary: one entry per scan, each scan has a list of one image per ROI
        self.motorDict  = {}  # additional dictionary : for every scan a dictionary with values associated to motor names
        self.recenterings = recenterings
                
    def loadscan_2Dimages(self, scannumbers,scantype='sty', isolateSpot = 0):

        print( " SONO IN loadscan_2Dimages monicolumns ", self.monicolumn)
        edfmat =  self.read_just_first_scanimage(scannumbers[0])
        if edfmat.shape[1] != self.DET_PIXEL_NUMx :
            self.DET_PIXEL_NUMy ,self.DET_PIXEL_NUMx  = edfmat.shape

        # check if there are ROIs defined
        if not self.roi_obj:
            print( 'Please define some ROIs first.')
            return
        # make sure scannumbers are iterable (list)
        if not isinstance(scannumbers,list):
            scannums = []
            scannums.append(scannumbers)
        else:
            scannums = scannumbers 
            maxvalue = 0

        for number in scannums:    # go through all scans
            scanname = 'Scan%03d' % number
            data, motors, counters, edfmats = self.readscan(number)
            if self.monicolumn in counters:
                intensities = counters[self.monicolumn]
                np.multiply( edfmats, (self.monitor_divider/np.array(intensities))[:,None,None],edfmats  )
                monitor =  intensities
            else:
                monitor = np.ones( edfmats.shape[0],"f")
                self.monitor_divider = 1.0
            maxvalue = max(maxvalue, edfmats.max())
                                
            position = counters[scantype.lower()]
            self.twoDimages[scanname] = []

            #######
            fn = self.path + self.filename
            spc = SpecIO.Specfile(fn)
            spcscn = spc.select(str(number))

            motodict=dict( list(zip(spcscn.allmotors(),  list(map( spcscn.motorpos, spcscn.allmotors())))) )
            self.motorDict[scanname] = motodict
            #######
                        
            for ii in range(len(self.roi_obj.x_indices)): # go through all ROIs
                # construct the image
                if len( self.roi_obj.x_indices[ii]):
                                        
                    roixinds = self.roi_obj.x_indices[ii]
                    roiyinds = self.roi_obj.y_indices[ii]
                    minx = np.amin(roixinds)
                    miny = np.amin(roiyinds)

                    if self.recenterings is not None:
                            if self.recenterings[ii].shape ==(2,2):
                                    ## to be refined by comparaison : zb,xb is the gol for new baricenter
                                    [[zb, za],[xb,xa]] = self.recenterings[ii]
                                    sx = xb-xa
                                    sz = zb-za
                            else:
                                    ## good shift is known in advance
                                    sz,sx = self.recenterings[ii]

                            if sx != 0.0:
                                    print( " SHIFTO DI " , sx)
                                    isx = int(1000+sx)-1000
                                    dsx = sx-isx
                                    edfmats2use = (1.0-dsx) * np.roll(edfmats , isx , axis=-1)   + dsx *  np.roll(edfmats , isx+1 , axis=-1)
                            if sz != 0.0:
                                    isz = int(1000+sz)-1000
                                    dsz = sz-isz
                                    edfmats2use = (1.0-dsz) * np.roll(edfmats2use , isz , axis=-2)   + dsz *  np.roll(edfmats2use , isz+1 , axis=-1)

                            tmp = np.zeros_like(edfmats2use)
                            tmp[:,roixinds,roiyinds] = edfmats2use[:,roixinds,roiyinds]

                            print( " ENERGIA adesso "  , tmp.sum())
                            print( " ENERGIA prima  "  ,edfmats[:,roixinds,roiyinds] .sum())
                            if self.recenterings[ii].shape ==(2,2):

                                    dimz,dimy,dimx = tmp.shape
                                    BX = (tmp*np.arange(dimx)).sum()/tmp.sum()
                                    print( " original barix for zone ", ii, " : " , xa)
                                    print( " actual / goal  : ", BX, xb)
                                    sx = xb-xa + xb-BX
                                    if sx != 0.0:
                                            print( " SHIFTO DI " , sx)
                                            isx = int(1000+sx)-1000
                                            dsx = sx-isx
                                            edfmats2use = (1.0-dsx) * np.roll(edfmats , isx , axis=-1)   + dsx *  np.roll(edfmats , isx+1 , axis=-1)

                                    tmp = np.zeros_like(edfmats2use)
                                    tmp[:,roixinds,roiyinds] = edfmats2use[:,roixinds,roiyinds]

                                    print( "E ENERGIA adesso "  , tmp.sum())
                                    dimz,dimy,dimx = tmp.shape
                                    BX = (tmp*np.arange(dimx)).sum()/tmp.sum()

                                    print( " original barix for zone ", ii, " : " , xa)
                                    print( " NOW actual / goal  : ", BX, xb)

                                    #############################################################

                                    BY = (tmp*(np.arange(dimy)[:,None])).sum()/tmp.sum()
                                    print( " original bariy for zone ", ii, " : " , za  )
                                    print( " actual / goal  : ", BY, zb)
                                    sz = zb-za + zb-BY
                                    if sz != 0.0:
                                            print( " SHIFTO DI " , sz)
                                            isz = int(1000+sz)-1000
                                            dsz = sz-isz

                                            edfmats2use = (1.0-dsz) * np.roll(edfmats2use , isz , axis=-2)   + dsz *  np.roll(edfmats2use , isz+1 , axis=-2)

                                    shifts = np.array([sz,sx]) 
                                    self.recenterings[ii] = shifts

                                    tmp = np.zeros_like(edfmats2use)
                                    tmp[:,roixinds,roiyinds] = edfmats2use[:,roixinds,roiyinds]

                                    print( "E ENERGIA adesso "  , tmp.sum())
                                    dimz,dimy,dimx = tmp.shape
                                    BX = (tmp*np.arange(dimx)).sum()/tmp.sum()
                                    BY = (tmp*(np.arange(dimy)[:,None])).sum()/tmp.sum()

                                    print( " original barix for zone ", ii, " : " , xa, za)
                                    print( " NOW actual / goal  : ", BX, xb ," and ", BY, zb )

                                                        
                    else:
                        edfmats2use = edfmats
                                        
                    axesrange = [0,roiyinds[-1],position[-1],position[0]]
                    inset = np.zeros([edfmats2use.shape[0],
                                      np.amax(roixinds)+1-minx,
                                      np.amax(roiyinds)+1-miny  ])
                    inset [:,roixinds-minx, roiyinds-miny ] =edfmats2use[:, roixinds , roiyinds ]



                    if self.sumto1D:
                            imageMat  = (np.sum(inset,axis=1))
                    else:
                            if isolateSpot:
                                    imageLines = np.sum(inset,axis=1)
                                    imageLines =imageLines- scipy.ndimage.filters.gaussian_filter( imageLines  ,[0,isolateSpot],mode='constant',cval=0)
                                    poss       = np.argmax(imageLines,axis=1)
                                    for i in range(len(poss)):
                                            inset[i,:, : max(0,poss[i]-isolateSpot)      ]=0
                                            inset[i,:, poss[i]+isolateSpot   :  ]=0
                            imageMat = inset
                    imageInst = LRimage(imageMat,position,
                                                            roiyinds, 
                                                            cornerpos = (minx, miny),
                                                            monitor   = monitor
                                                    )
                    self.twoDimages[scanname].append(imageInst)
                else:
                    imageInst = LRimage([],[],[], monitor = monitor)
                                        
                    self.twoDimages[scanname].append(imageInst)
        return maxvalue



    def save_state_hdf5(self, filename, groupname, comment="", myrank=0, factor = 1.0,  save_also_roi = False):
        import h5py
        h5 = h5py.File(filename,"a")

        if (myrank==0):
            if groupname in h5:
                    print( " levo " , groupname)
                    del h5[groupname]
                    h5.flush()
                    h5.close()
                    h5 = h5py.File(filename,"a")


        h5.require_group(groupname)
        h5group =  h5[groupname]

        if(myrank==0):
                h5group["README"]        = " Some metadata concerning the way data have been extracted \n"
                "comment : contains the containt of input file give to XRS_swissknife\n"
                "path  : the path to the aquisition ditrectory \n"
                "filename : the spec file \n"
                "encolumn : the spec column with the scan energy\n"
                "monicolumn : the column containing the beam reference (not used if not found) \n"
                "edfName    : the prefix of edf files contained in edf directory ( the aquired images)\n"

                h5group["path"]          = self.path
                h5group["filename"]      = self.filename
                h5group["edfName"]       = self.edfName
                h5group["monicolumn"]    = self.monicolumn
                h5group["encolumn"]      = self.encolumn
                h5group["comment"]       = comment

        h5group.require_group("scans")
        h5group = h5group["scans"]

        if save_also_roi and self.roi_obj is not None:
            xrs_rois.write_rois_toh5(h5group  ,  self.roi_obj.red_rois  )

        for scanname,  imageInst in  self.twoDimages.items():
            h5group.require_group(scanname)
            h5group_scan = h5group[scanname]
            if scanname in self.motorDict:
                    h5_md = h5group_scan.require_group("motorDict")
                    for mn, mv in self.motorDict[scanname].items():
                            h5_md[mn] = mv


            
            for ii, ima in enumerate(imageInst):
                print(" ii ima ", ii, ima)
                if ima is not None and ima.matrix !=[]:
                    h5group_scan.require_group( str(ii) )
                    h5group_scan_ii = h5group_scan[   str(ii)      ]
                    h5group_scan_ii["README"]=("barix : the x(fast index) component of the baricenter. May be used to change displacement relative to the roi. \n"
                    "bariy : the y component of the baricenter. May be used to change displacement relative to the roi. \n"
                    "monitor : by what the data have been divided. The denominator. If the monitor name was not given correctly it is a list of ones\n"
                    "monitor_divider :  a renormalisation for the monitor. Can be used to avoid numerically extreme situations where the monitor is huge or too small.\n"
                    "matrix : the stack of image roi data.\n"
                    "cornerpos :  where in the image the roi begins.\n"
                    "xscale : may have several meaning. In general the variable which changes during the scan\n"
                    "yscal  : a variable which changes inside each image of the scan. For example the pixel horizontal position\n")

                    h5group_scan_ii["name"]   = ima.name
                    h5group_scan_ii["matrix"] = ima.matrix*factor
                    h5group_scan_ii["xscale"] = ima.xscale
                    h5group_scan_ii["yscale"] = ima.yscale
                    h5group_scan_ii["tth"]    = ima.tth
                    h5group_scan_ii["cornerpos"]    = ima.cornerpos
                    barix = ( (  ima.matrix*np.arange(ima.matrix.shape[-1] )    ).sum()       )/(  ima.matrix.sum()    )
                    h5group_scan_ii["barix"]  = barix

                    bariy = ( (  ima.matrix*  ( np.arange(ima.matrix.shape[-2])[:,None]  )    ).sum()       )/(  ima.matrix.sum()    )
                    h5group_scan_ii["bariy"]  = bariy
                    h5group_scan_ii["monitor_divider"]    = self.monitor_divider
                    h5group_scan_ii["monitor"]    = ima.monitor


        h5.flush()
        h5.close()

    def load_state_hdf5(self, filename, groupname):
        import h5py

        h5 = h5py.File(filename,"r")

        h5group =  h5[groupname]

        self.path          = h5group["path"          ]          
        self.filename      = h5group["filename"      ]      
        self.edfName       = h5group["edfName"       ]       
        self.monicolumn = h5group["monicolumn" ] 
        self.monitor_divider = h5group["monitor_divider"]  
        self.encolumn  = h5group["encolumn"  ]  

        h5group =  h5[groupname+"/scans"]

        self.twoDimages = {}
        self.motorDict  = {}
        for scanname, group in h5group.items():
            if "motorDict" in group:
                    h5group_scan = group["motorDict"]
                    md = {}
                    for mn, mv in h5group_scan.items():
                            md[mn] = mv
                    self.motorDict[ scanname ] = md

                    
            iis = np.array(list(map( int , list(group.keys()))))
            N = iis.max()+1
            self.twoDimages[scanname]=N*["This items must be filled by load_state_hdf5 " +__file__    ] 
            
            for i  in iis:
                group_i = group[str(i)]
                name    = np.array(group_i["name"])
                matrix  = np.array(group_i["matrix"])
                xscale  = np.array(group_i["xscale"])
                yscale  = np.array(group_i["yscale"])
                tth     = np.array(group_i["tth"])
                cornerpos     = np.array(group_i["cornerpos"])

                ima = image(matrix, xscale, yscale, )
                ima.name = name 
                ima.tth  = tth
                ima.cornerpos=cornerpos
                            
                self.twoDimages[scanname][i] = ima
 
        h5.close()




class imageset:
    """
    class to make SR-images from list of LR-images
    """
    def __init__(self):
        self.list_of_images = []
        self.xshifts  = []
        self.yshifts  = []
        self.shifts   = []
        self.srimage  = []
        self.srxscale = []
        self.sryscale = []
        self.refimagenum = []

    def estimate_xshifts(self,whichimage=None):
        if not whichimage:
            ind = 0 # first image in list_of_images is the reference image
        else:
            ind = whichimage
        origx   = self.list_of_images[ind].xscale
        origy   = self.list_of_images[ind].yscale
        origim  = self.list_of_images[ind].matrix
        xshifts = []        
        for image in self.list_of_images:
            newx  = image.xscale
            newy  = image.yscale
            newim = image.matrix
            xshifts.append(estimate_xshift(origx,origy,origim,newx,newy,newim))
        self.refimagenum = ind
        self.xshifts = xshifts

    def estimate_yshifts(self,whichimage=None):
        if not whichimage:
            ind = 0 # first image in list_of_images is the reference image
        else:
            ind = whichimage
        origx   = self.list_of_images[ind].xscale
        origy   = self.list_of_images[ind].yscale
        origim  = self.list_of_images[ind].matrix
        yshifts = []        
        for image in self.list_of_images:
            newx  = image.xscale
            newy  = image.yscale
            newim = image.matrix
            yshifts.append(estimate_yshift(origx,origy,origim,newx,newy,newim))
        self.refimagenum = ind
        self.yshifts = yshifts

    def estimate_shifts(self,whichimage=None):
        if not whichimage:
            ind = 0 # first image in list_of_images is the reference image
        else:
            ind = whichimage
        origx   = self.list_of_images[ind].xscale
        origy   = self.list_of_images[ind].yscale
        origim  = self.list_of_images[ind].matrix
        shifts = []        
        for image in self.list_of_images:
            newx  = image.xscale
            newy  = image.yscale
            newim = image.matrix
            shifts.append(estimate_shift(origx,origy,origim,newx,newy,newim))
        self.refimagenum = ind
        self.shifts = shifts

    def interpolate_xshift_images(self,scaling,whichimages=None):
        if not whichimages:
            inds = list(range(len(self.list_of_images)))
        elif not isinstance(whichimages,list):
            inds = []
            inds.append(whichimages)
        else:
            inds = whichimages

        newim = np.zeros((len(inds),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling,np.shape(self.list_of_images[inds[0]].matrix)[1]))
        newx  = np.linspace(self.list_of_images[inds[0]].xscale[0]-self.xshifts[inds[0]],self.list_of_images[inds[0]].xscale[-1]-self.xshifts[inds[0]],len(self.list_of_images[inds[0]].xscale)*scaling)
        
        newy  = self.list_of_images[inds[0]].yscale

        for n in range(len(inds)):
            print( self.xshifts[inds[n]])
            oldim = self.list_of_images[inds[n]].matrix
            oldx  = self.list_of_images[inds[n]].xscale-self.xshifts[inds[n]]
            oldy  = self.list_of_images[inds[n]].yscale
            newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

        self.srimage = np.sum(newim,axis=0)
        self.srxscale = newx
        self.sryscale = newy

    def interpolate_yshift_images(self,scaling,whichimages=None):
        if not whichimages:
            inds = list(range(len(self.list_of_images)))
        elif not isinstance(whichimages,list):
            inds = []
            inds.append(whichimages)
        else:
            inds = whichimages

        newim = np.zeros((len(inds),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling,np.shape(self.list_of_images[inds[0]].matrix)[1]))
        newx  = self.list_of_images[0].xscale
        newy  = np.linspace(self.list_of_images[inds[0]].yscale[0]-self.yshifts[inds[0]],self.list_of_images[inds[0]].yscale[-1]-self.yshifts[inds[0]],len(self.list_of_images[inds[0]].yscale)*scaling)
        for n in range(len(inds)):
            oldim = self.list_of_images[inds[n]].matrix
            oldx  = self.list_of_images[inds[n]].xscale
            oldy  = self.list_of_images[inds[n]].yscale-self.yshifts[inds[n]]
            newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

        self.srimage = np.sum(newim,axis=0)
        self.srxscale = newx
        self.sryscale = newy

    def interpolate_shift_images(self,scaling,whichimages=None):
        if not whichimages:
            inds = list(range(len(self.list_of_images)))
        elif not isinstance(whichimages,list):
            inds = []
            inds.append(whichimages)
        else:
            inds = whichimages

        if len(scaling)<2:
            scaling = [scaling, scaling]
        print( inds, self.list_of_images[inds[0]].xscale[0], self.shifts[inds[0]], self.list_of_images[inds[0]].xscale[-1])
        newim = np.zeros((len(self.list_of_images),np.shape(self.list_of_images[inds[0]].matrix)[0]*scaling[0],np.shape(self.list_of_images[inds[0]].matrix)[1]*scaling[1]))
        newx  = np.linspace(self.list_of_images[inds[0]].xscale[0]-self.shifts[inds[0]][0],self.list_of_images[inds[0]].xscale[-1]-self.shifts[inds[0]][0],len(self.list_of_images[inds[0]].xscale)*scaling[0])
        newy  = np.linspace(self.list_of_images[inds[0]].yscale[0]-self.shifts[inds[0]][1],self.list_of_images[inds[0]].yscale[-1]-self.shifts[inds[0]][1],len(self.list_of_images[inds[0]].yscale)*scaling[1])

        for n in range(len(inds)):
            oldim = self.list_of_images[inds[n]].matrix
            oldx  = self.list_of_images[inds[n]].xscale-self.shifts[inds[n]][0]
            oldy  = self.list_of_images[inds[n]].yscale-self.shifts[inds[n]][1]
            newim[n,:,:] = interpolate_image(oldx,oldy,oldim,newx,newy)

        self.srimage = np.sum(newim,axis=0)
        self.srxscale = newx
        self.sryscale = newy

    def plotSR(self):
        X,Y = pylab.meshgrid(self.srxscale,self.sryscale)
        pylab.pcolor(X,Y,np.transpose(self.srimage))
        pylab.show(block=False)

    def plotLR(self,whichimage):
        if not isinstance(whichimage,list):
            inds = []
            inds.append(whichimage)
        else:
            inds = list(whichimage)

        for ind in inds:
            X,Y = pylab.meshgrid(self.list_of_images[ind].xscale,self.list_of_images[ind].yscale)
            pylab.figure(ind)
            pylab.pcolor(X,Y,np.transpose(self.list_of_images[ind].matrix))
            pylab.show(block=False)

    def save(self):
        pass

    def load(self):
        pass

    def loadkimberlite(self,matfilename):
        data_dict   = io.loadmat(matfilename)
        sorted_keys = sorted(data_dict.keys())
        sy    = data_dict['sy'][0]
        allsx = []
        for key in sorted_keys[3:12]:
            allsx.append(data_dict[key][0])

        allmats = []
        for key in sorted_keys[13:]:
            allmats.append(data_dict[key])

        alllengths = []
        for sx in allsx:
            alllengths.append(len(sx))
        
        ind = np.where(alllengths == np.max(alllengths))[0][0]
        # spline everything onto longest sx-scale
        for n in range(len(allmats)):
            print( np.shape(allsx[n]), np.shape(sy), np.shape(allmats[n]))
            ip         = RectBivariateSpline(allsx[n],sy,allmats[n])
            allmats[n] = ip(allsx[ind],sy)
            allsx[n]   = allsx[ind]

        allimages = []
        for n in range(len(allmats)):
            allimages.append(image(allmats[n],allsx[n],sy))
        self.list_of_images = allimages

    def loadhe3070(self,matfilename):
        data_dict  = io.loadmat(matfilename)
        sy    = data_dict['det'][0][0]['sy'][0][0][0]
        allsx = []
        allmats = []
        for n in range(9):
            allsx.append(np.reshape(data_dict['det'][0][n]['sx'][0][0]-data_dict['det'][0][n]['sx'][0][0][0],len(data_dict['det'][0][n]['sx'][0][0],)))
            allmats.append(data_dict['det'][0][n]['img'][0][0])

        alllengths = []
        for sx in allsx:
            alllengths.append(len(sx))

        ind = np.where(alllengths == np.max(alllengths))[0][0]
        for n in range(len(allmats)):
            print( np.shape(allsx[n]), np.shape(sy), np.shape(np.transpose(allmats[n])))
            ip         = RectBivariateSpline(allsx[n],sy,np.transpose(allmats[n]))
            allmats[n] = ip(allsx[ind],sy)
            allsx[n]   = allsx[ind]

        allimages = []
        for n in range(len(allmats)):
            allimages.append(image(allmats[n],allsx[n],sy))
        self.list_of_images = allimages

class image:
    """
    Container class to hold info of a single LR-image to be put togther in a SR-image by the imageset class
    """
    def __init__(self,matrix,xscale,yscale):
        self.name    = []
        self.matrix    = matrix
        self.xscale    = xscale
        self.yscale    = yscale
        self.tth    = []

    def plotimage(self):
        pass

    def shiftx(self):
        pass

    def shifty(self):
        pass

    def save(self):
        pass

    def load(self):
        pass


def interpolate_image(oldx,oldy,oldIM,newx,newy):
    """
    2d interpolation
    """
    interp = RectBivariateSpline(oldx,oldy,oldIM)
    return interp(newx,newy)

def estimate_xshift(x1,y1,im1,x2,y2,im2):
    """
    estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
    """
    funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2+a,y2))**2.0)
    res = leastsq(funct,0.0)
    return res[0]

def estimate_yshift(x1,y1,im1,x2,y2,im2):
    """
    estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
    """
    funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2,y2+a))**2.0)
    res = leastsq(funct,0.0)
    return res[0]

def estimate_shift(x1,y1,im1,x2,y2,im2):
    """
    estimate shift in x-direction only by stepwise shifting im2 by precision and thus minimising the sum of the difference between im1 and im2
    """
    funct = lambda a: np.sum((interpolate_image(x1,y1,im1,x1,y1) - interpolate_image(x2,y2,im2,x2+a[0],y2+a[1]))**2.0)
    res = fmin(funct,[0.0,0.0],disp=0)
    return res


class LRimage:
    """
    container class to hold info of a single LR-image to be put togther in a SR-image by the imageset class
    """
    def __init__(self,matrix,xscale,yscale,  cornerpos = [0,0], monitor = 1.0):
        self.name    = []
        self.matrix    = matrix
        self.xscale    = xscale
        self.yscale    = yscale
        self.tth    = []
        self.cornerpos = cornerpos
        self.monitor = monitor
    def plotimage(self):
        pass

    def shiftx(self):
        pass

    def shifty(self):
        pass

    def save(self):
        pass

    def load(self):
        pass








