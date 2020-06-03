from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import six
from six.moves import range
#!/usr/bin/python
# Filename: xrs_rois.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

import numpy as np
import copy
import h5py
import os
import matplotlib.pyplot as plt

# commented the *import because otherwise sphinx documents all the symbol of other packages 
# from xrs_utilities import *
# from math_functions import *
from . import xrs_utilities
from . import math_functions
###########################################


from matplotlib.widgets import Cursor, Button
from scipy.ndimage import measurements
from scipy import signal


def h5_assign_force(h5group, name, item):
    if name in h5group:
        del h5group[name]
    h5group[name] = item

class container:
    """
    Random container class to hold values
    """
    def __init__(self):
        pass




sl={}
sl[0  ] = slice(0  ,256)
sl[256] = slice(256,512)
sl[512] = slice(512,768)
sl[3*256] = slice(3*256    ,3*256 + 256)
sl[4*256] = slice(4*256    ,4*256 + 256)


V147 = [1,4,7,10,2,5,8,11,3,6,9,12]
V1296 = [12,9,6,3,11,8,5,2,10,7,4,1]
V1074 = [10,7,4,1,11,8,5,2,12,9,6,3     ]
V369  = [3,6,9,12,2,5,8,11,1,4,7,10             ]
V1234 =  [1,2,3,4]


def order_generator_ascending(a,b,c,d):
    return [ a,b,c,d,
             a+1,b+1,c+1,d+1,
             a+2,b+2,c+2,d+2]
def order_generator_descending(a,b,c,d):
    return [ a,b,c,d,
             a-1,b-1,c-1,d-1,
             a-2,b-2,c-2,d-2]

OG_inc = order_generator_ascending
OG_dec = order_generator_descending

geo_informations = {(256,768,None): { "DET_PIXEL_NUM":256, "geo":[256,768], "nofrois":36,
                                 "subnames": ["RD" ,"LU","B"],
                                 "subgeos" : [(sl[0] ,sl[0]),
                                              (sl[0],sl[256]),
                                              (sl[0],sl[512])]   ,
                                 "analyser_nIDs": {"LU":{"3x4":OG_inc(10,7,4,1),"Vertical": OG_inc(1,4,7,10)},  
                                                   "RD":{"3x4":OG_inc(1,4,7,10),"Vertical": OG_inc(1,4,7,10)},
                                                   "B": {"3x4":OG_inc(1,4,7,10),"Vertical": OG_inc(1,4,7,10)}  
                                                   }
                                 },
                    (512,768,None) : { "DET_PIXEL_NUM":256, "geo":[512,768],"nofrois":72,
                                  "subnames":["VD" , "VU","VB","HR" ,"HL","HB", ] ,
                                  "subgeos"  :[(sl[0],sl[0]     ),
                                               (sl[0],sl[256]   ),
                                               (sl[0],sl[512]   ),
                                               (sl[256],sl[0]   ),
                                               (sl[256],sl[256] ),
                                               (sl[256],sl[512] )],
                                  "analyser_nIDs": {"VD":{"3x4": OG_dec(12,9,6,3)  ,"Vertical":   OG_dec(12,9,6,3)  },
                                                    "VU":{"3x4": OG_dec(3,6,9,12)  ,"Vertical":   OG_dec(12,9,6,3)  },
                                                    "VB":{"3x4": OG_dec(3,6,9,12)  ,"Vertical":   OG_dec(12,9,6,3) },
                                                    "HR":{"3x4": OG_inc(1,4,7,10)  ,"Vertical":   OG_inc(1,4,7,10) },
                                                    "HL":{"3x4": OG_inc(10,7,4,1)  ,"Vertical":   OG_inc(1,4,7,10) },
                                                    "HB":{"3x4": OG_inc(10,7,4,1)  ,"Vertical":   OG_inc(1,4,7,10) },
                                                    }
                                  },
                    (256,256,None):{"DET_PIXEL_NUM":256, "geo":[256,256],"nofrois":1,"subnames":["DETECTOR"],"subgeos" : [(sl[0] ,sl[0])],
                               "analyser_nIDs": {"DETECTOR":{"3x4":V147,"Vertical": V147}}
                                },
                    (256,256,"1X1-4"):{"DET_PIXEL_NUM":256, "geo":[256,256],"nofrois":4,"subnames":["DETECTOR"],"subgeos" : [(sl[0] ,sl[0])],
                               "analyser_nIDs": {"DETECTOR":{ "Vertical": V1234} }
                                }
                }


geo_informations[(256,768,"1X3-12")]=geo_informations[(256,768,None)]
geo_informations[(512,768,"2X3-12")]=geo_informations[(512,768,None)]
geo_informations[(256,256,"1X1-12")]=geo_informations[(256,256,None)]


for NBCU in [255,256]:
        geo_informations[(NBCU,5*259,"1X5-1")] = { "DET_PIXEL_NUM":NBCU, "geo":[NBCU,5*259],"nofrois":5,
                                                   "subnames":["A" ,"B","C","D" , "E" ] ,
                                                   "subgeos"  :[(slice(0,NBCU),slice(259*0,259*1)  ),
                                                                (slice(0,NBCU),slice(259*1,259*2)  ),
                                                                (slice(0,NBCU),slice(259*2,259*3)  ),
                                                                (slice(0,NBCU),slice(259*3,259*4)   ),
                                                                (slice(0,NBCU),slice(259*4,259*5)   ),
                                                        ],
                                                   "analyser_nIDs": {"A":{"Vertical":[1]},
                                                                     "B":{"Vertical":[1]},
                                                                     "C":{"Vertical":[1]},
                                                                     "D":{"Vertical":[1]},
                                                                     "E":{"Vertical":[1]},
                                                             }
                                           }
        
        
        geo_informations[(NBCU,5*259+1,"1X5-1")]  = geo_informations[(NBCU,5*259,"1X5-1")] 

def get_geo_informations(shape):
    return geo_informations[shape]


class roi_object:
    """
    Container class to hold all relevant information about given ROIs.
    """
    def __init__(self):
        self.roi_matrix     = np.array([]) # single matrix of zeros, ones, twos, ... , n's (where n is the number of ROIs defined)
        self.red_rois       = {}           # dictionary, one entry for each ROI, each ROI has an origin and a rectangular box of ones and zeros defining the ROI
        self.indices        = [] # list of list of tuples (one list of tuples for each ROI)
        self.number_of_rois = 0  # number of ROIs defined
        self.kind           = [] # keyword (e.g. 'zoom', 'line', 'auto', etc.), certain features (esp. in imaging) are only available for certain kinds of ROIs
        self.x_indices      = [] # list of numpy arrays of x-indices (for each ROI)
        self.y_indices      = [] # list of numpy arrays of y-indices (for each ROI)
        self.masks          = [] # 3D numpy array with slices of zeros and ones (same size as detector image) for each roi
        self.input_image    = [] # 2D imput image that was used to define the ROIs
        
    def load_rois_fromMasksDict(self, masksDict, newshape=None, kind="zoom"):
        self.kind=kind
        self.red_rois = masksDict
        if newshape is not None:
            self.roi_matrix = np.zeros(newshape)
        self.roi_matrix = convert_redmatrix_to_matrix( masksDict,self.roi_matrix , offsetX=0, offsetY=0)

        self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)
        self.indices        = convert_matrix_rois_to_inds(self.roi_matrix)
        self.number_of_rois = int(np.amax(self.roi_matrix))
        self.x_indices      = convert_inds_to_xinds(self.indices)
        self.y_indices      = convert_inds_to_yinds(self.indices)

    def writeH5(self,fname):
        """ **writeH5**
        Creates an HDF5 file and writes the ROIs into it.

        Args:
          * fname (str) : Full path and filename for the HDF5 file to be created.

        """
        if self.indices:
            # check if file already exists
            if os.path.isfile(fname):
                os.remove(fname)
            f = h5py.File(fname, "w")
            f.require_group("rois_definition")
            f["rois_definition"]["image"] = self.input_image
            f["rois_definition"].require_group("rois_dict")
            write_rois_toh5(f["rois_definition"]["rois_dict"],self.red_rois)
            f.close()
        else:
            print('There are no ROIs to save.')

    def loadH5(self,fname):
        """ **loadH5**
        Loads ROIs from an HDF5 file written by the self.writeH5() method.

        Args:
          * fname (str) : Full path and filename for the HDF5 file to be read.
        """
        f   = h5py.File(fname, "r")
        self.input_image = f["rois_definition"]["image"][:]
        self.red_rois    = {}
        shape = load_rois_fromh5(f,self.red_rois)

        if 1:
            self.load_rois_fromMasksDict(self.red_rois ,  newshape = shape, kind="zoom")
        else:
    
            self.roi_matrix     = convert_redmatrix_to_matrix( self.red_rois, np.zeros_like(self.input_image), offsetX=0, offsetY=0)
            self.indices        = convert_matrix_rois_to_inds(self.roi_matrix)
            self.number_of_rois = int(np.amax(self.roi_matrix))
            self.x_indices      = convert_inds_to_xinds(self.indices)
            self.y_indices      = convert_inds_to_yinds(self.indices)
            self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)

    def load_shadok_h5(self, fname, group_name1, group_name2='ROI_AS_SELECTED' ):
        f   = h5py.File(fname, "r")
        self.input_image = f[group_name1][group_name2]["rois_definition"]["image"][:]
        self.red_rois    = {}
        load_rois_fromh5(f[group_name1][group_name2],self.red_rois)
        self.roi_matrix     = convert_redmatrix_to_matrix( self.red_rois, np.zeros_like(self.input_image), offsetX=0, offsetY=0)
        self.indices        = convert_matrix_rois_to_inds(self.roi_matrix)
        self.number_of_rois = int(np.amax(self.roi_matrix))
        self.x_indices      = convert_inds_to_xinds(self.indices)
        self.y_indices      = convert_inds_to_yinds(self.indices)
        self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)
                
    def append(self,roi_object):
        orig_length = len(self.red_rois)
        self.indices.extend(roi_object.indices) # list of list of tuples (one list of tuples for each ROI)
        self.number_of_rois =+ roi_object.number_of_rois  # number of ROIs defined
        self.x_indices.extend(roi_object.x_indices) # list of numpy arrays of x-indices (for each ROI)
        self.y_indices.extend(roi_object.y_indices) # list of numpy arrays of y-indices (for each ROI)
        #self.masks          = [] # 3D numpy array with slices of zeros and ones (same size as detector image) for each roi
        #self.input_image    += [] # 2D imput image that was used to define the ROIs
        roi_object.roi_matrix[roi_object.roi_matrix>0] += orig_length
        self.roi_matrix     += roi_object.roi_matrix  # single matrix of zeros, ones, twos, ... , n's (where n is the number of ROIs defined)
        
        for ii,key in enumerate(sorted(roi_object.red_rois)):
            new_key = 'ROI%02d'%(ii+orig_length)
            self.red_rois[new_key] = roi_object.red_rois[key]
            self.red_rois[new_key][1][ self.red_rois[new_key][1]>0 ] += orig_length

    def get_number_of_rois(self):
        return self.number_of_rois

    def get_indices(self):
        return self.indices

    def get_x_indices(self):
        return self.x_indices

    def get_y_indices(self):
        return self.y_indices

    def get_bounding_boxes(self):
        return self.bounding_boxes

    def get_masks(self):
        return self.masks

    def get_copy(self):
        """
        **get_copy**
        Returns a deep copy of self.
        """
        return copy.deepcopy(self)

    def strip_rois(self):
        """ **strip_rois**
        Strips extra zeros out of ROIs.

        """
        pass

    def delete_empty_rois(self):
        """ **delete_empty_rois**
        """
        pass
    
    def shift_rois(self,shiftVal,direction='horiz',whichroi=None):
        """
        **shift_rois**
        Displaces the defined ROIs by the provided value.

        Args
          * shiftVal : int
            Value by which the ROIs should be shifted.
          * direction : string
            Description of which direction to shit by.
          * whichroi : sequence
            Sequence (iterable) for which ROIs should be shifted.
        """
        the_indices = []

        if not whichroi:
            inds = list(range(len(self.indices)))
        else:
            inds = whichroi

        if direction == 'vert':
            for roi in self.indices:
                oneroi = []
                for pixel in roi:
                    oneroi.append( (pixel[0]+shiftVal,pixel[1]) )
                the_indices.append(oneroi)

        if direction == 'horiz':
            for roi in self.indices:
                oneroi = []
                for pixel in roi:
                    oneroi.append( (pixel[0], pixel[1]+shiftVal) )
                the_indices.append(oneroi)

        self.indices = the_indices

        self.roi_matrix     = convert_inds_to_matrix(self.indices,self.input_image.shape)
        self.red_rois       = convert_matrix_to_redmatrix(self.roi_matrix)
        self.x_indices      = convert_inds_to_xinds(self.indices)
        self.y_indices      = convert_inds_to_yinds(self.indices)
        self.masks          = convert_roi_matrix_to_masks(self.roi_matrix)

    def show_rois(self,colormap='jet',interpolation='nearest'):
        """ **show_rois**
        Creates a figure with the defined ROIs as numbered boxes on it.
        """
        roi_matrix = self.roi_matrix

        # check if there are ROIs defined
        if not np.any(roi_matrix):
            print( 'Please select some rois first.')

        # make a figure
        else:
            plt.ion()
            plt.cla()
            figure_obj = plt.imshow(roi_matrix,interpolation=interpolation)
            figure_obj.set_cmap(colormap)
            plt.xlabel('x-axis [pixel]')
            plt.ylabel('y-axis [pixel]')
            plt.show()
            # add a label with the number to the center of each ROI
            for ii in range(int(np.amax(roi_matrix))):
                # find center of the ROI and write the label
                inds    = np.where(roi_matrix[:,:] == ii+1)
                xcenter = np.mean(inds[1])
                ycenter = np.mean(inds[0])
                string  = '%02d' % (ii+1)
                plt.text(xcenter,ycenter,string)
                plt.draw()

def convert_redmatrix_to_matrix( masksDict,mask, offsetX=0, offsetY=0):
    for key, (pos,M)  in six.iteritems(masksDict):
        num=int("".join([c for c in key if c.isdigit()]))
        S = M.shape

        inset =    (slice(offsetY+pos[0]  , offsetY+pos[0]+S[0]   ), slice(  offsetX+pos[1]  , offsetX+pos[1]+S[1] ) )
        mask[  inset   ][M>0] =  num+1

    return mask

def convert_redmatrix_to_matrix_my( masksDict, mask, offsetX=0, offsetY=0):
    for key in masksDict:
        num = int("".join([c for c in key if c.isdigit()]))
        origin = masksDict[key][0]
        data   = masksDict[key][1]
        for xx in range(len(data[:,0])):
            for yy in range(len(data[0,:])):
                if data[xx,yy] >= 1.0:
                    mask[origin[0]+xx,origin[1]+yy] = num+1
    return mask


def convert_inds_to_matrix(ind_rois,image_shape):
    """
    Converts a ROI defined by a list of lists of tuples into a ROI
    that is defined by an array containing zeros, ones, twos, ..., n's, 
    where n is the number of ROIs.
    ind_rois    = list of lists with pairs of pixel indices
    image_shape = touple defining the shape of the matrix for which the ROIs are valid 
    """
    roi_matrix = np.zeros(image_shape)
    counter = 1
    for pixel in ind_rois:
        for xyind in pixel:
            roi_matrix[int(xyind[0]),int(xyind[1])] = counter
        counter += 1
    return roi_matrix

def convert_matrix_to_redmatrix(matrix_rois, labelformat= 'ROI%02d'):
    """
    Converts a ROI defined by an array containing zeros, ones, twos, ..., n's, 
    where n is the number of ROIs, into a dictionary with keys 'ROI00',
    'ROI01', ..., 'ROInn'. Each entry of the dictionary is a list containing a
    tuple with the origin and the reduced ROI.
    matrix_roi = numpy array
    """
    redmatrix = {}
    for ind in range(int(np.amax(matrix_rois))):
        try:
            indices   = np.where(matrix_rois == ind+1)
            origin    = (np.amin(indices[0]), np.amin(indices[1]))
            bound_box = matrix_rois[np.amin(indices[0]):np.amax(indices[0])+1,np.amin(indices[1]):np.amax(indices[1])+1]
            thekey    =labelformat % ind
            redmatrix[thekey] = [origin,bound_box]
        except: # handle empty ROIs
            indices   = np.where(matrix_rois == ind+1)
            origin    = (0, 0)
            bound_box = matrix_rois[0:0,0:0]
            thekey    =labelformat % ind
            redmatrix[thekey] = [origin,bound_box]
    return redmatrix

def convert_inds_to_xinds(roi_inds):
    """
    Converts ROIs defined in lists of lists of x-y-coordinate tuples into
    a list of x-coordinates only.
    """
    xind_rois = []
    for roi in roi_inds:
        xinds = []
        for pair in roi:
            xinds.append(pair[0])
        xind_rois.append(xinds)
    return xind_rois

def convert_inds_to_yinds(roi_inds):
    """
    Converts ROIs defined in lists of lists of x-y-coordinate tuples into
    a list of y-coordinates only.
    """
    yind_rois = []
    for roi in roi_inds:
        yinds = []
        for pair in roi:
            yinds.append(pair[1])
        yind_rois.append(yinds)
    return yind_rois

def convert_roi_matrix_to_masks(roi_matrix):
    """
    Converts a 2D ROI matrix with zeros, ones, twos, ..., n's (where n is the number of ROIs) to
    a 3D matrix with one slice of zeros and ones per ROI.
    """
    # get the shape

    roi_masks = np.zeros((int(np.amax(roi_matrix)),roi_matrix.shape[0],roi_matrix.shape[1]))

    for ii in range(int(np.amax(roi_matrix))):
        inds = np.where(roi_matrix[:,:] == ii+1)
        for jj in range(len(inds[0])):
            roi_masks[ii,inds[0][jj],inds[1][jj]] = ii+1
    return roi_masks

def convert_matrix_rois_to_inds(roi_matrix):
    """
    Converts a 2D ROI matrix with zeros, ones, twos, ..., n's (where n is the number of ROIs) to
    a list of lists each of which has tuples with coordinates for each pixel in each roi.
    """
    rois = []
    number_of_rois = int(np.amax(roi_matrix))
    for ii in range(int(number_of_rois)):
        inds = np.where(roi_matrix[:,:] == ii+1)
        oneroi = []
        for i in range(len(inds[0])):
            oneroi.append( (inds[0][i],inds[1][i]) )        
        rois.append(oneroi)
    return rois

def break_down_det_image(image,pixel_num):
    """
    Desomposes a Detector image into subimages. Returns a 3D matrix. 
    """
    # check that there are integer multiples of pixel_num in the big image
    if not image.shape[0] % pixel_num == 0 or not image.shape[1] % pixel_num == 0:
        print( 'There must be an integer number of \'pixel_num\' in the large image.')
        return

    detnum_row  = image.shape[0]//pixel_num
    detnum_col  = image.shape[1]//pixel_num
    num_of_dets = detnum_row * detnum_col
    det_images  = np.zeros((num_of_dets,pixel_num,pixel_num))
    x_ranges = []
    for ii in range(detnum_col):
        x_ranges.append((ii*pixel_num, (ii+1)*pixel_num))

    y_ranges = []
    for ii in range(detnum_row):
        y_ranges.append((ii*pixel_num, (ii+1)*pixel_num))

    counter = 0
    offsets = []
    for i in x_ranges:
        for j in y_ranges:
            det_images[counter,:,:]=image[j[0]:j[1],i[0]:i[1]]
            offsets.append((j[0],i[0]))
            counter += 1

    return det_images, offsets

def shift_roi_indices(indices,shift):
    """
    Applies a given shift (xshift,yshift) to given indices. \
    indices = list of (x,y)-tuples
    shift   = (xshift,yshift) tuple
    """
    for ind in indices:
        ind[0] += shift[0]
        ind[1] += shift[1]
    return indices

def merge_roi_objects_by_matrix(list_of_roi_objects,large_image_shape,offsets,pixel_num):
    """
    Merges several roi_objects into one big one using the roi_matrices.
    """
    # prepare a big roi matrix
    roi_matrix = np.zeros(large_image_shape)
    counter = 0
    for roi_obj in list_of_roi_objects:
        max_number = int(np.amax(roi_matrix))
        single_matrix = roi_obj.roi_obj.roi_matrix
        inds = single_matrix != 0
        single_matrix[inds] += max_number
        roi_matrix[offsets[counter][0]:offsets[counter][0]+pixel_num,offsets[counter][1]:offsets[counter][1]+pixel_num] = single_matrix
        #roi_obj.roi_obj.roi_matrix + max_number
        counter += 1

    # make a new object to return
    merged_obj = roi_object()
    merged_obj.roi_matrix     = roi_matrix
    merged_obj.indices        = convert_matrix_rois_to_inds(roi_matrix)
    merged_obj.red_rois       = convert_matrix_to_redmatrix(roi_matrix)
    merged_obj.number_of_rois = int(np.amax(roi_matrix))
    merged_obj.kind           = list_of_roi_objects[0].roi_obj.kind
    merged_obj.x_indices      = convert_inds_to_xinds(merged_obj.indices)
    merged_obj.y_indices      = convert_inds_to_yinds(merged_obj.indices)
    merged_obj.masks          = convert_roi_matrix_to_masks(roi_matrix)

    return merged_obj

def swap_indices_old_rois(old_indices):
    """
    Swappes x- and y-indices from indices ROIs.
    """
    new_indices = []
    for roi in old_indices:
        one_new_roi = []
        for point in roi:
            one_new_roi.append((point[1],point[0]))
        new_indices.append(one_new_roi)
    return new_indices



def  load_rois_fromh5_address(address):
    filename, groupname = xrs_utilities.split_hdf5_address(address)
    h5file  = h5py.File(filename, "r")
    h5group = h5file[groupname]
    masks={}
    newshape, imagesum = load_rois_fromh5(h5group,masks, retrieveImage=True)
    myroi = roi_object()
    myroi.load_rois_fromMasksDict(masks, newshape=newshape)
    return myroi


def  load_rois_fromh5(h5group_tot,md, retrieveImage=False, metadata = None):
    h5group = h5group_tot["rois_definition/rois_dict"]
    for key in h5group.keys():
        md[key]=[]
        md[key].append(h5group[key]["origin"][:])
        md[key].append(h5group[key]["mask"][:])

        if metadata is not None:
            newmeta = {}
            print("metadata", h5group[key], key)
            if "metadata" in h5group[key]:
                for kk in h5group[key]["metadata"]:
                    newmeta[kk] = h5group[key]["metadata"][kk].value
            metadata[key]=newmeta

    h5data = h5group_tot["rois_definition/image"]
    shape = h5data.shape
    if retrieveImage:
        return shape, np.array(h5data[:])
    else:
        return shape

def  write_rois_toh5(h5group,md, filterMask=None, metadata=None):
    for key in md.keys():
        if key in h5group:
                del h5group[key]
                
        h5group.require_group(key)

        README="""origin : the Y and X coordinates of the bottom left corner of the mask (putting the origin at the bottom left)
mask : the mask , 1 for considered pixel, zero for discarded ones.
"""
            
        h5_assign_force(h5group[key], "README" , README   )


        
        h5group[key]["origin"]=md[key][0]
        if filterMask is None:
                h5group[key]["mask"]=md[key][1]
        else:
                ori = md[key][0]
                sh  = md[key][1].shape
                Mfilter =  filterMask[ori[0]:ori[0]+sh[0], ori[1]:ori[1]+sh[1]    ]
                h5group[key]["mask"]=md[key][1]*Mfilter
                
        if metadata is not None:
            if key in metadata:
                metagroup = h5group[key].require_group("metadata")

                for kk,mdata in metadata[key].items():
                    metagroup[kk] = mdata







