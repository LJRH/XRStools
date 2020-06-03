from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import map
from six.moves import range
from six.moves import zip
#!/usr/bin/python
# Filename: xrs_utilities.py

__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

import numpy as np
import array as arr
import collections

# # try to import the fast PyMCA parsers
# try:
#     import PyMca.EdfFile as EdfIO
#     # use as: data = EdfIO.EdfFile(fname,"r").GetData(0)
#     import PyMca.specfilewrapper as SpecIO
#     # use as:
#     # sf = SpecIO.Specfile(filename)
#     # sf.scanno() # returns the number of scans in the SPEC file
#     # scan = sf.select('1') # returns the first scan in the specfile
#     # scan.data() # np.array of the data
#     # scan.alllabels() # Python list of all counter names
#     use_PyMca = True
# except:
#     use_PyMca = False


try:
    import PyMca5.PyMcaIO.EdfFile as EdfIO
    import PyMca5.PyMcaIO.specfilewrapper as SpecIO
    import fabio
    from silx.io.specfile import SpecFile

    use_PyMca = True
except:
    use_PyMca = False

print( " >>>>>>>>  use_PyMca " , use_PyMca)

__metaclass__ = type # new style classes



SHOW_LOADED_FILES = False  # this is useful when you want to take home a limited number
                           # of files for an experiment : you run the analysis at the lab :
                           # you look a the filenames that are printed and the you tar gz
                           # that ones for home. This without having to bring home the 
                           # whole shift

def SpecRead(filename,nscan):
    """Parses a SPEC file and returns a specified scan.

    Args:
      * filename (string): SPEC file name (inlc. path)
      * nscan (int): Number of the desired scan.
 
    Returns:
      * data (np.array): array of the data from the specified scan.
      * motors (list): list of all motor positions from the header of the specified scan.
      * counters (dict): all counters in a dictionary with the counter names as keys.

    """
    scannid   = '#S'
    countid   = '#L'
    motorid   = '#P'
    data      = []
    motors    = []
    counterss = []    
    f = open(filename,'r')
    while True:
        line = f.readline()
        if not line: break
        if line[0:2] == scannid:
            if int(line.split()[1]) == nscan:
                line = '##'+line
                while line and line[0:2]!='#S':
                    line = f.readline() 
                    if not line:
                        break                    
                    if line[0:2] == countid:
                        cline = '  '+line[2:]
                        counterss = [n.strip() for n in [_f for _f in cline.split('  ')[1:] if _f]]
                    if line[0:2] == motorid:
                        motors.append([float(n) for n in line.strip().split()[1:]])                    
                    if line[0] != '#':
                        data.append([float(n) for n in line.strip().split()])
    data.pop(-1) # the trailing empty line                    
    f.close()
    # put the data into a dictionary with entries from the counterss
    counters = {}
    for n in range(len(counterss)):
        counters[counterss[n].lower()] = [row[n] for row in data] # data[:,n]
    return data, motors, counters

def myEdfRead(filename):
    """
    Returns EDF-data, if PyMCA is not installed (this is slow).
    """    
    # get some info from header
    f = open(filename,'rb').readlines()
    counter = 0
    predata = []
    for entry in f:
        counter += 1
        if entry.strip().split()[0] == '}':
            break
    for entry in f[:counter]:
        if entry.strip().split()[0] == 'Dim_1':
            dim1 = int(entry.strip().split()[2])
        if entry.strip().split()[0] == 'Dim_2':
            dim2 = int(entry.strip().split()[2])
        if entry.strip().split()[0] == 'Size':
            size = int(entry.strip().split()[2])
    length = 0
    for line in f:
        length += len(line)
    headerlength = (length-size)//2
    # get the data
    f = open(filename,'rb')
    predata = arr.array('H')
    predata.fromfile(f,(headerlength+dim1*dim2))
    data = np.reshape(predata[headerlength:],(dim2,dim1))
    f.close()
    return data

def PyMcaSpecRead_my(filename,nscan):
    """
    Returns data, counter-names, and EDF-files using PyMCA.
    """
    sf   = SpecIO.Specfile(filename)
    scan = sf.select(str(nscan))
    data = scan.data()
    lables = scan.alllabels()
    counters = {}
    cou = 0
    for lable in lables:
        counters[lable.lower()] = data[cou,:]
        cou += 1
    #motors = [np.array(ii.split()[1::]).astype(np.float) for ii in scan.header('P')]
    #motors = collections.OrderedDict( zip(scan.allmotors(),  map( scan.motorpos, scan.allmotors())) )
    motors = dict( list(zip(scan.allmotors(),  list(map( scan.motorpos, scan.allmotors())))) )
    return data.T, motors, counters, lables

def SilxSpecRead( filename, nscan ):
    """
    Returns data, motors, counter-names, and labels using Silx.
    """
    nscan -= 1
    sf     = SpecFile( filename )
    data   = sf.data( nscan )
    lables = sf.labels( nscan )
    counters = {}
    cou = 0
    for lable in lables:
        counters[lable.lower()] = data[:,cou]
        cou += 1
    motors = dict( list(zip(sf.motor_names(nscan),
            list( map( sf.motor_position_by_name,
                        [nscan]*len(sf.motor_names(nscan)),
                           sf.motor_names( nscan ))))) )
    return data, motors, counters, lables

def PyMcaSpecRead(filename,nscan):
    """
    Returns data, counter-names, and EDF-files using PyMCA.
    """
    sf   = SpecIO.Specfile(filename)
    scan = sf.select(str(nscan))
    data = scan.data()
    lables = scan.alllabels()
    counters = {}
    cou = 0
    for lable in lables:
        counters[lable.lower()] = data[cou,:]
        cou += 1
    #motors = [np.array(ii.split()[1::]).astype(np.float) for ii in scan.header('P')]
    #motors = collections.OrderedDict( zip(scan.allmotors(),  map( scan.motorpos, scan.allmotors())) )
    motors = dict( list(zip(scan.allmotors(),  list(map( scan.motorpos, scan.allmotors())))) )
    return data.T, motors, counters

def PyMcaEdfRead(fname):
    """
    Returns the EDF-data using PyMCA.
    """
    data = EdfIO.EdfFile(fname,"r").GetData(0)
    return data

def FabioEdfRead(fname):
    """
    Returns the EDF-data using FabIO.
    """
    return fabio.edfimage.EdfImage().read(fname).data

def ReadScanFromFile(fname):
    """
    Returns a scan stored in a Numpy archive.
    """
    try:
        print( 'Trying to load scan from file.')
        scanname = 'Scan%03d' % scannumber
        scan     = np.load(fname)
        data     = list(scan['data'])
        motors   = list(scan['motors'])
        counters = scan['counters'].item()
        edfmats  = scan['edfmats']
        return data, motors, counters, edfmats
    except:
        print( 'Failed loading scan from file, will read EDF- and SPEC-file.')
        pass

def WriteScanToFile(fname,data,motors,counters,edfmats):
    """
    Writes a scan into a Numpy archive.
    """
    try:
        np.savez(fname, data=data, motors=motors, counters=counters, edfmats=edfmats)
        print( 'trying to save file in numpy-archive.')
    except:
        print( 'Could not write ' + fname + '.')
        pass

def PrepareEdfMatrix_TwoImages(scan_length,num_pix_x,num_pix_y):
    """
    Returns np.zeros for old data (horizontal and vertical Maxipix
    images in different files).
    """
    edfmats_h = np.zeros((scan_length,num_pix_y//2,num_pix_x))
    edfmats_v = np.zeros((scan_length,num_pix_y//2,num_pix_x))
    edfmats   = np.zeros((scan_length,num_pix_y,num_pix_x))
    return edfmats_h, edfmats_v, edfmats

def PrepareEdfMatrix(scan_length,num_pix_x,num_pix_y):
    """
    Returns np.zeros of the shape of the detector.
    """
    edfmats   = np.zeros((scan_length,num_pix_y,num_pix_x))
    return edfmats

def ReadEdfImages_TwoImages(ccdcounter,num_pix_x,num_pix_y,path,EdfPrefix_h,EdfPrefix_v, EdfNmae, EdfPostfix):
    """
    Reads a series of EDF-images and returs them in a 3D Numpy array
    (horizontal and vertical Maxipix images in different files).
    """
    edfmats_h, edfmats_v, edfmats = PrepareEdfMatrix_TwoImages(len(ccdcounter),num_pix_x,num_pix_y)
    for m in range(len(ccdcounter)):
        ccdnumber = ccdcounter[m]
        edfnameh   = path + EDF_PREFIXh + edfName + '_' + "%04d" % ccdnumber + EDF_POSTFIX
        edfnamev   = path + EDF_PREFIXv + edfName + '_' + "%04d" % ccdnumber + EDF_POSTFIX
        if use_PyMca == True:
            edfmatsh[m,:,:] = PyMcaEdfRead(edfnameh)
            edfmatsv[m,:,:] = PyMcaEdfRead(edfnamev)
            edfmats[m,0:num_pix_y//2,:] = edfmatsv[m,:,:]
            edfmats[m,num_pix_y//2:,:]  = edfmatsh[m,:,:]
        else:
            edfmatsh[m,:,:] = myEdfRead(edfnameh)
            edfmatsv[m,:,:] = myEdfRead(edfnamev)
            edfmats[m,0:num_pix_y//2,:] = edfmatsv[m,:,:]
            edfmats[m,num_pix_y//2:,:]  = edfmatsh[m,:,:]
    return edfmats

def EdfRead(fname):
    if use_PyMca == True:
        # print " using pymca " 
        edfmat = PyMcaEdfRead(fname)
    else:
        # print "NOT using pymca " 
        edfmat = myEdfRead(fname)
    return edfmat

def ReadEdf_justFirstImage(ccdcounter,  path, EdfPrefix, EdfName, EdfPostfix):

    m=0
    ccdnumber = ccdcounter[m]
    fname   = path + EdfPrefix + EdfName + '_' + "%04d" % ccdnumber + EdfPostfix
    if use_PyMca == True:
        print( " using pymca " )
        edfmat = PyMcaEdfRead(fname)
    else:
        print( "NOT using pymca " )
        edfmat = myEdfRead(fname)
    return edfmat

def ReadEdfImages(ccdcounter, num_pix_x, num_pix_y, path, EdfPrefix, EdfName, EdfPostfix):
    """
    Reads a series of EDF-images and returs them in a 3D Numpy array
    (horizontal and vertical Maxipix images in different files).
    """
    edfmats = PrepareEdfMatrix(len(ccdcounter),num_pix_x,num_pix_y)
    for m in range(len(ccdcounter)):
        ccdnumber = ccdcounter[m]
        fname   = path + EdfPrefix + EdfName + '_' + "%04d" % ccdnumber + EdfPostfix
        print( " for file ", fname)
        if use_PyMca == True:
            print( " using pymca " )
            edfmats[m,:,:] = PyMcaEdfRead(fname)
        else:
            print( "NOT using pymca " )
            edfmats[m,:,:] = myEdfRead(fname)
    return edfmats

def ReadEdfImages_my(ccdcounter, path, EdfPrefix, EdfName, EdfPostfix):
    """
    Reads a series of EDF-images and returs them in a 3D Numpy array
    (horizontal and vertical Maxipix images in different files).
    """

    fname   = path + EdfPrefix + EdfName + "%04d" % ccdcounter[0] + EdfPostfix
    xyShape = np.shape(myEdfRead(fname))
    edfmats = PrepareEdfMatrix(len(ccdcounter),xyShape[1],xyShape[0])
    for m in range(len(ccdcounter)):
        ccdnumber = ccdcounter[m]
        fname   = path + EdfPrefix + EdfName + "%04d" % ccdnumber + EdfPostfix
        if SHOW_LOADED_FILES : print( " LEGGO ", fname)
        edfmats[m,:,:] = myEdfRead(fname)

    return edfmats

def ReadEdfImages_PyMca(ccdcounter, path, EdfPrefix, EdfName, EdfPostfix):
    """
    Reads a series of EDF-images and returs them in a 3D Numpy array
    (horizontal and vertical Maxipix images in different files).
    """
    fname   = path + EdfPrefix + EdfName + "%04d" % ccdcounter[0] + EdfPostfix

    xyShape = np.shape(PyMcaEdfRead(fname))
    edfmats = PrepareEdfMatrix(len(ccdcounter),xyShape[1],xyShape[0])
    for m in range(len(ccdcounter)):
        ccdnumber = ccdcounter[m]
        fname   = path + EdfPrefix + EdfName + "%04d" % ccdnumber + EdfPostfix
        if SHOW_LOADED_FILES :  print( " LEGGO ", fname)
        edfmats[m,:,:] = PyMcaEdfRead(fname)
    return edfmats

def readbiggsdata(filename,element):
    """
    Reads Hartree-Fock Profile of element 'element' from values tabulated 
    by Biggs et al. (Atomic Data and Nuclear Data Tables 16, 201-309 (1975))
    as provided by the DABAX library (http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/ComptonProfiles.dat).
    input:
    filename = path to the ComptonProfiles.dat file (the file should be distributed with this package)
    element  = string of element name
    returns:
    data     = the data for the according element as in the file:
        #UD  Columns: 
        #UD  col1: pz in atomic units 
        #UD  col2: Total compton profile (sum over the atomic electrons
        #UD  col3,...coln: Compton profile for the individual sub-shells
    occupation = occupation number of the according shells
    bindingen  = binding energies of the accorting shells
    colnames   = strings of column names as used in the file
    """
    elementid = '#S'
    sizeid    = '#N'
    occid     = '#UOCCUP'
    bindingid = '#UBIND'
    colnameid = '#L'
    data = []
    f = open(filename,'r')
    istrue = True
    while istrue:
        line = f.readline()
        if line[0:2] == elementid:
            if line.split()[-1] == element:
                line = f.readline()
                while line[0:2] != elementid:
                    if line[0:2] == sizeid:
                        arraysize = int(line.split()[-1])
                        line = f.readline()
                    if line[0:7] == occid:
                        occupation = line.split()[1:]
                        line = f.readline()
                    if line[0:6] == bindingid:
                        bindingen = line.split()[1:]    
                        line = f.readline()
                    if line[0:2] == colnameid:
                        colnames = line.split()[1:]
                        line = f.readline()
                    if line[0]== ' ':
                        data.append([float(n) for n in line.strip().split()])
                        #data = np.zeros((31,arraysize))
                        line = f.readline()
                break
    length = len(data)
    data = (np.reshape(np.array(data),(length,arraysize)))
    return data, occupation, bindingen, colnames










