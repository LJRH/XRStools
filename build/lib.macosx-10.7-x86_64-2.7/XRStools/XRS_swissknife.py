from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import scipy
### __doc__="""


generality_doc = 1
"""
Documentation of XRS_swissknife
-------------------------------

The script is called in this way ::

  XRS_swissknife yourinput.yaml

The input file is in *yaml* format. In *yaml* format each line introduces an item
and the indentation expresses the hierarchy.
An example is ::

  Fancy_Reduction :
        parameterTom : 3.14
        parameterJerry : False
        choicesBob     : [1,2,3]

In this example we create an item called *Fancy_Reduction* which contains three subitems.

The *XRS_swissknife* expects that for each operation that you want to apply, you provide an item
named as the operation, and the associated subitems that provide that values for the parameters.

*XRS_swissknife* will do what you want provided you respect the proper indentation. A thing which helps
is using emacs and activating the *python mode*, because python uses the same indentation principle
to structure the code. 

Each processing item has the additional, optional, key **active**.
This key can be set to **0** or **1** to desactivate or not (default is **1**, active) the processing.
Here a desactivation example ::

  Fancy_Reduction :
        active : 0
        parameterTom : 3.14
        parameterJerry : False
        choicesBob     : [1,2,3]

The following documentation has been created automatically, for each functionality, based on the documentation string
written in the code for the fonctionality itself.  You can also write mathematical expressions :  :math:`\\int{ x dx} = \\frac { x^2 }{ 2 }`
and even include graphics.

"""
import collections

try:
    from mayavi import mlab
except:
    print( " WAS not able to load mayavi, some feature might be missing ")

import string

import numpy as np
import math

import yaml
from yaml import load, dump

import numbers
import re
import yaml
import yaml.resolver
import PyMca5.PyMcaIO.specfilewrapper as SpecIO
import fabio
from six import u


from silx.gui  import qt as Qt
## from  PyQt4 import Qt, QtCore
from XRStools.roiNmaSelectionGui   import roiNmaSelectionWidget
from XRStools import roiSelectionWidget

import h5py
import sys

yaml.resolver.Resolver
Resolver = yaml.resolver.Resolver
Resolver.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u(r"""^(?:[-+]?(?:[0-9][0-9_]*)(\.[0-9_]*)?(?:[eE][-+]?[0-9]+)?
                    |\.[0-9_]+(?:[eE][-+][0-9]+)?
                    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*
                    |[-+]?\.(?:inf|Inf|INF)
                    |\.(?:nan|NaN|NAN))$"""), re.X),
        list(u'-+0123456789.'))

try:
    from mpi4py import MPI
    myrank = MPI.COMM_WORLD.Get_rank()
    nprocs = MPI.COMM_WORLD.Get_size()
    procnm = MPI.Get_processor_name()
    comm = MPI.COMM_WORLD
    print( "MPI LOADED , nprocs = ", nprocs)
except:
    class FakeComm:
        def Barrier(self):
            pass
        def allreduce(self,number, operation):
            assert ( isinstance(number, numbers.Number)    )
            return number
        def Get_size(self):
            return 1

    myrank=0
    nprocs = 1
    comm = FakeComm()
    print( "no MPI LOADED , nprocs = ", nprocs)



def stripROI(t):
    if "ROI" in t:
        return t[3:]
    else:
        return t

def filterRoiList(l, strip=False, prefix="ROI"):
    
    if strip:
        result =  [ str(int(stripROI(t))) for t in l if ( t not in [ "motorDict", "response"] and t[:len(prefix)]==prefix )]
    else:
        result =  [ t for t in l  if ( t not in [ "motorDict", "response"] and  t[:len(prefix)]==prefix ) ]
        
    result.sort()
    return result

def filterScanList(l, prefix="Scan"):
    result =  [ t for t in l  if ( t not in [ "motorDict", "response"] and t[:len(prefix)]==prefix ) ]
    result.sort()
    return result
    

def checkNoParallel( routineName):
    if nprocs>1:
        msg = " ERROR : %s feature not yet parallel : relaunch it with 1 process only "
        print( msg)
        raise Exception( msg)

# try:
#     from yaml import CLoader as Loader, CDumper as Dumper
# except ImportError:
#     from yaml import Loader, Dumper            




nprocs = comm.Get_size()

if nprocs>1:
    # circumvent issus with mpi4py not stopping mpirun
    def excepthook(type, value, traceback):
        res = sys.__excepthook__(type, value, traceback)
        comm.Abort(1)
        return res
    sys.excepthook=excepthook



import os
from XRStools import xrs_rois
from XRStools import roifinder_and_gui
from XRStools import xrs_scans
from XRStools import xrs_read
from XRStools import rixs_read
from XRStools import theory
from XRStools import extraction

from XRStools import xrs_prediction
from XRStools import xrs_imaging


from XRStools import xrs_imaging
from XRStools import superr

#################################################################
##  THIS redefinition of yaml is used to keep the entry ordering
## when accessing yamlData keys
##
#yaml_anydict.py
import yaml
from yaml.representer import Representer
from yaml.constructor import Constructor, MappingNode, ConstructorError

from XRStools import fit_spectra
from XRStools import reponse_percussionelle

    
def dump_anydict_as_map( anydict):
    yaml.add_representer( anydict, _represent_dictorder)
def _represent_dictorder( self, data):
    return self.represent_mapping('tag:yaml.org,2002:map', data.items() )

class Loader_map_as_anydict( object):
    'inherit + Loader'
    anydict = None      #override
    @classmethod        #and call this
    def load_map_as_anydict( klas):
        yaml.add_constructor( 'tag:yaml.org,2002:map', klas.construct_yaml_map)

    'copied from constructor.BaseConstructor, replacing {} with self.anydict()'
    def construct_mapping(self, node, deep=False):
        if not isinstance(node, MappingNode):
            raise ConstructorError(None, None,
                    "expected a mapping node, but found %s" % node.id,
                    node.start_mark)
        mapping = self.anydict()
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            try:
                hash(key)
            except TypeError as exc:
                raise ConstructorError("while constructing a mapping", node.start_mark,
                        "found unacceptable key (%s)" % exc, key_node.start_mark)
            value = self.construct_object(value_node, deep=deep)
            mapping[key] = value
        return mapping

    def construct_yaml_map( self, node):
        data = self.anydict()
        yield data
        value = self.construct_mapping(node)
        data.update(value)


class myOrderedDict (collections.OrderedDict):
    def __setitem__(self,a,b):
        ## print "cucu",a,b
        if type(a)==type("") and a in self:
            self[a+"_tagkajs"]=b 
        else:
            ## print super(myOrderedDict, self)
            super(myOrderedDict, self).__setitem__(a,b )


def cleaned(key):
    while key[-8:]=="_tagkajs":
        key=key[:-8]
    return key
  

            
dictOrder = myOrderedDict

class Loader( Loader_map_as_anydict, yaml.Loader):
    anydict = dictOrder
Loader.load_map_as_anydict()
dump_anydict_as_map( dictOrder)


##
## END of yaml redefinition
###############################

def check_libre( h5 , groupname   ) :

    print("check_libre DESACTIVATED. RESULTS CAN BE OVERWRITTEN")
    return

    if type(h5)==type(""):
        h5f = h5py.File(h5, "r"  )
        h5f.close()
        if groupname in h5f:
            msg=(("ERROR: %s key already present in the hdf5 file %s. Erase it before if you dont need it.\n" % (groupname, h5))*10 )
            print( msg)
            raise Exception( msg)
        else:
            if groupname in h5:
                msg = (("ERROR: %s key already present in the hdf5 file. Erase it before if you dont need it.\n"%groupname)*10 )
                print( msg)
                raise Exception( msg)



inputtext=""

def main():
    global  inputtext
    filename = sys.argv[1]

    inputfile = open(filename,"r")
    yamlData = load(inputfile, Loader=Loader)
    inputtext = open(filename,"r").read()

    for key in list(yamlData.keys()):
      
        mydata = yamlData[key]
        if isinstance(mydata,dict) and    "active"  in mydata :
            if mydata["active"]==0:
                # print " continuo " 
                continue

        if key != "help":
            swissknife_operations[cleaned(key)](  mydata )
        else:

            if cleaned(key)  not in parallelised_operations:
                checkNoParallel( cleaned(key)  )
            
            swissknife_operations[cleaned(key)]( yamlData  )
    
    # workbench_file = yamlData["workbench_file"]


def help(yamlData):
    """
    **help**

    Displays doc on the operations. In the input file ::

       help :

    will trigger printing of all the available operation names ::

       help :
           create_rois
           load_scans

    will print the help on *create_rois* and the help about *load_scans*.
    By the way, it is the same that you can read here because the *sphinx* doc-generation tool
    reads the same docstrings contained in the code.
     
    """
    print( " HELP " *15)
    if yamlData ["help"] is None:
        print( """
              Printing all the function names
              To get help on a specific function:

              help : "functionName"
        """)
        for key,func in swissknife_operations.items():
            print( " FUNCTION : " , key)

    else:
        func = swissknife_operations[ yamlData ["help"]]
        print( "---------------------------------------")
        print( func.__doc__)


def Extraction(mydata):
    """ 
    **Extraction**

    Function to extract the interesting signal after removal of Compton profile,
    linear baselines,  Pearson profile
    Example ::

     Extraction :
         active : 1
         dataadress   : "pippo.h5:/ROI_A/loaded_datas"         # where load_scans wrote data
         HFaddress    : "pippo.h5:/ROI_A/loaded_datas/HF_O"    # where compton profiles have been calculated
         # prenormrange : [ 5 , .inf ]        

         analyzerAverage :                                     #  averaging over analysers
             active : 1
             which : [0,11  , 36,59   ]
             errorweighing  : False

         removeLinearAv :                                      #  fit a linear baseline and remove it
             active  : 1
             region1 :  [520.0,532.0]   
             region2 :  None 
             ewindow : 100 
             scale : 1

         removePearsonAv:                                      # fit a Pearson and remove it
             active  : 0
             region1 :  [520.0,532.0]   
             region2 :  None  
             guess :
                 Peak_position : 600.0
                 FWHM          :  10
                 Shape         : "Lorentzian" 
                 Peak_intensity: 100.0
                 linear_slope  : 1
                 linear_background : 0
                 scaling_factor : 1

         view   :   0
         target :   "myextraction"                            # path relative to dataadress where extracted signal will be written

    """

    reader , filename, groupname= read_reader(mydata, name="dataadress")
    HF     = read_HF(mydata, name="hfspectrum_address")

    extr  = extraction.extraction(reader , HF)

    if  ("analyzerAverage") in mydata :
        aa_data = mydata["analyzerAverage"]
        if gvord(aa_data,"active",True):
            which = aa_data["which"]
            errorweighing  = gvord(aa_data,"errorweighing",False)
            extr .analyzerAverage(which,errorweighing=errorweighing)

    if  ("removeLinearAv") in mydata :
        rla_data = mydata["removeLinearAv"]
        if gvord(rla_data,"active",True):
            region1 = rla_data["region1"]
            region2 = gvord( rla_data,"region2",None)
            ewindow = gvord( rla_data,"ewindow",100)
            scale = gvord( rla_data,"scale",100)
            extr .removeLinearAv(region1, region2=region2,ewindow=ewindow, 
                                 scale=scale, view = gvord(mydata,"view",False),                                 
                             ) 
    print( gvord(mydata,"view",False))


    groupname =  groupname+"/"+ mydata["target"]
    check_libre( filename , groupname   ) 

    extr.save_state_hdf5( filename,groupname, comment = inputtext )
   


def split_hdf5_address(dataadress):
 
    pos = dataadress.rfind(":")
    if ( pos==-1):
        raise Exception( """
roiaddress   must be given in the form  roiaddress : "myfile.hdf5:/path/to/hdf5/group"
but : was not found
""")
    filename, groupname = dataadress[:pos], dataadress[pos+1:]
    return filename, groupname 



def read_HF(mydata, name="hfspectrum_address"):

    dataadress = mydata[name]

    filename, groupname = split_hdf5_address(dataadress)

    HF = theory.HFspectrum(None,None,None, initialise=False)
    HF.load_state_hdf5( filename, groupname)
    return HF


def  HFspectrum(mydata):
    """
    **HFspectrum**

    function for building S(q,w) from tabulated Hartree-Fock Compton profiles to use in
    the extraction algorithm.

    EXAMPLE ::

      dataadress : "hdf5filename:full_nameof_signals_group"   # where load_scans wrote data
      formulas   :  ['O']     # list of strings of chemical sum formulas of which the sample is made up
      concentrations : [1.0]  # list of concentrations of how the different chemical formulas are mixed (sum should be 1)
      correctasym    : [[0.0,0.0,0.0]]  #  single value or list of scaling values for the HR-correction to 
                                # the 1s, 2s, and 2p shells. one value per element in the list of formulas
      hfspectrum_address : "nameofgroup" # Target group for writing Relative to dataadress (and in the same file)!!!!

 
    """


    reader , filename, groupname= read_reader(mydata, name="dataadress")

    hf   = theory.HFspectrum(reader ,   
                             mydata["formulas"]     ,
                             mydata["concentrations"]     ,
                             mydata["correctasym"]  
                             )

    groupname = groupname+"/"+ mydata["hfspectrum_address"]
    check_libre( filename , groupname   ) 

    hf.save_state_hdf5( filename,groupname , comment = inputtext )
    

def load_scans(mydata):
    """
    **load_scans**

    This command harvest the selected signals.
    the instructions on the scans to be taken must be in the form( as submembers ofload_scans ) ::


     load_scans :
         roiaddress :  "hdf5filename:nameofroigroup"  # the same given in create_rois
         expdata    :  "absolutepathtoaspecfile"  # this points to a spec file

         elastic_scans    : [623]
         fine_scans       : [626,630,634,638,642]
         n_loop           : 4
         long_scan        : 624

         signaladdress : "nameofsignalgroup"  # Target group for writing Relative to ROI (and in the same file)!!!!

         #############################################################
         # OPTIONALS
         #
         order : [0,1,2,3,4,5]  #  list of integers (0-5) which describes the order of modules in which the 
                                #    ROIs were defined (default is VD, VU, VB, HR, HL, HB; i.e. [0,1,2,3,4,5])

         rvd : -41              # mean tth angle of HL module (default is 0.0)
         rvu : 85               # mean tth angle of HR module (default is 0.0)
         rvb : 121.8            # mean tth angle of HB module (default is 0.0)
         rhl : 41.0             # mean tth angle of VD module (default is 0.0)
         rhr : 41.0             # mean tth angle of VU module (default is 0.0)
         rhb : 121.8            # mean tth angle of VB module (default is 0.0)


    #
    """

    roiaddress=None
    roiaddress = mydata["roiaddress"]

    filename, groupname = split_hdf5_address (roiaddress)
    file= h5py.File(filename,"r")
    rois = {}
    shape=xrs_rois.load_rois_fromh5(file[groupname],rois)
    file.close()

    roiob = xrs_rois.roi_object()
    roiob.load_rois_fromMasksDict(rois ,  newshape = shape, kind="zoom")
    
    reader = xrs_read.read_id20(mydata["expdata"] , monitorcolumn='kapraman')


    reader.set_roiObj(roiob)

    elastic_scans   =  mydata["elastic_scans"][:] 
    fine_scans      =  mydata["fine_scans"][:] 
    n_loop          =  mydata["n_loop"] 
    long_scan       =  mydata["long_scan"]

    reader.loadelasticdirect(elastic_scans)
    reader.loadloopdirect(fine_scans,n_loop)
    print( " LUNGO " )
    reader.loadlongdirect(long_scan)

    

    reader.getspectrum()

    reader.geteloss()
    reader.gettths(
        rvd   = gvord(mydata,"rvd",0.0) ,
        rvu   = gvord(mydata,"rvu",0.0) ,
        rvb   = gvord(mydata,"rvb",0.0) ,
        rhl   = gvord(mydata,"rhl",0.0) ,
        rhr   = gvord(mydata,"rhr",0.0) ,
        rhb   = gvord(mydata,"rhb",0.0) ,
        order = gvord(mydata,"order",        [0,1,2,3,4,5])
        )

    groupname = groupname+"/"+ mydata["signaladdress"]
    check_libre( filename , groupname   ) 

    reader.save_state_hdf5( filename, groupname , comment = inputtext )





def volume_from_2Dimages(mydata):
    """
    imagesaddress :  "test_imaging.hdf5:/ROI_A/images"  # where the data have been saved

    scan_interval    :  [372,375]                    # OPTIONAL : can be shorter then the scans effectively present in the file
    roi_n            : 0           # OPTIONAL. if not given, the first non empty found roi. Starts from 0

    imagesaddress : "myfile.hdf5:/path/to/hdf5/data"  # OPTIONAL. the target destination for volume. if not given mayavi is launched on the fly.

    """


    reader = xrs_imaging.oneD_imaging( "bidon"  , "bidon",  "bidon"  , "bidon")


    imagesaddress = mydata["imagesaddress"]
    filename, groupname = split_hdf5_address(imagesaddress)

    reader.load_state_hdf5( filename, groupname)

    scan_names = list( reader.twoDimages.keys() )
    scan_ids = map(int, [name[4:] for name in scan_names  ] )
    order = np.argsort(scan_ids)

    if not  ('scan_interval') in mydata :
        scan_names = [  scan_names[id] for id in order ] 
    else:
        scan_interval = mydata['scan_interval']

        print( order)
        print( scan_names)
        print( scan_interval)
        scan_names = [  scan_names[id] for id in order if scan_ids[id]>=scan_interval[0] and scan_ids[id]<scan_interval[1]  ] 
    
    first_name = scan_names[0]
    roi_n=0
    
    if not  ('roi_n' in mydata ):
        while(1):
            shape  = reader.twoDimages[first_name][roi_n].matrix.shape
            if shape != (0,) : 
                break
            roi_n+=1
    else:
        roi_n = mydata["roi_n"]
        shape  = reader.twoDimages[first_name][roi_n].matrix.shape
            
    Volume  = np.zeros((     shape[0],  shape[1]   ,     len(scan_names)    ))

    for i,scanname  in enumerate(scan_names):
        Volume[:,:,i] = reader.twoDimages[scanname][roi_n].matrix

    if   ('volumeaddress') in mydata :
        filename, groupname = split_hdf5_address( mydata['volumeaddress'] )

        h5=h5py.File(filename,'a')
        check_libre( h5 , groupname   ) 
        h5[groupname] = Volume
        h5.close()
        h5=None
    else:
       view_Volume_myavi_(Volume) 



def view_Volume_myavi(mydata):
    """
    volume_address : "myfile.hdf5:/path/to/hdf5/group"  # the target destination for volume. 
    """


    filename, groupname = split_hdf5_address( mydata['volume_address'] )
    h5=h5py.File(filename,'r')
    Volume = h5[groupname] [:]
    h5.close()
    h5=None


    isolevel = mydata['isolevel']
    opacity = mydata['opacity']

    
    view_Volume_myavi_(Volume,  isolevel, opacity)
    


def view_Volume_myavi_(V,  isolevel, opacity) :
    print( " IN view ")
    src = mlab.pipeline.scalar_field(V)
    mlab.pipeline.iso_surface(src, contours=[V.min()+isolevel*V.ptp(), ], opacity=opacity)
    mlab.show()
    src = mlab.pipeline.scalar_field(V)
    mlab.pipeline.volume(src,vmin=1000.0, vmax=2000.0)
    mlab.show()


def calculate_recenterings(mydata):
    """
    **calculate_recenterings**

    calculates offsets to go from baricenter A to baricenter B, for all the rois
    in 
 

      calculate_recenterings:
         bariA : "demo_rois.h5:/ROI_AS_SELECTED/images_broad/scans/scan342"
         bariB : "demo_rois.h5:/ROI_AS_SELECTED/energy_scan/scans/scan237"
         target: "recenterings.h5:/recenterings4rois" 
    #
    """
    
    bariA = mydata["bariA"]
    bariA_filename, bariA_groupname = split_hdf5_address( bariA )
    
    bariB = mydata["bariB"]
    bariB_filename, bariB_groupname = split_hdf5_address( bariB )

    target = mydata["target"]
    target_filename, target_groupname = split_hdf5_address( target )

    print( bariA_filename, bariA_groupname)
    print( " OPENIN FILE FOR RECENTERING ")
    h5A_f = h5py.File(bariA_filename,"r")
    h5A = h5A_f [bariA_groupname]
    if bariB_filename == bariA_filename :
        h5B_f =  h5A_f
    else:
        h5B_f=h5py.File(bariB_filename,"r")
        
    h5B = h5B_f[bariB_groupname]
    offs = {}
    chiavi = filterRoiList(h5A.keys())
    for c in chiavi:
        bA = h5A[c]["barix"].value + h5A[c]["cornerpos"][:][1]
        bB = h5B[c]["barix"].value + h5B[c]["cornerpos"][:][1]
        bAy = h5A[c]["bariy"].value + h5A[c]["cornerpos"][:][0]
        bBy = h5B[c]["bariy"].value + h5B[c]["cornerpos"][:][0]
        offs[c] = [[bBy, bAy ],[bB,bA]]
        
    if h5B_f is not h5A_f:
        h5B_f.close()
    h5A_f.close()


    if os.path.exists(target_filename):
        # check_libre( target_filename ,  target_groupname  )
        h5f = h5py.File(target_filename,"a")
        if target_groupname in h5f:
            del h5f[target_groupname]

        
    else:
        h5f = h5py.File(target_filename,"w")
    h5f.require_group(  target_groupname  )
    h5 = h5f[target_groupname]
    for c in chiavi :
       h5[c] = np.array(  offs[c]    )
    h5f.flush()
    h5f.close()
    h5f = None


def     sum_scans2maps(mydata):

    roiaddress=None
    roiaddress = mydata["mask_file"]
  
    filename, groupname = split_hdf5_address( roiaddress)

    file= h5py.File(filename,"r")
    rois = {}
    shape=xrs_rois.load_rois_fromh5(file[groupname],rois)
    file.close()
 
    specfile_name = mydata["spec_file"]
    Scan_Variable = mydata["Scan_Variable"]
    Motor_Variable = mydata["Motor_Variable"]

    
    specfile = SpecIO.Specfile( specfile_name )


    dirname  = os.path.dirname(   specfile_name )   
    basename = os.path.basename(   specfile_name )  

    
    scans_infos = []
    signals    = []

    s1 = int(mydata["first_scan"])
    s2 = int(mydata["last_scan"])

    roi_names = list(rois.keys())
    roi_list = [  rois[k] for k in roi_names  ]
    
    for i in range(s1,s2+1):
        # print " SCAN lettura " , i
        scan = specfile.select(str(i))
        scan_data = scan.data()
        
        scan_themotor   =  scan.motorpos( Motor_Variable  )
        scan_othermotors = [ scan.motorpos( name   )  for name in  scan.allmotors()  if name != Motor_Variable   ]

        othermotorsname = [name for name in  scan.allmotors()  if name != Motor_Variable]
        
        labels = scan.alllabels()
        scan_variable = scan_data[ labels.index( Scan_Variable )   , :]
        scan_ccdnos  = scan_data[ labels.index( "ccdno" )   , :].astype("i")
        signal = []
        for no in scan_ccdnos:
            print( " opening image ", os.path.join(   dirname   ,  "edf", basename+"_"+str(no)+".edf"))
            data = fabio.open(  os.path.join(   dirname   ,  "edf", basename+"_"+str(no)+".edf" )   ).data 
            tok = [ (data[corner[0]:corner[0]+mask.shape[0],   corner[1]:corner[1]+mask.shape[1]]*mask).sum() for corner, mask in roi_list  ]
            signal.append(tok)
            # print " OK "


        # print signal
        # print " Appendo " , signal

    
        signals.append(np.array(signal))
        # print "OK " 
        scans_infos.append( [scan_themotor, scan_othermotors,  scan_variable  ]  )

        # print " DONE scan " , i
        
    done = np.zeros( len(scans_infos) ,"i")
    DONE={}
    synthes = {}
    for kscan in range(len(scans_infos) ):
        # print " kscan " , kscan
        DONE[kscan]=0
        if done[kscan]:
            continue
        else:
            res  = np.array(signals[kscan])
            
            done[kscan]=1
            kinfos = scans_infos[kscan]
            kM, kOM,  kV = kinfos
            
            for oscan in  range(len(scans_infos)) :
                # print " oscan " , oscan
                if done[oscan]:
                    continue
                else:
                    oinfos = scans_infos[oscan]
                    
                    oM, oOM,  oV = oinfos

                    if kM==oM:
                        if True or (np.abs(np.array(kOM)-np.array(oOM)).sum()==0.0):
                            print( " SONO UGUALI " )
                            if len(oV)== len(kV) :
#                                   if  np.abs(kV-oV).sum()==0.0:
                                    print( " AGGIUNGO " )
                                    res = res+np.array(signals[oscan])
                                    done[oscan]=1


            print( " AGGIUNGO " , kM, len(kV), len(kOM))
            synthes[kscan] = [  kM, kV, kOM, res  ]
            done[kscan]    = 1


    # print " QUI " 

    target = mydata["scans_file"]
    target_filename, target_groupname = split_hdf5_address( target)
    if os.path.exists(target_filename):
        h5f = h5py.File(target_filename,"a")
    else:
        h5f = h5py.File(target_filename,"w")
    h5 = h5f.require_group(target_groupname)
    target_subdir = "collection_"+str(s1)+"_"+str(s2)
    if target_subdir in  h5:
        del h5[ target_subdir]
        
    h5 =  h5.require_group(target_subdir)
            
    for kscan in    list(synthes.keys() ):
        if DONE[kscan] :
            continue
        
        kM, kV, kOM, kdata  = synthes[kscan] 

        # print " kdata.shape" , kdata.shape
         
        
        myscans = [kdata]
        myM     = [kM    ]
        myV     = [kV    ]
        DONE[kscan] = 1

        for oscan in    list(synthes.keys() ):
            if DONE[oscan]:
                continue
            else:
                oM, oV, oOM, odata  = synthes[oscan]
                # print " COMPARO "
                # diff =  np.array(kOM)-np.array(oOM)
                # for nn ,dd in zip(othermotorsname, diff    ):
                #     if dd >0.0:
                #         print nn, dd
                if True or (np.abs(np.array(kOM)-np.array(oOM)).sum()==0.0):
                    print( " OM equal " ,  len(kV) , len(oV))
                    if( len(kV) == len(oV)):
                      # print np.abs(kV-oV)
                      # if np.abs(kV-oV).sum()==0.0:

                        assert( kM!=oM  )
                        DONE[oscan]=1

                        # print " odata.shape " ,  odata.shape
                        
                        myscans.append(odata)

                        # print "  myscans.shape " , np.array(myscans).shape
                        
                        myM.append(oM)
                        myV.append(oV)
        DONE[kscan]=1 
        myscans = np.array(myscans)



        if len(myM)>1:
            order = np.argsort(myM)
            myM     = np.array(myM)[order]
            myscans = np.array(myscans)[order]
            myV     = np.array(myV)
        else:
            myM     = np.array(myM)
            myV     = np.array(myV)
            
        ts = "scan_"+str(int(s1)+kscan)
        h5t = h5.require_group(ts)

        h5t["Variable"] = myV
        h5t["Motor"]    = myM

        s="f, axarr = pylab.plt.subplots(%d, sharex=True)\n"%len(roi_names)

        ly = myM[0]
        hy = myM[-1]

        if (  np.abs(myV[-1,:]-myV[0,:]).sum()==0.0 ):
              xlabel =  "An. Energy"
              lx = myV[0,0]
              hx = myV[0,-1]
        else:
              xlabel =  "Energy Loss"
              diff = myM[0]- myV[0,:]
              lx =  diff[0]
              hx =  diff[-1]


              
        for i,rn in enumerate(roi_names):
            tok = myscans[:,:,i]
            h5t["signal_"+rn] = tok
            
            s=s+"axarr[%d].imshow( %s ,extent=( %e ,%e ,%e ,%e) ) \n"%(i, "self.signal_"+rn, lx,hx,ly,hy)
            s=s+"axarr[%d].set_title('%s')\n"%(i,rn)
            if i==0:
                s=s+"axarr[%d].set_xlabel('%s')\n"%(i,xlabel)
            
        s=s+"pylab.plt.show()\n"
        h5t["python_plot"]=s
                                                             
 



            
    h5f.flush()
    h5f.close()

        
        
"""
     sum_scans2maps :
 spec_file : /mntdirect/_data_visitor/hc2892/id20/run7_hc2892/rixs
 first_scan : 127
 last_scan : 136
 Scan_Variable : Anal Energy
 Motor_Variable : energy
 scans_file : /scisoft/users/mirone/WORKS/matlabID20/rawdata/scansfile.h5:SCAN
"""  

def loadscan_2Dimages(mydata):
    """
    **loadscan_2Dimages**

    This command harvest the selected signals.
    the instructions on the scans to be taken must be in the form( as submembers ofload_scans ) ::


      loadscan_2Dimages :
         roiaddress :  "hdf5filename:nameofroigroup"  # the same given in create_rois
         expdata    :  "absolutepathtoaspecfile"  # this points to a spec file

         scan_interval    :  [372,423]   # from 372 to 422 included

         signaladdress : "nameofsignalgroup"  # Target group for writing Relative to nameofroigroup/ (and in the same file)!!!!
                                              # unless it is in the format filename:groupname

         energycolumn : 'sty'        # OPTIONAL, if not give defaults to sty
         monitorcolumn : 'kapraman'  # OPTIONAL , default is kapraman. If the key is not found in spec file, then no normalisation will be done
                                     # You can also write kapraman/1000 in this case the monito will be divided by 1000
                                     # (or the other number that you write)
         edfName    :  'edfprefix'   # OPTIONAL, if not given autonmatically determined

         sumto1D  : 1  # OPTIONAL, default 0

         isolateSpot : 0   # is different from zero selects on each image  ( when sumto1d=0 ) ROI  the spot region and sets to zero outside a radius =  isolateSpot

         # the following defaults to None
         recenterings : "recenterings.h5:/recenterings4rois"
 
    #
    """

    roiaddress=None
    roiaddress = mydata["roiaddress"]
  
    filename, groupname = split_hdf5_address( roiaddress)

    file= h5py.File(filename,"r")
    rois = {}
    shape=xrs_rois.load_rois_fromh5(file[groupname],rois)
    file.close()

    print( " carico maschere ")
    roiob = xrs_rois.roi_object()
    roiob.load_rois_fromMasksDict(rois ,  newshape = shape, kind="zoom")

    monitor_divider = 1.0
    
    if  ('monitorcolumn') in mydata :
        monitorcolumn = mydata['monitorcolumn']
    else:
        monitorcolumn = 'kapraman'

    pos = monitorcolumn.find("/")
    if pos != -1:
        monitor_divider = float(monitorcolumn[pos+1:])
        monitorcolumn = monitorcolumn[:pos]

    print( "monitorcolumn : " , monitorcolumn)
    print( "monitor_divider : " , monitor_divider)
        
    is_by_refinement =False
    if  ('recenterings') in mydata :
        recenterings = mydata['recenterings']
        recenterings_filename, recenterings_groupname = split_hdf5_address( recenterings )
        h5f = h5py.File(recenterings_filename,"r")
        h5 = h5f[recenterings_groupname]
        recenterings= {}
        chiavi = filterRoiList(h5.keys())
        for c in chiavi:
            recenterings[int(c)]= h5[c][:]
            if recenterings[int(c)].shape == (2,2):
                if nprocs>1:
                    raise Exception("When using recentering with refinement parallelism cannote be used")
                is_by_refinement = True
        h5f.close()
        if is_by_refinement:
            recenterings_confirmed = mydata['recenterings_confirmed']
            recenterings_confirmed_filename, recenterings_confirmed_groupname = split_hdf5_address( recenterings_confirmed )
           
    else:
        recenterings = None

        
    if  ('energycolumn') in mydata :
       energycolumn  = mydata['energycolumn']
    else:
       energycolumn  = 'sty'

    if  ('edfName') in mydata :
        edfName = mydata['edfName']
    else:
       edfName  = None

    if  ('sumto1D') in mydata :
        sumto1D = mydata['sumto1D']
    else:
        sumto1D  = 0
        
    if  ('isolateSpot') in mydata :
        isolateSpot = mydata['isolateSpot']
    else:
        isolateSpot  = 0

    print( " creo oggetto ", energycolumn)
    print( " DIVIDER   ", monitor_divider)
    
    reader = xrs_imaging.oneD_imaging( mydata["expdata"]  ,monitorcolumn = monitorcolumn , monitor_divider=monitor_divider,
                                       energycolumn = energycolumn  , edfName = edfName, sumto1D = sumto1D,
                                       recenterings=recenterings)
    reader.set_roiObj(roiob)

    scan_interval = mydata["scan_interval"]
    print( " LOAD ")
    todo_list = []
    ninterval = len(scan_interval)//2
    for i in range(ninterval):
        todo_list = todo_list + list(range(    int(scan_interval[2*i]),  int(scan_interval[2*i+1])) ) #   *scan_interval[2*i :2*i+2])

    mytodo = np.array_split(todo_list, nprocs) [myrank]
    print( mytodo)
    maxvalue = reader.loadscan_2Dimages( list(mytodo) ,scantype=energycolumn, isolateSpot = isolateSpot)
    if nprocs>1:
        maxvalue = comm.allreduce(maxvalue, op=MPI.MAX)

    if is_by_refinement :
        if nprocs>1:
            raise Exception("When using recentering with refinement parallelism cannote be used")
        if os.path.exists(recenterings_confirmed_filename):
            check_libre(  recenterings_confirmed_filename    , recenterings_confirmed_groupname       )
            print( " APRO IN MODO a ", recenterings_confirmed_filename)
            h5f = h5py.File(recenterings_confirmed_filename,"a")
        else:
            h5f = h5py.File(recenterings_confirmed_filename,"w")
        h5 = h5f.require_group(recenterings_confirmed_groupname)
        for c in chiavi:
            h5[c]  =  reader.recenterings[int(c)]
        h5f.flush()
        h5f.close()
        h5f = None

    signaladdress = mydata["signaladdress"]
    if ":" not in signaladdress:
        groupname = groupname+"/"+ mydata["signaladdress"]+"/"
        check_libre( filename , groupname   )
        save_also_roi = False
        
    else:
        filename , groupname = split_hdf5_address(signaladdress )
        save_also_roi = True

    

    for iproc in range(nprocs):
        comm.Barrier()
        if iproc==myrank:
            reader.save_state_hdf5( filename, groupname, comment = inputtext, myrank = myrank, save_also_roi = save_also_roi  )# factor = 16000.0/maxvalue)

 

def loadscan_2Dimages_galaxies(mydata):
    """
    **loadscan_2Dimages_galaxies**

    This command harvest the selected signals.
    the instructions on the scans to be taken must be in the form( as submembers of load_scans ) ::


     loadscan_2Dimages_galaxies :
         roiaddress :  "hdf5filename:nameofroigroup"  # the same given in create_rois
         expdata    :  "kapton_%05d_01.nxs:/root_spyc_config1d_RIXS_00001/scan_data/data_07"  

         scan_interval    :  [1,2]   # from 1 to 1 included ( kapton_00001_01.nxs)

         Ydim : 16
         Zdim : 16
         Edim : 19

         signalfile : "filename"  # Target file for the signals

 
    #
    """

    roiaddress=None
    roiaddress = mydata["roiaddress"]
  
    filename, groupname = split_hdf5_address( roiaddress)

    file= h5py.File(filename,"r")
    rois = {}
    shape, image=xrs_rois.load_rois_fromh5(file[groupname],rois, retrieveImage = True)
    file.close()

    print( " carico maschere ")
    roiob = xrs_rois.roi_object()
    roiob.load_rois_fromMasksDict(rois ,  newshape = shape, kind="zoom")
    roiob.input_image = image
    
    monitor_divider = 1.0
 

    scan_interval = mydata["scan_interval"]
    todo_list = []
    ninterval = len(scan_interval)//2
    for i in range(ninterval):
        todo_list = todo_list + list(range(    int(scan_interval[2*i]),  int(scan_interval[2*i+1])) ) #   *scan_interval[2*i :2*i+2])


    Ydim = mydata['Ydim']
    Zdim = mydata['Zdim']
    Edim = mydata['Zdim']

    hf = h5py.File(mydata["signalfile"],"w")
    for iE in range(Edim):
        egroup = "E%d/"%iE
        hf.require_group(egroup)
        hf[  egroup+ "image"  ]  = roiob.input_image
        
        for roikey, (origin, box) in roiob.red_rois.items():
            roigroup=roikey+"/"
            hf.require_group(  egroup+roigroup   )
            hf[  egroup+roigroup +"origin"  ] =  origin
            hf[  egroup+roigroup +"mask"    ] =  box
        for iZ in range(    Zdim   ):
            scangroup = "Scan%d/"% iZ
            hf.require_group(egroup+scangroup)
            for roikey, (origin, box) in roiob.red_rois.items():
                roigroup=roikey[3:]+"/"  # with "ROI" removed, nly numerical part

                hf.require_group(egroup+scangroup+roigroup)
                hf.require_dataset(egroup+scangroup+roigroup+"matrix",  [Ydim, box.shape[0], box.shape[1]]    , "f",exact="True")
                hf[egroup+scangroup+roigroup+"cornerpos"] = origin
           
    for iscan in todo_list:
        iZ = (iscan-scan_interval[0]) %  Zdim
        iY = (iscan-scan_interval[0]) // Zdim

        filename, dataname = split_hdf5_address( mydata["expdata"]   % iscan      )
        print(" working on ", filename, dataname )

        
        data =  np.array(h5py.File(filename,"r")[dataname][:])
        
        for roikey, (origin, box) in roiob.red_rois.items():
            
            ##roigroup=roikey+"/"
            roigroup=roikey[3:]+"/"  # with "ROI" removed, nly numerical part

            sliced = data[:, origin[0]:origin[0]+box.shape[0], origin[1]:origin[1]+box.shape[1]] * box

            for iE in range(Edim):
                egroup = "E%d/"%iE
                scangroup = "Scan%d/"% iZ
                print (  " iE  iY, shape", iE, iY, sliced.shape,  hf[ egroup+scangroup+roigroup+"matrix" ].shape )
                hf[ egroup+scangroup+roigroup+"matrix" ][ iY ] = sliced[iE]
    hf.close()



def loadscan_2Dimages_galaxies_foilscan(mydata):
    """
    **loadscan_2Dimages_galaxies_foilscan**

    This command harvest the selected signals.
    the instructions on the scans to be taken must be in the form( as submembers ofload_scans ) ::


      loadscan_2Dimages_galaxies_foilscan :
         roiaddress :  "hdf5filename:nameofroigroup"  # the same given in create_rois
         expdata    :  "kapton_00001_01.nxs:/root_spyc_config1d_RIXS_00001/scan_data/data_07"  

         signalfile : "filename"  # Target file for the signals

         isolateSpot : 0   # if different from zero selects on each image  ( when sumto1d=0 ) ROI  the spot region and sets to zero outside a radius =  isolateSpot
 
    #
    """

    roiaddress=None
    roiaddress = mydata["roiaddress"]
  
    filename, groupname = split_hdf5_address( roiaddress)

    file= h5py.File(filename,"r")
    rois = {}
    shape, image=xrs_rois.load_rois_fromh5(file[groupname],rois,retrieveImage = True)
    file.close()

    print( " carico maschere ")
    roiob = xrs_rois.roi_object()
    roiob.load_rois_fromMasksDict(rois ,  newshape = shape, kind="zoom")
    roiob.input_image = image

    monitor_divider = 1.0
 
        
    if  ('isolateSpot') in mydata :
        isolateSpot = mydata['isolateSpot']
    else:
        isolateSpot  = 0

    hf = h5py.File(mydata["signalfile"],"w")

    hf[   "image"  ]  = roiob.input_image

    for roikey, (origin, box) in roiob.red_rois.items():
        roigroup=roikey+"/"
        hf.require_group(  roigroup   )
        hf[  roigroup +"origin"  ] =  origin
        hf[  roigroup +"mask"    ] =  box
       
        

    filename, dataname = split_hdf5_address( mydata["expdata"]       )

    data =  np.array(h5py.File(filename,"r")[dataname][:])

    for roikey, (origin, box) in roiob.red_rois.items():

        roigroup=roikey[3:]+"/"  # with "ROI" removed, nly numerical part
        sliced = data[:, origin[0]:origin[0]+box.shape[0], origin[1]:origin[1]+box.shape[1]] * box

        if isolateSpot:
            imageLines = np.sum(sliced,axis=1)
            imageLines =imageLines- scipy.ndimage.filters.gaussian_filter( imageLines  ,[0,isolateSpot],mode='constant',cval=0)
            poss       = np.argmax(imageLines,axis=1)
            for i in range(len(poss)):
                sliced[i,:, : max(0,poss[i]-isolateSpot)      ]=0
                sliced[i,:, poss[i]+isolateSpot   :  ]=0
            
        scangroup = "Scan%d/"% 0
            
        hf[ scangroup+roigroup+"matrix" ] = sliced
        hf[ scangroup+roigroup+"cornerpos" ] = origin

    hf.close()



def extract_spectra(mydata):
    """
    **extract_spectra**


    parameters ::
    
      extract_spectra :
        reference_address : "demo_rois.h5:/ROI_AS_SELECTED/energy_scanb/scans/Scan237/"
        sample_address    : "demo_rois.h5:/ROI_AS_SELECTED/images2/scans/"
        roiaddress      : "demo_rois.h5:/ROI_AS_SELECTED/"
        reference_scan             : 237
        scan_interval             : [342,343]
        DE                : 5
        zmargin           : 4
        niterLip          : 100
        niter             : 500
        beta              : 0.0
        target            : "extracted_spectra.h5:/spectra_scan_342"
        final_plot        : "PLOT"  # or "NOPLOT"
    """

    target_filename , target_groupname  = split_hdf5_address( mydata["target"])    


    sample_address = mydata["sample_address"]
    sample_file, sample_groupname  = split_hdf5_address(sample_address)

    roiaddress = mydata["roiaddress"]
    rois_file, rois_groupname  = split_hdf5_address(roiaddress)

    scan_interval = mydata["scan_interval"]
    scans = []
    extratags = []
    
    for i in range(    len(scan_interval)//2     ):
        tok = list(range( int(scan_interval[2*i]),  scan_interval[2*i+1] ))
        scans = scans + tok
        if( isinstance(    scan_interval[2*i]  , int ) )  :
            extratags = extratags +  [0]*len(tok)
        else:
            extratags = extratags + tok 
        
    if "DE" in mydata:    
        DE  = mydata["DE"]
    else:
        DE = 0


    if "discard_threshold" in mydata:
        discard_threshold = mydata["discard_threshold"]
    else:
        discard_threshold =0
        
    if "threshold_fraction" in mydata:
        threshold_fraction = mydata["threshold_fraction"]
    else:
        threshold_fraction =0

    if "zmargin" in mydata:
        zmargin  = mydata["zmargin"]
    else:
        zmargin = 0

    h5frois = h5py.File(rois_file,"r" )
    h5rois  = h5frois[rois_groupname]["rois_definition/rois_dict"]    
    rois_keys_orig = filterRoiList(h5rois.keys(), strip=True)


    references = {}
    niter    = 0
    niterLip = 0
    beta     = 0
    
    slopeInfos = {}
    
    if "flatTriggerer" in mydata:
        rois_keys = rois_keys_orig
        for k,c in enumerate(  rois_keys_orig   ):
            slopeInfos[c] = {           "slope":0.0 ,    "zrate":1.0e38, "estep":1.0 ,   "mask":   h5rois["ROI%02d"%int(c)]["mask"][:] }
            references[c] = None
    elif "slope" in mydata:
        rois_keys = rois_keys_orig
        for k,c in enumerate(  rois_keys_orig   ):
            slopeInfos[c] = {           "slope":mydata["slope"][k]  ,    "zrate":mydata["dh4estep"][k] , "estep":mydata["estep"],   "mask":   h5rois["ROI%02d"%int(c)]["mask"][:] }
            references[c] = None
    else:
        for k,c in enumerate(  rois_keys_orig   ):
            slopeInfos[c] = None
                    
        if "fitted_response" in mydata:
            deconvolve = False
            chiave = "fitted_response"
            optical_response_name = "data"
        else:
            deconvolve = True
            chiave = "reference_address"
            optical_response_name = "optical_response"
        
            if "niter" in mydata:
                niter    = mydata["niter"]
                niterLip = mydata["niterLip"]
                beta     = mydata["beta"]


        if "weight_by_response" in mydata:
            weight_by_response = mydata["weight_by_response"]
        else:
            weight_by_response = 1

                    
        reference_address = mydata[chiave]
        reference_file, reference_groupname  = split_hdf5_address(reference_address)
        if deconvolve:
            nrefscan = mydata["reference_scan"]
            reference_groupname = reference_groupname +"/scans/Scan%03d"%nrefscan

        references = {}
        h5f = h5py.File(reference_file,"r" )
        h5  = h5f[reference_groupname]
        print( " FILTRO ", list(h5.keys())  ) 
        rois_keys = filterRoiList(h5.keys(),prefix="")
        print(reference_file, reference_groupname, list(h5.keys()  )  )
        print("CONFRONTO roiskeys ", rois_keys, rois_keys_orig)
        rois_keys = list(set.intersection(  set(rois_keys), set(rois_keys_orig)       ) )
        

        
        incidentE = None
        if "motorDict/energy" in h5:
            incidentE =  h5["motorDict/energy"][()]
        for k in rois_keys:
            if deconvolve:
                mm = h5[k]["matrix"][:]
            else:
                mm = None
            print(reference_file, reference_groupname, list(h5[k].keys()  )  )
            zscale = h5[k]["xscale"].value*1000.0
            mask = h5rois["ROI%02d"%int(k)]["mask"][:]
            cummask = np.cumsum(mask,axis=0)
            mask[cummask<=zmargin]=0
            mask[(cummask.max(axis=0) -cummask)<zmargin]=0
            oggetto_line = None
            if "line" in h5[k]:
                nref = h5[k]["nref"].value
                rebinned =  h5[k][ optical_response_name  ].value
                s0,s1 = rebinned.shape
                rebinned = np.reshape(rebinned, [  s0//nref, nref,  s1//nref,  nref  ] )//nref
                rebinned = (rebinned.sum(axis=3)).sum(axis=1)
                dico_line = {"line"       : h5[k]["line"].value * np.array([1.0/nref,1.0]),
                             "Xintercept" : h5[k]["Xintercept"].value,
                             "Yintercept" : h5[k]["Yintercept"].value,
                             "Xslope"     : h5[k]["Xslope"].value,
                             "Yslope"     : h5[k]["Yslope"].value,
                             "nref"       : h5[k]["nref"].value,
                             "optical_response" : rebinned,
                             "weight_by_response" : weight_by_response
                         }

                oggetto_line=  type('MyObjectPourLeLignes', (object,), dico_line)
            if deconvolve:
                dico = {"mm":(mm).astype("f")*(mask).astype("f"),"zscale":zscale, "mask":mask,"line_infos":oggetto_line,   "incidentE" : incidentE }
            else:
                dico = {"mm":None,"zscale":zscale, "mask":mask,"line_infos":oggetto_line,   "incidentE" : incidentE }
            references[k] =   type('MyObjectPourDecrireLesReferences', (object,), dico)
        h5f.close()
    
    sample_s = {}
    ene_s    = {}
    h5f = h5py.File(sample_file,"r" )
    h5_sample_group  = h5f[sample_groupname]
    #### rois_keys = list(h5.keys())
    for scan_num, extra  in zip(scans, extratags) :
        sample = {}
        scan_name = "scans/Scan%03d"%scan_num
        print( " FILE was ", sample_file)
        print( "  sample_groupname " , sample_groupname)
        print( "  scan_name " , scan_name)
        
        h5 = h5_sample_group[scan_name]
        scan_energy_0 = h5["motorDict/energy"].value
        print(" ROISKEYS ",  rois_keys ) 
        denominator =  h5[ rois_keys[0]   ]["monitor"].value/(float(h5[   rois_keys[0]    ]["monitor_divider"].value))
        
        for k in rois_keys:
            print( " KKKKK  " , k)
            mm = h5[k]["matrix"][:]
            
            zscale = h5[k]["xscale"][:]*1000
            
            mask = h5rois["ROI%02d"%int(k)]["mask"][:]
            cummask = np.cumsum(mask,axis=0)
            mask[cummask<=zmargin]=0
            mask[(cummask.max(axis=0) -cummask)<zmargin]=0
            
            sample[k] = MyObject = type('MyObject', (object,), {"mm":mm*mask,"zscale":zscale,"denominator": denominator})
            
        sample_s[scan_num] = sample
        ene_s[ (scan_energy_0, extra) ] = ene_s.get((scan_energy_0,extra) ,[]) +[ scan_num ]
        

    myres={}
    enbyscan = {}
    for (ENE0, extra) , myscans in  ene_s.items():
        
        ascan = myscans[0]
        enbyscan[ascan] = ENE0
        mysample_s = { mykey: sample_s[mykey] for mykey in myscans }
        myres[ascan ] = fit_spectra.fit_spectra_main( references  , mysample_s , DE , beta,   niter, niterLip , slopeInfos,  discard_threshold = discard_threshold , threshold_fraction = threshold_fraction  )
        
    filenametxt = target_filename.replace(".h5","") +"_"+  target_groupname.replace("/","_")    +".txt"

    mat4txt = []
    mat4txtheader=""
    
    h5f=h5py.File(target_filename,'a')
    h5f.require_group(target_groupname)
    h5 = h5f [target_groupname]
    listakeys = list(myres[ascan].keys())
    listakeys.sort()

    
    for k in listakeys:
        
        # if str(k) in h5:         piuttosto aggiungo incrementalmente nuovi spettri
        #     del h5[str(k)]
            
        h5k = h5.require_group(str(k))

        plot_string = "fig,ax=pylab.plt.subplots()\n"
        plotfit_string = "fig,ax=pylab.plt.subplots()\n"
        plotconv_string = "fig,ax=pylab.plt.subplots()\n"
        README=""
        for ene0_key in myres.keys():
            # if len(mat4txt)==0:
            #     mat4txt = [ myres[ene0_key][k]["energies"]]
            #     mat4txtheader="# energies "
            # mat4txtheader +=  "roi"+k+"_sp " +    "roi"+k+"_err "
            # mat4txt.append(myres[ene0_key][k]["spectraByLine"])
            # mat4txt.append(myres[ene0_key][k]["errors"])

            envar = "energies_"+str(ene0_key)
            spvar = "spectraByLine_"+str(ene0_key)
            ervar = "errors_"+str(ene0_key)
            E0var = "E0_"+str(ene0_key)
            README += ("""energies_{key} : the energies fo scan {key} 
spectraByLine_{key} : the spectra for scan {key}
error_{key} : the errors fo scan {key}
E0_{key} : the monocromator energy for scan {key}
            """).format(key=str(ene0_key))

            for lab in [envar,spvar,ervar,E0var , "python_plot_spectra_byline" , "README", "python_plot_spectra_byfit" , "python_plot_convergency"  ]:
                if lab in h5k:
                    del h5k[lab]
            

            
            h5k[envar] = myres[ene0_key][k]["energies"]
            h5k[spvar] = myres[ene0_key][k]["spectraByLine"]
            h5k[ervar] = myres[ene0_key][k]["errors"]
            h5k[E0var] = enbyscan[ene0_key]*1000

            plot_string +=("ax.plot(self.%s  - self.%s  , self.spectraByLine_%s,label=\"spectra %f\")\n"
                           "ax.plot(self.%s  - self.%s  , 3*self.%s, label = \"3*sigma %f\")\n"
                           "ax.legend(loc=\"best\")\n") %( envar, E0var , ene0_key, enbyscan[ene0_key],
                                                           envar, E0var , ervar   , enbyscan[ene0_key])


            if niter:
                h5k["spectraByFit_"+str(ene0_key)] = myres[ene0_key][k]["spectraByFit"]
                h5k["fit_errList_"+str(ene0_key)] = myres[ene0_key][k]["fit_errList"]
                h5k["sintesi_"+str(ene0_key)] = myres[ene0_key][k]["sintesi"]

                plotfit_string+="ax.plot(self.energies_%s  -self.E0_%s,self.spectraByFit_%s, label=\"spectra %f\")\n"%(ene0_key,ene0_key,ene0_key,enbyscan[ene0_key])
                plotconv_string+="ax.plot(self.fit_errList%s)\n"%(ene0_key)
                

        plot_string     += "pylab.plt.show()\n"
        plotfit_string  += "pylab.plt.show()\n"
        plotconv_string += "pylab.plt.show()\n"
        

        h5k["python_plot_spectra_byline"] = plot_string
        README+="python_plot_spectra_byline : double click on this to run the code to plot the data \n"
        h5k["README"]=README
        if niter:
            h5k["python_plot_spectra_byfit" ]  = plotfit_string
            h5k["python_plot_convergency"   ]    = plotconv_string
        
        
    h5f.flush()
    h5f.close()
    h5frois.close()

    # mat4txt=np.array(mat4txt).T
    # np.savetxt(filenametxt, mat4txt, header= mat4txtheader)
    
            
    # if "final_plot" in mydata:
    #     if mydata["final_plot"]=="PLOT":
    #         import pylab
    #         mat4txt=np.array(mat4txt).T
    #         nc, npts = mat4txt.shape
    #         npls = (nc-1)/2
            
   
    #         RATIO=1

    #         nlato = int(math.sqrt(npls/1.0/RATIO))
    #         if nlato*nlato*RATIO<npls :
    #             nlato+=1

    #         pylab.plt.figure()

    #         for i in range( npls ) :
    #             x = mat4txt[0]
    #             y = mat4txt[1+2*i]
    #             ye = mat4txt[1+2*i+1]
    #             ax = pylab.plt.subplot(nlato, nlato, i + 1)
                
    #             ax.plot(x,y,label="spectra")
    #             ax.plot(x,3*ye, label = "3*sigma")
    #             ax.legend(loc="best")
    #         pylab.plt.show()
             
    
#det value or default
def gvord(yamldata, key, default=None):
    if  (key) in yamldata :
        return yamldata[key]
    else:
        return default


def h5_assign_force(h5group, name, item):
    if name in h5group:
        del h5group[name]
    h5group[name] = item

def hdf5path_exists(filename, groupname):    
    dn = os.path.dirname(filename)
    if dn=="":
        dn="./"
    os.stat( dn  )
    if not  os.path.exists(filename):
        return False
    if not  h5py.is_hdf5(filename):
        return False
    
    os.system("touch %s" % filename)
    
    h5 = h5py.File(filename,"r" )
    res = groupname in h5
    h5.close()
    return res
    
def create_rois(mydata):
    """
    **create_rois**

    This launches the roi selector.
    If yamlData contains instructions on the scan to be taken , the roi 
    selector will pop up with the image ready to be used. Otherwise you have to
    browse some file around.

    At the end you have the possibility to write the rois on a  hdf5 file.

    In the extreme case when you give no arguments ( parameters) beside the name of 
    the method, the input file looks like this ::

        create_rois :

    you must select an image from the menu or, from the Local menu, 
    select the experimental scans from which an image is summed up.
    The selection will be operated on this image. 

    The experimental scans can be given in input using the following structure ::

      create_rois :

        expdata :  "absolutepathtoaspecfile"  # this points to a spec file
        scans   : [623,624]                   # a list containing one or more scans 
                                              # for the elastic part.
                                              # They will be summed to an image
                                              # and rois will then be extracted from this image.
        roiaddress : "myfile.hdf5:/path/to/hdf5/group"  # the target destination for rois

    The selection will be written if you confirm before exiting ( File menu).
    If roiaddress is not given,  before quitting, and after the validation of the selection,
    you have the possibility to write the rois into a hdf5 file.
    You do this by choosing a filename for the hdf5 file, and a name for the roi selection.
    An entry, with this  name, will be created into the hdf5 file ( all subsequent
    treatements, done by other methods than this, which uses this roi selection 
    ( called by its name )  will be reported into subentries of such entry)
    
    The followings are  optionals arguments.
    First the filter_path ::

        filter_path : filename.h5:/groupname

    It is formed by  hdf5 file AND path  where the filter is stored. It must be a matrix  dataset with 1 for good pixels , 0 for bad ones.
    
    This parameter is used to create a filter mask by hand ::

        masktype : filter

    then the target will be written with a mask of zero and one.
    The mask will have the same size as the image and will be zero for points to be discarded.

    """

    image4roi = None
    layout = None

    if mydata is not None and  ("expdata" in mydata)  :
        repertorio = mydata["expdata"]
        scans = mydata["scans"]
        experiment = xrs_read.read_id20( repertorio ,monitorcolumn='kapraman')
        edfmat =  experiment.read_just_first_scanimage(scans[0])

        if edfmat.shape == (512,768) :
            layout="2X3-12"
        elif edfmat.shape in  [(255,259*5),(255,259*5+1),(255+1,259*5+1) ]  :
            experiment = rixs_read.read_id20(repertorio  ,energyColumn='Anal Energy',monitorColumn='kap4dio', edf_shape=edfmat.shape)
            layout="1X5-1"
        else:
            raise Exception(" Image format not recognised yet. It is : %s (unknow till now)"% str(edfmat.shape))
         
        image4roi =  experiment.SumDirect( scans )

    roiaddress=None
    if mydata is not None and  ("roiaddress" in mydata)  :
        roiaddress = mydata["roiaddress"]
        roiaddress = split_hdf5_address(roiaddress)
    filterMask=None
    if mydata is not None and  ("filter_path" in mydata ) :
        filter_path = mydata["filter_path"]
        if filter_path is not None and len(filter_path):
            filename, dataname = split_hdf5_address(filter_path)
            h5f = h5py.File( filename,"r"   )
            filterMask =  h5f[dataname][:]
    masktype = "ROI"
    if mydata is not None and  ("masktype" in mydata)  :
        masktype = mydata["masktype"]
    app=Qt.QApplication([])
    w4r = roiSelectionWidget.mainwindow(layout=layout)
    if image4roi is not None:
        if filterMask is not None:
            print( image4roi)
            print( filterMask)
            image4roi = image4roi * filterMask
        w4r.showImage( image4roi , xrs_rois.get_geo_informations(  image4roi.shape +(layout,) ))

        if roiaddress is not None:
            if hdf5path_exists(roiaddress[0], roiaddress[1]):
                try:
                    w4r.load_maskDict_from_givenhdf5andgroup(roiaddress[0], roiaddress[1] )
                except:
                    print(  "Problem trying to reload old mask  " , file=sys.stderr  )
            
    w4r.show()
    app.exec_()
    print(" USCITA ", w4r.isOK)
    if not  w4r.isOK:
        sys.stderr.write('ROI CREATION SEEMS TO HAVE BEEN STOPPED BY USER')
        sys.exit(1)
    if w4r.isOK:
        if roiaddress is  None:
            roiaddress =  Qt.QFileDialog.getSaveFileName()

        
            if isinstance(roiaddress, tuple):
                roiaddress  = roiaddress[0]
            
            if roiaddress is not None :
                roiaddress=str(roiaddress)
                if len(roiaddress)==0: roiaddress=None
            if roiaddress is not None :
                roiaddress = split_hdf5_address(roiaddress)
        if roiaddress is not None:
            w4r.saveMaskDictOnH5( roiaddress, masktype , filterMask   ) 



def create_rois_galaxies(mydata):
    """
    **create_rois_galaxies**

    This launches the roi selector.
 
    At the end you have the possibility to write the rois on a  hdf5 file.

    The experimental scans can be given in input using the following structure ::

      create_rois_galaxies :

        expdata :  "myfile.hdf5:/path/to/hdf5/group/dataset"  # this pointsto, inside a hdf5 file, a stack of images z,y,x
        roiaddress : "myfile.hdf5:/path/to/hdf5/group"  # the target destination for rois

    The selection will be written if you confirm before exiting ( File menu).
    If roiaddress is not given,  before quitting, and after the validation of the selection,
    you have the possibility to write the rois into a hdf5 file.
    You do this by choosing a filename for the hdf5 file, and a name for the roi selection.
    An entry, with this  name, will be created into the hdf5 file ( all subsequent
    treatements, done by other methods than this, which uses this roi selection 
    ( called by its name )  will be reported into subentries of such entry)
    
    The followings is   optionals arguments ::

        filter_path : filename.h5:/groupname/dataset

    It is formed by  hdf5 file AND path  where the filter is stored. It must be a matrix  dataset with 1 for good pixels , 0 for bad ones.

    The way to create a good mask file may vary greatly depending on the situation.
    Here an example ( for low signal high noise) ::

                  import h5py
                  import numpy as np
                  fh5 = h5py.File( "kapton_00001_01.nxs" ,"r")
                  mystack = fh5["root_spyc_config1d_RIXS_00001/scan_data/data_07"][:]
                  myimage = mystack.sum(axis=0)
                  mymask = np.less(myimage, 50).astype("f")
                  h5py.File( "mymask.h5" ,"w")["mymask"] = mymask

    """

    image4roi = None
    layout = None


    repertorio = mydata["expdata"]

    filename, dataname = split_hdf5_address(repertorio)


    
    image4roi =  (np.array(h5py.File(filename,"r")[dataname])).sum(axis=0)
    

    if image4roi.shape == (512,768) :
        layout="2X3-12"
    elif image4roi.shape in  [(255,259*5),(255,259*5+1),(255+1,259*5+1) ]  :
        layout="1X5-1"
    elif image4roi.shape in  [(256,256) ]  :
        layout="1X1-4"
    else:
        raise Exception(" Image format not recognised yet. It is : %s (unknow till now)"% str(edfmat.shape))
    
    roiaddress = mydata["roiaddress"]
    roiaddress = split_hdf5_address(roiaddress)
    
    filterMask=None
    if mydata is not None and  ("filter_path" in mydata ) :
        filter_path = mydata["filter_path"]
        if filter_path is not None and len(filter_path):
            filename, dataname = split_hdf5_address(filter_path)
            h5f = h5py.File( filename,"r"   )
            filterMask =  h5f[dataname][:]
            
    masktype = "ROI"
    if mydata is not None and  ("masktype" in mydata)  :
        masktype = mydata["masktype"]
        
    app=Qt.QApplication([])
    w4r = roiSelectionWidget.mainwindow(layout=layout)
    if image4roi is not None:
        if filterMask is not None:
            image4roi = image4roi * filterMask
        w4r.showImage( image4roi , xrs_rois.get_geo_informations(  image4roi.shape +(layout,) ))

        if roiaddress is not None:
            if hdf5path_exists(roiaddress[0], roiaddress[1]):
                try:
                    w4r.load_maskDict_from_givenhdf5andgroup(roiaddress[0], roiaddress[1] )
                except:
                    print(  "Problem trying to reload old mask  " , file=sys.stderr  )
            
    w4r.show()
    app.exec_()
    print(" USCITA ", w4r.isOK)
    if not  w4r.isOK:
        sys.stderr.write('ROI CREATION SEEMS TO HAVE BEEN STOPPED BY USER')
        sys.exit(1)
    if w4r.isOK:
        if roiaddress is  None:
            roiaddress =  Qt.QFileDialog.getSaveFileName()
            roiaddress=str(roiaddress)
            if len(roiaddress)==0: roiaddress=None
            if roiaddress is not None :
                roiaddress = split_hdf5_address(roiaddress)
        if roiaddress is not None:
            w4r.saveMaskDictOnH5( roiaddress, masktype  , filterMask ) 



            
def create_spectral_rois(mydata):
    """
    **create_spectral_rois**

    This launches the spectral roi selector.
    If yamlData contains instructions on the scan to be taken , the roi 
    selector will pop up with the images ready to be used and the ROIs preselection.
    Otherwise you have to select files and scans ("Load Data" menu).

    At the end you have the possibility to write the new rois on a  hdf5 file.

    In the extreme case when you give no argument ( parameters) beside the name of the method

      create_spectral_rois :

        expdata :  /data/id20/inhouse/data/run5_17/run7_ihr/hydra'  # this points to a directory containing a spec file or to the specfile itself
                                                                    # You can simply pass the directory. By default the file hydra will be selected
        scans   : [613,617,621,625]                  # a list containing one or more scans 
                                              # They will be summed up  to form an image to display 
                                              # and analysed in the depth direction for the spectral part to refine
                                              # the rois by NNMA decomposition
        spatial_roiaddress  : "myfile.hdf5:/path/to/hdf5/group"  # a previous ROIs spatial selection
        spectral_roiaddress : "spectral_myfile.hdf5:/path/to/hdf5/group"  # the target destination for the new refined  rois


    If not parameters are given in the input file they can be selected through the "Load Files" menu.

    """
    sf, fn , ns_s  = None, None, None
    spectral_roiaddress = None
    
    if mydata is not None and  ("expdata" in mydata)  :
        sf = mydata["expdata"]
        ns_s = mydata["scans"]
        fn   = mydata["spatial_roiaddress"]

 
    app=Qt.QApplication([])
 
    spectral_w4r = roiNmaSelectionWidget.mainwindow()
    if sf is not None:
        spectral_w4r.load_rois(sf, fn , ns_s )
    spectral_w4r.show()
    
    ####### TO RELOAD A PREVIOUS NNMA SELECTION inside the widget  ########################
    if mydata is not None and  ("spectral_roiaddress" in mydata)  :
        spectral_roiaddress = mydata["spectral_roiaddress"]

        if os.path.exists(spectral_roiaddress) and os.path.isfile(spectral_roiaddress):
            spectral_w4r.loadMaskDictFromH5( mydata["spectral_roiaddress"] )
        
    ###################################################################
    app   .exec_()
    if spectral_w4r.isOK : 
        if spectral_roiaddress is  None:
             spectral_roiaddress =  Qt.QFileDialog.getSaveFileName()
             
             if isinstance(spectral_roiaddress, tuple):
                   spectral_roiaddress = spectral_roiaddress[0]

             spectral_roiaddress=str( spectral_roiaddress)
             if len( spectral_roiaddress)==0:
                 spectral_roiaddress=None
        if  spectral_roiaddress is not None :
            spectral_w4r.saveMaskDictOnH5(spectral_roiaddress  )
        else:
            sys.stderr.write('NO OUTPUT FILE WAS SPECIFIED FOR WRITING SPECTRAL ROI CREATION')
            sys.exit(1)
                        
    else:
        sys.stderr.write('SPECTRAL ROI CREATION SEEMS TO HAVE BEEN STOPPED BY USER')
        sys.exit(1)
            
def Fourc_extraction( yamlData ):
    """ **Fourc_extraction**

    Launches the data extraction from the FOURC spectrometer.

    example ::

        Fourc_extracion :
            active : 1

            data :
                path (str): Absolute path to directory holding the data.
                SPECfname   (str): Name of the SPEC-file ('rixs' is the default).
                EDFprefix   (str): Prefix for the EDF-files ('/edf/' is the default).
                EDFname     (str): Filename of the EDF-files ('rixs_' is the default).
                EDFpostfix  (str): Postfix for the EDF-files ('.edf' is the default).
                en1_column  (str): Counter mnemonic for the energy motor ('anal energy' is the default).
                en2_column  (str): Counter mnemonic for the energy motor ('anal energy' is the default).
                moni_column (str): Mnemonic for the monitor counter ('izero' is the default).
                EinCoor    (list): Coordinates, where to find the incident energy value in the 
                                        SPEC-file (default is [10,1])

            rois :
                scan_number  (int): Scan number (as in SPEC file) of scan to be used for ROI definition.
                zoom_roi (boolean):    Keyword if ROIs should be defined by zooming (default is True).
                poly_roi (boolean): Keyword if ROIs should be defined by selecting polygons (default is False).
                auto_roi (boolean): Keyword if ROIs should be defined automatically (default is False).
                load_address (str): Absolute filename for loading ROIs (HDF5 format).   
                save_address (str): Absolute filename for saving ROIs (HDF5 format).

            elastic :
                scan_number   (int): Scan number (as in SPEC file).
                line_comp (boolean): Keyword, if line-by-line compensation should be used (default is True).
                roi_number    (int): Number of ROI for which to find the compensation factor (default is 0).

            scans :
                scan_numbers (int or list): Integer or list of scan numbers that should be integrated.
                scan_type            (str): Type of scans to be loaded (default is 'inelastic').
                comp_factor        (float): Optional compensation factor to be used (default is None).
                rot_angles          (list): Optional rotation angles (in degrees) for tilted ROIs (default is None).
            saving :
                path     (str): Absolute path to directory where the data should be saved as ASCII files.
                f_name   (str): Base file name for the files to be created.
                post_fix (str): Post-fix for the files to be created (default is '.dat').
                scan_numbers (int or list): Integer or list of integers of scans to be written into the files. 
                header   (str): Optional string defining a header-line for the files (default is '').

    """
    mydata = yamlData #["Fourc_extraction"]
    if mydata is not None and  ("active" in mydata ):
        if mydata["active"]==0:
            return

    if mydata is not None:
        # manage file locations and naming conventions
        try:
            data = mydata["data"]

        except:
            data  = {}

        path        = gvord(data,"path",None)
        SPECfname   = gvord(data,"SPECfname", "rixs")
        EDFprefix   = gvord(data,"EDFprefix", "/edf/")
        EDFname     = gvord(data,"EDFname", "rixs_")
        EDFpostfix  = gvord(data,"EDFpostfix", ".edf")
        moni_column = gvord(data,"moni_column", "izero")
        EinCoor     = gvord(data,"EinCoor", [10,1])

        Fourc_obj = xrs_read.Fourc(path, SPECfname, EDFprefix, EDFname, EDFpostfix, moni_column, EinCoor )

        # manage ROIs
        try:
            rois = mydata["rois"]

        except:
            rois  = {}

        # manage ROIs
        try:
            rois = mydata["rois"]
        except:
            rois  = {}

        scan_number     = gvord( rois, "scan_number", None)
        roi_type        = gvord( rois, "roi_type", 'zoom')
        load_address    = gvord( rois, "load_address", None)
        save_address    = gvord( rois, "save_address", None)
        refinement_key  = gvord( rois, "refinement_key", None)
        refinement_scan = gvord( rois, "refinement_scan", None)

        if load_address:
            roifinder_obj = roifinder_and_gui.roi_finder()
            roifinder_obj.roi_obj.loadH5(load_address)
        else:
            roifinder_obj = roifinder_and_gui.roi_finder()
            image4rois    = Fourc_obj.SumDirect( scan_number )
            if   roi_type == 'zoom':
                roifinder_obj.get_zoom_rois( image4rois )
            elif roi_type == 'poly':
                roifinder_obj.get_poly_rois( image4rois )
            elif roi_type == 'auto':
                roifinder_obj.get_auto_rois( image4rois )
            elif roi_type == 'line':
                roifinder_obj.get_linear_rois( image4rois )
            else:
                print( 'ROI type ' + roi_type + ' not supported. Will end here.' )
                return

        Fourc_obj.set_roiObj( roifinder_obj.roi_obj )

        if refinement_key and refinement_scan:
            print( '>>>>>>> refining ROIs.' )
            # make refinement_key iterable:
            ref_keys = []
            if not isinstance(refinement_key,list):
                ref_keys.append( refinement_key )
            else:
                ref_keys = refinement_key

            # go through keys and do refinements
            for ref_key in ref_keys:
                if ref_key == 'NNMF':
                    roifinder_obj.refine_rois_MF( Fourc_obj, refinement_scan )
                elif ref_key == 'pw':
                    roifinder_obj.refine_rois_PW( Fourc_obj, refinement_scan )
                elif ref_key == 'cw':
                    roifinder_obj.refine_rois_CW( Fourc_obj, refinement_scan )
                else:
                    print ('Unknown refinement keyword, will end here.')
                    return
                Fourc_obj.set_roiObj( roifinder_obj.roi_obj )

        if save_address:
            roifinder_obj.roi_obj.writeH5(save_address)

        # manage elastic line scan
        try:
            elastic = mydata["elastic"]
            print( '>>>>>>> Integrating elastic line scan.' )
            scan_number = gvord( elastic, "scan_number", None )
            roi_number  = gvord( elastic, "roi_number", 0 )
            method      = gvord( elastic, "method", 'sum' )
            Fourc_obj.get_compensation_factor( scan_number, method=method, roi_number=roi_number )
        except:
            elastic  = {}

        # manage scans to be read
        try:
            scans = mydata["scans"]
        except:
            scans  = {}

        print( '>>>>>>> Reading scans.' )
        scan_numbers = gvord(scans, "scan_numbers",None)
        scan_type    = gvord(scans, "scan_type", 'inelastic')
        comp_factor  = gvord(scans, "comp_factor", None)
        rot_angles   = gvord(scans, "rot_angles", None)
        if comp_factor:
            Fourc_obj.load_scan(scan_numbers, direct=True, comp_factor=comp_factor, scan_type=scan_type, rot_angles=rot_angles)
        else:
            Fourc_obj.load_scan(scan_numbers, direct=True, scan_type=scan_type, rot_angles=rot_angles)

        # should there be a spectrum constructed?
        try:
            spectrum = mydata["spectrum"]
            xes = gvord( spectrum, "xes", False )
            if xes:
                Fourc_obj.get_XES_spectrum()
            if rixs:
                print (' NOT implemented yet!' )
        except:
            spectrum = {}

        # manage saving of scans
        try:
            print( '>>>>>>> Saving scans.' )
            saving = mydata["saving"]
            path         = gvord(saving, "path",None)
            f_name       = gvord(saving, "f_name", None)
            post_fix     = gvord(saving, "post_fix", ".dat")
            scan_numbers = gvord(saving, "scan_numbers", None)
            print( '>>>>>>>', scan_numbers, type(scan_numbers))
            header       = gvord(saving, "header", "")
            Fourc_obj.dump_scans_ascii( scan_numbers, path, f_name, post_fix, header )
        except:
            saving  = {}

        print( '>>>>>>> All finished.' )


def Hydra_extraction( yamlData ):
    """ **Hydra_extraction**

    Launches the data extraction from the FOURC spectrometer.

    example ::

        Hydra_extracion :
            active : 1

            data :
                path (str): Absolute path to directory holding the data.
                SPECfname   (str): Name of the SPEC-file ('rixs' is the default).
                EDFprefix   (str): Prefix for the EDF-files ('/edf/' is the default).
                EDFname     (str): Filename of the EDF-files ('rixs_' is the default).
                EDFpostfix  (str): Postfix for the EDF-files ('.edf' is the default).
                en_column   (str): Counter mnemonic for the energy motor ('energy' is the default).
                moni_column (str): Mnemonic for the monitor counter ('izero' is the default).

            rois :
                scan_number  (int): Scan number (as in SPEC file) of scan to be used for ROI definition.
                zoom_roi (boolean):    Keyword if ROIs should be defined by zooming (default is True).
                poly_roi (boolean): Keyword if ROIs should be defined by selecting polygons (default is False).
                auto_roi (boolean): Keyword if ROIs should be defined automatically (default is False).
                load_address (str): Absolute filename for loading ROIs (HDF5 format).   
                save_address (str): Absolute filename for saving ROIs (HDF5 format).

            elastic :
                scan_number   (int): Scan number (as in SPEC file).
                line_comp (boolean): Keyword, if line-by-line compensation should be used (default is True).
                roi_number    (int): Number of ROI for which to find the compensation factor (default is 0).

            scans :
                scan_numbers (int or list): Integer or list of scan numbers that should be integrated.
                scan_type            (str): Type of scans to be loaded (default is 'inelastic').
                comp_factor        (float): Optional compensation factor to be used (default is None).
                
            saving :
                path     (str): Absolute path to directory where the data should be saved as ASCII files.
                f_name   (str): Base file name for the files to be created.
                post_fix (str): Post-fix for the files to be created (default is '.dat').
                scan_numbers (int or list): Integer or list of integers of scans to be written into the files. 
                header   (str): Optional string defining a header-line for the files (default is '').

    """
    mydata = yamlData#["Hydra_extraction"]
    if mydata is not None and  ("active" in mydata) :
        if mydata["active"]==0:
            return

    if mydata is not None:
        # manage file locations and naming conventions
        try:
            data = mydata["data"]
        except:
            data  = {}

        path        = gvord(data,"path",None)
        SPECfname   = gvord(data,"SPECfname", "hydra")
        EDFprefix   = gvord(data,"EDFprefix", "/edf/")
        EDFname     = gvord(data,"EDFname", "hydra_")
        EDFpostfix  = gvord(data,"EDFpostfix", ".edf")
        en_column   = gvord(data,"en_column", "energy")
        moni_column = gvord(data,"moni_column", "izero")

        Hydra_obj = xrs_read.Hydra(path, SPECfname, EDFprefix, EDFname, EDFpostfix, en_column, moni_column )


        # manage ROIs
        try:
            rois = mydata["rois"]
        except:
            rois  = {}

        scan_number  = gvord(rois,"scan_number",None)
        roi_type     = gvord(rois,"roi_type",'zoom')
        load_address = gvord(rois,"load_address",None)
        save_address = gvord(rois,"save_address",None)
        refinement_key = gvord(rois,"refinement_key",None)
        refinement_scan= gvord(rois,"refinement_scan",None)

        if load_address:
            roifinder_obj = roifinder_and_gui.roi_finder()
            roifinder_obj.roi_obj.loadH5(load_address)
        else:
            roifinder_obj = roifinder_and_gui.roi_finder()
            image4rois    = Hydra_obj.SumDirect(scan_number)
            if   roi_type == 'zoom':
                roifinder_obj.get_zoom_rois(image4rois)
            elif roi_type == 'poly':
                roifinder_obj.get_polygon_rois(image4rois)
            elif roi_type == 'auto':
                roifinder_obj.get_auto_rois(image4rois)
            elif roi_type == 'line':
                roifinder_obj.get_linear_rois(image4rois)
            else:
                print( 'ROI type ' + roi_type + ' not supported. Will end here.' )
                return

        Hydra_obj.set_roiObj(roifinder_obj.roi_obj)

        if refinement_key is not None and refinement_scan is not None:
            print( '>>>>>>> refining ROIs.' )
            # make refinement_key iterable:
            ref_keys = []
            if not isinstance(refinement_key,list):
                ref_keys.append(refinement_key)
            else:
                ref_keys = refinement_key

            # go through keys and do refinements
            for ref_key in ref_keys:
                if ref_key == 'NNMF':
                    roifinder_obj.refine_rois_MF( Hydra_obj, refinement_scan )
                elif ref_key == 'pw':
                    roifinder_obj.refine_rois_PW( Hydra_obj, refinement_scan )
                elif ref_key == 'cw':
                    roifinder_obj.refine_rois_CW( Hydra_obj, refinement_scan )
                else:
                    print ('Unknown refinement keyword, will end here.')
                    return
                Hydra_obj.set_roiObj(roifinder_obj.roi_obj)

        if save_address:
            roifinder_obj.roi_obj.writeH5(save_address)

        # manage scans to be read
        try:
            scans = mydata["scans"]
        except:
            scans  = {}

        scan_meth  = gvord(scans, "method", 0)
        print ('>>>>>>>>>>>>>>>> ', type(scan_meth))
        if scan_meth == 0:
            scan_method = 'sum'
        elif scan_meth == 1:
            scan_method = 'pixel'
        elif scan_method == 2:
            scan_mehtod = 'row'
        direct       = gvord(scans, "direct", True)
        scaling      = gvord(scans, "scaling", None)
        scan_numbers = gvord(scans, "scan_numbers", None)
        scan_types   = gvord(scans, "scan_types", 'generic')
        comp_factor  = gvord(scans, "comp_factor", None)
        include_elastic = gvord(scans, "include_elastic", True)

        for scan_num, scan_type in zip(scan_numbers, scan_types):
            print( '>>>>>>>>>>HHHHHH>>>>>>> ', scan_type, type(scan_type))
            Hydra_obj.load_scan(scan_num, scan_type=scan_type, direct=direct, scaling=scaling, method=scan_method )

        Hydra_obj.get_spectrum_new( method=scan_method, include_elastic=include_elastic, abs_counts=False, interpolation=False )

        # manage saving of scans
        try:
            saving = mydata["output"]
        except:
            saving  = {}

        if saving:
            print( '>>>>>>> saving results.' )
            save_format  = gvord(saving, "format", 'ascii')
            file_name    = gvord(saving, "file_name", 'default.dat')
            group_name   = gvord(saving, "group_name", 'spectrum')
            comment      = gvord(saving, "comment", '')

            if save_format == 'ascii':
                Hydra_obj.dump_spectrum_ascii( file_name )
            elif save_format == 'hdf5':
                Hydra_obj.dump_spectrum_hdf5( file_name, group_name, comment=comment )
            else:
                print( 'Format not supported. Will end here.' )
                return

        print( '>>>>>>> All finished.' )


def superR_getVolume_fullfit(mydata):
    """
    Here an example of the input file dedicated section ::

      superR_getVolume_fullfit :
         sample_address : "demo_imaging.hdf5:ROI_B_FIT8/images/scans/"
         delta_address  : "demo_imaging.hdf5:ROI_B_FIT8/scanXX/scans/Scan273/"


         scalprods_address         : "scalprods.hdf5:scal_prods/"
         target_filename            : "volume.hdf5"

         niter                     :  20
         beta                      :  1.0e-8
         eps                       : 0.000002

         ###################################
         # optional

         nbin           : 5                               # defaults to 1
                                      # it will  bin 5 xpixels in one

         roi_keys       :  [60, 64, 35, 69, 34, 24, 5, 6, 71, 70, 39, 58, 56, 33]
         # roi_keys default to all the keys present in delta_address
         optional_solution : /mntdirect/_scisoft/users/mirone/WORKS/Christoph/XRSTools/volumes_gasket.h5:/ref408_bis_423_548/Volume0p03
         ## a     solution with dimensions  [ZDIM,YDIM,XDIM] 
         ## If given, will be used to balance analyzer factors

    The volume will be written in file target_filename( which must not exist already),
    in the datagroup Volume.
    The parameter debin dafaults to [1,1]
    It is used to increase a dimension Z,Y or both , to make it match with X

    The parameters for the Fista optimisation cicle are :
       - niter : the number of fista cycles
       - beta : the factor of the Total Variation  penalisation term
       - eps  : a parameter for the convergence of the Chambolle-Pock TV denoising phase

    """
    

    
    delta_address = mydata["delta_address"]
    delta_filename, delta_groupname = split_hdf5_address(delta_address)
    h5f = h5py.File(delta_filename,"r")
    h5  = h5f[delta_groupname]
    if not  ('nbin' in mydata) :
        width = 1
    else:
        width  = mydata['nbin']
        
    if not  ('optional_solution' in mydata ):
        solution   = None
    else:
        solution_address = str(mydata["optional_solution"])
        print( "solution_address  " , solution_address)
        if solution_address=="None" or solution_address is None or solution_address.strip()=="":
            solution = None
        else:
            solution_filename, solution_groupname = split_hdf5_address(solution_address)    
            solution  =        h5py.File(solution_filename,"r")
            solution  = solution[solution_groupname][:]
    roi_keys  = filterRoiList(h5.keys(),prefix="")
    roi_keys = [str(t) for t in roi_keys]
    if ('roi_keys' in mydata) :
        u_roi_keys  = mydata['roi_keys']
        u_roi_keys = [str(t) for t in roi_keys]
        roi_keys = [ t for t in u_roi_keys if t in roi_keys]
    roi_keys = [str(t) for t in roi_keys]


    ## ##########################
    ## DELTA >>>>>>>>>>>>>>
    sonde = {}  ## The response is added to sonde (probe)
    XDIM = None
    for t in roi_keys:
        m = np.array(h5[t+"/matrix"][:],"d")
        if width != 1 :   # REBINNING
            nbin = width
            assert(nbin>1)
            m=m[:(m.shape[0]//nbin)*nbin].reshape(-1, nbin, m.shape[1],m.shape[2]).sum(axis=1)/nbin
        sonde [t]  =  m   ##  PROBE STUFF GOES HERE
        if XDIM is None:
            XDIM = m.shape[0]
        else:
            assert(XDIM==m.shape[0])
    ##  DELTA <<<<<<<<<<<<<<<<<<<<<
    ## #############################
    h5f.close()

    ## ##########################
    ## >>>>>>>>>>>>>>>>>>> SAMPLE
    ##
    sample_address = mydata["sample_address"]
    sample_filename, sample_groupname = split_hdf5_address(sample_address)
    h5f = h5py.File(sample_filename,"r")
    h5  = h5f[sample_groupname]
    if not  ('scan_interval' in mydata) :
        zscan_keys =  sorted(  filterScanList(h5.keys())     , key = lambda x:  int(''.join(filter(str.isdigit, str(x) ))) )
    else:
        zscan_keys =[ "Scan%03d"%i for i in range(*list(mydata['scan_interval']))]
    ZDIM = len(zscan_keys)

    ## in case of parallelism the workload is splitted between processes
    ## myrange = np.array_split( np.arange( ZDIM ), nprocs   )[myrank]
    
    myrange = np.arange( ZDIM )
    myZDIM =   len(myrange)


    YDIM = None
    rois_to_be_removed=[]
    for ro in roi_keys:
        if ro in h5[zscan_keys[0]]:
            m = h5[zscan_keys[0]][ro]["matrix"][:]        
            if YDIM is not None:
                assert(YDIM == m.shape[0]) ## we take the Y lenght from the first roi of the first scan :
                                           ## this lenght is supposed to be the same are supposed to be the same for all scans
            else:
                YDIM = m.shape[0]
        else:
            rois_to_be_removed.append(ro)
            del sonde[ro]
            del integrated_images[ro]
    for ro in rois_to_be_removed:
        print ( " RIMUOVO ", ro )
        roi_keys.remove(ro)
            
    if YDIM is None:
        return 

    fattori = {} ## This is used for balancing. Will stay to 1.0 for all rois if solution is not given in input
    for i,rk in enumerate(roi_keys):
        fattori[rk] = 1.0
    ## IF solution is given then balancing factors are calculated     
    if solution is not None:
        for rk in roi_keys:
            scal_dd=np.array([0.0],"d")    
            scal_ds = np.array([0.0],"d")
            scal_ss = np.array([0.0],"d")            
            probes      = sonde [rk]
            SS  =  np.tensordot( probes, probes, axes = [  [1,2], [1,2] ] ) 
            for iz in range(myZDIM):
                zkey = zscan_keys[myrange[iz] ]
                m = np.array(h5[ zkey ][ rk ]["matrix"][:],"d")
                msum = m.sum(axis=0)
                probes = sonde [rk]
                assert( probes.shape[1:] == m.shape[1:])
                assert( XDIM == probes.shape[0] )
                assert( YDIM == m.shape[0]      )
                plane_contrib  = np.tensordot( m, probes, axes = [  [1,2], [1,2] ] ) 
                scal_dd     += (m*m).sum()
                
                print( "  zscan_keys " , zscan_keys)
                keypos = zscan_keys.index( zkey)
                scal_ds[:] = scal_ds +  np.tensordot( plane_contrib , solution[keypos], axes = [  [0,1], [0,1] ] )
                scal_ss[:] = scal_ss +  np.tensordot( np.tensordot(SS,solution[keypos],axes=[[1],[1]]) ,   solution[keypos],axes=[[0,1],[1,0]])    
            fattori[rk] = scal_ds/scal_ss
    ### Renormalising the overall strenght of all the factors
    sum = 0.0
    for rk in roi_keys:
        sum = sum + fattori[rk]*fattori[rk]
    for rk in roi_keys:
        fattori[rk] = fattori[rk]/np.sqrt( sum/len(roi_keys) )
    
    scalprods_address = mydata["scalprods_address"]
    scalprods_filename, scalprods_groupname = split_hdf5_address(scalprods_address)
    
    target_address = mydata["target_address"]
    target_filename, target_groupname = split_hdf5_address(target_address)

    niter = mydata["niter"]
    beta = mydata["beta"]
    eps = mydata["eps"]

    if not  ('debin' in mydata) :
        debin = [1,1]
    else:
        debin  = mydata['debin']

    ## moltiplica per fattori sonde
    for rk in roi_keys:
        sonde[rk][:] = sonde[rk] * fattori[rk]

    h5f = h5py.File(scalprods_filename, "r")
    h5_sampledata  = h5f   [scalprods_groupname]
    
    ## modo di impiego
    #   per ogni rk in roi_keys
    #       sonde[rk] e' una serie di risposte Nech * Ny * Nx
    #   per ogni zk in zscan_keys
    #       m = np.array(h5[ zkey ][ rk ]["matrix"][:],"d")
    #          ydim = m.shape[0] 
    #          probes.shape[1:] == m.shape[1:]  ( le immagini)
    #
    #         Zdim = len(  zscan_keys  )
    #         Ydim = m.shape[0]
    #         Xdim = probes.shape[0]
    
    Volume = superr.superr_fullfit(  ZDIM, YDIM, XDIM,
                                     roi_keys, sonde,
                                     zscan_keys,
                                     h5_sampledata, 
                                     niter=niter, beta=beta)
    
    if os.path.exists(target_filename):
        h5 = h5py.File(target_filename,"a")
    else:
        h5 = h5py.File(target_filename,"w")

    print( h5.keys())
    print( target_groupname)
    print( target_groupname in h5)
    if target_groupname in h5:
        del h5[target_groupname]
        
    h5[target_groupname] = Volume
    h5.close()


def superR_getVolume(mydata):
    """
    Here an example of the input file dedicated section ::

      superR_getVolume :
         scalprods_address         : "scalprods.hdf5:scal_prods/"
         target_filename            : "volume.hdf5"
         debin:                    : [2,1]

         niter                     :  20
         beta                      :  1.0e-8
         eps                       : 0.000002

    scalprods_address points to the scalDS, scalDD, scalSS scalar products
    calculated by superR_scal_deltaXimages.
    The volume will be written in file target_filename( which must not exist already),
    in the datagroup Volume.
    The parameter debin dafaults to [1,1]
    It is used to increase a dimension Z,Y or both , to make it match with X

    The parameters for the Fista optimisation cicle are :
       - niter : the number of fista cycles
       - beta : the factor of the Total Variation  penalisation term
       - eps  : a parameter for the convergence of the Chambolle-Pock TV denoising phase

    """

    scalprods_address = mydata["scalprods_address"]
    scalprods_filename, scalprods_groupname = split_hdf5_address(scalprods_address)
    
    target_address = mydata["target_address"]
    target_filename, target_groupname = split_hdf5_address(target_address)

    niter = mydata["niter"]
    beta = mydata["beta"]
    eps = mydata["eps"]

    if not  ('debin' in mydata) :
        debin = [1,1]
    else:
        debin  = mydata['debin']
    
    h5f = h5py.File(scalprods_filename, "r")
    h5  = h5f   [scalprods_groupname]
    
    scalDS = h5["scalDS"][:]
    scalDD = h5["scalDD"][:]
    scalSS = h5["scalSS"][:]
    h5f.close()

    DIMZ,DIMY,DIMX = scalDS.shape

    if debin != [1,1]:
        
        nuovoDS = np.zeros([DIMZ, debin[0], DIMY, debin[1], DIMX ],  "d")
        scalDS.shape =  DIMZ,1,DIMY,1,DIMX
        
        nuovoDS[:] = scalDS
        scalDS = nuovoDS
        scalDS.shape = DIMZ*debin[0], DIMY*debin[1], DIMX
    
    Volume = superr.superr( scalDD, scalDS, scalSS, niter=niter, beta=beta)

    if os.path.exists(target_filename):
        h5 = h5py.File(target_filename,"a")
    else:
        h5 = h5py.File(target_filename,"w")

    print( h5.keys())
    print( target_groupname)
    print( target_groupname in h5)
    if target_groupname in h5:
        del h5[target_groupname]

        
    h5[target_groupname] = Volume
    h5.close()

    
    
def superR_scal_deltaXimages(mydata):
    """ This step supposes that you  have:

      -  already extracted 2D images with the **loadscan_2Dimages** command.
      The **loadscan_2Dimages** has then already accomplished the following requirements
      which are listed below for informative purposes :

          - these images must reside at *sample_address*
          - Under *sample_address* there must be a a set of datagroups 
          with name *ScanZ* where Z is an integer. The number of these datagroups 
          will be called ZDIM
          - Inside each *ScanZ* there must be a a set of datagroup
          with name N where N is the ROI number. 
          - inside each roi datagroup there is the dataset *matrix*.
          This is a three_dimensional array :
             - first dimension is YDIM : the number of steps in the Y direction
             - the other two dimensions are the dimensions of the ROI

      - Obtained the optical PSF of all 
      desired analyzers, and the maxipix response function. This can be done
      with the **iPSF** commands which will have provided the responses
      for a dirac Delta placed at different positions along X direction.
      The **iPSF** has then already  taken care of placing  in the 
      *delta_address* data_group the following(listed for informational purposes):

          - a list of datagroup with name N, N being the number of the ROI.
          - Inside each datagroup a dataset called *matrix* exists
              - the matrix has 3 Dimensions
              - The first dimension is the number for steps done 
              with the thin foil in the X direction to get super-resolution.
              This will be called XDIM
              - The other two dimensions are the dimensions of the ROI.
              They must be equal to those appearing in the the sample datas describe 
              informatively above.
   
    Here an example of the input file dedicated section ::

      superR_scal_deltaXimages :
         sample_address : "demo_imaging.hdf5:ROI_B_FIT8/images/scans/"
         delta_address  : "demo_imaging.hdf5:ROI_B_FIT8/scanXX/scans/Scan273/"
         target_address         : "scalprods.hdf5:scal_prods/"

         ###################################
         # optional

         nbin           : 5                               # defaults to 1
                                      # it will  bin 5 xpixels in one

         roi_keys       :  [60, 64, 35, 69, 34, 24, 5, 6, 71, 70, 39, 58, 56, 33]
         # roi_keys default to all the keys present in delta_address

         orig_delta_address  : "demo_imaging.hdf5:ROI_B/foil_scanXX/scans/Scan273/"
         # defaults to None. If given the integrated image and the average line will be written
         # to check the superposability between the foil scans and the sample scans

         ###
         ## optional
         optional_solution : /mntdirect/_scisoft/users/mirone/WORKS/Christoph/XRSTools/volumes_gasket.h5:/ref408_bis_423_548/Volume0p03
         ## a     solution with dimensions  [ZDIM,YDIM,XDIM] 
         ## If given, will be used to balance analyzer factors

    If nbin is given the dimensios of the superresolution axis, will be reduced or increased,
    by binning together the foil PSFs.
    What the program will produce, under *target_address* datagroup, is 

         -    scalDS  which is  an array  [ZDIM,YDIM,XDIM]  , type "d" .
         -    scalDD  which is the total sum of the squared datas.
         -    scalSS  which is an array [XDIM,XDIM]  , type "d" .

    From these three quantities the volume can be reconstructed with iterative procedure
    in subsequent steps.

    Here what they are :

        - scalSS  is a 2D matrix, one of its elements is  the scalar product of the response function
          for a given position of the foil, along X, with the response function for another position
          of the foil. The sum over ROIS is implicitely done.
        - scalDS is a 3D array. One of its element is the scalar product of the sample image
          for a given Z,Y position of the sample, with the reponse function for a given X position
          of the foil.  The sum over the ROIs is implicitely done.



"""
    
    delta_address = mydata["delta_address"]
    delta_filename, delta_groupname = split_hdf5_address(delta_address)

    if  ('orig_delta_address' in mydata) :
        orig_delta_address = mydata["orig_delta_address"]
        orig_delta_filename, orig_delta_groupname = split_hdf5_address(orig_delta_address)
    else:
        orig_delta_filename, orig_delta_groupname = None, None

    h5f = h5py.File(delta_filename,"r")
    h5  = h5f[delta_groupname]
       
    if not  ('nbin' in mydata) :
        width = 1
    else:
        width  = mydata['nbin']

    if not  ('optional_solution' in mydata ):
        solution   = None
    else:
        solution_address = str(mydata["optional_solution"])
        print( "solution_address  " , solution_address)
        if solution_address=="None" or solution_address is None or solution_address.strip()=="":
            solution = None
        else:
            solution_filename, solution_groupname = split_hdf5_address(solution_address)
            
            solution  =        h5py.File(solution_filename,"r")
            solution  = solution[solution_groupname][:]


    roi_keys  = filterRoiList(h5.keys(),prefix="")
    roi_keys = [str(t) for t in roi_keys]
        
    if ('roi_keys' in mydata) :
        u_roi_keys  = mydata['roi_keys']
        u_roi_keys = [str(t) for t in roi_keys]
        roi_keys = [ t for t in u_roi_keys if t in roi_keys]

    roi_keys = [str(t) for t in roi_keys]

    ## ##########################
    ## DELTA >>>>>>>>>>>>>>
    sonde = {}  ## The response is added to sonde (probe)
    XDIM = None

    ## integrated_images  ::
    ###  much of the code is just for checking the trajectory and compare between the sample measurement
    ###  and the foil measurement.
    ###  integrated_images will be used only if is given in input.
    ###  If you are interested only in the building of the scalar products then just skip everything related to integrated_image
    integrated_images={} ## to monitor if the foil move along a line which superimposes unto the sample integrated image

    for t in roi_keys:
        m = np.array(h5[t+"/matrix"][:],"d")

        integrated_images [t]=[  m.sum(axis=0) ,0, 0, None, None]  # the second entry is left free for the sample integrated image, the third for the orig
        # the 4th the cornerpos, to be initialised by sample data, the fifth cornerpos by origdelta images ( if origdelta is read)

        if width != 1 :   # REBINNING
            nbin = width
            assert(nbin>1)
            m=m[:(m.shape[0]//nbin)*nbin].reshape(-1, nbin, m.shape[1],m.shape[2]).sum(axis=1)/nbin
            
        sonde [t]  =  m   ##  PROBE STUFF GOES HERE
        
        if XDIM is None:
            XDIM = m.shape[0]
        else:
            assert(XDIM==m.shape[0])
    ##  DELTA <<<<<<<<<<<<<<<<<<<<<
    ## #############################
    h5f.close()

    ## ###############################################
    ## >>>>>>>>>>>>>> ORIG  : checking related thing.
    ## The idea is that, originally, the delta files was a large scan resintetised
    ## from a small reference scan. So in this section, that you can skip,
    ## on is taking information about the original scan. If delta and original delta  are the same it's OK
    ## In most cases , if you want to do check, the original delta file is the same file as the delta file
    ## and the interesting part will be comparing the integrated line of the reference and the sample
    if orig_delta_filename is not None:
        h5f = h5py.File(orig_delta_filename,"r")
        h5  = h5f[orig_delta_groupname]
        for t in roi_keys:
            m = np.array(h5[t+"/matrix"][:],"d")
            integrated_images [t][2] += m.sum(axis=0)
            cornerpos = np.array(h5[t+"/cornerpos"][:])
            integrated_images [t][4] =  cornerpos
            
        h5f.close()
    ## ORIG <<<<<<<<<<<<<<<<<<<<<<<<<
    ## #############################
    
    ## ##########################
    ## >>>>>>>>>>>>>>>>>>> SAMPLE
    ##

    sample_address = mydata["sample_address"]
    sample_filename, sample_groupname = split_hdf5_address(sample_address)
    h5f = h5py.File(sample_filename,"r")
    h5  = h5f[sample_groupname]

    if not  ('scan_interval' in mydata) :
        zscan_keys =  sorted(  filterScanList(h5.keys())     , key = lambda x:  int(''.join(filter(str.isdigit, str(x) ))) )
    else:
        zscan_keys =[ "Scan%03d"%i for i in range(*list(mydata['scan_interval']))]
    
    ZDIM = len(zscan_keys)

    ## in case of parallelism the workload is splitted between processes
    myrange = np.array_split( np.arange( ZDIM ), nprocs   )[myrank]

    myZDIM =   len(myrange)


    YDIM = None
    rois_to_be_removed=[]
    for ro in roi_keys:
        if ro in h5[zscan_keys[0]]:
            m = h5[zscan_keys[0]][ro]["matrix"][:]
            
            if YDIM is not None:
                assert(YDIM == m.shape[0]) ## we take the Y lenght from the first roi of the first scan :
                                           ## this lenght is supposed to be the same are supposed to be the same for all scans
            else:
                YDIM = m.shape[0]
        else:
            rois_to_be_removed.append(ro)
            del sonde[ro]
            del integrated_images[ro]

            
    for ro in rois_to_be_removed:
        print ( " RIMUOVO ", ro )
        roi_keys.remove(ro)
            
    if YDIM is None:
        return 

    fattori = {} ## This is used for balancing. Will stay to 1.0 for all rois if solution is not given in input
    for i,rk in enumerate(roi_keys):
        fattori[rk] = 1.0

    ## IF solution is given then balancing factors are calculated     
    if solution is not None:
        for rk in roi_keys:
            ## This are scalars, so why am I using a np.array?
            ##  I am doing that in prevision of mpi AllReduce with mpi4py
            scal_dd=np.array([0.0],"d")    
            scal_ds = np.array([0.0],"d")
            scal_ss = np.array([0.0],"d")
            
            probes      = sonde [rk]
            SS  =  np.tensordot( probes, probes, axes = [  [1,2], [1,2] ] ) 

            for iz in range(myZDIM):
                zkey = zscan_keys[myrange[iz] ]

                m = np.array(h5[ zkey ][ rk ]["matrix"][:],"d")
                msum = m.sum(axis=0)


                probes = sonde [rk]
                assert( probes.shape[1:] == m.shape[1:])
                assert( XDIM == probes.shape[0] )
                assert( YDIM == m.shape[0]      )

                plane_contrib  = np.tensordot( m, probes, axes = [  [1,2], [1,2] ] ) 

                scal_dd     += (m*m).sum()

                print( "  zscan_keys " , zscan_keys)
                keypos = zscan_keys.index( zkey)
                
                scal_ds[:] = scal_ds +  np.tensordot( plane_contrib , solution[keypos], axes = [  [0,1], [0,1] ] )
                scal_ss[:] = scal_ss +  np.tensordot( np.tensordot(SS,solution[keypos],axes=[[1],[1]]) ,   solution[keypos],axes=[[0,1],[1,0]]) 

            if nprocs>1:
                comm.Allreduce([np.array(scal_ss), MPI.DOUBLE],
                            [scal_ss, MPI.DOUBLE],
                            op=MPI.SUM)
                comm.Allreduce([np.array(scal_dd), MPI.DOUBLE],
                            [scal_dd, MPI.DOUBLE],
                            op=MPI.SUM)
                comm.Allreduce([np.array(scal_ds), MPI.DOUBLE],
                               [scal_ds, MPI.DOUBLE],
                               op=MPI.SUM)
                
                comm.Barrier()
            
            fattori[rk] = scal_ds/scal_ss

    ### Renormalising the overall strenght of all the factors        
    sum = 0.0
    for rk in roi_keys:
        sum = sum + fattori[rk]*fattori[rk]
    for rk in roi_keys:
        fattori[rk] = fattori[rk]/np.sqrt( sum/len(roi_keys) )

    ## These arrays below will contain, after summation and for each process,
    ## the contribution of the process's workload
    ## They will be mpi Reduced to get the final result
    
    scalDS = np.zeros( [myZDIM,YDIM,XDIM]  ,"d" )
    scalDD = 0.0
    scalSS = np.zeros( [XDIM,XDIM]  ,"d" )

    for i,rk in enumerate(roi_keys):
        if i%nprocs == myrank:
            probes      = sonde [rk]
            ## Consider that, below, factor is a factor which is applied to the probe to better adapt it to the sample strenght
            ## variations from roi to roi.
            scalSS[:]  =  scalSS[:] +   np.tensordot( probes, probes, axes = [  [1,2], [1,2] ] ) *fattori[rk]*fattori[rk]


    doppio_filtro_done=0
    for iz in range(myZDIM):
        
        ## each scan is at fixed Z and contains many ys. So we proceed along z
        zkey = zscan_keys[myrange[iz] ]
        print( " process %d analyzing scan : "%myrank , zkey)


        my_roi_keys  = filterRoiList(h5[ zkey ].keys(),prefix="")
        my_roi_keys = [str(t) for t in my_roi_keys]


        ## The following piece of code makes sure that our roi_keys list
        ## contains rois which can be found in the datas
        if not doppio_filtro_done:
            new_roi_keys=[]
            for ok in roi_keys:
                if ok not in my_roi_keys:
                    print(" MANCA ", ok , end= {False:"", True:"\n"} [ok==roi_keys[-1]] )
                    del integrated_images[ok]
                else:
                    new_roi_keys.append(ok)
            doppio_filtro_done=1
            roi_keys = new_roi_keys
        
        
        for rk in roi_keys:
            m = np.array(h5[ zkey ][ rk ]["matrix"][:],"d")
            msum = m.sum(axis=0)
            if iz:
                if msum.shape != integrated_images[rk][1].shape:
                    msg =  " ERROR : the yscan elements have different shapes.\n selects homogeneous scans."
                    print( msg)
                    raise Exception( msg)

            integrated_images[rk][1] = integrated_images[rk][1]+msum

            cornerpos = np.array(h5 [ zkey ][ rk ]["cornerpos"][:])
            integrated_images [rk][3] =  cornerpos


            probes = sonde [rk]
            assert( probes.shape[1:] == m.shape[1:])
            assert( XDIM == probes.shape[0] )
            assert( YDIM == m.shape[0]      )


            ## At the end all boils down to these three lines of code. Note that they sum-up
            ## contributions from all the rois alltogether
            plane_contrib  = np.tensordot( m, probes, axes = [  [1,2], [1,2] ] ) 
            scalDS[iz] =  scalDS[iz]+   plane_contrib*fattori[rk]

            scalDD     += (m*m).sum()

    ## sum-Reducing the final result 
    if nprocs>1:
        comm.Reduce([np.array(scalSS), MPI.DOUBLE],
                    [scalSS, MPI.DOUBLE],
                    op=MPI.SUM, root=0)

        comm.Barrier()


                
    h5f.close()

    ##  
    ## ######################


    ## stuffs for checking (integrated_images)
    if nprocs>1:    
        for n in  list(integrated_images.keys()):
            if myrank:
                print( myrank, "A ",  integrated_images[n][0].dtype, integrated_images[n][0].shape)
                comm.Reduce([integrated_images[n][0], MPI.DOUBLE], None, op=MPI.SUM, root=0)
                print( myrank, "B ", integrated_images[n][1].dtype, integrated_images[n][1].shape)
                comm.Reduce([integrated_images[n][1], MPI.DOUBLE], None, op=MPI.SUM, root=0)
            else:
                print( myrank, " A " , integrated_images[n][0].dtype, integrated_images[n][0].shape)
                comm.Reduce(  [integrated_images[n][1], MPI.DOUBLE] , [integrated_images[n][0], MPI.DOUBLE],  op=MPI.SUM, root=0)
                print( myrank, " B " , integrated_images[n][1].dtype, integrated_images[n][1].shape)
                comm.Reduce(  [np.array(integrated_images[n][1]), MPI.DOUBLE], [integrated_images[n][1], MPI.DOUBLE],  op=MPI.SUM, root=0)

                
    ## All the remaining part is just writing            
    target_address = mydata["target_address"]
    target_filename, target_groupname = split_hdf5_address(target_address)

    for iproc in range(nprocs):
        comm.barrier()
        if iproc != myrank:
            continue

        h5f = h5py.File(target_filename,"a")
        
        if myrank==0:
            if h5f.__contains__( target_groupname ):
                del h5f[target_groupname]
                h5f.close()
                h5f = h5py.File(target_filename,"a")
                
            h5f.require_group(target_groupname )
            h5  = h5f[target_groupname]
            h5["scalSS"] = scalSS

            h5.create_dataset("scalDS", ( ZDIM  ,  YDIM  ,   XDIM ), dtype='d')
            h5.create_dataset("scalDD", ( 1, ), dtype='d')
            h5["scalDS"][:]=0
            h5["scalDD"][:]=0
            
            for n in  list(integrated_images.keys()):
                print (" in key " , n)
                B=integrated_images[n][1]
                A=integrated_images[n][0]
                # B=B.sum(axis=0)
                pesiA = A.sum(axis=0)
                pesiB = B.sum(axis=0)
                ## print(" pesi ", pesiA, pesiB)
                medieA = (np.arange(A.shape[0])[:,None]*A).sum(axis=0)/pesiA
                medieB = (np.arange(B.shape[0])[:,None]*B).sum(axis=0)/pesiB

                h5.require_group(n)
                h5n=h5[n]
                h5n["delta_poss"] = medieA
                h5n["sample_poss"] = medieB

                h5n["delta_integrated"  ] = integrated_images[n][0]
                h5n["sample_integrated" ] = integrated_images[n][1]           
                h5n["sample_integrated_weight"   ] = pesiB
                if orig_delta_filename is not None:
                    corner_C = np.array(integrated_images[n][4])
                    corner_B = np.array(integrated_images[n][3])
                    diff     = corner_C-corner_B

                    C     = integrated_images[n][2]
                    pesiC = C.sum(axis=0)
                    medieC = (np.arange(C.shape[0])[:,None]*C).sum(axis=0)/pesiC
                    coords  = np.arange(len( medieC  )) + diff[1]
                    h5n["orig_delta_poss"             ] = np.array(  medieC+diff[0] )
                    h5n["orig_delta_poss_coord"       ] = np.array(  coords )

                    inset = integrated_images[n][2]  
                    tmp = np.zeros_like( integrated_images[n][1]     )
                    target = tmp [ diff[0]:diff[0]+  inset.shape[0], diff[1]:diff[1]+  inset.shape[1]]  
                    target[:] = inset[  :target.shape[0], :target.shape[1] ]
                    h5n["orig_delta_integrated" ] = tmp

        h5f.require_group(target_groupname )
        h5  = h5f[target_groupname]
        
        h5["scalDD"][:]                           +=  scalDD
        h5["scalDS"][myrange[0]:myrange[-1]+1]    +=  scalDS

        h5.require_group("Mean_Poss")
        h5=h5["Mean_Poss"]

        h5f.flush()
        h5f.close()

def XRSprediction( yamlData ):
    """ **prediction**

    This launches the XRS prediction routines.

    If yamlData contains information about: the sample, the incident beam, 
    the analyzer, the detector, the polarization, and the HF compton profiles,
    this will create the desired predicted XRS data.

    At the end you have the possibility to write the predicted profiles into a container hdf5 file.

    In the extreme case when you give no argument ( parameters) ::

        xrs_prediction :

    The following canonical example will be run.

    example ::
    
        xrs_prediction :
            active : 1

            sample :
                chem_formulas : ['C']    # list of strings of chemical sum formulas
                concentrations : [1.0]    # list of concentrations, should contain values between 0.0 and 1.0
                densities : [2.266]        # list of densities of the constituents [g/cm^3]
                angle_tth : 35.0        # scattering angle [deg]
                sample_thickness : 0.1    # sample thickness/diameter in [cm]
                angle_in : None            # incident beam angle in [deg] relative to sample surface normal
                angle_out : None        # beam exit angle in [deg] relatice to sample surface normal
                                        # (negative for transmission geometry)
                shape : 'sphere'        # keyword, can be 'slab' or 'sphere'
                molar_masses : [12.0]    # list of molar masses of all constituents

            incident_beam :
                i0_intensity : 1e13        #  # number of incident photons [1/sec]
                beam_height : 10.0        # in micron
                beam_width : 20.0        # in micron
                
            analyzer :
                material : 'Si'            # analyzer material (e.g. 'Si', 'Ge')
                hkl : [6,6,0]            # [hkl] indices of reflection used
                mask_d : 60.0            # analyzer mask diameter in [mm]
                bend_r : 1.0            # bending radius of the crystal [mm]
                energy_resolution : 0.5 # energy resolution [eV]
                diced : False            # boolean (True or False) if a diced crystal is used or not (defalt is False)
                thickness : 500.0        # thickness of the analyzer crystal
                database_dir : installation_dir

            compton_profiles :
                eloss_range : np.arange(0.0,1000.0,0.1)
                E0 : 9.7

            detector :
                energy : 9.7            # analyzer energy [keV]
                thickness : 500.0        # thickness of the active material [microns]
                material : 'Si'            # detector active material

            thomson :
                scattering_plane : 'vertical'    # keyword to indicate scattering plane relative to lab frame ('vertical' or 'horizontal')
                polarization : 0.99                # degree of polarization (close to 1.0 for undulator radiation)

            saveaddress : "myfile.hdf5:/path/to/hdf5/group"  # the target destination, if data should be saved
   
    """
    mydata = yamlData#["XRSprediction"]

    if mydata is not None and  ("active" in mydata ):
        if mydata["active"]==0:
            return

    if mydata is not None:
        try:
            sample = mydata["sample"]
        except:
            sample = {}
        chem_formulas    = gvord(sample,"chem_formulas",["C"])
        concentrations   = gvord(sample,"concentrations", [1.0])
        densities        = gvord(sample,"densities", [2.266])
        angle_tth        = gvord(sample,"angle_tth", 35.0)
        sample_thickness = gvord(sample,"sample_thickness", 0.1)
        angle_in         = gvord(sample,"angle_in", None)
        angle_out        = gvord(sample,"angle_out", None)
        shape            = gvord(sample,"shape", 'sphere')
        molar_masses     = gvord(sample,"molar_masses", [12.0])
        sample_obj = xrs_prediction.sample(chem_formulas, concentrations, densities, angle_tth, sample_thickness, angle_in, angle_out, shape, molar_masses)
    
        try:
                    incident_beam = mydata["incident_beam"]
        except:
            incident_beam = {}
        i0_intensity  = gvord(incident_beam,"i0_intensity", 1e13)
        beam_height   = gvord(incident_beam,"beam_height", 10.0)
        beam_width    = gvord(incident_beam,"beam_width", 20.0)
        beam_obj = xrs_prediction.beam(i0_intensity, beam_height, beam_width, 0.0)

        try:
            analyzer = mydata["analyzer"]
        except:
            analyzer = {}
        material          = gvord(analyzer,"material", 'Si')
        hkl               = gvord(analyzer,"hkl", [6,6,0])
        mask_d            = gvord(analyzer,"mask_d", 60.0)
        bend_r            = gvord(analyzer,"bend_r", 1.0)
        energy_resolution = gvord(analyzer,"energy_resolution", 0.5)
        diced             = gvord(analyzer,"diced", False)
        thickness         = gvord(analyzer,"thickness", 500.0)


        datadir_default = os.path.join(   os.path.dirname(   __file__ )   ,  "data"     )


        database_dir      = gvord(analyzer,"database_dir",  datadir_default)
        analyzer_obj = xrs_prediction.analyzer(material, hkl, mask_d, bend_r, energy_resolution, diced, thickness, database_dir)

        try:
            detector = mydata["detector"]
        except:
            detector = {}
        energy    = gvord(detector,"energy", 9.7)
        thickness = gvord(detector,"thickness", 500.0)
        material  = gvord(detector,"material", 'Si')
        detector_obj  = xrs_prediction.detector(energy, thickness, material, [256,768])

        try:
            compton_profile = mydata["compton_profile"]
        except:
            compton_profile = {}
        eloss_range = gvord(compton_profile,"eloss_range", np.arange(0.0,1500.0,0.1))
        E0          = gvord(compton_profile,"E0", 9.7)
        compton_profile_obj = xrs_prediction.compton_profiles(sample_obj, eloss_range, E0)

        try:
            thomson = mydata["polarization"]
        except:
            thomson = {}
        scattering_plane = gvord(thomson,"scattering_plane", 'vertical')
        polarization     = gvord(thomson,"polarization", 0.99)
        thomson_obj = xrs_prediction.thomson(compton_profile_obj.get_energy_in_keV(),compton_profile_obj.get_E0(),compton_profile_obj.get_tth())

        abs_cross_section_obj = xrs_prediction.absolute_cross_section(beam_obj, sample_obj, analyzer_obj, detector_obj, thomson_obj, compton_profile_obj)
        abs_cross_section_obj.plot_abs_cross_section()

def XRS_matrix_elements( yamlData ):

    """ **XRS_matrix_elements**

    This calculates transition matrix elements and plots them vs. q.

    The following canonical example will be run.

    example ::
    
        XRS_matrix_elements :
            active : 1

            atom :
                Z : 47                    # atomic number.
                initial_n    : 1        # initial state main quantum number.
                initial_l    : 0        # initial state orbital quantum number.
                final_n      : 2        # final state wave function.
                final_l      : 1        # final state orbital quantum number.
                k            : [1,3,5]    # calculate dipole, octupole, and triacontadipole transition matrix elements.

            plotting : 1

            saving :
                ascii :
                fname :
   
    """

    mydata = yamlData
    if mydata is not None and  ("active" in mydata) :
        if mydata["active"]==0:
            return

    if mydata is not None:
        try:
            data = mydata["atom"]
        except:
            data  = {}

        Z   = gvord(data, "atom",47)
        n_i = gvord(data, "initial_n", 3 )
        l_i = gvord(data, "initial_l", 0 )
        n_f = gvord(data, "final_n", 3 )
        l_f = gvord(data, "final_f", 1 )
        k   = gvord(data, "k", [1,3,5] )

        R1 = xrs_prediction.radial_wave_function()
        R1.load_from_sympy( Z, n_i, l_i )
        R2 = xrs_prediction.radial_wave_function()
        R2.load_from_sympy( Z, n_f, l_f )

        Mel = xrs_prediction.matrix_element(R1, R2)
        Mel.compute(k)
        
    if mydata is not None:
        if mydata["plotting"] == 1:
            import matplotlib.pyplot as plt
            plt.plot(Mel.q/0.5291, Mel.Mel**2)
            plt.xlabel('q [A$^{-1}$]')
            plt.ylabel('squared matrix elements')
            plt.legend(['k = '+str(ii) for ii in k])
            plt.show()

    if mydata is not None:
        try:
            data = mydata["saving"]
        except:
            data  = {}

        ascii = gvord(data, "ascii", False)
        fname = gvord(data, "fname", None)

        if ascii:
            Mel.write_ascii()
        else:
            Mel.write_H5()


def read_reader(mydata, name="dataadress"):

    dataadress = mydata[name]

    filename, groupname =   split_hdf5_address(dataadress)
    reader = xrs_read.read_id20(None)
    reader.load_state_hdf5( filename, groupname)
    return reader, filename, groupname


def superR_recreate_rois(mydata):
    """ 
    This command extend the rois and creates an extrapolated foil scan
    
    The parameters are as following ::
    
       superR_recreate_rois :


         ### we have calculated the responses in responsefilename
         ### and we want to enlarge the scan  by a margin of 3 times
         ### the original scan on the right and on the left 
         ###  ( so for a total of a 7 expansion factor )

         responsefilename :  "responses.h5:/fit"
         nex : 3

         ## the old scan covered by the old rois
         old_scan_address : "../nonregressions/demo_imaging.hdf5:ROI_B/foil_scanXX/scans/Scan273/"


         ## where new rois and bnew scan are written
         target_filename : "newrois.h5:ROI_B_FIT8/"
         filter_rois      : 1
    """
    foil_scan_address = mydata["old_scan_address"]
    foil_filename ,foil_groupname   = split_hdf5_address(foil_scan_address)

    roisgroupname = foil_groupname
    newscanstarget = ""
    for i in range(3):
        pos = roisgroupname.rfind("/")
        newscanstarget = roisgroupname[pos:]+ newscanstarget
        roisgroupname=roisgroupname[:pos]

    responsefilename= mydata["responsefilename"]
    responsefilename, responsepath = split_hdf5_address( responsefilename)
    
    nex = mydata["nex"]

    target_filename , roisgroupname_target= split_hdf5_address( mydata["target_filename"])
    
    newscanstarget = newscanstarget[1:]

    
    
    # if os.path.exists(target_filename):
    #     sys.stdout.write("Error : file %s exists already. Remove it yourself\n"%target_filename)
    #     sys.stderr.write("Error : file %s exists already. Remove it yourself\n"%target_filename)
    #     sys.exit(1)
        
    if  ("filter_rois" in mydata)  :
        filter_rois  =  mydata["filter_rois"]
    else:
        filter_rois      = 1
    


    print( "LANCIO ")
    dic = {"filename" : foil_filename  , "groupname" : foil_groupname,
                                  "roisgroupname" : roisgroupname,
                                  "target_filename":target_filename, 
                                  "roisgroupname_target" : roisgroupname_target ,
                                  "newscanstarget"    : newscanstarget,
                                  "responsefilename" :  responsefilename,
                                  "responsepath"     : responsepath, 
                                  "nex" : nex, "filter_rois":filter_rois}
    print( dic)

    if  "recenterings_refined" in mydata :
        recenterings_refined = mydata["recenterings_refined"]

        recenterings_filename, recenterings_groupname = split_hdf5_address( recenterings_refined )
        h5f = h5py.File(recenterings_filename,"r")
        h5 = h5f[recenterings_groupname]
        recenterings= {}
        chiavi = filterRoiList(h5.keys())
        for c in chiavi:
            recenterings[int(c)]= h5[c][:]
            assert(recenterings[int(c)].shape == (2,))
        h5f.close()
    else:
        recenterings= None



    filterMask=None
    if mydata is not None and  ("filter_path" in mydata ) :
        filter_path = mydata["filter_path"]
        if filter_path is not None and len(filter_path):
            filename, dataname = split_hdf5_address(filter_path)
            h5f = h5py.File( filename,"r"   )
            filterMask =  h5f[dataname][:]


        
    
    reponse_percussionelle.DOROIS(filename = foil_filename  , groupname = foil_groupname,
                                  roisgroupname = roisgroupname,
                                  target_filename=target_filename, 
                                  roisgroupname_target = roisgroupname_target ,
                                  newscanstarget    = newscanstarget,
                                  responsefilename =  responsefilename,
                                  responsepath     = responsepath, 
                                  nex = nex, filter_rois=filter_rois, recenterings=  recenterings ,
                                  filterMask = filterMask)


    
def superR_fit_responses(mydata):
    """
    superR_fit_responses :
       foil_scan_address : "demo_foilroi.h5:/ROI_FOIL/foil_scan/scans/Scan273"
       nref : 5                 # the number of subdivision per pixel dimension used to 
                                # represent the optical response function at higher resolution
       niter_optical  :  100    # the number of iterations used in the optimisation of the optical
                                # response
       beta_optical  :  0.1     # The L1 norm factor in the regularisation 
                                #  term for the optical functions
       pixel_dim : 6            # The pixel response function is represented with a 
                                #  pixel_dim**2 array
       niter_pixel : 100        # The number of iterations in the pixel response optimisation
                                # phase. A negative number stands for ISTA, positive for FISTA
       beta_pixel  :  1000.0    # L1 factor for the pixel response regularisation

       ## The used trajectories are always written whith the calculated response 
       ## They can be reloaded and used as initialization(and freezed with do_refine_trajectory : 0 )
       ## Uncomment the following line if you want to reload a set of trajectories
       ## without this options trajectories are initialised from the spots drifts
       ##
       #   reload_trajectories_file : "response.h5"

       ######
       ## The method first find an estimation of the foil scan trajectory on each roi
       ## then, based on this, obtain a fit of the optical response function
       ## assuming a flat pixel response. Finally the pixel response is optimised
       ##
       ## There is a final phase where a global optimisation
       ## is done in niter_global steps.
       ##
       ## Each step is composed of optical response fit, followed by a pixel response fit.
       ## If do_refine_trajectory is different from zero, the trajectory is reoptimised at each step
       ## 
       niter_global  :  20

       ## if do_refine_trajectory=1 the start and end point of the trajectory are free
       ##  if =2 then the start and end point are forced to a trajectory which is obtained
       ##  from a reference scan : the foil scan may be short, then one can use the scan of
       ##   an object to get another one : key *trajectory_reference_scan_address*
       ##

       do_refine_trajectory : 2

       ## optional: only if do_refine_trajectory = 2

       trajectory_reference_scansequence_address : "demo_newrois.h5:/ROI_FOIL/images/scans/"
       trajectory_threshold   : 0.1

       ## if the pixel response function is forced to be symmetrical 

       simmetrizza : 1

       ## where the found responses are written

       target_file : "demo_responses_bis.h5"

    """ 

    if "foil_scan_address" in mydata:
        foil_scan_address = mydata["foil_scan_address"]
    else:
        ref_scan_number =  mydata["ref_scan_number"]
        foil_scan_address = mydata["response_scan_address"]+"/scans/Scan"+ ("%03d"% ref_scan_number)
        
    foil_filename ,foil_groupname   = split_hdf5_address(foil_scan_address)

    nref = mydata["nref"]
    niter_optical = mydata["niter_optical"]
    beta_optical = mydata["beta_optical"]
    beta_pixel   = mydata["beta_pixel"]
    niter_pixel = mydata["niter_pixel"]
    niter_global = mydata["niter_global"]
    pixel_dim = mydata["pixel_dim"]
    simmetrizza = mydata["simmetrizza"]
    do_refine_trajectory = mydata["do_refine_trajectory"]
    target_file = mydata["target_file"]

    target_file, target_groupname = split_hdf5_address( target_file )

    
    # if os.path.exists(target_file):
    #     sys.stdout.write("Error : file %s exists already. Remove it yourself\n"%target_file)
    #     sys.stderr.write("Error : file %s exists already. Remove it yourself\n"%target_file)
    #     sys.exit(1)
        
    
    if  do_refine_trajectory==2 :
        ## optional: only if do_refine_trajectory = 2
        trajectory_reference_scansequence_address = mydata["trajectory_reference_scansequence_address"]
        trajectory_reference_scansequence_filename , trajectory_reference_scansequence_groupname = split_hdf5_address(trajectory_reference_scansequence_address)
        trajectory_threshold = mydata["trajectory_threshold"]
        
    else:
        trajectory_reference_scansequence_filename , trajectory_reference_scansequence_groupname  = None, None
        trajectory_threshold =0


    if  ("reload_trajectories_file" in mydata)  :
        trajectory_file =  mydata["reload_trajectories_file"]
    else:
        trajectory_file       = None

    if  ("filter_rois" in mydata ) :
        filter_rois  =  mydata["filter_rois"]
    else:
        filter_rois      = 1

    if  ("fit_lines" in mydata)  :
        fit_lines  =  mydata["fit_lines"]
    else:
        fit_lines  = 0



        
    reponse_percussionelle.DOFIT(filename=foil_filename, groupname=foil_groupname, nref=nref, niter_optical=niter_optical, beta_optical=beta_optical ,
                                 beta_pixel=beta_pixel, niter_pixel = niter_pixel,
                                 niter_global  = niter_global, pixel_dim=pixel_dim, simmetrizza=simmetrizza,
                                 do_refine_trajectory=do_refine_trajectory, target_file=target_file, target_groupname = target_groupname, 
                                 trajectory_reference_scansequence_filename =  trajectory_reference_scansequence_filename ,
                                 trajectory_reference_scansequence_groupname = trajectory_reference_scansequence_groupname ,
                                 trajectory_threshold = trajectory_threshold, trajectory_file = trajectory_file, filter_rois=filter_rois, fit_lines = fit_lines)
swissknife_operations={

    "help"                           :  help,
    "create_rois"                    :  create_rois,
    "create_rois_galaxies"                    :  create_rois_galaxies,
    "loadscan_2Dimages_galaxies" : loadscan_2Dimages_galaxies,
    "loadscan_2Dimages_galaxies_foilscan" : loadscan_2Dimages_galaxies_foilscan,
    "create_spectral_rois"           :  create_spectral_rois,
    "load_scans"                     :  load_scans,
    "HFspectrum"                     :  HFspectrum,
    "Extraction"                     :  Extraction,
    "loadscan_2Dimages"              :  loadscan_2Dimages,
    "volume_from_2Dimages"           :  volume_from_2Dimages,
    "view_Volume_myavi"              :  view_Volume_myavi,
    "superR_scal_deltaXimages"       :  superR_scal_deltaXimages,
    "superR_fit_responses"           :  superR_fit_responses,
    "superR_recreate_rois"           :  superR_recreate_rois,
    "superR_getVolume"               :  superR_getVolume,
    "calculate_recenterings"         :  calculate_recenterings,
    "extract_spectra"                :  extract_spectra,
    "sum_scans2maps"                 :  sum_scans2maps,
    "XRSprediction" : XRSprediction,
    "Fourc_extraction"     : Fourc_extraction,
    "Hydra_extraction"     : Hydra_extraction,
    "XRS_matrix_elements"  : XRS_matrix_elements
}

parallelised_operations = [ "loadscan_2Dimages" , "superR_scal_deltaXimages" , "superR_fit_responses" ,  "superR_recreate_rois"     ]


if __name__=="__main__":
    main()








