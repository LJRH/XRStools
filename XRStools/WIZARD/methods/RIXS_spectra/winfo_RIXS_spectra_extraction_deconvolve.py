from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import collections



######## Managing versioning postpending to project name ################
# Get the name of the launching script,
# the launching script (an principal module)
# is post pended by version name
# We use this to retrieve the version
import inspect
import os
import importlib

curframe = inspect.currentframe()
calframe = inspect.getouterframes(curframe, 2)


global version
caller   = calframe[-1][1]
if caller[-3:]==".py":
    print( " caller :" , caller , " //  " ,  os.path.dirname(os.path.dirname(caller)))
    # imported by Wizard.py to see what are the subtasks
    version = os.path.basename(os.path.dirname(os.path.dirname(caller)))[len("XRStools"):]
else:
    # called by the XRS_wizard script
    caller   = os.path.basename(caller)
    version = caller[len("XRS_wizard"):]
print( " Caller version==", version)


XRStools = importlib.import_module( "XRStools"+version    )

exec("from XRStools%s.WIZARD.Wizard import *"%version)


###################################################################


def simplecopy(par):
    if hasattr(par,"value" ):
        res = copy.deepcopy(par.value)
    else:
        res = copy.deepcopy(par)
    return res

class functor_sigterm_handler:
    def __init__(self, p):
        self.p = p
    def __call__(self, *x):
        print( x)
        print( " KILLO " , self.p.pid)
        os.kill(self.p.pid, signal.SIGKILL)

def swissknife_runner( yamltext, where ):
    inputname = os.path.join(where,"input.yaml")
    stdname = os.path.join(where,"stdout.txt")
    errname = os.path.join(where,"stderr.txt")
    open(inputname,"w").write(yamltext)
    # os.system("XRS_swissknife %s 1> %s  2>%s"%(inputname,stdname,  errname   ) )
    #p = subprocess.call(  string.split("XRS_swissknife %s 1> %s  2> %s"%(inputname,stdname,  errname ))    , shell=False, stdout= subprocess .PIPE ) #  preexec_fn=os.setsid
    # p= subprocess.Popen(string.split( "XRS_swissknife %s 1> %s  2>%s"%(inputname,stdname,  errname ))    , shell=False,  stdout= subprocess .PIPE ,  stderr= subprocess .PIPE  )
    stdout = open(stdname,"w")
    stderr = open(errname,"w")
    p= subprocess.Popen(( "XRS_swissknife%s %s 1> %s  2>%s"%(version,inputname,stdname,  errname )).split()    , shell=False,  stdout= stdout ,  stderr= stderr )
    # p= subprocess.Popen( "XRS_swissknife %s 1> %s  2>%s"%(inputname,stdname,  errname )    , shell=True )
    signal.signal(signal.SIGTERM, functor_sigterm_handler(p))
    p.communicate()
    # import time
    # time.sleep(100)
    # p.communicate()

class Functor_adapt:
    def __init__(self, dico):
        self.dico=dico
    def __call__(self):
        chiavi = list(self.dico.keys())[1:]
        for c in chiavi :

            subdic = self.dico[c]
            for name,par in subdic.items():
                if isinstance(par,Parametro):
                    if par.defaults2 is not None:
                        source, method = par.defaults2

                        if hasattr(source,"value" ):
                            par.value = method(source.value)
                        else:

                            par.value = method(source)

        return self.dico["CREATION_PARS"]["reference_address"].value+"__"+   self.dico["CREATION_PARS"]["scans"] .value

def dic2yaml_extraction_scan_with_recentering(dicodic):
    s = "loadscan_2Dimages :\n"
    s = s+ "  expdata : %s \n" % dicodic["expdata"].render() 
    s = s+ "  roiaddress : %s \n" % dicodic["roiaddress"].render() 
    s = s+ "  scan_interval : %s \n" % dicodic["scans"].render() 
    s = s+ "  monitor_column : %s \n" % dicodic["monitor_column"].render() 
    ## s = s+ "  recenterings  : %s \n" % dicodic["recenterings"].render()
    s = s+ "  signaladdress : %s \n" % dicodic["signaladdress"].render()
    s = s+"""
  sumto1D  : 0
  energycolumn : 'Anal Energy'
  monitorcolumn : %s
""" %  dicodic["monitor_column"].render()
    return s

def dic2yaml_generic(dicodic):
    s = "%s :\n"%dicodic["method_name"]
    for n,t in dicodic.items():
        if isinstance(t, Parametro):
            s = s +" %s : %s\n"%(n,t.render())
    return s

        
def getMethod(options):
    
    metodo = collections.OrderedDict()

    c_p = collections.OrderedDict()
    metodo["CREATION_PARS"] = c_p

    c_p["spec_file"]      = FilePath(doc="The path to the specfile\n file must exist", fileme=1)

    # c_p["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
    #                                     " where the ROI for calibration extraction can be found  \n"+
    #                                     "  Format toto.hdf5:datagroupname  ",
    #                                     fileme=1, gnameme=1)

    c_p["scans"]   = many_Number4SpecScan(  c_p["spec_file"] ,isinterval=1, doc=
                                            "This is the list of two number (**Use Python syntax [start,end+1]**)\n"+
                                            " defining the interval scan numbers containing data.\n"+
                                            "You can also give more intervals in the same list.\n"+
                                            # " Must containing two numbers : [ start,end+1] \n"+
                                            "",   # "Will also be used for ROI selection",
                                            nmin=2,nmax=200)

    if 0:
        c_p["recenterings"]      = Hdf5FilePath(doc=("The  hdf5 file where recentering informations can be found\n"
                                                     "Better if you give a refined one so that it will be the same\n"
                                                     "for all spectra (will not be refined each time" ),
                                                fileme=1, gnameme=1)

    c_p["signaladdress"]   =  Hdf5FilePath(doc=("The  hdf5 file where extracted data will be written\n"
                                                     "Must be a free name\n"), fileme=0.5, gnameme=0)

    c_p["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                  initial_value = "kapraman/1000"
                                              )

    c_p["reference_address"] =Hdf5FilePath(doc="The full path where the references can be found\n It is a groupname ending by ScanXXX contained in \n"+
                                           " group called scans",
                                           fileme=1, gnameme=1)
    # c_p["reference_scan"]   = aNumber(   doc="This is the scan number containing references.\n") 

    c_p["target"]            =  Hdf5FilePath(doc="The hdf5 full path on which spectra will be written. Must NOT exist   ",
                                             fileme=0.5, gnameme=0)

    c_p["signaladdress"]    .isresult = 1
    c_p["target"]    .isresult = 1


    c_p["define_new"] = Functor_adapt( metodo )



    ## ==============================================================================
    es_moved  = collections.OrderedDict()
    metodo["EXTRACTION_DATA_SAMPLE_CENTERED"] = es_moved


    es_moved["getyaml"]  =   dic2yaml_extraction_scan_with_recentering

    es_moved["HELP"] = " Recollect the sample data, with\n the rois centered on calibration scan\n" ## , and with the recentering  displacement"

    es_moved["method_name"]  =  "loadscan_2Dimages"
    es_moved["expdata"]      =  FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )

    es_moved["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                             " where the ROI for calibration extraction can be found  \n"+
                                             "  Format toto.hdf5:datagroupname  ",
                                             defaults2 = ( c_p["reference_address"] , Functor_add_datagroup("", remplace_last= 3)  ), 
                                             fileme=1, gnameme=1)

    es_moved["scans"]   = many_Number4SpecScan(  es_moved["expdata"]  ,isinterval=1, doc=
                                                 "This is the list of two number (**Use Python syntax [start,end+1]**)\n"+
                                                 " defining the interval scan numbers containing data.\n"+
                                                 "You can also give more intervals in the same list.\n"+
                                                 # " Must containing two numbers : [ start,end+1] \n"+
                                                 "",   # "Will also be used for ROI selection",
                                                 nmin=2,nmax=200,
                                                 defaults2 = ( c_p["scans"]   , simplecopy))

    if 0:
        es_moved["recenterings"]      = Hdf5FilePath(doc=("The  hdf5 file where recentering informations can be found\n"
                                                          "Better if you give a refined one so that it will be the same\n"
                                                          "for all spectra (will not be refined each time" ),
                                                     defaults2 = ( c_p["recenterings"], simplecopy), fileme=1, gnameme=1)

    es_moved["signaladdress"]   =  Hdf5FilePath(doc=("The  hdf5 file where extracted data will be written\n"
                                                     "Must be a free name\n"),
                                                defaults2 = (  ( c_p["signaladdress"],
                                                                 es_moved["scans"]   ) ,
                                                               Functor_append_numinterval_to_name() )
                                                , fileme=0.5, gnameme=0)

    es_moved["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                       defaults2=(c_p["monitor_column"], simplecopy)
                                                   )

    es_moved["signaladdress"]  .isresult = 1
    es_moved["swissknife_runner"]  =  swissknife_runner  



    # ==============================================================================
    extract_spectra  = collections.OrderedDict()
    metodo["EXTRACT_SPECTRA"] = extract_spectra


    extract_spectra["getyaml"]  =   dic2yaml_generic

    extract_spectra["HELP"] = " Extract the spectra by fitting\nthe data with superpositions\ of references"

    extract_spectra["method_name"]  =  "extract_spectra"


    extract_spectra["reference_address"] =Hdf5FilePath(doc="Where the references can be foun.d\nIn general a subgroup of ROI_something\n"+
                                                       "containing a group called scans.\n The precise scan will be indicated by reference_scan below",
                                                       defaults2 = (  c_p["reference_address"]   , Functor_add_datagroup("", remplace_last= 2)  ),
                                                       fileme=1, gnameme=1)


    extract_spectra["sample_address"] =Hdf5FilePath(doc="Where the sample data  can be found.\nIn general a subgroup of ROI_something\n"+
                                                       "containing a goup called scans",
                                                       defaults2 = ( es_moved["signaladdress"],simplecopy  ),
                                                       fileme=1, gnameme=1)

    extract_spectra["roiaddress"] =Hdf5FilePath(doc="Where the roi can be found\nA datagroup containing rois_definition",
                                                defaults2 =  (es_moved["roiaddress"], simplecopy),
                                                fileme=1, gnameme=1)


    extract_spectra["reference_scan"]   = aNumber(   doc="This is the scan number containing references.\n",
                                                     defaults2 = (c_p["reference_address"], Functor_add_datagroup("", extract_numeric= 1)   )) 

    extract_spectra["scan_interval"]   = many_Number4SpecScan(  es_moved["expdata"]  , isinterval=1,  doc=
                                                        "This is the list of two number (**Use Python syntax [start,end+1 ,start2,end2+1  ]**)\n"+
                                                                " defining the interval scan numbers containing data.\n"+
                                                                "Will also be used for ROI selection\n",
                                                 nmin=2,nmax=200,
                                                 defaults2 = ( es_moved["scans"]   , simplecopy))

    extract_spectra["DE"]   = aNumberFloat(  doc="The energy stem in eV for the resulting spectra",
                                             defaults2 = ("0.02", simplecopy) )

    extract_spectra["zmargin"]   = aNumber(  doc="The used ROIs will be reduced by this number of pixels vertically, at both top and bottom of the ROI",
                                             defaults2 = ("4", simplecopy) )

    extract_spectra["niterLip"]   = aNumber(  doc="The number of power method iterations used for estimation of the Lipschitz factor",
                                              defaults2 = ("20", simplecopy) )

    extract_spectra["niter"]   = aNumber(  doc="The number of Fista iterations",
                                              defaults2 = ("50", simplecopy) )

    extract_spectra["beta"]   = aNumberFloat(  doc="Optionally you can use this L1 penalization factor if you need to regularize the problem",
                                             defaults2 = ("0.0", simplecopy) )

    extract_spectra["target"]   =  Hdf5FilePath(doc="The hdf5 full path on which spectra will be written. Must NOT exist   ",

                                         defaults2 = (  ("extracted_spectra.h5:/spectra_scan", extract_spectra["scan_interval"]   ) , Functor_append_numinterval_to_name() ),


                                         fileme=0.5, gnameme=0)

    extract_spectra["target"].isresult = 1
    extract_spectra["swissknife_runner"]  =  swissknife_runner  

    return metodo 







