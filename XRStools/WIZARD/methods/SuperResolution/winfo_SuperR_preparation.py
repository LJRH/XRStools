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


import yaml 


def dic2yaml_generic(dicodic):
    s = "%s :\n"%dicodic["method_name"]
    for n,t in dicodic.items():
        if isinstance(t, Parametro):
            s = s +" %s : %s\n"%(n,t.render())
    return s

class functor_sigterm_handler:
    def __init__(self, p):
        self.p = p
    def __call__(self, *x):
        print( x)
        print( " KILLO " , self.p.pid)
        os.kill(self.p.pid, signal.SIGKILL)

def swissknife_runner( yamltext, where, ret_dico ):

    mydata =  yaml.load(yamltext)
    mname, mydata =  list(mydata.items())[0]
    print( mydata)
    if "MPI_N_PROCS" in mydata:
        MPI_N_PROCS = mydata["MPI_N_PROCS"]
    else:
        MPI_N_PROCS = 1
    print( " MPI_N_PROCS  " ,  MPI_N_PROCS)
        
    inputname = os.path.join(where,"input.yaml")
    stdname = os.path.join(where,"stdout.txt")
    errname = os.path.join(where,"stderr.txt")

    ret_dico["input"] =inputname 
    ret_dico["stdout"]=stdname   
    ret_dico["stderr"]=errname   


    
    open(inputname,"w").write(yamltext)
    stdout = open(stdname,"w")
    stderr = open(errname,"w")
    if MPI_N_PROCS==1:
        p= subprocess.Popen(( "XRS_swissknife%s %s 1> %s  2>%s"%(version,inputname,stdname,  errname )).split()    , shell=False,  stdout= stdout ,  stderr= stderr )
    else:
        p= subprocess.Popen(( "mpirun -n %d XRS_swissknife%s %s 1> %s  2>%s"%(version, MPI_N_PROCS  , inputname,stdname,  errname )).split()    , shell=False,  stdout= stdout ,  stderr= stderr )
    signal.signal(signal.SIGTERM, functor_sigterm_handler(p))
    p.communicate()
    ret_dico["return_code"]=p.returncode
    
    
class Functor_adapt:
    def __init__(self, dico):
        self.dico=dico
    def __call__(self):
        chiavi = list(self.dico.keys())[1:]
        for c in chiavi :
            
            subdic = self.dico[c]
            if not isinstance(subdic,dict):
                continue
            for name,par in subdic.items():
                if isinstance(par,Parametro):
                    if par.defaults2 is not None:
                        source, method = par.defaults2
                        print( source)
                        par.value = method(source)
                   
        return self.dico["CREATION_PARS"]["result_file"].value


def dic2yaml_create_rois(dicodic):
    s = "create_rois:\n"
    s = s+ "   expdata : %s \n" % dicodic["expdata"].render()
    print( dicodic["scans"].render)
    print( dicodic["scans"].rendering)
    print( dicodic["scans"].getValue())
    print( dicodic["scans"].rendering(dicodic["scans"].getValue()))
    s = s+ "   scans : %s \n" % dicodic["scans"].render ()
    s = s+ "   roiaddress : %s \n" % dicodic["result_file"].render()
    tok = dicodic["filter_path"].render()
    if tok !="None":
        s = s+ "   filter_path : %s \n" % tok
    
    # s = s+ "   masktype   : filter \n"
    return s


def dic2yaml_extraction_scan(dicodic):
    s = "loadscan_2Dimages :\n"
    s = s+ "  expdata : %s \n" % dicodic["expdata"].render() 
    s = s+ "  roiaddress : %s \n" % dicodic["roiaddress"].render() 
    s = s+ "  monitor_column : %s \n" % dicodic["monitor_column"].render() 
    s = s+ "  scan_interval : %s \n" % dicodic["scans"].render() 
    s = s+ "  signaladdress : %s \n" % dicodic["signaladdress"].render()
    s = s+ "  isolateSpot : %s \n" % dicodic["isolateSpot"].render()
    s = s+"""
  sumto1D  : 0
  energycolumn : 'stx'
  monitorcolumn : %s

""" % dicodic["monitor_column"].render()
    
    return s  

class Functor_Compose_A_Scan_Address:
    def __init__(self, take_first=0):
        self.take_first = take_first
    def __call__(self,par):
        a,b = par
        print( " ====================== " )
        print( a,b)

        va=a.getFullValue()
        vb=b.getValue()
        print( va, vb)
        if isinstance(vb,int):
            return va+"/scans/Scan%03d/"%vb
        elif isinstance(vb,list):
            return va+"/scans/Scan%03d/"%vb[0]
        else:
            raise Exception(" unprepared for this instance : "+str(vb) +" from "+ str(b.value))
        
class Functor_Compose_An_Absolute_Address:
    def __init__(self, take_first=0, is_retrieved_scan=0):
        self.take_first = take_first
        self.is_retrieved_scan = is_retrieved_scan
    def __call__(self,par):

        if self.is_retrieved_scan:
            par, N = par
            if hasattr(N,"value"):
                N=N.value
        
        if isinstance(par, hdf5_relative_path):
            b = par
            a = b.base_h5file
        else:
            a,b = par
            
        print( " ====================== " )
        print( a,b)

        va=a.getFullValue()
        vb=b.getValue()
        print( va, vb)
        res = va+"/"+vb
        if self.is_retrieved_scan:
            res = res+"/scans/Scan%03d"%int(N)
        return res
  


def dic2yaml_response_fit(dicodic):
    print( list(dicodic.keys()))
    s = "superR_fit_responses:\n"
    s = s+ "   foil_scan_address : %s \n" % dicodic["response_scan_address"].render() 
    s = s+ "   nref : %s \n" % dicodic["nref"].render ()
    s = s+ "   niter_optical : %s \n" % dicodic["niter_optical"].render ()
    s = s+ "   beta_optical : %s \n" % dicodic["beta"].render ()
    s = s+ "   niter_global : %s \n" % dicodic["niter_global"].render ()
    s = s+ "   pixel_dim  : 1 \n" 
    s = s+ "   niter_pixel : 0 \n" 
    s = s+ "   beta_pixel :  0 \n" 
    s = s+ "   do_refine_trajectory  : 1 \n"
    s = s+ "   simmetrizza : 1\n"
    s = s+ "   filter_rois : 0\n"
    s = s+ "   target_file : %s \n" % dicodic["result_file"].render()
    s = s+ "   MPI_N_PROCS : %s \n" % dicodic["MPI_N_PROCS"].render ()
    s = s+ "   fit_lines   : 1 \n"
    return s

class Functor_set_new_root:
    def __init__(self):
        pass
    def __call__(self,par):
        a,b = par
        va=a.getFullValue()
        print( va, b)
        return va+"/"+b

def dic2yaml_resynt_fit(dicodic):
    s = "superR_recreate_rois :\n"
    s = s+ "   responsefilename : %s \n" % dicodic["responsefilename"].render() 
    s = s+ "   nex : 0 \n"
    s = s+ "   old_scan_address : %s \n" % dicodic["response_scan_address"].render ()

    if "recenterings_refined" in dicodic:
        s = s+ "   recenterings_refined : %s \n" % dicodic["recenterings_refined"].render ()

    s = s+ "   filter_rois : 0\n"
    s = s+ "   target_filename : %s \n" % dicodic["result_file"].render() 
    return s



def getMethod(options):
    
    metodo = collections.OrderedDict()

    c_p = collections.OrderedDict()
    metodo["CREATION_PARS"] = c_p

    c_p["spec_file"]      = FilePath(doc="The path to the specfile\n file must exist", fileme=1)
    c_p["reference_scan"] = aNumber4SpecScan(c_p["spec_file"]  ,"This is the scan number containing references.\n Will also be used for ROI selection")


    c_p["filter_path"]    = Hdf5FilePath(doc="OPTIONAL(you can let it to None)\nThe default hdf5 file AND path  where the filter is stored\n must be a matrix  dataset with 1 for good pixels , 0 for bad ones.\n the format is filename:groupname. ", fileme=1, gnameme=1, canbeNone=True)

    

    c_p["result_file"]    = Hdf5FilePath(doc="The default hdf5 file where results will be stored\n must be a new file.\n the format is filename:groupname. The datagroup name must be a new one for the entry to be ready", fileme=0.5, gnameme=0)


    c_p["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                  " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                  " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                  initial_value = "kapraman/1000"
                                              )

    c_p["result_file"]    .isresult = 1

    c_p["define_new"] = Functor_adapt( metodo )
    c_r = collections.OrderedDict()
    ## ======================================================================
    metodo["CREATION_ROIS"] = c_r
    c_r["HELP"]  =  """ Create the ROIS for sample data"""

    c_r["getyaml"]  =  dic2yaml_create_rois  
    c_r["swissknife_runner"]  =  swissknife_runner  
    c_r["method_name"] = "create_rois"
    c_r["expdata"] = FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )
    c_r["scans"]   = aNumber4SpecScan( c_r["expdata"] ,  doc="This is the scan number containing references.\n"+
                                       " Will  be used here for ROI selection",
                                      defaults2 = (c_p["reference_scan"], simplecopy),
                                      rendering = render_numero_in_una_lista )


    c_r["filter_path"]    = Hdf5FilePath(doc="OPTIONAL(you can let it to None)\nThe default hdf5 file AND path  where the filter is stored\n must be a matrix  dataset with 1 for good pixels , 0 for bad ones.\n the format is filename:groupname. ", fileme=1, gnameme=1, canbeNone=True,
                                        defaults2 = (c_p["filter_path"], simplecopy))

    
    
    c_r["result_file"]      = Hdf5FilePath(doc="The  hdf5 file and datagroup where results will be stored. Format toto.hdf5:datagroupname",
                                       defaults2 = (c_p["result_file"], Functor_add_datagroup("ROI_AS_SELECTED") ), fileme=0.5, gnameme=0)

    c_r["result_file"].isresult = 1

    ## ==============================================================================

    ed_calib  = collections.OrderedDict()
    ed_calib["HELP"]  =  """ Collects the calibration scans """

    ed_calib["getyaml"]  =    dic2yaml_extraction_scan

    metodo["EXTRACTION_DATA_CALIBRATION"] = ed_calib
    ed_calib["method_name"] = "loadscan_2Dimages"
    ed_calib["expdata"] = FilePath(doc="The path to the specfile", defaults2 = (c_r["expdata"], simplecopy)    )
    ed_calib["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                            " where the ROI for calibration extraction can be found  \n"+
                                             "  Format toto.hdf5:datagroupname  ",
                                             defaults2 = (c_r["result_file"], simplecopy  ), 
                                             fileme=1, gnameme=1)
    ed_calib["scans"]   = aNumber4SpecScan( ed_calib["expdata"] ,  doc="This is the scan number containing references.\n",
                                            defaults2 = (c_r["scans"], simplecopy),
                                            rendering = render_numero_in_una_lista_piu_uno )


    ed_calib["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                       defaults2=(c_p["monitor_column"], simplecopy)
                                                   )


    ed_calib["isolateSpot"]   = aNumber( doc="If different from zero :\n"+
                                       " clean the area outside a radius=isolateSpot  from the maxium",
                                        defaults2 = ("7", simplecopy),nmin=0, nmax=40
                            )



    
    ed_calib["signaladdress"]   =  hdf5_relative_path( ed_calib["roiaddress"] ,
                                                       doc="The datagroup path, relative to roiaddress,\n"+
                                                       "where extracted data will be written. Must be a new name",
                                                       defaults2 = ( "calibration_scan" ,  simplecopy ) ,
                                                       gnameme = 0)

    ed_calib["swissknife_runner"]  =  swissknife_runner  

    ed_calib["signaladdress"]    .isresult = 1

    ## ==============================================================================

    return metodo







