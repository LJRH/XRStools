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



print( " Caller version", version)

XRStools = importlib.import_module( "XRStools"+version    )
exec("from XRStools%s.WIZARD.Wizard import *"%version)
###################################################################

import yaml 
import PyMca5.PyMcaIO.specfilewrapper as SpecIO


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
        p= subprocess.Popen(( "mpirun -n %d XRS_swissknife%s %s 1> %s  2>%s"%(version,MPI_N_PROCS  , inputname,stdname,  errname )) .split()   , shell=False,  stdout= stdout ,  stderr= stderr )
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
            for name,par in subdic.items():
                if isinstance(par,Parametro):
                    if par.defaults2 is not None:
                        source, method = par.defaults2
                        print( source)
                        par.value = method(source)
                   
        return self.dico["CREATION_PARS"]["result_directory"].value


def dic2yaml_create_rois(dicodic):
    s = "create_rois:\n"
    s = s+ "   expdata : %s \n" % dicodic["spec_file"].render()

    s = s+ "   scans : %s \n" % dicodic["elastic_scan"].render ()
    s = s+ "   roiaddress : %s \n" % dicodic["mask_file"].render()
    # tok = dicodic["filter_path"].render()
    # if tok !="None":
    #     s = s+ "   filter_path : %s \n" % tok
    
    # # s = s+ "   masktype   : filter \n"
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
    c_p["elastic_scan"] = aNumber4SpecScan(c_p["spec_file"]  ,"This is the scan number containing the elastic scan.\n Will also be for ROI selection")
    c_p["first_scan"]     = aNumber4SpecScan(c_p["spec_file"]  ,"This is the first scan of the to-be-summed interval.")
    c_p["last_scan"]     = aNumber4SpecScan(c_p["spec_file"]  ,"This is the first scan of the to-be-summed interval.")
    c_p["result_directory"]      = FilePath(doc="The path to the directory\n where results will be written", fileme=1, isadir=1)


    c_p["define_new"] = Functor_adapt( metodo )
    
    c_r = collections.OrderedDict()
    ## ======================================================================
    metodo["CREATION_ROIS"] = c_r
    c_r["HELP"]  =  """ Create the ROIS for sample data"""

    c_r["getyaml"]  =  dic2yaml_create_rois  
    c_r["swissknife_runner"]  =  swissknife_runner  
    c_r["method_name"] = "create_rois"
    c_r["spec_file"] = FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )

    
    c_r["elastic_scan"]   = aNumber4SpecScan( c_r["spec_file"] ,  doc="This is the scan number containing the elastic scan.\n"+
                                              " Will  be used here for ROI selection",
                                              defaults2 = (c_p["elastic_scan"], simplecopy),
                                              rendering = render_numero_in_una_lista )

    
    c_r["mask_file"] = Hdf5FilePath(doc="The  hdf5 file and datagroup where results the MASK will be WRITTEN. Format toto.hdf5:datagroupname",
                                    defaults2 = (c_p["result_directory"], Functor_add_datagroup("ROI", add_basename="roifile.h5") ), fileme=0.5, gnameme=0.5)


    c_r["mask_file"].master= True

    
    ########################## c_r["mask_file"].automatic_forward = 2
    c_r["mask_file"].isresult = 1

    # ## ==============================================================================


    s_s = collections.OrderedDict()
    ## ======================================================================
    metodo["SUM_SCANS"] = s_s
    s_s["HELP"]  =  """ Do the sum over scans providing two dimensional maps
    A label is selected for the Y axis.
    A motor is selected for the X axis.
    Scans having all equals motors and equal range for Y are summed together,
    otherwise they are summed into new maps
 """

    s_s["getyaml"]  =  dic2yaml_generic  
    s_s["swissknife_runner"]  =  swissknife_runner  
    s_s["method_name"] = "sum_scans2maps"


    s_s_spec_file = FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )
    

    s_s["first_scan"]    = aNumber4SpecScan(s_s_spec_file  ,"This is the first scan of the to-be-summed interval.", defaults2=(c_p["first_scan"], simplecopy))
    s_s["last_scan"]     = aNumber4SpecScan(s_s_spec_file  ,"This is the first scan of the to-be-summed interval.", defaults2=(c_p["last_scan"] , simplecopy))

    s_s["first_scan"].master= True
    s_s["last_scan"].master= True
    
    s_s["Scan_Variable"]           =   Parameter_Choices(  "Variable moving along the scan. To be chosed between the scans labels" ,
                                                           choices_functor = Functor_choicesFromSpecLabels( s_s_spec_file , s_s["first_scan"] ,  s_s["last_scan"]  ),
                                                           defaults2=( "Anal Energy" , simplecopy)
                                                       )

    s_s["Motor_Variable"]           =   Parameter_Choices(  "Motor changing from scan to scan. To be chosed between the specfile motors" ,
                                                           choices_functor = Functor_choicesFromSpecLabels( s_s_spec_file , s_s["first_scan"] ,  s_s["last_scan"] , labels=0 ),
                                                           defaults2=( "energy" , simplecopy)
                                                       )
    s_s["Scan_Variable"] .master= True
    s_s["Motor_Variable"].master= True

    s_s["spec_file"] = s_s_spec_file
    s_s["mask_file"] = Hdf5FilePath(doc="The  hdf5 file and datagroup where results the MASK will be WRITTEN. Format toto.hdf5:datagroupname",
                                    defaults2 = (c_r["mask_file"], simplecopy ), fileme=0.5, gnameme=0.5)

    s_s["scans_file"] = Hdf5FilePath(doc="The  hdf5 file and datagroup where the SUMs  will be WRITTEN. Format toto.hdf5:datagroupname",
                                     defaults2 = ((c_p["result_directory"],  c_r["mask_file"]  ), Functor_add_datagroup("SCAN",  add_basename="scansfile.h5") ), fileme=0.5, gnameme=0.5 )

    
    s_s["scans_file"] .always_visible = True
    
    s_s["scans_file"].isresult = 1    

    

    

    # ## ==============================================================================

    return metodo







