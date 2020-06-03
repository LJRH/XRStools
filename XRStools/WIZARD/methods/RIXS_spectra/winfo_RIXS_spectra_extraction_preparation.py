from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import collections
import os
import inspect


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


version = caller[len("XRS_wizard"):]
print( " Caller version", version)

XRStools = importlib.import_module( "XRStools"+version    )
exec("from XRStools%s.WIZARD.Wizard import *"%version)


###################################################################


import yaml 


" modification to test branching "

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
    if not sys.platform=="win32":
        if MPI_N_PROCS==1:
            p= subprocess.Popen(( "XRS_swissknife%s %s 1> %s  2>%s"%(version, inputname,stdname,  errname )).split()    , shell=False,  stdout= stdout ,  stderr= stderr )
        else:
            p= subprocess.Popen(( "mpirun -n %d XRS_swissknife%s %s 1> %s  2>%s"%(version, MPI_N_PROCS  , inputname,stdname,  errname )) .split()   , shell=False,  stdout= stdout ,  stderr= stderr )
    else:
        packages = (os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(XRStools.__file__ ) ))))
        script = os.path.join(packages,"scripts","XRS_swissknife%s.bat"%version )
        if MPI_N_PROCS==1:
            p= subprocess.Popen(( "%s %s"%(script, inputname)) .split()   , shell=False,  stdout= stdout ,  stderr= stderr )
        else:
            p= subprocess.Popen(( "mpiexec -n %d   %s %s"%(MPI_N_PROCS  ,script, inputname) ) .split()  , shell=False,  stdout= stdout ,  stderr= stderr )


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
                   
        return self.dico["CREATION_PARS"]["result_file"].value


def dic2yaml_create_rois(dicodic):
    s = "create_rois:\n"
    s = s+ "   expdata : %s \n" % dicodic["expdata"].render()
    
    tok = dicodic["filter_path"].render()
    if tok not in ["None",None]:
        s = s+ "   filter_path : %s \n" % tok

    s = s+ "   scans : %s \n" % dicodic["scans"].render ()
    s = s+ "   roiaddress : %s \n" % dicodic["result_file"].render() 
    return s


def dic2yaml_extraction_scan(dicodic):
    s = "loadscan_2Dimages :\n"
    s = s+ "  expdata : %s \n" % dicodic["expdata"].render() 
    s = s+ "  roiaddress : %s \n" % dicodic["roiaddress"].render() 
    s = s+ "  monitor_column : %s \n" % dicodic["monitor_column"].render() 
    s = s+ "  scan_interval : %s \n" % dicodic["scans"].render() 
    s = s+ "  signaladdress : %s \n" % dicodic["signaladdress"].render()
    s = s+"""
  sumto1D  : 0
  energycolumn : 'Anal Energy'
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
  

class Functor_compose_recentering_name:
    def __init__(self,expd, scans):
        self.expd = expd
        self.scans = scans
    def __call__(self, fn):
        if isinstance(fn, Parametro):
            fn = fn.render()
        # a,b = string.split(fn,"/")
        expd = self.expd.render()
        expd=expd.replace("/","_")
        expd=expd.replace("\\","_")

        scans = self.scans.getValue()
        i,f = scans[:2]
        # name = a+expd+"/"+b+"_"+str(i)+"_"+str(f)
        name =  fn +"/recentering"+"_"+str(i)+"_"+str(f)
        return name
        

    
def dic2yaml_calculate_recenterings(dicodic):
    s = "calculate_recenterings :\n"
    s = s+ "   bariA : %s \n" % dicodic["bariA"].render() 
    s = s+ "   bariB : %s \n" % dicodic["bariB"].render() 
    s = s+ "   target : %s \n" % dicodic["recenterings"].render() 
    return s


def dic2yaml_extraction_scan_with_recentering(dicodic):
    s = "loadscan_2Dimages :\n"
    s = s+ "  expdata : %s \n" % dicodic["expdata"].render() 
    s = s+ "  roiaddress : %s \n" % dicodic["roiaddress"].render() 
    s = s+ "  scan_interval : %s \n" % dicodic["scans"].render() 
    s = s+ "  monitor_column : %s \n" % dicodic["monitor_column"].render() 
    s = s+ "  recenterings  : %s \n" % dicodic["recenterings"].render()
    s = s+ "  recenterings_confirmed : %s \n" % dicodic["recenterings_refined"].render()
    s = s+ "  signaladdress : %s \n" % dicodic["signaladdress"].render()
    s = s+"""
  sumto1D  : 0
  energycolumn : 'Anal Energy'
  monitorcolumn : %s
""" %  dicodic["monitor_column"].render()
    return s


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

    if "recenterings_refined" in dicodic and dicodic["recenterings_refined"].render () not in ["", None]:
        s = s+ "   recenterings_refined : %s \n" % dicodic["recenterings_refined"].render ()

    s = s+ "   filter_rois : 0\n"
    s = s+ "   target_filename : %s \n" % dicodic["result_file"].render() 
    return s



def getMethod(options):

    
    shift_the_reference = True
    if "shift_the_reference" in options:
        shift_the_reference = options["shift_the_reference"]


    
    metodo = collections.OrderedDict()


    c_p = collections.OrderedDict()
    metodo["CREATION_PARS"] = c_p

    c_p["spec_file"]      = FilePath(doc="The path to the specfile\n file must exist", fileme=1)
    
    
    c_p["reference_scan"] = aNumber4SpecScan(c_p["spec_file"]  ,"This is the scan number containing references.\n Will also be used for ROI selection")

    c_p["sample_scans"]   = many_Number4SpecScan(c_p["spec_file"]  ,
                                                 "This defines the list(**Use Python syntax [1,2,3]**)\n"+
                                                 "of scan numbers containing data.\n"+
                                                 # " Must containing two numbers : [ start,end+1] \n"+
                                                 "Will also be used for ROI selection",
                                                 nmin=1,nmax=100)


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

    c_r["getyaml"]  =  dic2yaml_create_rois  , {"Reference_scan":"scans","SpecExpFile":"expdata"}
    c_r["swissknife_runner"]  =  swissknife_runner  
    c_r["method_name"] = "create_rois"
    c_r["SpecExpFile"] = FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )

    c_r["SpecExpFile"].master = 1
    
    c_r["Reference_scan"]   = aNumber4SpecScan( c_r["SpecExpFile"] ,  doc="This is the scan number containing references.\n"+
                                       " Will  be used here for ROI selection",
                                      defaults2 = (c_p["reference_scan"], simplecopy),
                                      rendering = render_numero_in_una_lista )


    c_r["filter_path"]    = Hdf5FilePath(doc="OPTIONAL(you can let it to None)\nThe default hdf5 file AND path  where the filter is stored\n must be a matrix  dataset with 1 for good pixels , 0 for bad ones.\n the format is filename:groupname. ", fileme=1, gnameme=1, canbeNone=True,
                                        defaults2 = (c_p["filter_path"], simplecopy))
   
    
    c_r["Reference_scan"].master = 1
    
    c_r["result_file"]      = Hdf5FilePath(doc="The  hdf5 file and datagroup where results will be stored. Format toto.hdf5:datagroupname",
                                       defaults2 = (c_p["result_file"], Functor_add_datagroup("ROI_AS_SELECTED") ), fileme=0.5, gnameme=0)

    c_r["result_file"].isresult = 1
    c_r["result_file"].always_visible = True

    ## ==============================================================================
    c_r_sample = collections.OrderedDict()
    c_r_sample["HELP"]  =  """ Create the ROIS for sample data"""

    c_r_sample["getyaml"]  =  dic2yaml_create_rois  

    metodo["CREATION_ROIS_SAMPLE"] = c_r_sample
    c_r_sample["method_name"] = "create_rois"
    c_r_sample["expdata"] = FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )
    c_r_sample["scans"]   = many_Number4SpecScan(  c_r_sample["expdata"]  ,
                                                   "This is the list(**Use Python syntax [1,2,3]**)\n"+
                                                   "of scan numbers containing data.\n"+
                                                   "Will  be used here for ROI selection\n",
                                                   # " Must containing two numbers : [ start,end+1] \n",
                                                   nmin=1,nmax=100,
                                                   defaults2 = (c_p["sample_scans"], simplecopy))
    c_r_sample["scans"]  .master = True


    c_r_sample["filter_path"]    = Hdf5FilePath(doc="OPTIONAL(you can let it to None)\nThe default hdf5 file AND path  where the filter is stored\n must be a matrix  dataset with 1 for good pixels , 0 for bad ones.\n the format is filename:groupname. ", fileme=1, gnameme=1, canbeNone=True,
                                                defaults2 = (c_p["filter_path"], simplecopy))
    
    
    
    
    c_r_sample["result_file"]      = Hdf5FilePath(doc="The  hdf5 file and datagroup where results will be stored. Format toto.hdf5:datagroupname",
                                                  defaults2 = (c_p["result_file"], Functor_add_datagroup("ROI_FROMSAMPLE") ), fileme=0.5, gnameme=0)

    c_r_sample["result_file"] .always_visible = True
    
    c_r_sample["result_file"].isresult = 1
    c_r_sample["swissknife_runner"]  =  swissknife_runner  



    ## ==============================================================================
    ed_calib  = collections.OrderedDict()
    ed_calib["HELP"]  =  """ Collects the calibration scans """

    ed_calib["getyaml"]  =    dic2yaml_extraction_scan

    metodo["EXTRACTION_DATA_CALIBRATION"] = ed_calib
    ed_calib["method_name"] = "loadscan_2Dimages"
    ed_calib["expdata"] = FilePath(doc="The path to the specfile", defaults2 = (c_r["SpecExpFile"], simplecopy)    )
    ed_calib["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                            " where the ROI for calibration extraction can be found  \n"+
                                             "  Format toto.hdf5:datagroupname  ",
                                             defaults2 = (c_r["result_file"], simplecopy  ), 
                                             fileme=1, gnameme=1)
    ed_calib["scans"]   = aNumber4SpecScan( ed_calib["expdata"] ,  doc="This is the scan number containing references.\n",
                                            defaults2 = (c_r["Reference_scan"], simplecopy),
                                            rendering = render_numero_in_una_lista_piu_uno )
    
    ed_calib["scans"] .master = True

    ed_calib["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                       defaults2=(c_p["monitor_column"], simplecopy)
                                                   )

    ed_calib["monitor_column"]  .master = True


    ed_calib["signaladdress"]   =  hdf5_relative_path( ed_calib["roiaddress"] ,
                                                       doc="The datagroup path, relative to roiaddress,\n"+
                                                       "where extracted data will be written. Must be a new name",
                                                       defaults2 = ( "calibration_scan" ,  simplecopy ) ,
                                                       gnameme = 0)

    ed_calib["signaladdress"] .always_visible = True

    ed_calib["swissknife_runner"]  =  swissknife_runner  

    ed_calib["signaladdress"]    .isresult = 1


    ## ==============================================================================
    es_nonmoved  = collections.OrderedDict()
    es_nonmoved["HELP"]  =  """ Collects the scans data for the sample, based on the sample ROIS"""

    es_nonmoved["getyaml"] = dic2yaml_extraction_scan

    metodo["EXTRACTION_DATA_SAMPLE_UNCENTERED"] = es_nonmoved
    es_nonmoved["method_name"] = "loadscan_2Dimages"
    es_nonmoved["expdata"] = FilePath(doc="The path to the specfile", defaults2 = (ed_calib["expdata"], simplecopy)  )
    es_nonmoved["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                                " where the ROI for calibration extraction can be found  \n"+
                                                "  Format toto.hdf5:datagroupname  ",
                                                defaults2 = (c_r_sample["result_file"] , simplecopy  ), 
                                                fileme=1, gnameme=1)
    
    es_nonmoved["roiaddress"] .always_visible = True

    es_nonmoved["scans"]   = many_Number4SpecScan(  es_nonmoved["expdata"]  ,isinterval=1, doc=
                                                    "This is the list of two number (**Use Python syntax [start,end+1]**)\n"+
                                                   " defining the intervals scan numbers containing data.\n"+
                                                    "Will also be used for ROI selection\n"+
                                                    "You can also give more intervals :  [start,end+1,start2,end2+1] in the same list   \n"
                                                    "EXAMPLE  [263,264] will read only scan 263\n"
                                                    "         [263,264, 270, 276]  will read 263, 270,271,272,273,274,275\n"
                                                    # " Must containing two numbers : [ start,end+1] \n"+
                                                    "",   # "Will also be used for ROI selection",
                                                    nmin=2,nmax=200,
                                                    defaults2 = (c_r_sample["scans"], makeintervalfromfirst))

    # es_nonmoved["scans"].master = True
    

    es_nonmoved["signaladdress"]   =  hdf5_relative_path( c_r_sample["result_file"] ,
                                                          doc="The datagroup path, relative to roiaddress,\n"+
                                                          "where extracted data will be written. Must be a new name",
                                                          defaults2 = ( "uncentered_sample_scans" ,  simplecopy ) ,
                                                          gnameme = 0)

    es_nonmoved["signaladdress"].always_visible = True

    
    es_nonmoved["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                       defaults2=(c_p["monitor_column"], simplecopy)
                                                   )


    es_nonmoved["signaladdress"]    .isresult = 1
 
    
    es_nonmoved["swissknife_runner"]  =  swissknife_runner  


    ## ==============================================================================
    calc_rec  = collections.OrderedDict()
    metodo["CALCULATE_RECENTERINGS"] = calc_rec
    calc_rec["HELP"] = "Calculates the horizontal shifts to go from baricenters A to baricenters B"

    calc_rec["getyaml"]  =    dic2yaml_calculate_recenterings
    calc_rec["bariA"]   =  Hdf5FilePath(doc=" The scan on which baricenters are calculated for sample.  (hdf5 full path. must exist)  ",
                                        defaults2 = (   (es_nonmoved["signaladdress"],  es_nonmoved["scans"]  ) , Functor_Compose_A_Scan_Address(take_first=1) ),
                                        fileme=1, gnameme=1)

    calc_rec["bariB"]   =  Hdf5FilePath(doc=" The scan on which baricenters are calculated for reference. (hdf5 full path. must exist) ",
                                        defaults2 = (   (ed_calib["signaladdress"],  ed_calib["scans"]  ) , Functor_Compose_A_Scan_Address() ),
                                        fileme=1, gnameme=1)

    calc_rec["recenterings"]   =  Hdf5FilePath(doc=("The hdf5 full path on which shifts will be written. Must NOT exist\n"
                                                    "The shifts informations will consist in the baricenter goal B,\n"
                                                    " and the starting baricenter B"),
                                               defaults2 = ( c_p["result_file"]   , Functor_compose_recentering_name( c_r["SpecExpFile"] ,es_nonmoved["scans"]  ) ),
                                               fileme=0.5, gnameme=0)


    calc_rec["recenterings"]  .always_visible = True

    calc_rec["recenterings"].isresult = 1
    calc_rec["swissknife_runner"]  =  swissknife_runner  

    ## ==============================================================================
    es_moved  = collections.OrderedDict()
    metodo["EXTRACTION_DATA_SAMPLE_CENTERED"] = es_moved

    es_moved["getyaml"]  =   dic2yaml_extraction_scan_with_recentering

    es_moved["HELP"] = " Recollect the sample data, with\n the rois centered on calibration scan\n, and with the recentering  displacement"

    es_moved["method_name"]  =  "loadscan_2Dimages"
    es_moved["expdata"]      =  FilePath(doc="The path to the specfile", defaults2 = (c_r["SpecExpFile"], simplecopy)  )
    es_moved["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                                " where the ROI for calibration extraction can be found  \n"+
                                                "  Format toto.hdf5:datagroupname  ",
                                                defaults2 = (c_r["result_file"], simplecopy  ), 
                                                fileme=1, gnameme=1)

   
    es_moved["scans"]   = many_Number4SpecScan(  es_moved["expdata"]  ,isinterval=1, doc=
                                                 "This is the list of two number (**Use Python syntax [start,end+1]**)\n"+
                                                 " defining the interval scan numbers containing data.\n"+
                                                 "You can also give more intervals in the same list.\n"+
                                                 # " Must containing two numbers : [ start,end+1] \n"+
                                                 "",   # "Will also be used for ROI selection",
                                                 nmin=2,nmax=200,
                                                 defaults2 = ( es_nonmoved["scans"]   , simplecopy))
    # es_moved["scans"].master = True

    es_moved["recenterings"]      = Hdf5FilePath(doc="The  hdf5 file where recentering informations can be found",
                                                  defaults2 = ( calc_rec["recenterings"], simplecopy), fileme=1, gnameme=1)

    es_moved["recenterings_refined"]   =  Hdf5FilePath(doc=("The hdf5 full path on which shifts will be written. Must NOT exist\n"
                                                            "The shift informations will consist of the refined shift"),
                                                       defaults2 = (  es_moved["recenterings"]   ,
                                                                      Functor_add_datagroup("confirmed_shift", remplace_last= True, keep_last_numbers=True)),
                                                       fileme=0.5, gnameme=0)
    
    es_moved["recenterings_refined"] .always_visible = True

    es_moved["signaladdress"]   =  hdf5_relative_path( es_moved["roiaddress"] ,
                                                       doc="The datagroup path, relative to roiaddress,\n"+
                                                       "where extracted data will be written. Must be a new name",
                                                       defaults2 = ( "Centered_sample_scans" ,  simplecopy ) ,
                                                       gnameme = 0)

    es_moved["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                       defaults2=(es_nonmoved["monitor_column"], simplecopy)
                                                   )


    es_moved["signaladdress"]  .isresult = 1
    es_moved["recenterings_refined"] .isresult = 1
    es_moved["swissknife_runner"]  =  swissknife_runner  

    # =============================================================================================================
    f_r = collections.OrderedDict()
    metodo["RESPONSE_FIT"] = f_r

    f_r["HELP"]  =  """    The method first find an estimation of the response scan trajectory on each roi
    then, based on this, obtain a fit of the  response function
    assuming a flat pixel response. 
    There is a final phase where a global optimisation
    is done in niter_global steps.
    Each step is composed of  response functionfit, followed by a trajectory fit.
    """
    f_r["getyaml"]  =  dic2yaml_response_fit
    f_r["swissknife_runner"]  =  swissknife_runner  
    f_r["method_name"] = "response_fit"

    f_r["response_scan_address"]      = Hdf5FilePath(doc="The full address of the scan which measures the response",
                                                     defaults2 = ((ed_calib["signaladdress"], ed_calib["scans"] ) ,
                                                                  Functor_Compose_An_Absolute_Address( is_retrieved_scan=1)),
                                                     fileme=1, gnameme=1)
    f_r["nref"]      =  aNumber(  doc="The number of subdivision per pixel dimension used to represent\n"+
                                  "the optical response function at higher resolution",
                                  defaults2 = ("2", simplecopy) )

    f_r["nref"].master = True
    
    f_r["niter_optical"]      =  aNumber(  doc="the number of iterations used in the \n"+
                                           " optimisation of the optical response",
                                           defaults2 = ("50", simplecopy) )

    f_r["niter_optical"]   .master = True
    
    f_r["beta"] =  aNumberFloat(  doc="The L1 norm factor in the regularisation\n"+
                                          "term for the response functions",
                                          defaults2 = ("0.1", simplecopy) )

    f_r["beta"].master = True
    
    f_r["niter_global"] = aNumber(  doc="Each step is composed of optical response fit, followed \n"+
                                    " trajectory reoptimisation. All this niter_global times ",
                                    defaults2 = ("2", simplecopy) )

    f_r["niter_global"].master = True
    
    f_r["MPI_N_PROCS"] = aNumber(  doc="If bigger than 1, run it with mpi.\n"
                                   " useful if you have several ROIs \n"
                                   " It is not useful to haveMPI_N_PROCS > n rois  \n",
                                    defaults2 = ("1", simplecopy) )

    
    f_r["MPI_N_PROCS"].master = True

    f_r["result_file"]    = Hdf5FilePath(doc=("The default hdf5 file where results will be stored\n"
                                              " must be a new address.\n"
                                              " the format is filename:groupname.\n"
                                              " The datagroup name must be a new one for the entry to be ready")
                                         , defaults2 = ( f_r["response_scan_address"]  ,
                                                         Functor_add_datagroup("ScanFittedResponse", transform_last= 1)),
                                         fileme=0.5, gnameme=0)
    f_r["result_file"].isresult = 1
    f_r["result_file"].always_visible = True

    ## ======================================================================

    r_r = collections.OrderedDict()
    metodo["RESYNTHETISE_FIT"] = r_r
    r_r["HELP"]  =  """   Recreate the reference scan
    """

    r_r["getyaml"]  =  dic2yaml_resynt_fit, {"ResynthScan":"result_file"}
    r_r["swissknife_runner"]  =  swissknife_runner  
    r_r["method_name"] = "Scan Resynth."

    r_r["responsefilename"]      = Hdf5FilePath(doc="The  address on which the fitted responsed was written",
                                                     defaults2 = (f_r["result_file"]  , simplecopy),
                                                     fileme=1, gnameme=1)

    r_r["response_scan_address"]      = Hdf5FilePath(doc="The original address of the scan which measures the response",
                                                     defaults2 = (f_r["response_scan_address"] , simplecopy),
                                                     fileme=1, gnameme=1)


    if shift_the_reference:
        r_r["recenterings_refined"]    =  Hdf5FilePath(doc=("The hdf5 full path containing the refined shift"),
                                                       defaults2 = ( es_moved["recenterings_refined"] , simplecopy ),
                                                       fileme=1, gnameme=1)



    r_r["ResynthScan"]    = Hdf5FilePath(doc="The default hdf5 file where results will be stored\n"+
                                         "must be a new file.\n"+
                                         " The format is filename:groupname. \n"+
                                         "The datagroup name must be a new one\n"+
                                         "for the entry to be ready",
                                         defaults2 = ( f_r["result_file"], Functor_add_datagroup("ScanReSynth", transform_last= 1)      ),
                                         fileme=0.5, gnameme=0)

    r_r["ResynthScan"].isresult = 1
    r_r["ResynthScan"].always_visible = True


    return metodo







