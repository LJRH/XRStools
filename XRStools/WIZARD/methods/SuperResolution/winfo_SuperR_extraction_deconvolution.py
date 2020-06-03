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

    if not sys.platform=="win32":
        if MPI_N_PROCS==1:
            p= subprocess.Popen(( "XRS_swissknife%s %s 1> %s  2>%s"%(version, inputname,stdname,  errname )) .split()   , shell=False,  stdout= stdout ,  stderr= stderr )
        else:
            p= subprocess.Popen(( "mpirun -n %d XRS_swissknife%s %s 1> %s  2>%s"%(version, MPI_N_PROCS  , inputname,stdname,  errname )) .split()   , shell=False,  stdout= stdout ,  stderr= stderr )
    else:
        packages = (os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(XRStools.__file__ ) ))))
        script = os.path.join(packages,"scripts","XRS_swissknife%s.bat"%version )
        if MPI_N_PROCS==1:
            p= subprocess.Popen(( "%s %s"%(script, inputname)) .split()    , shell=False,  stdout= stdout ,  stderr= stderr )
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

                        if hasattr(source,"value" ):
                            par.value = method(source.value)
                        else:

                            par.value = method(source)

                   
        return self.dico["CREATION_PARS"]["target_address"].value


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

class functor_finalise_last_part_of_the_path:
    def __init__(self, cp_num):
        self.cp_num=cp_num


    def __call__(self, cp_add ) :
        cp_num  = self.cp_num
        path = cp_add
        num = cp_num.getValue()
        # path_items = path.split("/")[-3:]
        
        res = path + "/scans/Scan%03d" % num
        
        # res = path + "/" + path_items[0] + "/" + path_items[1] +"/" + "Scan%03d" % num
        # res = path + "/" + path_items[0] + "/" 
        print( " RITORNO " , res)
        return res




def dic2yaml_extraction_scan(dicodic):
    s = "loadscan_2Dimages :\n"
    s = s+ "  expdata : %s \n" % dicodic["expdata"].render() 
    s = s+ "  roiaddress : %s \n" % dicodic["roiaddress"].render() 
    s = s+ "  monitor_column : %s \n" % dicodic["monitor_column"].render() 
    s = s+ "  scan_interval : %s \n" % dicodic["scans"].render()
    s = s+ "  energy_column : %s \n" % 'sty'
    s = s+ "  signaladdress : %s \n" % dicodic["signaladdress"].render()
    s = s+"""
  sumto1D  : 0
  monitorcolumn : %s

""" % dicodic["monitor_column"].render()
    
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

    c_p["scans"]   = many_Number4SpecScan(  c_p["spec_file"] ,isinterval=1, doc=
                                            "This is the list of two number (**Use Python syntax [start,end+1]**)\n"+
                                            " defining the interval scan numbers containing data.\n"+
                                            "You can also give more intervals in the same list.\n"+
                                            # " Must containing two numbers : [ start,end+1] \n"+
                                            "",   # "Will also be used for ROI selection",
                                            nmin=2,nmax=200)

    
    c_p["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                  " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                  " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                  initial_value = "kapraman/1000"
                                              )
    
    c_p["reference_address"] =Hdf5FilePath(doc=
                                           "The full path where the references can be found\n"
                                           "It is  the hdf5 file plus the path pointing to a datagroup\n"
                                           " this datagroup must contain a datagroup named 'scans' \n"
                                           " which at its turn contains ScanXXX for each reference scan " ,
                                           fileme=1, gnameme=1)

    
    c_p["reference_scan"]   = aNumber(   doc=
                                         "This is the  Reference scan Number\n"
                                         " It is used to point to ScanXXX and select the scan number "
                                    ) 

    c_p["signaladdress"]   =  Hdf5FilePath(doc=("The  hdf5 file where extracted data will be written\n"
                                                     "Must be a free name\n"), fileme=0.5, gnameme=0)

    
 

    c_p["target_address"]            =  Hdf5FilePath(doc="The hdf5 full path on which Volumes will be written. Must NOT exist   ",
                                             fileme=0.5, gnameme=0)

    c_p["signaladdress"]    .isresult = 1
    c_p["target_address"]    .isresult = 1


    c_p["define_new"] = Functor_adapt( metodo )

    ## ==============================================================================
    es  = collections.OrderedDict()
    metodo["DATA_EXTRACTION"] = es

    es["getyaml"]  =   dic2yaml_extraction_scan 

    es["HELP"] = " Collect the sample data\n"

    es["method_name"]  =  "loadscan_2Dimages"
    es["expdata"]      =  FilePath(doc="The path to the specfile", defaults2 = (c_p["spec_file"], simplecopy)  )


    es["roiaddress"]   =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                       " where the ROI for calibration extraction can be found  \n"+
                                       "  Format toto.hdf5:datagroupname  ",
                                       defaults2 = (  c_p["reference_address"],
                                                      Functor_add_datagroup(  "", remplace_last=1) ), 
                                             fileme=1, gnameme=1)



    
    es["scans"]   = many_Number4SpecScan(  es["expdata"]  ,isinterval=1, doc=
                                                 "This is the list of two number (**Use Python syntax [start,end+1]**)\n"+
                                                 " defining the interval scan numbers containing data.\n"+
                                                 "You can also give more intervals in the same list.\n"+
                                                 # " Must containing two numbers : [ start,end+1] \n"+
                                                 "",   # "Will also be used for ROI selection",
                                                 nmin=2,nmax=200,
                                                 defaults2 = ( c_p["scans"]   , simplecopy))

    
    
    es["signaladdress"]   =  Hdf5FilePath(doc=("The  hdf5 file where extracted data will be written\n"
                                                     "Must be a free name\n"),
                                                defaults2 = (  ( c_p["signaladdress"],
                                                                 es["scans"]   ) ,
                                                            Functor_append_numinterval_to_name() )
                                        , fileme=0.5, gnameme=0)

    es["monitor_column"]    = NameWithNormaliser(doc="The name of the spec monitor, used to renormalise the intensities\n"+
                                                       " this name maybe absent from the spec file, in this case no normalisation will be done\n"+
                                                       " You can also write something like kapraman/100000 to renormalise too big monitor intensities ",
                                                       defaults2=(c_p["monitor_column"], simplecopy)
                                                   )

    es["signaladdress"]  .isresult = 1
    es["swissknife_runner"]  =  swissknife_runner  

    # ==============================================================================

    
    scal  = collections.OrderedDict()
    metodo["SCALAR_PRODUCTS"] = scal
    scal["getyaml"]  =   dic2yaml_generic 
    scal["method_name"]  =  "superR_scal_deltaXimages"

    scal["sample_address"] =  Hdf5FilePath(doc=" The  hdf5 file and datagroup \n"+
                                           " where the extracted scans   can be found  \n"+
                                           " It contains ScanXXX for every extracted scan ",
                                           defaults2 = ( es["signaladdress"] ,
                                                         Functor_add_datagroup(  "scans") ), 
                                           fileme=1, gnameme=1)

 
    scal["delta_address"] =  Hdf5FilePath(doc="Where the references can be found\nIn general a subgroup of ROI_something\n"+
                                                       "containing a group called scans.\n The precise scan will be indicated by reference_scan below",
                                                       defaults2 = (   c_p["reference_address"]      , functor_finalise_last_part_of_the_path(c_p["reference_scan"] ) ),
                                                       fileme=1, gnameme=1)



    scal["nbin"]   = aNumber( doc="If different from 1 :\n"+
                              " the reference scan positions will be binned",
                              defaults2 = ("1", simplecopy),nmin=1, nmax=40
                          )

    scal["optional_solution"]   =  Hdf5FilePath(doc=("OPTIONAL : The  hdf5 file where scalar volumes has been written\n"
                                                   "If given, used for balancing analyzer factors\n"),
                                              defaults2 = ( "", simplecopy  )
                                              , fileme=0.5, gnameme=0.5)


    



    scal["target_address"]   =  Hdf5FilePath(doc=("The  hdf5 file where scalar products  will be written\n"
                                          "Must be a free name\n"),
                                             defaults2 = (  (  c_p["target_address"],
                                                               es["scans"]   ) ,
                                                            Functor_append_numinterval_to_name( plus="scal_prods") )
                                             , fileme=0.5, gnameme=0)
    


    
    scal["HELP"] = " Scalar products between sample data ROIs and reference scan\n"
    
    scal["target_address"]  .isresult = 1
    scal["swissknife_runner"]  =  swissknife_runner  

    # ==============================================================================

    
    gv  = collections.OrderedDict()
    metodo["GET_VOLUME"] = gv
    gv["getyaml"]  =   dic2yaml_generic 
    gv["method_name"]  =  "superR_getVolume"

    
    
    gv["scalprods_address"] =  Hdf5FilePath(doc=" The  hdf5 file and datagroup where\n"+
                                            " scalar products can be found",
                                            defaults2 = (scal["target_address"]  ,simplecopy), 
                                            fileme=1, gnameme=1)
    

    gv["target_address"]   =  Hdf5FilePath(doc=("The  hdf5 file where volumes will be written\n"
                                                "Must be a free name\n"),
                                           defaults2 = (  (  c_p["target_address"],
                                                             es["scans"]   ) ,
                                                          Functor_append_numinterval_to_name( plus="Volume") )
                                           , fileme=0.5, gnameme=0)

    
    gv["niter"]  =    aNumber( doc="The number of iterations",
                               defaults2 = ("100", simplecopy),nmin=1, nmax=40000000
                           )
    
    gv["beta"]  =    aNumber( doc="the TV preconditioning factor", float_also=True,
                              defaults2 = ("0.06", simplecopy),nmin=0, nmax=1.0e38
                          )
    
    gv["eps"]  =    aNumber( doc="a parameter for the convergence of the Chambolle-Pock TV denoising phase", float_also=True,
                              defaults2 = (" 0.000002", simplecopy),nmin=0, nmax=1.0e38
                          )
    
    gv["debin"]  =  many_Number(   doc="The parameter debin dafaults to [1,1]\n"
                                   "It is used to increase a dimension Z,Y or both , to make it match with X",
                                    defaults2 = ("[1,1]", simplecopy),
                                   nmin=2, nmax=2)

    gv["HELP"] = "Obtaining the volume \n"
    
    gv["target_address"]  .isresult = 1
    gv["swissknife_runner"]  =  swissknife_runner  
    # ==============================================================================

    
    vis  = collections.OrderedDict()
    metodo["VIS_VOLUME"] = vis
    vis["getyaml"]  =   dic2yaml_generic 
    vis["method_name"]  =  "view_Volume_myavi"

    vis["volume_address"] =  Hdf5FilePath(doc=" The  hdf5 file and datagroup where volume is \n",
                                          defaults2 = ( gv["target_address"]  ,simplecopy), 
                                          fileme=1, gnameme=1)
    
    
    vis["isolevel"]  =    aNumber( doc="The isolevel : a number  between 0 and 1 ( 0 represents minimum, 1, maximum)", float_also=True,
                              defaults2 = ("0.1", simplecopy),nmin=0, nmax=1.0
                          )

    vis["opacity"]  =    aNumber( doc="The opacity : a number  between 0 and 1 ", float_also=True,
                                  defaults2 = ("0.6", simplecopy),nmin=0, nmax=1.0)
           

    vis["HELP"] = "Visualising teh volume \n"

    vis["swissknife_runner"]  =  swissknife_runner  

    vis["opacity"]  .master = True
    vis["isolevel"] .master = True
    gv["beta"]      .master = True
    gv["niter"]     .master = True
    gv["target_address"] . master = True
    return metodo


    







