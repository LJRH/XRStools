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

class functor_sigterm_handler:
    def __init__(self, p):
        self.p = p
    def __call__(self, *x):
        print( x)
        print( " KILLO " , self.p.pid)
        os.kill(self.p.pid, signal.SIGKILL)

def swissknife_runner( yamltext, where ):

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
    open(inputname,"w").write(yamltext)
    stdout = open(stdname,"w")
    stderr = open(errname,"w")
    if MPI_N_PROCS==1:
        p= subprocess.Popen(( "XRS_swissknife%s %s 1> %s  2>%s"%(version, inputname,stdname,  errname )) .split()   , shell=False,  stdout= stdout ,  stderr= stderr )
    else:
        p= subprocess.Popen(( "mpirun -n %d XRS_swissknife%s %s 1> %s  2>%s"%(version, MPI_N_PROCS  , inputname,stdname,  errname )).split()    , shell=False,  stdout= stdout ,  stderr= stderr )
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
                                   
        return self.dico["CREATION_PARS"]["result_file"].value

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
    s = s+ "   filter_rois : 0\n"
    s = s+ "   target_filename : %s \n" % dicodic["result_file"].render() 
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


def getMethod(options):


    metodo = collections.OrderedDict()

    c_p = collections.OrderedDict()
    metodo["CREATION_PARS"] = c_p

    c_p["response_scan_address"]      = Hdf5FilePath(doc="The full address of the scan which measures the response", fileme=1, gnameme=1)

    c_p["result_file"]    = Hdf5FilePath(doc="The default hdf5 file where results will be stored\n must be a new file.\n the format is filename:groupname. The datagroup name must be a new one for the entry to be ready", fileme=0.5, gnameme=0)


    c_p["result_file"]    .isresult = 1






    c_p["define_new"] = Functor_adapt( metodo)


    f_r = collections.OrderedDict()
    ## ======================================================================
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
                                                     defaults2 = (c_p["response_scan_address"], copy.deepcopy),
                                                     fileme=1, gnameme=1)
    f_r["nref"]      =  aNumber(  doc="The number of subdivision per pixel dimension used to represent\n"+
                                  "the optical response function at higher resolution",
                                  defaults2 = ("2", copy.deepcopy) )
    f_r["niter_optical"]      =  aNumber(  doc="the number of iterations used in the \n"+
                                           " optimisation of the optical response",
                                           defaults2 = ("50", copy.deepcopy) )
    f_r["beta"] =  aNumberFloat(  doc="The L1 norm factor in the regularisation\n"+
                                          "term for the response functions",
                                          defaults2 = ("0.1", copy.deepcopy) )
    f_r["niter_global"] = aNumber(  doc="Each step is composed of optical response fit, followed \n"+
                                    " trajectory reoptimisation. All this niter_global times ",
                                    defaults2 = ("10", copy.deepcopy) )
    f_r["MPI_N_PROCS"] = aNumber(  doc="If bigger than 1, run it with mpi \n",
                                    defaults2 = ("1", copy.deepcopy) )




    f_r["result_file"]    = Hdf5FilePath(doc="The default hdf5 file where results will be stored\n must be a new file.\n the format is filename:groupname. The datagroup name must be a new one for the entry to be ready", defaults2 = (c_p["result_file"], copy.deepcopy),
                                         fileme=0.5, gnameme=0)
    f_r["result_file"].isresult = 1



    r_r = collections.OrderedDict()
    ## ======================================================================
    metodo["RESYNTHETISE_FIT"] = r_r
    r_r["HELP"]  =  """   Recreate the reference scan
    """

    r_r["getyaml"]  =  dic2yaml_resynt_fit
    r_r["swissknife_runner"]  =  swissknife_runner  
    r_r["method_name"] = "Scan Resynth."

    r_r["responsefilename"]      = Hdf5FilePath(doc="The  address on which the fitted responsed was written",
                                                     defaults2 = (f_r["result_file"]  , copy.deepcopy),
                                                     fileme=1, gnameme=1)

    r_r["response_scan_address"]      = Hdf5FilePath(doc="The original address of the scan which measures the response",
                                                     defaults2 = (f_r["response_scan_address"] , copy.deepcopy),
                                                     fileme=1, gnameme=1)

    r_r["result_file"]    = Hdf5FilePath(doc="The default hdf5 file where results will be stored\n"+
                                         "must be a new file.\n"+
                                         " The format is filename:groupname. \n"+
                                         "The datagroup name must be a new one\n"+
                                         "for the entry to be ready",
                                         defaults2 = ( (f_r["result_file"],"syntesis") , Functor_set_new_root()),
                                         fileme=0.5, gnameme=0)

    r_r["result_file"].isresult = 1
    
    return metodo 









