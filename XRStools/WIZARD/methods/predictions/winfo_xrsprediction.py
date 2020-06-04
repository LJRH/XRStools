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
from six.moves import map
from six.moves import zip
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
xrs_utilities = importlib.import_module( "XRStools"+version +".xrs_utilities"   )


###################################################################

import yaml
import glob


class functor_sigterm_handler:
    def __init__(self, p):
        self.p = p
    def __call__(self, *x):
        print( x)
        print( " KILLO " , self.p.pid)
        os.kill(self.p.pid, signal.SIGKILL)

def swissknife_runner( yamltext, where , ret_dico):
    inputname = os.path.join(where,"input.yaml")
    stdname = os.path.join(where,"stdout.txt")
    errname = os.path.join(where,"stderr.txt")
    
    ret_dico["stdout"] = stdname
    ret_dico["stderr"] = errname
    ret_dico["input"] = inputname
    
    
    open(inputname,"w").write(yamltext)
    # os.system("XRS_swissknife %s 1> %s  2>%s"%(inputname,stdname,  errname   ) )
    #p = subprocess.call(  string.split("XRS_swissknife %s 1> %s  2> %s"%(inputname,stdname,  errname ))    , shell=False, stdout= subprocess .PIPE ) #  preexec_fn=os.setsid
    # p= subprocess.Popen(string.split( "XRS_swissknife %s 1> %s  2>%s"%(inputname,stdname,  errname ))    , shell=False,  stdout= subprocess .PIPE ,  stderr= subprocess .PIPE  )
    stdout = open(stdname,"w")
    stderr = open(errname,"w")

    if not sys.platform=="win32":
        p= subprocess.Popen( ("XRS_swissknife%s %s 1> %s  2>%s"%(version,inputname,stdname,  errname )).split()    , shell=False,  stdout= stdout ,  stderr= stderr )
    else:
        packages = (os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(XRStools.__file__ ) ))))
        script = os.path.join(packages,"scripts","XRS_swissknife%s.bat"%version )
        print( " SCRIPT " , script)
        p= subprocess.Popen(( "%s %s"%(script, inputname)).split()    , shell=False,  stdout= stdout ,  stderr= stderr )

    # p= subprocess.Popen(string.split( "XRS_swissknife %s 1> %s  2>  %s"%(inputname,stdname,  errname ))    , shell=False,  stdout= stdout ,  stderr= stderr )

    # p= subprocess.Popen(string.split( "c:\\users\aless\\packaes\\scripts\\XRS_swissknife")    , shell=False,  stdout= stdout ,  stderr= stderr )
    # p= subprocess.Popen( "XRS_swissknife %s 1> %s  2>%s"%(inputname,stdname,  errname )    , shell=True )
    signal.signal(signal.SIGTERM, functor_sigterm_handler(p))
    p.communicate()
    ret_dico["return_code"]=p.returncode
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
                            
        return self.dico["CREATION_PARS"]["chem_formulas"].value+"__"



class many_ChemFormula(Parametro):
    def __init__(self, doc="" , nmin=1, nmax=100, defaults2=None, rendering = str, just_one = False):
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = None
        self.nmin = nmin
        self.nmax = nmax
        self.just_one = just_one

        data_dir = os.path.join(   os.path.dirname(   XRStools.__file__ )  ,  "data"     )
        afile = os.path.join(data_dir,"AtomicData.yaml")
        self.adata = yaml.load(open(afile,"r").read())


    def prevalid(self,value):
        value=value.strip()
        if not self.just_one:
            if value[0]!="[" or value[-1]!="]":
                return None
            value = value[1:-1]
            return value.split(",")
        else:
            return [value]
        
        
    def isReady(self):
        return self.isValid()

    def getValue(self):
        res = self.prevalid(self.value)
        return res
    def render(self):
        if not self.just_one:
            return self.rendering( self.getValue())
        else:
            return self.rendering( self.getValue()[0])
       
    def isValid(self):
        res  = self.prevalid( self.value)
        if res is None:
            return 0
        
        test = self.getValue()

        if (not self.just_one) and (not isinstance(test,list)):
            return 0
        
        n=len(test)
        
        if(not self.just_one) and ( n<self.nmin or n>self.nmax):
            return 0
        good = 1
        for t in test:
            elements,stoichiometries = xrs_utilities.parseformula(t)
            for el in elements:
                if el not in self.adata:
                    good=0
        return good
    
class Functor_list_from_list:
    def __init__(self,method):
        self.method=method
    def __call__(self, source):
        print( source)
        vals = source.getValue()
        if hasattr(self.method,"__call__"):
            print( vals)
            res = [self.method(t) for t in  vals    ]
        else:
            print( " num ", vals)
            res = [self.method for t in vals    ]
        return str(res)

class Functor_retrieveDensity:
    def __init__(self, fobj):
        self.fobj = fobj
        
    def __call__(self, t):
        
        elements,stoichiometries = xrs_utilities.parseformula(t)
        
        if len(elements)>1:
            return -1.0
        else:
            el = elements[0]
        return self.fobj.adata[el]["Density_g__ccm_"]
        
class Functor_retrieveMolarMasses:
    def __init__(self, fobj):
        self.fobj = fobj
    def __call__(self, t):
        elements,stoichiometries = xrs_utilities.parseformula(t)

        massa=0.0        
        for el,st in zip(elements,stoichiometries):
            massa += st*self.fobj.adata[el]["AtomicMass"]
        return massa


def dic2yaml_prediction(dicodic):
    s =     "XRSprediction :\n"
    s = s + "     sample :\n"
    s = s + "         chem_formulas  :   %s   \n"   % dicodic["sample"]["chem_formulas"].render()
    s = s + "         concentrations  :   %s   \n"   % dicodic["sample"]["concentrations"].render()
    s = s + "         densities  :   %s   \n"   % dicodic["sample"]["densities"].render()
    s = s + "         angle_tth  :   %s   \n"   % dicodic["sample"]["angle_tth"].render()
    s = s + "         sample_thickness  :   %s   \n"   % dicodic["sample"]["sample_thickness"].render()
    s = s + "         shape  :   %s   \n"   % dicodic["sample"]["shape"].render()
    if dicodic["sample"]["shape"].render()=="sphere":
        s = s + "         angle_in  :   %s   \n"   % None
        s = s + "         angle_out  :   %s   \n"   % None
    else:
        s = s + "         angle_in  :   %s   \n"   % dicodic["sample"]["angle_in"].render()
        s = s + "         angle_out  :   %s   \n"   % dicodic["sample"]["angle_out"].render()
    s = s + "         molar_masses   :   %s   \n"   % dicodic["sample"]["molar_masses"].render()

    s = s + "     incident_beam :\n"
    s = s + "         i0_intensity  :   %s   \n"   % dicodic["incident_beam"]["i0_intensity"].render()
    s = s + "         beam_height  :   %s   \n"   % dicodic["incident_beam"]["beam_height"].render()
    s = s + "         beam_width  :   %s   \n"   % dicodic["incident_beam"]["beam_width"].render()
   
    tavola = ((dicodic["analyzer"]["crystal_reflection"]).render()[len("chitable_"):])[:-len(".dat")]
    elemento=tavola
    numpart=""
    while(elemento[-1].isdigit()):
        numpart=elemento[-1]+numpart
        elemento = elemento[:-1]
    splitting_scheme = {3:[1,1,1], 4:[2,1,1],5:[2,2,1],6:[2,2,2]}[len(numpart)]
    ss = splitting_scheme
    indexes = [int(numpart[0:ss[0]]), int(numpart[ss[0]: ss[0]+ss[1]])   , int(numpart[ ss[0]+ss[1]: ss[0]+ss[1]+ss[2]]) ]

    elemento = elemento[0].upper()+elemento[1:]
    
    s = s + "     analyzer :\n"
     
    s = s + "         material  :   %s   \n"   %  elemento
    s = s + "         hkl       :   %s   \n"   %  indexes
    s = s + "         mask_d :   %s   \n"   % dicodic["analyzer"]["mask_d"].render()
    s = s + "         bend_r :   %s   \n"   % dicodic["analyzer"]["bend_r"].render()
    s = s + "         energy_resolution :   %s   \n"   % dicodic["analyzer"]["energy_resolution"].render()
    s = s + "         diced :   %s   \n"   % dicodic["analyzer"]["diced"].render()
    s = s + "         thickness :   %s   \n"   % dicodic["analyzer"]["thickness"].render()

    datadir = os.path.join(   os.path.dirname(   XRStools.__file__ )  ,  "data"     )
    
    s = s + "         database_dir :   %s   \n"   %  datadir

    
    
    s = s + "     compton_profiles  :\n"
    s = s + "         eloss_range  :   np.arange(%s,%s,%s)   \n"   %  (  dicodic["compton_profiles"]["eloss_start"].render(),
                                                                         dicodic["compton_profiles"]["eloss_end"].render(),
                                                                         dicodic["compton_profiles"]["eloss_step"].render())
    s = s + "         E0  :  %s  \n"   %  dicodic["compton_profiles"]["E0"] .render()
    
    s = s + "     detector  :\n"
    s = s + "         energy  :   %s   \n"      %   dicodic["detector"]["energy"] .render()
    s = s + "         thickness  :   %s   \n"   %   dicodic["detector"]["thickness"] .render()
    s = s + "         material  :   %s   \n"    %   dicodic["detector"]["material"] .render()

    s = s + "     thomson  :\n"
    s = s + "         scattering_plane  :   %s   \n"      %      dicodic["thomson"]["scattering_plane"] .render()
    s = s + "         polarization  :   %s   \n"   %   dicodic["thomson"]["polarization"] .render()

    saveaddress = dicodic["saveaddress"].render()
    if isinstance(saveaddress,str):
        if len(saveaddress):
            s = s + "     saveaddress : %s  \n" % saveaddress
    return s        

    
def getMethod(options):
    
    if "shift_the_reference" in options:
        shift_the_reference = options["shift_the_reference"]

    if "do_deconvolution" in options:
        do_deconvolution = options["do_deconvolution"]


    
    metodo = collections.OrderedDict()

    c_p = collections.OrderedDict()
    metodo["CREATION_PARS"] = c_p

    c_p["chem_formulas"]   = many_ChemFormula(   doc=" a list of formulas like [SiO2,C] ",
                                                 nmin=1,nmax=200,
                                                 defaults2=("[C]",simplecopy))

 
    c_p["define_new"] = Functor_adapt( metodo )



    ######################################################3
    xrs_p = collections.OrderedDict()
    metodo["xrs_prediction"] = xrs_p

    xrs_p["swissknife_runner"]  =  swissknife_runner  

    
    xrs_p["getyaml"]  = dic2yaml_prediction
    
    xrs_p["HELP"] = """
    This launches the XRS prediction routines.

    If yamlData contains information about: the sample, the incident beam, 
    the analyzer, the detector, the polarization, and the HF compton profiles,
    this will create the desired predicted XRS data.

    At the end you have the possibility to write the predicted profiles into a container hdf5 file.
    """

    
    sample =  collections.OrderedDict()
    metodo["xrs_prediction"]["sample"] = sample


    sample["chem_formulas"] =  many_ChemFormula(   doc=" a list of formulas like [SiO2,C] ",
                                                   nmin=1,nmax=200,
                                                   defaults2 = ( c_p["chem_formulas"], simplecopy  )
                                               )
    
    sample["concentrations"]  = many_Number(  "list of concentrations, should contain values between 0.0 and 1.0"   ,
                                              nmin=1, nmax=200, float_also = 1,
                                              defaults2=( sample["chem_formulas"] , Functor_list_from_list(  1.0   ))
                                          )
    
    
    
    sample["densities"]  = many_Number(  "list of densities of the constituents [g/cm^3]"   ,
                                         nmin=1, nmax=200, float_also = 1,
                                         defaults2=( sample["chem_formulas"] ,     Functor_list_from_list(  Functor_retrieveDensity(sample["chem_formulas"])   ) )
    )




    if 0:
        sample["angle_tth"]       = aNumber(doc = "scattering angle [deg]",
                                            defaults2=( "35.0" , simplecopy) ,
                                            float_also = True)
    else:

        sample["angle_tth"]       = many_Number(doc = "scattering angles [deg1,deg2, deg2], dont foret the square brakets or try if it works",
                                            defaults2=( "[35.0]" , simplecopy) ,
                                            float_also = True)
    
    
    sample["sample_thickness"]       = aNumber(doc = "sample thickness/diameter in [cm]",
                                        defaults2=( "0.1" , simplecopy) ,
                                        float_also = True)


    sample["shape"]           =   Parameter_Choices(  " keyword, can be 'slab' or 'sphere' " ,
                                                      ["sphere", "slab" ],
                                                      defaults2=( "sphere" , simplecopy)
                                                  )
    
    sample["angle_in"]       = aNumber(doc = "beam exit angle in [deg] relative to sample surface normal",
                                       defaults2=( "45.0" , simplecopy) ,
                                       float_also = True)

    sample["angle_out"]       = aNumber(doc = "incident beam angle in [deg] relative to sample surface normal",
                                        defaults2=( "45.0" , simplecopy) ,
                                        float_also = True)
    
    sample["shape"]   . visibility_depends_on = {("slab", ):[sample["angle_in"],  sample["angle_out"] ]   }
    
    sample["molar_masses"]  = many_Number(  "list of densities of the constituents [g/cm^3]"   ,
                                            nmin=1, nmax=200, float_also = 1,
                                        defaults2=( sample["chem_formulas"] ,     Functor_list_from_list(  Functor_retrieveMolarMasses(sample["chem_formulas"])   ) )
    )
    
  
      
    incident_beam =  collections.OrderedDict()
    metodo["xrs_prediction"]["incident_beam"] = incident_beam


    
    incident_beam["i0_intensity"]       = aNumber(doc = " number of incident photons [1/sec]",
                                                  defaults2=( "1.0e+13" , simplecopy) ,
                                                  float_also = True)
    incident_beam["beam_height"]       = aNumber(doc = "in micron",
                                                 defaults2=( "10.0" , simplecopy) ,
                                                 float_also = True)
    
    incident_beam["beam_width"]       = aNumber(doc = "in micron ",
                                                defaults2=( "20.0" , simplecopy) ,
                                                float_also = True)
    
    
    analyzer =  collections.OrderedDict()
    metodo["xrs_prediction"]["analyzer"] = analyzer

    #####
    datadir = os.path.join(   os.path.dirname(   XRStools.__file__ )  ,  "data"     )
    
    chitables = glob.glob(datadir+"/chitable_*.dat")
    chitables=list(map(os.path.basename,chitables))
    ####
    analyzer["crystal_reflection"]           =   Parameter_Choices(  "analyzer material (e.g. 'Si', 'Ge') and reflection.\n"
                                                                     " Possibilities are given by the chitables in data directory : " +datadir,
                                                                     chitables,
                                                                     defaults2=( "chitable_si660.dat" , simplecopy)
                                                                 )
    analyzer["mask_d"]       = aNumber(doc = "analyzer mask diameter in [mm] ",
                                        defaults2=( "60.0" , simplecopy) ,
                                        float_also = True)
    
    analyzer["bend_r"]       = aNumber(doc = "bending radius of the crystal [mm] ",
                                        defaults2=( "1.0" , simplecopy) ,
                                        float_also = True)
    
    analyzer["energy_resolution"]       = aNumber(doc = "energy resolution [eV] ",
                                        defaults2=( "0.5 " , simplecopy) ,
                                        float_also = True)
    
    analyzer["diced"]           =   Parameter_Choices(  "boolean (True or False) if a diced crystal is used or not (defalt is False)" ,
                                                       ["False","True"],
                                                       defaults2=( "False" , simplecopy))
                                                        
    analyzer["thickness"]       = aNumber(doc = " thickness of the analyzer crystal",
                                          defaults2=( "500.0" , simplecopy) ,
                                          float_also = True)
                                                                        

    compton_profiles =  collections.OrderedDict()
    metodo["xrs_prediction"]["compton_profiles"] = compton_profiles
    
    
    compton_profiles["eloss_start"]   = aNumber(doc = "eloss range start eV",
                                        defaults2=( "0.0" , simplecopy) ,
                                        float_also = True)


    compton_profiles["eloss_end"]     = aNumber(doc = " eloss range end eV",
                                        defaults2=( "1000.0" , simplecopy) ,
                                        float_also = True)


    compton_profiles["eloss_step"]     = aNumber(doc = " eloss range step eV ",
                                        defaults2=( "0.1" , simplecopy) ,
                                        float_also = True)


    compton_profiles["E0"]       = aNumber(doc = "E0 KeV ",
                                        defaults2=( "9.7" , simplecopy) ,
                                        float_also = True)



    detector =  collections.OrderedDict()
    metodo["xrs_prediction"]["detector"] = detector

    detector["energy"]       = aNumber(doc = "analyzer energy [keV] ",
                                        defaults2=( "9.7" , simplecopy) ,
                                        float_also = True)
    
    
   
    detector["thickness"]       = aNumber(doc = "thickness of the active material [microns] ",
                                        defaults2=( " 500.0" , simplecopy) ,
                                        float_also = True)
        
   
    detector["material"]           =   many_ChemFormula( doc="detector active material " , nmin=1, nmax=100, defaults2=( "Si" , simplecopy), rendering = str, just_one = True)

    thomson =  collections.OrderedDict()
    metodo["xrs_prediction"]["thomson"] = thomson

    
   
    thomson["scattering_plane"]           =   Parameter_Choices(  "keyword to indicate scattering plane relative to lab frame ('vertical' or 'horizontal')detector active material",
                                                                  ["vertical","horizontal"],
                                                                  defaults2=( "vertical" , simplecopy)
                                                              )


    
    thomson["polarization"]       = aNumber(doc = "degree of polarization (close to 1.0 for undulator radiation)",
                                            defaults2=( "0.99" , simplecopy) ,
                                            float_also = True)

    
    metodo["xrs_prediction"]["saveaddress"]   =  Hdf5FilePath(doc="the target destination, if data should be saved\n"
                                                              "let a blanck line if not",
                                                              defaults2 = (  "myfile.hdf5:/path/to/hdf5/group" , simplecopy ),
                                                              fileme=0.5, gnameme=0)



    return metodo







