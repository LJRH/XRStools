from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six import u

version=""

##########################################################################
# Copyright (C) 2016 European Synchrotron Radiation Facility
#  
#  European Synchrotron Radiation Facility, Grenoble,France
#  WIZARD.py
#  author : Alessandro Mirone 
#
# This program is free software; you can redistribute it and/or modify it 
# under the terms of BSD License
########*/
import os 
import sys

DEBUG = False

import tokenize
import imp

import re
import yaml
import yaml.resolver

yaml.resolver.Resolver
Resolver = yaml.resolver.Resolver
Resolver.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u(r'''^(?:[-+]?(?:[0-9][0-9_]*)(\.[0-9_]*)?(?:[eE][-+]?[0-9]+)?
                    |\.[0-9_]+(?:[eE][-+][0-9]+)?
                    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*
                    |[-+]?\.(?:inf|Inf|INF)
                    |\.(?:nan|NaN|NAN))$'''), re.X),
        list(u'-+0123456789.'))

os.environ['QT_API'] = 'pyqt'
timerPila = [0]

try:
    import sip
    sip.setapi("QString", 2)
    sip.setapi("QVariant", 2)
except:
    print( " SIP NO ")
    pass


import copy
import collections
from  PyQt4 import Qt, QtCore
import PyQt4
import PyMca5.PyMcaIO.specfilewrapper as SpecIO
import h5py
import pickle
from multiprocessing import Process, Manager
# manager_context=Manager()  poneva dei problemi in fase di import per via del forking
import subprocess
import os
import time
import tempfile
## import thread
import threading
import string
import signal
import sys, os
import traceback
from PyQt4 import QtGui
import numpy as np
import math
import matplotlib
matplotlib.use('Qt4Agg')
import pylab
from . import Wizard_safe_module


from XRStools.installation_dir import installation_dir
my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]

global animation_file 





global dico_copy_buffer
dico_copy_buffer = {}





def split_hdf5_address(dataadress):
    pos = dataadress.rfind(":")
    if ( pos==-1):
        return None
    filename, groupname = dataadress[:pos], dataadress[pos+1:]
    return filename, groupname 



def initialise_specmem( spec_mem,  specfile_name  ):
    new_size = os.path.getsize(specfile_name)
    spec_mem["sf"   ]    = None
    spec_mem["sf"   ]    = SpecIO.Specfile( specfile_name )
    spec_mem["size"   ]  = new_size
    scanlist = []
    for t in  spec_mem["sf"   ]:
        scanlist.append(t.number())
    spec_mem["scanlist"   ]  = scanlist
    
def  get_Specfile_Cached(  specfile_name, buffer ={}  ):

    if specfile_name in buffer:
        # print specfile_name , " eset in buffer " 
        spec_mem = buffer[specfile_name]

        old_size = spec_mem["size"]
        new_size = os.path.getsize(specfile_name)
        if new_size != old_size:            
            initialise_specmem( spec_mem,  specfile_name  )
        return spec_mem
    
    try:

        spec_mem = {}
        initialise_specmem( spec_mem,  specfile_name  )

        buffer[specfile_name] = spec_mem
        return spec_mem
    
    except:
        traceback.print_exc(file=sys.stdout)

        return None
    

def specScan_exists(specfile_name  , t  ):    
    # sf   = SpecIO.Specfile( specfile_name )

    sf = get_Specfile_Cached(  specfile_name  )

    if sf is None:
        return 0
    
    sf = sf["sf"]
    try:
        scan = sf.select(str(t))
        good=1
    except:
        good = 0

    return good


def render_numero_in_una_lista(value):
    return "[%s]" % value
def render_numero_in_una_lista_piu_uno(value):
    return "[%s,%d]"%(value, int(value)+1)

    
def RepresentsFloat(s):
    try:
        s="test : "+s
        test = yaml.load(s)["test"]
        
        if isinstance(test,float):
            return True
        else:
            return False
    except :
        return False


    
def RepresentsInt(s):
    try:
        s="test : "+s
        test = yaml.load(s)["test"]
        if isinstance(test,int):
            return True 
        else:
            return False
    except :
        return False

def reduce_to_value(par)   :
    if hasattr(par,"value" ):
        res = par.value
    else:
        res = par
    return res

class Functor_append_numinterval_to_name:
    def __init__(self, plus=""):
        self.plus = plus
        pass
    def __call__(self,arg):
        nome, intervallo = arg

        if hasattr(nome,"value"):
            nome=nome.value
        
        seque = intervallo.getValue()
        nome=nome+"_"
        for n,i in enumerate(seque):
            nome = nome+str(i)
            if n<len(seque)-1:
                nome=nome+"_"
                if n%2:
                    nome=nome+"_"
        if self.plus != "":
            nome = nome+"/"+self.plus            
        return nome


class Functor_add_datagroup:
    def __init__(self, datagroupname, remplace_last=False, transform_last = None ,  extract_numeric=None, keep_last_numbers=False, add_basename=None):
        self.datagroupname=datagroupname
        self.remplace_last =remplace_last
        self.transform_last = transform_last
        self.extract_numeric = extract_numeric
        self.keep_last_numbers  = keep_last_numbers
        self.add_basename = add_basename
        
    def __call__(self,arg):
        res=""
        datagroupname2=None

        if type(arg)==type(()):
            datagroupname2 = arg[1]
            arg = reduce_to_value(arg[0])
        else:
            arg = reduce_to_value(arg)
        
        arg = reduce_to_value(arg)

        datagroupname = self.datagroupname

        if datagroupname2 is not None:
            if type(datagroupname2) == type(""):
                s = datagroupname2
            else:
                s = datagroupname2.getValue()
                p=s.rfind("/")
                if p!=-1:
                    s=s[p+1:]
                p=s.rfind(":")
                if p!=-1:
                    s=s[p+1:]
                
            datagroupname =  datagroupname+"_"+ s

        if ":" in arg:
            
            last_numbers=""
            if self.keep_last_numbers :
                s = arg
                while len(s) and s[-1] in "_1234567890":
                    last_numbers =  s[-1] + last_numbers
                    s=s[:-1]
                datagroupname = datagroupname + last_numbers
                    
            
            tocut = self.remplace_last
            if self.transform_last is not None:
                tocut =  self.transform_last
                remainder = ""
            if self.extract_numeric is not None:
                tocut =  self.extract_numeric
                remainder = ""
            if tocut :
                for i in range(tocut):
                    pos = arg.rfind("/")
                    remainder = arg[pos+1:]
                    if pos>=0:
                        arg=arg[:pos]
                    else:
                        pos = arg.rfind(":")
                        arg=arg[:pos+1]


            if (self.extract_numeric is not None):
                s=""
                while(len(remainder)):
                    if remainder[-1].isdigit():
                        s=remainder[-1]+s
                        remainder = remainder[:-1]
                    else:
                        remainder = ""
                        return s
                return s
                        
            if len(datagroupname):
                if  self.transform_last is  None:
                    return   arg+"/"+datagroupname
                elif   (self.transform_last is not None)     :
                    s=""
                    while(len(remainder)):

                        if remainder[-1].isdigit():
                            s=remainder[-1]+s
                            remainder = remainder[:-1]
                        else:
                            remainder = ""
                    return   arg+"/"+datagroupname+"_"+s

                    
            else:
                return   arg
        else:
            if self.add_basename is not None:
                return   arg+"/"+ self.add_basename+  ":"+ datagroupname
            else:
                return   arg+":"+datagroupname
            

def simplecopy(par):
    if hasattr(par,"value" ):
        res = copy.deepcopy(par.value)
    else:
        res = copy.deepcopy(par)
    return res
    
class Parametro:
    tipo = "text"
    isresult=0
    mywarning=""
    enable=True
    automatic_forward = 0
    master = False
    always_visible = False
    # guiding = False
    
    def pump_default(self):
        if self.defaults2 is not None:
            source, method = self.defaults2
            self.value = method(source)

            
    def render(self):
        # print " A" , self
        # print self.rendering
        # print self.getValue()
        return self.rendering( self.getValue())
    def getValue(self):
        return self.value
    def getFullValue(self):
        return self.value
    def depends_from(self):
        return []

class   Functor_choicesFromSpecLabels:

    def __init__( self, spec_file_par , first_scan_par ,  last_scan_par  , labels=1):
        self.spec_file_par = spec_file_par
        self.first_scan_par = first_scan_par
        self.last_scan_par = last_scan_par
        self.labels = labels

    def __call__(self):
        filename = self.spec_file_par.render()
        s1 = int(self.first_scan_par.render() )
        s2 = int( self.last_scan_par.render() )
        sf = get_Specfile_Cached(  filename  )

        if sf is None:
            return []

        scan = sf["sf"].select(str(s1))

        if self.labels:    
            res = scan.alllabels()
        else:
            res = scan.allmotors()
        return res

    

    
class Parameter_Choices(Parametro):
    tipo ="multiplechoices"
    def __init__(self, doc="", choices = [], rendering = str, defaults2=None, choices_functor = None, initial_value = None):
        
        self.choices_functor = choices_functor
        self.rendering = rendering
        self.defaults2 = defaults2
        
        self.doc = doc
        self.value = initial_value
        self.choices = choices
        self.pump_default()
        
    def getIndex(self):
        if self.value not in self.choices:
            self.value = self.choices[0]
        index = self.choices.index(self.value)
            
        return index
        
    def isReady(self):
        return self.isValid()

    def isValid(self):
        return self.value in self.choices
   
    def getValue(self):
        return self.value


    

class FilePath(Parametro):
    tipo = "simplefile"
    isadir=0
    def __init__(self, doc="" , defaults2=None, rendering = str, fileme=1, isadir=0, initial_value = None):
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = initial_value
        self.fileme = fileme
        self.isadir=isadir
        
    def isReady(self):
        path=self.value
        if self.fileme==1:
            if self.isadir==0:
                return   self.isValid() and      os.path.exists(path ) and os.path.isfile(path)
            else:
                return   self.isValid() and      os.path.exists(path ) and os.path.isdir(path)
            
        elif self.fileme==0:
            return   self.isValid() and  not    os.path.exists(path )
        elif self.fileme==0.5:
            return   self.isValid()
    def isValid(self):
        path=self.value
        if  type(path)!=type("") or len(path)==0:
            return 0
        else:
            bd = os.path.dirname(path)
            if not os.path.exists(bd) or len(bd)==0:
                return 0
        return 1
    def getValue(self):
        return self.value


def hdf5path_exists(filename, groupname):
    # os.system("touch %s"%filename)
    # print " checco ", filename
    dn = os.path.dirname(filename)
    if dn=="":
        dn="./"
    os.stat( dn  )

    os.system("touch %s" % filename)
    
    
    h5 = h5py.File(filename,"r" )
    
    res = groupname in h5
    
    h5.close()
    
    return res

def is_hdf5(filename):
    try:
        h5 = h5py.File(filename,"r" )
        h5.close()
        return 1
    except:
        return 0

    
class Hdf5FilePath(FilePath):
    tipo = "simplefile hdf5"

    def __init__(self, doc="" , defaults2=None, rendering = str, fileme=1, gnameme=1, canbeNone=False, initial_value=None):
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = initial_value
        self.fileme = fileme
        self.gnameme = gnameme
        self.canbeNone = canbeNone



    def  hdf5path_exists(self):
        bf, bg = split_hdf5_address(self.getValue())

        if not    is_hdf5(bf):

            return 0
        if not hdf5path_exists(bf, bg):

            return 0

        return 1

    def fileexists(self):
        if self.canbeNone :
            if self.value=="None":
                return False
        
        if not self.isValid() :
            return 0
        filename, groupname = split_hdf5_address(self.value)
        return os.path.exists(filename ) 
        
        
    def isReady(self):
        if self.canbeNone :
            if self.value=="None":
                return True
        
        path=self.value
        if not self.isValid() :
            return 0
        filename, groupname = split_hdf5_address(path)
        (fileme, gnameme) = (self.fileme, self.gnameme)

        if fileme ==0.5:
            if  os.path.exists(filename ) and len(groupname):
                fileme = 1
            else:
                fileme=0
                gnameme=0
            
        if (fileme, gnameme)==(1,1):
            if not is_hdf5(filename):
                return 0
            if  os.path.exists(filename ):
                pass
            else:
                return 0
            return   (hdf5path_exists( filename, groupname  ))*1
        
        elif  (fileme, gnameme)==(0,0):
            
            if  os.path.exists(filename ):
                return 0
            else:
                return -1

        elif  (fileme, gnameme)==(1,0):
            if not is_hdf5(filename):
                return 0
            if  os.path.exists(filename ):
                return (not hdf5path_exists( filename, groupname  ))*(-1)
            else:
                return 0
        else:
            pass
            # raise Exception," invalid combination of fileme, gnameme in Hdf5FilePath definition "

    def corrected_dirname(self,path):
        pos = path.rfind(":")
        
        if pos==-1:
            return os.path.dirname(path)
        else:
            return os.path.dirname(path[:pos])
            
    def isValid(self):
        if self.canbeNone :
            if self.value=="None":
                return True
        
        path=self.value
        isvalid=1
        if type(path)!=type(""):
            return 0
        res=split_hdf5_address(path)
        if res is None:
            return 0

        path, groupname =res
        if  type(path)!=type("") or len(path)==0 or len(groupname)==0 :
            return 0
        else:
            bd = self.corrected_dirname(path)

            if not ( os.path.exists(bd) or len(bd)==0):
                return 0
            
        return 1



class NameWithNormaliser(Parametro):
    def __init__(self,  doc="" , defaults2=None, rendering = str, initial_value=""):
        self.doc = doc
        self.defaults2 = defaults2
        self.rendering = rendering
        self.pump_default()
        self.value=initial_value
        
    def getatypelikethis(self):
        return "s"
    
    def isReady(self):
        return self.isValid()

    def isValid(self):
        value=self.value
        if  type(value)!=type("") or len(value)==0:
            return 0
        else:
            if not len(value.split())==1:
                return 0
            if value.find("/")==-1:
                return 1
            pos =value.find("/")
            number = value[pos+1:]
            return RepresentsInt(number)
            
        
     
class hdf5_relative_path(Parametro):
    def __init__(self, base_h5file, doc="" , defaults2=None, rendering = str, gnameme=1, initial_value = None):
        self.gnameme = gnameme
        self.base_h5file  = base_h5file
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = initial_value
    def depends_from(self):
        return [ self.base_h5file ]
    def getatypelikethis(self):
        return "s"

    def fileexists(self):
        return self.base_h5file.fileexists()

    def getFullValue(self):
        res = self.base_h5file.getValue()+"/"+self.getValue()
        return res
    
    def isReady(self):
        path=self.value
        if not self.isValid():
            return 0
        bf, bg = split_hdf5_address(self.base_h5file.getValue())
        if not    is_hdf5(bf):
            return 0
        if not hdf5path_exists(bf, bg):
            return 0
        
        if self.gnameme==1:
            return  hdf5path_exists(bf, bg+"/"+path)
        else:

            if hdf5path_exists(bf, bg+"/"+path):
                return 0
            else:
                return -1
        
    def isValid(self):
        path=self.value
        if  type(path)!=type("") or len(path)==0:
            return 0
        else:
            return 1
    def getValue(self):
        return self.value


class aNumber(Parametro):
    def __init__(self, doc="", defaults2=None, rendering = str, float_also = False, nmin=-1.0e38,  nmax=1.0e38 , initial_value = None):
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = initial_value
        self.float_also = float_also
        self.nmin = nmin
        self.nmax = nmax
        
    def getatypelikethis(self):
        return "1"
    def isReady(self):
        return self.isValid()
    def isValid__(self):
        if self.float_also :
            return RepresentsInt(self.value) or  RepresentsFloat(self.value)
        else :
            return RepresentsInt(self.value) 

    def isValid(self):
        if not self.isValid__():
            return 0
        val = self.getValue()
        return (val>=self.nmin and val<=self.nmax)



    def getValue(self):
        if RepresentsInt(self.value):
            return int(self.value)
        else:
            return float(self.value)
        
class aNumberFloat(Parametro):
    def __init__(self, doc="", defaults2=None, rendering = str , initial_value = None):
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = initial_value
    def getatypelikethis(self):
        return "1.0"
    def isReady(self):
        return self.isValid()
    def isValid(self ):
        return RepresentsFloat(self.value)



    def getValue(self):
        return float(self.value)
    
def intervals2sequence(scan_interval):
    # print " scan_interval " , scan_interval
    todo_list = []
    ninterval = len(scan_interval)//2
    for i in range(ninterval):
        todo_list = todo_list + list(range( int(scan_interval[2*i]) , scan_interval[2*i+1]))
    return todo_list
  

class many_Number(Parametro):
    def __init__(self, doc="" , nmin=1, nmax=100, defaults2=None, rendering = str, isinterval = 0, float_also=0, initial_value = None):
        self.rendering = rendering
        self.defaults2 = defaults2
        self.doc = doc
        self.value = initial_value
        self.nmin = nmin
        self.nmax = nmax
        self.float_also = float_also
        self.isinterval=isinterval
        
    def getatypelikethis(self):
        return "[1,2]"
    
    def isReady(self):
        return self.isValid()

    def getValue(self):
        test=None
        try:
            ss="test="+self.value
            if Wizard_safe_module.python_code_is_safe(ss):
                resdic = {"math":math,"numpy":np,"np":np,"pylab":pylab}
                exec(ss,{"math":math,"numpy":np,"np":np,"pylab":pylab},resdic)
                test = resdic["test"]
            else:
                print( "cowardly refusing to execude suspicious code :", ss)
        except:
            pass
        return test
        
    def isValid(self):
        test = self.getValue()
        if not isinstance(test,list):
            return 0        
        n=len(test)
        if n<self.nmin or n>self.nmax:
            return 0
        good = 1
        for it,t in enumerate(test):
            if not ( isinstance(t, int) or  ( self.isinterval and (int(t)==t) and it%2==0  ) or (self.float_also and  isinstance(t, float)   ) ):
                good=0
        if good and self.isinterval:
            if ( n%2!=0):
                good = 0
        return good
    
def makeintervalfromfirst(value):
    value = reduce_to_value(value)
    test=None
    try:
        ss="test="+value
        if Wizard_safe_module.python_code_is_safe(ss):
            resdic={}
            exec(ss,{"math":math,"numpy":np,"np":np,"pylab":pylab},resdic)
            test = resdic["test"]
            val = test[0]
        else:    
            print( "cowardly refusing to execude suspicious code :", ss)
            return "unsafe code : "+value

        return "[%d,%d]"%(val,val+1)
    except:
        print( " failed in  makeintervalfromfirst with value = %s"% value)
        return value+" <<== PROBLEM "

    
class many_Number4SpecScan(many_Number):
    def __init__(self, specfile, doc="" , nmin=1, nmax=100, defaults2=None, rendering = str,isinterval=0, initial_value = None):
        self.rendering = rendering
        self.defaults2 = defaults2
        assert( isinstance(specfile, FilePath)    )
        self.doc = doc
        self.value = initial_value
        self.nmin = nmin
        self.nmax = nmax
        self.specfile=specfile
        self.isinterval=isinterval

                
    def isReady(self):
        specfile_name =self.specfile.value
        if not (  self.isValid() and  self.specfile.isReady(  ) ):
            return 0
        test = self.getValue()
        if self.isinterval:
            test=intervals2sequence(test)
        good = 1
        # print test
        for t in test[:4]+test[-4:] :
            if not specScan_exists(specfile_name  , t  ):
                good = 0
        return good
    
class aNumber4SpecScan(aNumber):
    def __init__(self, specfile, doc="", defaults2=None , rendering = str, float_also=0, nmin=-1.0e38,  nmax=1.0e38, initial_value = None):
        self.nmin = nmin
        self.nmax = nmax
        self.float_also = float_also
        self.rendering = rendering
        self.defaults2 = defaults2
        assert( isinstance(specfile, FilePath)    )
        self.doc = doc
        self.value = initial_value
        self.specfile=specfile

    def checkMyWarning(self):
        specfile_name =self.specfile.getValue()
        warning=""
        if not self.specfile.isReady(   ):
            if DEBUG : print(  "\nWARNING : cannot open file ", specfile_name)
            warning += "\nWARNING : cannot open file "
            self.mywarning = warning
            return
        if not self.isValid():
            if DEBUG :print( "\nWARNING : not a valid number")
            warning += "\nWARNING : not a valid number"
            # sf   = SpecIO.Specfile( specfile_name )

        if not self.isReady():
            if DEBUG :print( "\nWARNING : Scan Number is not in specfile")
            warning += "\nWARNING : Scan Number is not in specfile"

            
        sf_mem = get_Specfile_Cached(  specfile_name  )
        # print " sf_mem " , sf_mem
        if sf_mem is not None:
            sf = sf_mem["sf"]
            slist = sf_mem["scanlist"]
            warning += "\nNOTE : scans range from %d to %d\n" %(min(slist), max(slist))
            if DEBUG : print( warning)
            
        self.mywarning = warning

    def isReady(self):
        specfile_name =self.specfile.getValue()
        return  self.isValid() and  self.specfile.isReady(   )  and    specScan_exists(specfile_name  , self.getValue()   )


def DicRecursiveSearch(dico, key ):
    if key in dico:
        return dico[key]
    else:
        for k in dico.keys():
            d=dico[k]
            if isinstance(d,dict):
                res =  DicRecursiveSearch(d, key )
                if res is not None:
                    return res
    return None
    
def aggiorna_dics(reloaded , defined):
    for key, t in defined.items():
        if not ( key in reloaded):
            reloaded[key]=t
        else:
            r = reloaded[key] 
            if isinstance(t,dict):
                assert(    isinstance(r ,dict))
                aggiorna_dics(r , t)
            elif not isinstance(t, Parametro)  :
                reloaded[key]=t
            
class widget_base:
    def show_w(self):
        self.activity=1
        self.show()
        
    def hide_w(self):
        self.activity=0
        self.hide()

    def show(self):
        if self.activity==1:
            # print " SHOW ", self
            super(multichoices_widget, self).show()
    def hide(self):
        if self.activity==0:
            # print " HIDE ", self
            super(multichoices_widget, self).hide()

    def post_initialisation(self):
        
        if self.par.master:
            self.par.always_visible = True
            
        if self.par.master:
            # print " METTO BORDO OOOOOOOOOOOOOOOOO " 
            self.label.setStyleSheet("QLabel { border: 5px solid red ; }")
            
            # self.label.setStyleSheet("QLineEdit { background: rgb(200, 255, 200); selection-background-color: rgb(233, 99, 0); }")

    def textChangedFilter(self, nimportequoi):
        self.textChanged()
            
                
class multichoices_widget(Qt.QWidget, widget_base):
    def __getstate__(self):
        return {}
    
    def __init__(self, parent, name , par, dicodic, map_par2w=  {} ):
        self.map_par2w = map_par2w
        self.dicodic = dicodic
        # Qt.QDialog.__init__(self, parent)
        #Qt.QWidget.__init__(self, parent)

        super(multichoices_widget,self).__init__(parent)
        Qt.loadUi(  os.path.join(installation_dir  ,"resources" , my_relativ_path, "multichoices_widget.ui" ), self)

        self.name = name
        self.par  = par
        

        
        self.label.setText(name)
        
        self.comboBox.insertItems( 0,  par.choices )
        
        self.comboBox.connect(self.comboBox, QtCore.SIGNAL("activated(int)"), self.textChangedFilter)

        extraActions = {}
        self.comboBox.  contextMenuEvent  = Functor_contextMenuEvent(par, dicodic,self.comboBox , extraActions = extraActions )

        self.activity = 1
        if not hasattr(self.par,"visibility_depends_on"):
            self.par.visibility_depends_on = {}

        self.post_initialisation()
            
    def init_delayed(self):
        self.take_par_value( par)

    def take_par_value(self, par, no_automatic = 0 ):

        if self.par.choices_functor is not None:
            newchoices = self.par.choices_functor()
            if newchoices != self.par.choices:
                self.par.choices = newchoices
                self.comboBox.clear()
                self.comboBox.insertItems( 0, self.par.choices)
            
        
        if  str(self.comboBox.currentText()) != par.value:
            
            
            index = par.getIndex()
            
            self.comboBox.setCurrentIndex( index) 
        self.par.value = par.value
        self.textChanged(no_automatic)
        
    def textChanged(self, no_automatic = 0):
        s = str(self.comboBox.currentText())
        self.par.value=s
        valid = self.par.isValid()
        if valid:
            self.validentry_checkbox.setChecked(True)
        else:
            self.validentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()
        valid = self.par.isReady()
        if valid:
            self.readyentry_checkbox.setChecked(True)
        else:
            self.readyentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()

        value = self.par.value
        dico = self.par.visibility_depends_on
        for key, deps in dico.items():
                for d in deps:
                    if value in key:
                        if self.isVisible():
                            self.map_par2w[d].show_w()
                    else:
                        self.map_par2w[d].hide_w()
        
        # if self.par.automatic_forward and not no_automatic :
        #    self.lineEdit.  contextMenuEvent.propagateForward(self.par)


def dum():
    print (" gestate qlineedit")
    return {}


class filepath_widget(Qt.QWidget, widget_base):
    def __getstate__(self):
        return {}
    def __init__(self, parent, name , par, dicodic, map_par2w={}):
        self.map_par2w = map_par2w
        self.dicodic = dicodic
        #Qt.QWidget.__init__(self, parent)

        super(filepath_widget,self).__init__(parent)

        # Qt.QDialog.__init__(self, parent)
        Qt.loadUi(  os.path.join(installation_dir  ,"resources" , my_relativ_path, "filepath_widget.ui" ), self)

        self.name = name
        self.par  = par

        self.lineEdit.__getstate__ = dum
        self.label.setText(name)
        
        # self.browseButton.connect(self.browseButton, QtCore.SIGNAL("clicked()"), self.browse)
        # self.pymcaviewButton.connect(self.pymcaviewButton, QtCore.SIGNAL("clicked()"), self.pymcaview)
        self.lineEdit.connect(self.lineEdit, QtCore.SIGNAL("textChanged(const QString &)"), self.textChangedFilter)

        extraActions = collections.OrderedDict([("Browse",self.browse),("pymcaview" ,self.pymcaview)] )
        self.lineEdit.  contextMenuEvent  = Functor_contextMenuEvent(par, dicodic,self.lineEdit , extraActions = extraActions )
        self.activity = 1
        self.post_initialisation()

    def init_delayed(self):
        self.take_par_value( par)


        
    def take_par_value(self, par, no_automatic = 0 ):
        if  str(self.lineEdit.text()) != par.value:
            self.lineEdit.setText(par.value )
        self.par.value = par.value
        self.textChanged(no_automatic)
        
    def textChanged(self, no_automatic=0):
        
        diff = False
        
        if self.par.enable == False:
            s = self.par.value
            self.lineEdit.setText(s)
            #return
        else:
            s = str(self.lineEdit.text())
            diff = (s!= self.par.value)            
            self.par.value=s


        valid = self.par.isValid()
        if valid:
            self.validentry_checkbox.setChecked(True)
        else:
            self.validentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()
        valid = self.par.isReady()
        if valid:
            self.readyentry_checkbox.setChecked(True)
        else:
            self.readyentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()    

        if self.par.automatic_forward and not no_automatic :
            if diff :
                self.lineEdit.  contextMenuEvent.automatic_propagateForward(self.par)

        
    def pymcaview(self):
        filename = self.lineEdit.text()
        os.system("pymca -f %s &"%filename)
        
    def browse(self):

        if self.par.fileme==1:
            filename = self.lineEdit.text()
            if len(filename):
                if self.par.isadir :
                    filename = Qt.QFileDialog.getExistingDirectory(None,  "select" , filename)
                else:                    
                    filename =  Qt.QFileDialog.getOpenFileName(None, "select", filename)

                
            else:
                if self.par.isadir :
                    if DEBUG : print( " SETTO DirectoryOnly " )
                    filename = Qt.QFileDialog.getExistingDirectory(None,  "select" )
                else:
                    filename = Qt.QFileDialog.getOpenFileName(None, "select", )
                
        elif self.par.fileme==0:
            if self.par.isadir :
                filename =  Qt.QFileDialog.getExistingDirectory(None,  "select" )
            else:
                filename =  Qt.QFileDialog.getSaveFileName()
                
        elif self.par.fileme==0.5:
            if self.par.isadir :
                filename = Qt.QFileDialog.getExistingDirectory(None,  "select" )
            else:
                filename =  dialog.getSaveFileName()

        if isinstance(filename, tuple):
            filename = filename[0]


        filename=str(filename)
        if len(filename):
            self.lineEdit.setText(filename)

def     check_for_deps( papar, f):
    nm =  [a for a in dir(f) if not a.startswith('__')]
    tl = [ getattr(f,a ) for a in nm  ]
    return papar in tl

            
class Functor_contextMenuEvent:
    def __init__(self, par, dicodic, ql, extraActions):
        self.par = par
        self.dicodic = dicodic
        self.ql=ql
        self.extraActions = extraActions


    def automatic_propagateForward(self, papar=None):
        if not papar.isValid():
            return 
        if papar.automatic_forward==1:
            self.propagateForward( papar)
        elif papar.automatic_forward==2:
            self.propagateForwardRecursive( )
            
    def propagateForward(self, papar=None):
        if DEBUG : print( " propagateForward ", papar)
        timerS ,  global_dico   =  self.dicodic["refresher_stuff"]

        changed=[]
        if papar is  None:
            papar = self.par
        for nn,tt in global_dico.items():
            if nn == "CREATION_PARS":
                continue

            if isinstance(tt,collections.OrderedDict):
                stack=[iter(tt.items())]

                # for n,t in tt.items():
                while(len(stack)):
                    try:
                        name, t = next(stack[-1])
                        if isinstance(t, collections.OrderedDict):
                            stack.append( list(t.items()) )
                    except StopIteration:
                        stack = stack[:-1]

                    if isinstance(t,Parametro):
                        if hasattr(t,"defaults2"):
                            if t.defaults2 is not None:
                                aa,f = t.defaults2
                                if papar is aa or (  isinstance(aa,(list,tuple)) and papar in aa)  or papar in t.depends_from() or check_for_deps(papar , f):

                                    t. pump_default()
                                    changed.append(t)
        timerPila[0]+=1
        timerS.start()
        return changed



    
    def propagateForwardRecursive(self):
        changed_all=[]
        changed = self.propagateForward()
        while( len(set(changed)-set(changed_all))   ):
            changed_new = []
            for t in  (set(changed)-set(changed_all)):
                
                changed_new.extend   ( self.propagateForward(t) )
            changed_all.extend(changed )
            changed = changed_new

    def __call__(self, event):

        if not hasattr(self.ql,"createStandardContextMenu"):
            self.mem=[]

            testAction= Qt.QAction("Propagate Forward",None)
            self.ql.addAction(testAction);
            self.ql.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu )
            testAction.connect(testAction, Qt.SIGNAL("triggered()"),  self.propagateForward)
            self.mem.append(testAction)

            testAction= Qt.QAction("Propagate Forward Recursively",None)
            self.ql.addAction(testAction);
            self.ql.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu )
            testAction.connect(testAction, Qt.SIGNAL("triggered()"),  self.propagateForwardRecursive)
            self.mem.append(testAction)
            
        else:
            menu = self.ql.createStandardContextMenu()

            
            self.mem=[]
            menu.insertSeparator(menu.actions()[0])
            
            testAction = Qt.QAction("Propagate Forward",None);
            menu.insertAction(menu.actions()[0], testAction)
            testAction.connect(testAction, Qt.SIGNAL("triggered()"),  self.propagateForward)
            self.mem.append(testAction)
            
            testAction = Qt.QAction("Propagate Forward Recursively",None);
            menu.insertAction(menu.actions()[0], testAction)
            testAction.connect(testAction, Qt.SIGNAL("triggered()"),  self.propagateForwardRecursive)
            self.mem.append(testAction)
            
            menu.insertSeparator(menu.actions()[0])
            
            #menu.insertSeparator(menu.actions()[0])
            
            for na, aa in list(self.extraActions.items())[::-1]:
                testAction = Qt.QAction( na   ,None);
                self.mem.append(testAction)
                menu.insertAction(  menu.actions()[0]   , testAction)
                testAction.connect(testAction, Qt.SIGNAL("triggered()"),  aa )
            
                
            menu.exec_(event.globalPos());
            del menu



        
            

class hdf5filepath_widget(Qt.QWidget, widget_base):
    def __getstate__(self):
        return {}
    
    def __init__(self, parent, name , par, dicodic, map_par2w={}):
        self.map_par2w = map_par2w
        self.dicodic = dicodic

        super(hdf5filepath_widget,self).__init__(parent)

        # Qt.QWidget.__init__(self, parent)
        # Qt.QDialog.__init__(self, parent)
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path, "hdf5filepath_widget.ui" ), self)

        self.name = name
        self.par  = par

        self.label.setText(name)

        # self.browseButton.connect(self.browseButton, QtCore.SIGNAL("clicked()"), self.browse)
        # self.pymcaviewButton.connect(self.pymcaviewButton, QtCore.SIGNAL("clicked()"), self.pymcaview)
        # self.hdfviewButt<on.connect(self.hdfviewButton, QtCore.SIGNAL("clicked()"), self.hdfview)

        self.lineEdit.connect(self.lineEdit, QtCore.SIGNAL("textChanged(const QString &)"), self.textChangedFilter)
        # if par.fileme!=1 or par.gnameme==0:
        #     self.lineEdit.setStyleSheet("QLineEdit { background: rgb(255, 200, 200); selection-background-color: rgb(233, 99, 0); }")           
            

        extraActions = collections.OrderedDict([("Browse",self.browse),("pymcaview" ,self.pymcaview),("hdfview",self.hdfview ), ("To Console",self.toConsole) ] )
        self.lineEdit.  contextMenuEvent  = Functor_contextMenuEvent(par, dicodic,self.lineEdit , extraActions = extraActions )

        if False and self.par.enable == False:
            self.lineEdit.setEnabled(False)
            
        self.post_initialisation()

        
# void LineEdit::contextMenuEvent(QContextMenuEvent *event)
# {
#     QMenu *menu = createStandardContextMenu();
#     menu->addAction(tr("My Menu Item"));
#     //...
#     menu->exec(event->globalPos());
#     delete menu;
# }
                    
        # self.lineEdit.setFrameStyle(Qt.QFrame.Panel | Qt.QFrame.Sunken);
        self.activity = 1

    def init_delayed(self):
        self.take_par_value( par)


        
    def take_par_value(self, par, no_automatic = 0):
        if  str(self.lineEdit.text()) != par.value:
            self.lineEdit.setText(par.value )
        self.par.value = par.value
        self.textChanged(no_automatic)

    def pymcaview(self):
        filename = str(self.lineEdit.text())
        filename, groupname =  split_hdf5_address(filename)
        os.system("pymca -f %s &"%filename)
    def hdfview(self):
        filename = str(self.lineEdit.text())
        filename, groupname =  split_hdf5_address(filename)

        os.system("hdfview %s &"%filename)
        
    def browse(self):
        filename = str(self.lineEdit.text())
        fn, dg = hdf5_filedialog(filename, self.par.fileme)
        if fn is not None:
            filename=str(fn)
            if dg is not None:
                filename = filename +":"+str(dg)
            self.lineEdit.setText(filename)



    @staticmethod
    def myaction(fn, dg, varname="grabbed"):
        global MainWindow
        filename=str(fn)
        if dg is not None:
            groupname = str(dg)
            h5 = h5py.File(filename,"r")
            h5group = h5[groupname]

            runit, oggetto = Wizard_safe_module.make_obj( h5group ,h5,  groupname     )

            if runit is not None:
                runit(oggetto)

            update = {varname:oggetto}
            h5.close()
            if runit is None:
                if MainWindow.console is None:
                    msgBox = Qt.QMessageBox(None)
                    msgBox.setText("Launch a Console First")
                    msgBox.setInformativeText("You can launche the console from the mainwindow menu")
                    msgBox.setStandardButtons(  Qt. QMessageBox.Ok)
                    msgBox.setDefaultButton( Qt.QMessageBox.Ok)
                    ret = msgBox.exec_()
                    return
                MainWindow.console.updateNamespace(update)
                MainWindow.console.showMessage("# %s is now in namespace "% list(update.keys())[0])


            
    def toConsole(self):
 
        filename = str(self.lineEdit.text())
        fn, dg = hdf5_filedialog(filename, 1, self.myaction, modal=1)

        if fn is not None:
            self.myaction(fn, dg)


            
    def textChanged(self, no_automatic=0):
        diff = False
        if self.par.enable == False:
            s = self.par.value
            self.lineEdit.setText(s)
            # return
        else:
            s = str(self.lineEdit.text())
            diff = (s!= self.par.value)            

            self.par.value=s

            
        valid = self.par.isValid()
        if valid:
            self.validentry_checkbox.setChecked(True)
        else:
            self.validentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()
        valid = self.par.isReady()

        if valid:
            self.readyentry_checkbox.setChecked(True)
        else:
            self.readyentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()    

        if self.par.fileexists() and self.par.hdf5path_exists():
            warning = "\n STATUS : file plus  groupname point to an existing datagroup"
        else:
            if self.par.fileexists() :
                warning = "\n STATUS : file exists but  groupname dont"
            else:
                warning = "\n STATUS : file does not exists"
                
        self.lineEdit.setToolTip( self.par.doc+warning )

        if self.par.isresult:

            if self.par.fileexists() and self.par.hdf5path_exists():
                self.lineEdit.setStyleSheet("QLineEdit { background: rgb(200, 255, 200); selection-background-color: rgb(233, 99, 0); }")
            else:
                self.lineEdit.setStyleSheet("QLineEdit { background: rgb(255, 200, 200); selection-background-color: rgb(233, 99, 0); }")
            
            
            # if valid==1 or (not self.par.fileexists()) :
            #     self.lineEdit.setStyleSheet("QLineEdit { background: rgb(255, 200, 200); selection-background-color: rgb(233, 99, 0); }")
            #     if self.par.fileexists():
            #         warning = "\n STATUS : file seems  to  exists and groupname seems to be free  for writing over"
            #     else:
            #         warning = "\n STATUS : filename seems to be new"
            # elif (  self.par.fileexists()==1 ):
            #     self.lineEdit.setStyleSheet("QLineEdit { background: rgb(200, 255, 200); selection-background-color: rgb(233, 99, 0); }")
            #     warning = "\n STATUS : file plus  groupname point to an existing datagroup"
                


        if self.par.automatic_forward and not no_automatic :
            if diff :
                self.lineEdit.  contextMenuEvent.automatic_propagateForward(self.par)
                
        # lineEdit validentry_checkBox readyentry_checkBox  browseButton pymcaviewButton



class text_widget(Qt.QWidget, widget_base):
    def __getstate__(self):
        return {}
    
    def __getstate__(self):
        return {}
    
    def __init__(self, parent, name , par, dicodic , map_par2w={}):
        self.map_par2w = map_par2w
        self.dicodic = dicodic


        super(text_widget,self).__init__(parent)
        
        # Qt.QWidget.__init__(self, parent)
        # Qt.QDialog.__init__(self, parent)
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path, "text_widget.ui" ), self)

        self.name = name
        self.par  = par

        self.label.setText(name)
        self.lineEdit.connect(self.lineEdit, QtCore.SIGNAL("textChanged(const QString &)"), self.textChangedFilter)
 
        extraActions = collections.OrderedDict( )
        self.lineEdit.  contextMenuEvent  = Functor_contextMenuEvent(par, dicodic,self.lineEdit , extraActions = extraActions )
        self.activity = 1
        self.post_initialisation()
        
    def init_delayed(self):
        self.take_par_value( par)

    def take_par_value(self, par, no_automatic = 0):
        if  str(self.lineEdit.text()) != par.value:
            self.lineEdit.setText(par.value )
        self.par.value = par.value
        self.textChanged(no_automatic)
        
    def textChanged(self, no_automatic=0):
        diff = False
        if self.par.enable == False:
            # print " NOT ENABLED " 
            s = self.par.value
            self.lineEdit.setText(s)
            # return
        else:
            # print " IN TEXTCHANGED",  str(self.lineEdit.text()), self.par.value, self.par.value.__class__, self.par.automatic_forward, no_automatic
            s = str(self.lineEdit.text())
            
            diff = (s!= self.par.value)
            # print " diff " , diff
            self.par.value=s
            # if no_automatic not in [0,1,True, False]:
            #     raise


        valid = self.par.isValid()
        if valid:
            self.validentry_checkbox.setChecked(True)
        else:
            self.validentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()

        valid = self.par.isReady()
        if valid:
            self.readyentry_checkbox.setChecked(True)
        else:
            self.readyentry_checkbox.setChecked(False)
        self.readyentry_checkbox.show()    
        
        if isinstance(self.par, hdf5_relative_path):
            if self.par.isresult:
                if(valid==-1 or not self.par.fileexists() ):
                    self.lineEdit.setStyleSheet("QLineEdit { background: rgb(255, 200, 200); selection-background-color: rgb(233, 99, 0); }")
                    warning = "\nSTATUS :  root path seems to be already there and relative path seems to be free "   
                else:
                    if self.par.base_h5file.hdf5path_exists():
                        self.lineEdit.setStyleSheet("QLineEdit { background: rgb(200, 255, 200); selection-background-color: rgb(233, 99, 0); }")
                        warning = "\nSTATUS :  root path seems to be already there and relative path seems to be exist already "   
                    else:
                        self.lineEdit.setStyleSheet("QLineEdit { background: rgb(255, 20, 20); selection-background-color: rgb(233, 99, 0); }")
                        if self.par.base_h5file.fileexists():
                            warning = "\nSTATUS :  relative path seems to be non writable :  the file exists but the root path is not ready"
                        else:
                            warning = "\nSTATUS :  relative path seems to be non writable :  the file does not exist"
                            
                self.lineEdit.setToolTip( self.par.doc+warning )
        elif hasattr(self.par, "checkMyWarning"):
            self.par.checkMyWarning()
            # print " MYWARNING " , self.par.mywarning 
            if self.par.mywarning!="":
                # print self.par.doc+self.par.mywarning 
                self.lineEdit.setToolTip( self.par.doc+self.par.mywarning )

        if self.par.automatic_forward and not no_automatic :
            # print " forse propagate ", self.par, diff
            if diff :
                self.lineEdit.  contextMenuEvent.automatic_propagateForward(self.par)
                


class creationtab(Qt.QWidget):
    def __getstate__(self):
        return {}
    
    def __init__(self, parent):
        self.parent=parent
        # Qt.QWidget.__init__(self, parent)

        super(creationtab,self).__init__(parent)

        # Qt.QDialog.__init__(self, parent)
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path, "creationtab.ui" ), self)


        
class Functor_Copy_Dict:
    def __init__(self,dic):
        self.dic = dic
    def __call__(self):
        global dico_copy_buffer
        dico_copy_buffer = self.dic


def  GetYamlDicoTraduc( depth, subdic  ) :
    if depth==0:
        foo = subdic["getyaml"]
        if foo.__class__ == tuple:            
            foo, traduc = foo
            dico={}
            for key in subdic:
                dico[ traduc.get(key,key) ] = subdic[key]
        else:
            dico = subdic
        yamltext = foo(dico)
    else:
        yamltext=""
        for c,dico in subdic.items():
            if isinstance(dico, dict) and dico.has_key("getyaml"):
                yamltext=yamltext+ GetYamlDicoTraduc( 0,   dico  ) +"\n"
    return yamltext

    
class Functor_Show_Yaml:
    def __init__(self, subdic, depth=0):
        self.subdic=subdic
        self.depth=depth
    def __call__(self):


        yamltext =  GetYamlDicoTraduc( self.depth, self.subdic  ) 
        

        self.mle = Qt.QTextEdit(None)
        self.mle.setLineWrapMode(0)
        self.mle.setText(  yamltext  )
        self.mle.show()
        



class MyProcess(Process):
    def __init__(self, target,args, toimport):
        Process.__init__(self)

        self.toimport = toimport
        self.target = pickle.dumps( target)
        self.args = args
        

    def run(self):
        if DEBUG : print( 'run called...')
        if DEBUG : print( "importo ", self.toimport)
        module=imp.load_source(self.toimport[0], self.toimport[1])
        self.target = pickle.loads(self.target)
        self.target(*self.args)

        # print 'running...', module.step()


class Functor_Run_Yaml:
    Runthis=0
    Stopthis=1
    Viewthis=2
    Viewthiserr=3
    Center=4
    Pendings=5
    PlugOut=6
    ident =0 
    
    def __init__(self, subdic, depth=0, action = 0,sub_functors_list=[]  , identity = None ):
        
        self.plugmeout    = ( action == self.PlugOut  )
        self.run_pendings = ( action == self.Pendings )
        
        self.subdic=subdic
        self.depth=depth
        self.action = action
        self.sub_functors_list = sub_functors_list
        self.identity = identity
        self.subdic["multifunctors_isstopped"] = 0
        if identity is not None:
            self.identity_number, self.identity_widget  =  identity
        
        ## self.process={"mp":None, "dir":None}

        
    def viewthis(self, stderr=0):
        if len(self.sub_functors_list): return 
        
        shadok_movie, shadok_action, shadok_returncode,  run_informations, tb = self.subdic["shadok_data"]
        process_directory = run_informations["process_directory"]
        
        if process_directory is None:
            self.identity_widget.shadok_message_signal.emit("No information about previous runs")
            return

        
        p, where  = process_directory
        if stderr==0:
            text  = open(os.path.join(where,"stdout.txt") ,"r").read()
        else:
            text  = open(os.path.join(where,"stderr.txt") ,"r").read()
        if stderr==0:
            title = 'stdout at %s'%where
        else:
            title = 'stderr at %s'%where
        self.identity_widget.shadok_textshow_signal.emit( title   , text )

        
    def __call__(self):

        if self.action == self.Runthis:
            if len(self.sub_functors_list):


                t = threading.Thread(target=self.runthis, args=())
                # threads.append(t)
                t.start()
                
                ## thread.start_new_thread( self.runthis, () )
            else:
                self.runthis()
        elif self.action in [ self.PlugOut, self.Pendings ]:
             self.runthis()
        elif self.action == self.Center:
            self.doCenter()
        elif self.action == self.Stopthis:
            self.stopthis()
        elif self.action == self.Viewthis:
            self.viewthis(0)
        elif self.action == self.Viewthiserr:
            self.viewthis(1)


            
    def get_yaml(self):

        yamltext = GetYamlDicoTraduc( self.depth ,  self. subdic  ) 

        # if self.depth==0:
        #     yamltext = self.subdic["getyaml"](self.subdic)
        # else:
        #     yamltext=""
        #     for c,dico in self.subdic.items():
        #         if isinstance(dico, dict) and dico.has_key("getyaml"):
        #             yamltext=yamltext+dico["getyaml"](dico)+"\n"
        
        return yamltext


    
    def stopthis(self):
        if DEBUG : print( " IN STOP THIS ")
        if len(self.sub_functors_list):
            if DEBUG : print( " PER TUTTI ")
            self.subdic["multifunctors_isstopped"]=1
            msgBox = Qt.QMessageBox(None)
            msgBox.setText("Do you really want to stop all processes?")
            msgBox.setInformativeText("Do you really want to stop all processes?")
            msgBox.setStandardButtons( Qt.QMessageBox.Cancel | Qt. QMessageBox.Ok)
            msgBox.setDefaultButton( Qt.QMessageBox.Cancel)
            ret = msgBox.exec_()
            if not ret:
                return 
            for f in self.sub_functors_list:
                f.stopthis()
            return
        
        shadok_movie, shadok_action, shadok_returncode,  run_informations, tb = self.subdic["shadok_data"]
        process_directory = run_informations["process_directory"]
        if process_directory is not None:
            p, where = process_directory
            if p is not None:
                os.kill(p.pid, signal.SIGTERM)
                # p.terminate()


    def isrunning(self):
        
        shadok_main_movie, shadok_main_action, shadok_returncode,  run_informations, tb = self.subdic["shadok_data"]
        pd = run_informations["process_directory"]
        if pd is not None:
            p, where =pd
            if p is not None:
                return 1
        return 0

    def get_method_name(self, yamltext ):
        pos = yamltext.find(":")
        if pos !=-1:
            yamltext= yamltext[:pos]
        yamltext=yamltext.strip()
        return yamltext
    
    def runthis(self):
        global MainWindow
       
        if self.plugmeout:
            if DEBUG : print( " PLUGOUT")
            global MainWindow
            copia = copy.copy(self)
            copia.plugmeout=0
            update = {"grabbed":copia}
            copia.sysout = sys.stdout
            MainWindow.console.updateNamespace(update)
            MainWindow.console.showMessage("The new Functor object is available in the namespace with name grabbed")
            return
        
        if self.run_pendings:
            shadok_main_movie, shadok_main_action, shadok_returncode,  run_informations, tb = self.subdic["shadok_data"]
            sub_functors_list = run_informations["sub_functors_list"]
        else:
            sub_functors_list = self.sub_functors_list
 
        if len(sub_functors_list):
            self.subdic["multifunctors_isstopped"] = 0
            for f in sub_functors_list:
                if self.run_pendings:
                    
                    if f.identity < self.identity:
                        continue
                    todo = 0
                    for k,t in f.subdic.items():
                        if isinstance(t,Parametro):
                            if t.isresult:
                                nameval = t.getValue()
                                if t.isReady():
                                    todo=1
                    if not todo :
                        continue
                    
                if DEBUG :print (" runno " , f.identity, f)
                ## print " runno " , f, len(  self.sub_functors_list   )
                f.runthis()
                while(f.isrunning()):
                    # print " is running " 
                    time.sleep(0.2)
                    if self.subdic["multifunctors_isstopped"]:
                        # self.subdic["multifunctors_isstopped"]=0
                        return
            if DEBUG :print( " RETURN " )
            return

                   
        yamltext = self.get_yaml() ##
        methodname = self.get_method_name(yamltext)
        
        if self.isrunning():
            if methodname != "view_Volume_myavi" :
                self.identity_widget.shadok_message_signal.emit("This Job is already running %s"%methodname)
                return 


        if hasattr(self,"sysout"):
            oldsysout = sys.stdout
            sys.stdout = self.sysout
 
        
        where = os.path.join("wizardRUNS", methodname)
        if not os.path.exists(where):
            os.makedirs( where)
        where = tempfile.mkdtemp(suffix='', prefix='tmp', dir=where)

        self.manager_context=Manager()
        ret_dico = self.manager_context.dict()

        self.ret_dico=ret_dico

        self.identity_widget.shadok_returncode_signal.emit(self.identity_number, 0 )
        
        # p = Process(target=self.subdic["swissknife_runner"] , args=(yamltext,where, ret_dico))

        timerS ,  global_dico   =  self.subdic["refresher_stuff"]


        p = MyProcess(self.subdic["swissknife_runner"] , (yamltext,where, ret_dico),  global_dico[  "tobeimported"  ]  )
        p.start()

        if hasattr(self,"sysout"):
            sys.stdout = oldsysout
            
        
        shadok_main_movie, shadok_main_action, shadok_returncode,  run_informations, tb = self.subdic["shadok_data"]
        run_informations["process_directory"] = (p,where)

        t = threading.Thread(target=self.ProcessMonitor, args=())
        # threads.append(t)
        t.start()
        # thread.start_new_thread( self.ProcessMonitor, () )

        
        self.identity_widget.shadok_start_signal.emit(self.identity_number)

    def doCenter(self):
        if len(self.sub_functors_list):
            return


        self.identity_widget.shadok_center_signal.emit(self.identity_number)
        return 

        
    def ProcessMonitor(self):
        shadok_movie, shadok_action, shadok_returncode,  run_informations, tb = self.subdic["shadok_data"]
        timerS ,  global_dico   =  self.subdic["refresher_stuff"]
        p,where = run_informations["process_directory"]
        p.join()
        run_informations["process_directory"]=[None,where]
        if DEBUG :print( " JOINED ")
        if DEBUG :print (" STOPPO " , self.identity_number)
        self.identity_widget.shadok_stop_signal.emit(self.identity_number)

        return_code = self.ret_dico["return_code"]
        
        self.identity_widget.shadok_returncode_signal.emit(self.identity_number, return_code )

        if DEBUG :print (" il codice di ritorno est ", return_code)

        if return_code==0:
            if not MainWindow.actionOption.isChecked():
                if DEBUG :print (" cancello " , [self.ret_dico["input"], self.ret_dico["stdout"], self.ret_dico["stderr"]])
                for filename in [self.ret_dico["input"], self.ret_dico["stdout"], self.ret_dico["stderr"]]:
                    if DEBUG :print( filename)
                    if os.path.exists(filename):
                        os.remove(filename)    

        
        ## lanciare una thread per join. thread fa anche refresh e stop shadok dopo join
        ##  Metter anche in stop refresh e stop shadok. Guardare che se stop anche sia join
        

# def sleeper(name, seconds):
#    print 'starting child process with id: ', os.getpid()
#    print 'parent process:', os.getppid()
#    print 'sleeping for %s ' % seconds
#    time.sleep(seconds)
#    print "Done sleeping"

 
# if __name__ == '__main__':
#    print "in parent process (id %s)" % os.getpid()
#    p = Process(target=sleeper, args=('bob', 5))
#    p.start()
#    print dir(p)
#    print p.is_alive()
#    print p.terminate.__doc__
#    print p.terminate()
#    raw_input()
#    print p.is_alive()
#    raw_input()
#    print "in parent process after child process start"
#    print "parent process about to join child process"
#    p.join()
#    print p.is_alive()
#    print "in parent process after child process join" 
#    print "parent process exiting with id ", os.getpid()
#    print "The parent's parent process:", os.getppid()




        

class Functor_Run_Yaml_and_pendings:
    def __init__(self, subdic):
        self.subdic=subdic
    def __call__(self):

        self.mle = Qt.QTextEdit(None)
        self.mle.setText( str(self.subdic)   )
        self.mle.show()


def put_shadok_in_toolbar(toolBar):
    global animation_file 
    movie = Qt.QMovie(animation_file);
    processLabel = Qt.QLabel(toolBar)
    processLabel.setToolTip("Les Shadoks\n  at work \n  for you")
    processLabel.setMovie(movie)
    action = toolBar.addWidget(processLabel)

    
    movie.start()
    movie.stop()
    action.setVisible(False)


    processLabel = Qt.QLabel(toolBar)
    processLabel.setText("ERROR")
    # processLabel.setToolTip("No return code yet")
    returncode = toolBar.addWidget(processLabel)
    returncode.setVisible(False)
    
    
    return movie, action, (returncode, processLabel)


        

class creationtab2(Qt.QWidget):

    shadok_start_signal = QtCore.pyqtSignal(int)
    shadok_stop_signal = QtCore.pyqtSignal(int)
    shadok_returncode_signal = QtCore.pyqtSignal(int,int)
    shadok_message_signal = QtCore.pyqtSignal(str)
    shadok_textshow_signal = QtCore.pyqtSignal(str,str)
    shadok_center_signal = QtCore.pyqtSignal(int)
    
    def __getstate__(self):
        return {}

    def __init__(self, parent, dico = None):
        self.parent=parent


        super(creationtab2,self).__init__(parent)

        
        # Qt.QWidget.__init__(self, parent)
        # Qt.QDialog.__init__(self, parent)
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path, "creationtab2.ui" ), self)
        
        self.shadok_dictios=[]

        self.shadok_start_signal.connect(self.start_shadok)
        self.shadok_center_signal.connect(self.center_shadok)
        self.shadok_stop_signal.connect (self.stop_shadok )
        self.shadok_returncode_signal.connect (self.returncode_shadok )
        self.shadok_message_signal.connect (self.message_shadok )
        self.shadok_textshow_signal.connect (self.textshow_shadok )

        if dico is None:
            return
    def textshow_shadok(self,title, message):
        self.mle = Qt.QTextEdit(None)
        self.mle.setLineWrapMode(0)
        self.mle.setText( message   )
        self.mle.setWindowTitle(title)
        self.mle.show()                
        
    def message_shadok(self, message):
        msgBox = Qt.QMessageBox(None)
        msgBox.setText(message)
        msgBox.setInformativeText(message)
        msgBox.setStandardButtons(  Qt. QMessageBox.Ok)
        msgBox.setDefaultButton( Qt.QMessageBox.Ok)
        ret = msgBox.exec_()
        
        
    def returncode_shadok(self, identity,returncode):
        if DEBUG :print( "returncode SHADOK" , returncode)
        subdic = self.shadok_dictios[identity]
        for sd in ["shadok_data","shadok_data_bis" ] :
            shadok_movie, shadok_action, (shadok_returncode, labello) = subdic[sd][:3]
            timerS ,  global_dico   =  subdic["refresher_stuff"]
            if returncode:
                labello.setText("!!! attention : ERROR : %d "% returncode)
                shadok_returncode.setVisible(True)
            else:
                shadok_returncode.setVisible(False)

    def stop_shadok(self, identity):
        if DEBUG :print( "STOP SHADOK" , identity)
        subdic = self.shadok_dictios[identity]
        shadok_movie, shadok_action, shadok_returncode,  run_informations, tb  = subdic["shadok_data"]
        timerS ,  global_dico   =  subdic["refresher_stuff"]
        timerS.start()
        shadok_movie.stop()
        shadok_action.setVisible(False)
        shadok_movie, shadok_action, shadok_returncode,  tb  = subdic["shadok_data_bis"]
        shadok_movie.stop()
        shadok_action.setVisible(False)
       
    def center_shadok(self, identity):
        if DEBUG :print( "CENTER  SHADOK" , identity)
        subdic = self.shadok_dictios[identity]
        shadok_movie, shadok_action, shadok_returncode,  run_informations, tb = subdic["shadok_data"]
        self.scrollArea.ensureWidgetVisible(tb)
        
    def start_shadok(self, identity):
        if DEBUG :print( "START SHADOK" , identity)
        subdic = self.shadok_dictios[identity]
        shadok_movie, shadok_action, shadok_returncode,  run_informations, tb = subdic["shadok_data"]
        timerS ,  global_dico   =  subdic["refresher_stuff"]
        shadok_movie.start()
        shadok_action.setVisible(True)
        shadok_movie, shadok_action, shadok_returncode,   tb = subdic["shadok_data_bis"]
        shadok_movie.start()
        shadok_action.setVisible(True)

    def run_everything(self):
        if DEBUG :print( " run everything ")

    def refresh_everything(self):
        self.timerS.start()

    def refresh(self):
        if self.visibleRegion().isEmpty():
            return 
            # print " NO REFRESH non est visibile "
        else:
            pass
            # print " REFRESH "
        dicodic  = self.dico
        for edw in self.edws:
            if edw.par.value is None:
                edw.par.value = ""
            if type(edw.par.value) != type(""):
                if True or DEBUG :
                    print( " in refresh normale " , edw.par.value)
                edw.par.value="problem"
            else:
                edw.take_par_value(edw.par, no_automatic=0 )
                
    def refresh_from_forward(self):
        dicodic  = self.dico
        for edw in self.edws:
            if type(edw.par.value) != type(""):
                if True or DEBUG :
                    print( " in refresh forward ",  edw.par.value)
                edw.par.value="problem"
            else:
                edw.take_par_value(edw.par, no_automatic=1 )

    def save_parameters(self):
        filename = (Qt.QFileDialog.getSaveFileName(None, "select"))
        if isinstance(filename, tuple):
            filename = filename[0]
        filename = str(filename)
        if len(filename ):
            f=open(filename, "wb" )
            s = pickle.dumps(  self.dico)
            #s =  s.replace( "XRStools"+version+".", "XRStools." )
            f.write(s)
            # pickle.dump(  self.dico, f )
            f.close()

    def toggle_deps_visibility(self):
        if self.deps_vis_toggle_status:
            for w in self.w2toggle:
                w.hide_w()
        else:
            for w in self.w2toggle:
                w.show_w()
                
        self.deps_vis_toggle_status = 1- self.deps_vis_toggle_status
            
    def populate(self, dico, guided=False):

        self.guided=guided
        self.timer = Qt.QTimer()
        self.timerS = Qt.QTimer()
        
        self.dico=dico

        menu = Qt.QMenu()
        gAction = Qt.QAction("Run All", self);
        menu.addAction(gAction)
        
        gAction_stop = Qt.QAction("Stop All", self);
        menu.addAction(gAction_stop)
        
        gAction2 = Qt.QAction("Refresh", self);
        menu.addAction(gAction2)
        
        gAction3 = Qt.QAction("Display Yaml", self);
        menu.addAction(gAction3)

        gAction4 = Qt.QAction("Copy the dictionary", self);
        menu.addAction(gAction4)
        
        if guided:
            gAction_deps = Qt.QAction("Show/Hide depending parameters", self);
            menu.addAction(gAction_deps)
            self.w2toggle=[]
            self.deps_vis_toggle_status=0
            
        toolButton = Qt.QToolButton()
        toolButton.setMenu(menu)
        toolButton.setPopupMode(Qt.QToolButton.InstantPopup)
        toolButton.setToolButtonStyle(  2)   ## 1 per ToolButtonTextOnly
        toolButton.setText(  "Global Actions"   )
        toolButton.setIcon(  Qt.QIcon(Qt.QPixmap(bacchetta32_xpm))     )
        
        toolBar = Qt.QToolBar()
        toolBar.addWidget(toolButton)

        but =  Qt.QToolBar()
        toolBar.addWidget(but)

        shadok_main_movie, shadok_main_action, shadok_returncode  = put_shadok_in_toolbar(toolBar )
        self.verticalLayout_3.addWidget(toolBar )
        
        run_informations = {"ALL_dic": self.dico, "process_directory":None}
        self.dico["shadok_data"] = [shadok_main_movie, shadok_main_action, shadok_returncode,  run_informations, None ]
        identity = (len(self.shadok_dictios), self)
        self.shadok_dictios.append(self.dico)
        sub_functors_list = []
        gAction      .connect(gAction , Qt.SIGNAL("triggered()"),     Functor_Run_Yaml(self.dico , depth=1, sub_functors_list = sub_functors_list, action =  Functor_Run_Yaml.Runthis ,identity=identity ) )
        gAction_stop .connect(gAction_stop , Qt.SIGNAL("triggered()"),     Functor_Run_Yaml(self.dico , depth=1, sub_functors_list = sub_functors_list, action =  Functor_Run_Yaml.Stopthis,identity=identity ) )
        
        gAction2.connect(gAction2, Qt.SIGNAL("triggered()"),   self.refresh_everything  )
        gAction3.connect(gAction3, Qt.SIGNAL("triggered()"),   Functor_Show_Yaml( self.dico, depth=1))
        gAction4.connect(gAction4, Qt.SIGNAL("triggered()"),   Functor_Copy_Dict( self.dico))

        if self.guided:
            gAction_deps.connect(gAction_deps, Qt.SIGNAL("triggered()"), self.toggle_deps_visibility   )

        chiavi = list(self.dico.keys())[1:]
        if DEBUG :print( " CHIAVI ", chiavi)
        self.edws=[]
        for c in chiavi :

            if DEBUG :print( " CONSIDERO ", c)
            subdic = self.dico[c]
            if not isinstance(subdic, dict):
                continue
            #######################################
            
            identity = (len(self.shadok_dictios), self)
            
            for frame in range(2):
                menu = Qt.QMenu()
                testAction_center = Qt.QAction("Center", self);
                menu.addAction(testAction_center)
                
                testAction = Qt.QAction("display yaml inputs", self);
                menu.addAction(testAction)
                testAction2 = Qt.QAction("Run this", self);
                menu.addAction(testAction2)
                testAction3 = Qt.QAction("Stop this", self);
                menu.addAction(testAction3)
                testAction4 = Qt.QAction("View StdOutput", self);
                menu.addAction(testAction4)
                testAction5 = Qt.QAction("View StdErr", self);
                menu.addAction(testAction5)
                # testAction6 = Qt.QAction("Run this plus pending ones", self);
                # menu.addAction(testAction6)

                testAction7 = Qt.QAction("Plug out to Console", self);
                menu.addAction(testAction7)
                
                # toolButton = Qt.QToolButton()
                this = self
                class MyQToolButton(Qt.QToolButton):
                    mywidget = self
                    mysubdic = subdic
                    # def mouseMoveEvent(self, event):
                    #     self.mywidget.scrollArea.ensureWidgetVisible(self.mysubdic["shadok_data"][-1] )
                        
                toolButton = MyQToolButton()
                toolButton.setMenu(menu)
                toolButton.setPopupMode(Qt.QToolButton.InstantPopup)
                toolButton.setToolButtonStyle(  2)   ## 1 per ToolButtonTextOnly
                toolButton.setText(  c   )
                toolButton.setIcon(  Qt.QIcon(Qt.QPixmap(bacchetta32_xpm))     )
                toolButton.setToolTip( subdic["HELP"]   )
                toolBar = Qt.QToolBar(self)
                toolBar.addWidget(toolButton)
                if frame :
                    toolButton.setMouseTracking(True)
                if frame==0:
                    self.verticalLayout_5.addWidget( toolBar )
                    shadok_movie, shadok_action, shadok_returncode  = put_shadok_in_toolbar(toolBar )
                    subdic["shadok_data"] = [shadok_movie, shadok_action, shadok_returncode,  run_informations, toolBar ]
                else:
                    self.verticalLayout_50.addWidget( toolBar )
                    shadok_movie, shadok_action, shadok_returncode = put_shadok_in_toolbar(toolBar )
                    subdic["shadok_data_bis"] = [shadok_movie, shadok_action, shadok_returncode ,toolBar]

            
                testAction.connect(testAction, Qt.SIGNAL("triggered()"),  Functor_Show_Yaml(subdic    ))
                testAction.connect(testAction2, Qt.SIGNAL("triggered()"),  Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.Runthis     ,identity=identity  ))
                testAction.connect(testAction3, Qt.SIGNAL("triggered()"),  Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.Stopthis    ,identity=identity  ))
                testAction.connect(testAction4, Qt.SIGNAL("triggered()"),  Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.Viewthis    ,identity=identity  ))
                testAction.connect(testAction5, Qt.SIGNAL("triggered()"),  Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.Viewthiserr ,identity=identity  ))
                testAction.connect(testAction_center, Qt.SIGNAL("triggered()"),  Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.Center ,identity=identity  ))
                # testAction.connect(testAction6, Qt.SIGNAL("triggered()"), Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.Pendings ,identity=identity  ))
                testAction.connect(testAction7, Qt.SIGNAL("triggered()"), Functor_Run_Yaml(subdic  , action =  Functor_Run_Yaml.PlugOut ,identity=identity  ))
                
                
            run_informations = {"ALL_dic": self.dico, "process_directory":None, "sub_functors_list":sub_functors_list}


            self.shadok_dictios.append(subdic)

            sub_functors_list.append(Functor_Run_Yaml(subdic , action =  Functor_Run_Yaml.Runthis     ,identity=identity   ))
            
            ####################################
            subdic["refresher_stuff"] =  [ self.timerS ,  self.dico    ]

            lista = []
            class ToggleVisibility(Qt.QToolButton ): ##Qt.QPushButton):
                def __init__(self, lista, nome, parent):
                    self.lista = lista
                    self.visible = 1
                    # super(ToggleVisibility,self).__init__(nome, parent)
                    super(ToggleVisibility,self).__init__()
                    self.setToolButtonStyle(  1)
                    self.setText(  nome  )
                    
                    menu = Qt.QMenu()
                    testAction_toggle_visibility = Qt.QAction("toggle visibility", self);
                    testAction.connect(testAction_toggle_visibility, Qt.SIGNAL("triggered()"), self.toggle )
                    menu.addAction(testAction_toggle_visibility)
                    self.setMenu(menu)
                    self.setPopupMode(Qt.QToolButton.InstantPopup)



                    # butt.setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
                def toggle(self):
                    if self.visible:
                        self.visible=0
                        for w in self.lista:
                            w.hide()
                    else:
                        self.visible=1
                        for w in self.lista:
                            w.show()
                        

            # tbutton = ToggleVisibility(lista, "subparameters",self )
            # menu = Qt.QMenu()
            # testAction_toggle_visibility = Qt.QAction("toggle visibility", tbutton);
            # testAction.connect(testAction_toggle_visibility, Qt.SIGNAL("triggered()"), tbutton.toggle )
            # menu.addAction(testAction_toggle_visibility)
            # tbutton.setMenu(menu)
            # tbutton.setPopupMode(Qt.QToolButton.InstantPopup)
            # toolBar = Qt.QToolBar(self)
            # toolBar.addWidget(tbutton)
            # self.verticalLayout_5.addWidget(toolBar)

            map_par2w = collections.OrderedDict ()

            lista = None
            stack=[iter(subdic.items())]
            # for name,par in subdic.items():

            counttbutton=0
            while( len(stack) ):
                try:
                    name, par = next(stack[-1])
                    if isinstance(par, collections.OrderedDict):
                        stack.append( iter(par.items()) )
                    else:
                        if guided and  isinstance(par, Parametro):
                            if par.master:
                                par.enable = True
                                par.automatic_forward = 2
                                # par.guiding = True
                            else:
                                par.enable = False
                    
                except StopIteration:
                    stack = stack[:-1]
                    if(len(stack)):
                        if counttbutton : 
                            tbutton.toggle()
                        counttbutton+=1
                    continue
                    
                if isinstance(par, collections.OrderedDict):
                    lista = []
                    tbutton = ToggleVisibility(lista, name ,self )
                    toolBar = Qt.QToolBar(self)
                    toolBar.addWidget(tbutton)
                    self.verticalLayout_5.addWidget(toolBar)
                else:
                    if not isinstance(par, Parametro):
                        continue

                    
                    if DEBUG :print( " >>>>>>>>>> " , par.tipo, edition_widgets[par.tipo], par.doc)
                 
                    nw = edition_widgets[par.tipo](self,  name, par , subdic ,map_par2w= map_par2w )

                    if self.guided:
                        if not nw.par.always_visible:
                            self.w2toggle.append( nw)
                            nw.hide_w()

                    
                    # lista.append(nw)
                    self.edws.append(nw) 
                    nw.setToolTip( par.doc   )
                    self.verticalLayout_5.addWidget(nw)
                    map_par2w[par ] = nw
                    if lista is not None:
                        lista.append(nw)
        
            for par,nw in map_par2w.items():
                if par.value is None:
                    par.pump_default()
                if par.value is not None:
                    nw.take_par_value(par)


        spacer = Qt.QSpacerItem(40, 20,  Qt.QSizePolicy.Minimum, Qt.QSizePolicy.Expanding);
        self.verticalLayout_5.addSpacerItem(spacer);
                    
        self.dico["refresher_stuff"] =  [ self.timerS ,  self.dico    ]
        self.connect( self.timer, QtCore.SIGNAL("timeout()"),  self.refresh );
        self.connect( self.timerS, QtCore.SIGNAL("timeout()"),  self.refresh_from_forward );
        self.timerS.setSingleShot(True)
        self.timer.start( 2000)
        self.refresh()
 

class   MyTree( QtGui.QTreeView):
    def __init__(self, *args):
        QtGui.QTreeView.__init__(self, *args)
        self.clicked[QtCore.QModelIndex].connect(self.itemClicked)
        self.doubleClicked[QtCore.QModelIndex].connect(self.itemDoubleClicked)

    def set_other(self, t):
        self.other=t

    def itemClicked(self, modelIndex):

        m = self.other.selectionModel()
        if m is not None:
            m.clearSelection()
        
        event ="itemClicked"
        if DEBUG :print( event, modelIndex)
        
        index = self.selectedIndexes()[0]
        
        
        if False and not( sys.platform =="win32"):
            index = self.model().mapToSource(index )
        pippo = index.model().data(index)
            
        if DEBUG :print( pippo)
        
        pippo = index.model().filePath(index)
        print( pippo)
        
        pippo_name = os.path.basename(pippo)
        if pippo_name[-3:]==".py":
            pippo_name=pippo_name[:-3]
        
        foo = imp.load_source(pippo_name,pippo)
        parent = self.paparent
        for i in reversed(range(parent.scrolledGrid.count())): 
            widgetToRemove = parent.scrolledGrid.itemAt( i ).widget()
            parent.scrolledGrid.removeWidget( widgetToRemove )
            widgetToRemove.setParent( None )

        if hasattr(foo, "OPTIONS"):
            options  = foo.OPTIONS
            print( "  LE OPZIONI SONO " , options)
            i=0
            for chiave, (valore,aiuto) in options.items():
                rb = QtGui.QCheckBox(chiave)
                rb.setChecked(valore)
                
                parent.scrolledGrid.addWidget( rb, i+1, 0,)
                rb.setToolTip( aiuto  )
                i=i+2
            
            # parent.scrolledGrid.show()=
    def itemDoubleClicked(self, modelIndex):

            event ="itemDoubleClicked"
            if DEBUG :print( event, modelIndex)
            
            index = self.selectedIndexes()[0]


            if False and not( sys.platform =="win32"):
                index = self.model().mapToSource(index )
            pippo = index.model().data(index)
            
            if DEBUG :print( pippo)
            
            pippo = index.model().filePath(index)
            if DEBUG : print( str(pippo))
            pippo=str(pippo)
            
            if not( sys.platform =="win32"):
                os.system("emacs %s &"%pippo)
            else:
                comando = "emacs %s  "%pippo
                child_process = subprocess.Popen(comando.split())
                #os.system("emacs %s &"%pippo)
               
           

class MethodsProxyModel(Qt.QAbstractProxyModel):
    def __init__(self, parent=None, model=None):
        Qt.QAbstractProxyModel.__init__(self, model)
        self.model=model
    
    def columnCount( self, index ):
        return 4
    
    def rowCount( self, index ):
        sindex = self.mapToSource(index)
        res = self.model.rowCount(sindex)
        return res
        return 1
    
    def setRoots(self, rootIndexes):
        self.rootIndexes = rootIndexes
    
    def mapToSource(self, proxyIndex):
        if DEBUG :print( " in map to source ")
        index = proxyIndex
        coords = []
        while index != QtCore.QModelIndex():
            coords.append( [ index.row(), index.column()  ]  )
            index = index.parent()
            if DEBUG :print( index)
        coords.reverse()
        if len(coords)==0:
            if DEBUG :print( " RITORNO 0")
            return QtCore.QModelIndex()
        r,c = coords[0]
        if c>0 or r>= len(self.rootIndexes):
            if DEBUG :print( " RITORNO 0b")
            return QtCore.QModelIndex()
        index = self.rootIndexes[c]
        for r,c in coords[1:]:
            index = index.child(r,c)
        if DEBUG :print( " RITORNO ", index)
        return index
        
    def mapFromSource(self, sourceIndex):
        if DEBUG :print( " in map to proxy ")
        index = sourceIndex
        coords = []
        while index not in [ None, QtCore.QModelIndex()] and index not in self.rootIndexes  :
            coords.append( [ index.row(), index.column()  ]  )
            index = index.parent()
        if index in [None, QtCore.QModelIndex()]:
            return QtCore.QModelIndex()

        r,c =  self.rootIndexes.index(index)    , 0
        coords.append([r,c])
        coords.reverse()

        index = QtCore.QModelIndex()
        for r,c in coords:
            index = index.child(r,c)
        return index
        
    def index( self,row,  column,  parent = QtCore.QModelIndex()):
        if (parent.isValid()):
            if parent==QtCore.QModelIndex():
                r,c = row,  column
                if c>0 or r>= len(self.rootIndexes):
                    return QtCore.QModelIndex()
                index = self.rootIndexes[c]
                if DEBUG :print( " RITORNO UN O DEI PRIMI")
                return index
            if DEBUG :print( " DEVO MAPPAR " )
            sourceParent = self.mapToSource(parent)
            return sourceParent.child(row,column)
        else:
            return QtCore.QModelIndex()

  
class G_MethodsProxyModel(Qt.QSortFilterProxyModel):
    def __init__(self, parent=None, model=None, rootIndexes=[]):
        Qt.QSortFilterProxyModel.__init__(self, parent)
        self.model=model
        self.setSourceModel(model )
        self.rootIndexes = rootIndexes
    
    # def columnCount( self, index ):
    #     return 4
    
    # def rowCount( self, index ):
    #     sindex = self.mapToSource(index)
    #     res = self.model.rowCount(sindex)
    #     return res
    #     return 1
    
    def setRoots(self, rootIndexes):
        self.rootIndexes = rootIndexes
        path0 = rootIndexes[0]
        #proxyModel.setFilterRegExp(QRegExp(".png", Qt::CaseInsensitive,
        #                                    QRegExp::FixedString));

    def filterAcceptsRow(self , sourceRow, sourceParent):
        index0 = self.sourceModel().index(sourceRow, 0, sourceParent)
        filepath = index0.model().filePath(index0)
        accepted = False
        for root in self.rootIndexes:
            root =  root.model().filePath(root)
             

            if  False  and sys.platform != "win32":
                if     root[:len(filepath)] ==  filepath :
                    if ((len(root)> len(filepath)) and (root[len(filepath)] not in ["/","\\"])) and  filepath[-1] not in ["/","\\"] :
                        pass
                    else:
                        accepted = True
                        break

                if     filepath[:len(root)] ==   root : 

                    if ((len(filepath)> len(root)) and (filepath[len(root)] not in ["/","\\"]))  and  root[-1] not in ["/","\\"]:
                        pass
                    else:
                        accepted = True
                        break
            else:
                if     root[:len(filepath)].lower() ==  filepath.lower() :
                    if ((len(root)> len(filepath)) and (root[len(filepath)] not in ["/","\\"])) and  filepath[-1] not in ["/","\\"] :
                        pass
                    else:
                        accepted = True
                        break

                if     filepath[:len(root)].lower() ==   root.lower() : 

                    if ((len(filepath)> len(root)) and (filepath[len(root)] not in ["/","\\"]))  and  root[-1] not in ["/","\\"]:
                        pass
                    else:
                        accepted = True
                        break
        return accepted
        
    # def mapToSource(self, proxyIndex):
    #     print " in map to source "
    #     index = proxyIndex
    #     coords = []
    #     while index != QtCore.QModelIndex():
    #         coords.append( [ index.row(), index.column()  ]  )
    #         index = index.parent()
    #         print index
    #     coords.reverse()
        
    #     if len(coords)==0:
    #         print " RITORNO 0"
    #         return QtCore.QModelIndex()

        
    #     r,c = coords[0]
    #     if c>0 or r>= len(self.rootIndexes):
    #         print " RITORNO 0b"
    #         return QtCore.QModelIndex()
        
    #     index = self.rootIndexes[r]
    #     for r,c in coords[1:]:
    #         index = index.child(r,c)
    #     print " RITORNO ", index
    #     return index
    
    # def mapFromSource(self, sourceIndex):
    #     print " in map to proxy "
    #     index = sourceIndex
    #     coords = []
    #     while index not in [ None, QtCore.QModelIndex()] and index not in self.rootIndexes  :
    #         coords.append( [ index.row(), index.column()  ]  )
    #         index = index.parent()
    #     if index in [None, QtCore.QModelIndex()]:
    #         return QtCore.QModelIndex()

    #     r,c =  self.rootIndexes.index(index)    , 0
    #     coords.append([r,c])
    #     coords.reverse()

    #     index = self.index(0,0,       QtCore.QModelIndex()          )
    #     for r,c in coords:
    #         index = index.child(r,c)
    #     return index
        
    # def index( self,row,  column,  parent = QtCore.QModelIndex()):
    #     if (parent.isValid()):
    #         if parent==QtCore.QModelIndex():
    #             r,c = row,  column
    #             if c>0 or r>= len(self.rootIndexes):
    #                 return QtCore.QModelIndex()
    #             index = self.rootIndexes[c]
    #             print " RITORNO UN O DEI PRIMI"
    #             return index
    #         print " DEVO MAPPAR " 
    #         sourceParent = self.mapToSource(parent)
    #         return sourceParent.child(row,column)
    #     else:
    #         return QtCore.QModelIndex()

  

def exit_confirmation():
    msgBox = Qt.QMessageBox(None)
    msgBox.setText("The document has been modified.")
    msgBox.setInformativeText("DO you want to risk loosing unsaved work?")
    msgBox.setStandardButtons( Qt.QMessageBox.Cancel | Qt. QMessageBox.Ok)
    msgBox.setDefaultButton( Qt.QMessageBox.Cancel)
    ret = msgBox.exec_()
    return ret
        
edition_widgets =  { "simplefile": filepath_widget,
                     "simplefile hdf5": hdf5filepath_widget,
                     "text": text_widget,
                     "multiplechoices": multichoices_widget
                 }

class mainwindow(Qt.QMainWindow):
    def __getstate__(self):
        return {}
   
    def Exit(self):
        ret = exit_confirmation()
        if ret ==  Qt.QMessageBox.Discard:
            self.close()
    def LaunchConsole(self):
        from .console import  ExampleWidget
        namespace = {"Qt":Qt,"QtCore":QtCore, "QtGui":QtGui, "np":np, "pylab":pylab }
        self.namespace = namespace
        # self.console = LaunchAnotherConsole(namespace)
        self.consoleW = ExampleWidget()
        self.console = self.consoleW.console
        self.console.updateNamespace(self.namespace)
        self.consoleW.show()
    def __init__(self, InterfaceDescription ,  extrawpaths=[], options={},   parent=None):

        self.InterfaceDescription = copy.deepcopy(InterfaceDescription)
        self.console = None
        self.options=options
        Qt.QMainWindow.__init__(self, parent )
        Qt.loadUi(  os.path.join(installation_dir  ,"resources" , my_relativ_path, "mainwindow.ui" ), self)
       
        self.actionExit.triggered.connect(self.Exit)
        self.actionConsole_1.triggered.connect(self.LaunchConsole)
        self.setWindowIcon(Qt.QIcon(Qt.QPixmap(bacchetta32_xpm)))
        ## self.AAA.setIcon( Qt.QIcon(Qt.QPixmap(bacchetta32_xpm))    )
  
        
        # self.viewsTab.clear()
#        self.methodSelection = creationtab(self.viewsTab)

#         self.viewsTab.addTab(self.methodSelection , "Method Selection")

#         self.metodi_keys = list(self.InterfaceDescription.keys())
#         self.method_selection_list = Qt.QComboBox(self.viewsTab)
#         self.method_selection_list.addItems(self.metodi_keys )
#         self.methodSelection.verticalLayout_3.addWidget( self.method_selection_list )
   
#         self.methodSelection.Creation.setText("Create Selected Method Wizard")
#         self.methodSelection.Creation .connect(self.methodSelection.Creation , QtCore.SIGNAL("clicked()"), self.displayMethodCreator)
        
        self.create_new_calculation .connect(self.create_new_calculation , QtCore.SIGNAL("clicked()"), self.displayMethodCreatorFromTree)
        
        self.methods_model = QtGui.QFileSystemModel()
        self.methods_model.setNameFilters(["winfo_*.py"])
        self.methods_model.setNameFilterDisables(False)

        self.tree_view = MyTree(self )
        self.tree_view.paparent = self
        # @@@@@@@@@@
        self.vsplitter = QtGui.QSplitter()
        self.vsplitter.setOrientation(QtCore.Qt.Vertical) 
        self.verticalLayout_5_b_3.addWidget(self.vsplitter)

        self.optionScroll = QtGui.QScrollArea()
        self.vsplitter.addWidget( self.optionScroll)

        scrolledWidget = QtGui.QWidget()
        self.scrolledGrid   = QtGui.QGridLayout()
 

        self.scrolledGrid.addWidget( QtGui.QLabel("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"), 0, 0);
        self.scrolledGrid.addWidget( QtGui.QLabel("B"), 1, 0);
        self.scrolledGrid.addWidget( QtGui.QLabel("C"), 2, 0);
        self.scrolledGrid.addWidget( QtGui.QLabel("D"), 3, 0);
        self.scrolledGrid.addWidget( QtGui.QLabel("E"), 4, 0); 

        scrolledWidget.setLayout(self.scrolledGrid)
        self.optionScroll.setWidget(scrolledWidget)

        for i in reversed(range(self.scrolledGrid.count())): 
            widgetToRemove = self.scrolledGrid.itemAt( i ).widget()
            self.scrolledGrid.removeWidget( widgetToRemove )
            widgetToRemove.setParent( None )

        
        # self.verticalLayout_5_b_3.addWidget( self.tree_view)
        self.vsplitter.addWidget( self.tree_view)


        self.methods_model_bis = QtGui.QFileSystemModel()
        self.methods_model_bis.setNameFilters(["winfo_*.py"])
        self.methods_model_bis.setNameFilterDisables(False)

        
        self.tree_view_bis = MyTree(self )
        self.tree_view_bis.paparent = self
        self.tree_view_bis.set_other(self.tree_view  ) 
        self.tree_view.set_other(self.tree_view_bis  ) 


        # self.verticalLayout_5_b_3.addWidget( self.tree_view_bis)
        self.vsplitter.addWidget( self.tree_view_bis)
 
        folder = os.path.join(os.path.dirname(__file__),"methods")

        # self.methods_model2 = QtGui.QFileSystemModel()
        # self.methods_model.appendRow(  self.methods_model2     )

        if  False and (not (sys.platform=="win32")) :
  
            class MyLoop(QtCore.QEventLoop):
                def __init__(self):
                    QtCore.QEventLoop.__init__(self)
                    
                def quit(self,path):
                    if DEBUG :print( " loaded " , path)
                    QtCore.QEventLoop.quit(self)

            My_loop = MyLoop()
            self.methods_model.connect(self.methods_model,QtCore.SIGNAL("directoryLoaded(QString)"),  My_loop.quit )
            #My_loop.exec_()


            # self.methods_proxy.setSourceModel(self.methods_model)
            rindexes=[]
            wpaths =extrawpaths+ [folder]
            for p in wpaths:
                rindexes.append(self.methods_model.index(p))

            self.methods_model.setRootPath("/")
            self.methods_proxy = G_MethodsProxyModel( model=self.methods_model, rootIndexes = rindexes )
 
            # self.methods_model.setRootPath("\\")
            self.tree_view.setModel(self.methods_proxy)
            self.tree_view.expandAll()
        else:

            
            self.methods_model.setRootPath(folder)
            self.tree_view.setModel(self.methods_model)
            self.tree_view.setRootIndex(self.methods_model.index(folder))



            if len(extrawpaths):
                folder_bis = extrawpaths[0]
                self.methods_model_bis.setRootPath(folder_bis)
                self.tree_view_bis.setModel(self.methods_model_bis)
                self.tree_view_bis.setRootIndex(self.methods_model_bis.index(folder_bis))


            
           
        self.tree_view.header().setStretchLastSection(False);
        self.tree_view.header().setResizeMode(0, Qt.QHeaderView.Stretch);  
        self.tree_view_bis.header().setStretchLastSection(False);
        self.tree_view_bis.header().setResizeMode(0, Qt.QHeaderView.Stretch);  

        # view.setRootIndex(model.index("C:\\Users")) 
        # self.setCentralWidget(view)
    
    def displayMethodCreator(self, ):
        metodo = str(self.method_selection_list.currentText())
        interfaceDesc = copy.deepcopy( self.InterfaceDescription[metodo]  )        
        view = mainwindow2( interfaceDesc, self.viewsTab  )
        self.viewsTab.addTab(view, "Go create : %s"% metodo)


        
    def displayMethodCreatorFromTree(self, count=[0] ):
        count[0]+=1
        if len(  self.tree_view.selectedIndexes()  ):
            tree = self.tree_view
        elif len(  self.tree_view_bis.selectedIndexes()  ):
            tree = self.tree_view_bis
        else:
            return
            
        index = tree.selectedIndexes()[0]
        if False and not (sys.platform == "win32"):
            index = tree.model().mapToSource(index )

        pippo = index.model().data(index)

        # print pippo
    
        pippo = index.model().filePath(index)
        pippo = str(pippo)

        pippo_name = os.path.basename(pippo)
        if pippo_name[-3:]==".py":
            pippo_name=pippo_name[:-3]


        # modulename =     'pippo'+str(count[0])
        modulename =      pippo_name 
        
        foo = imp.load_source(modulename,pippo)

        options={}
        for i in reversed(range(self.scrolledGrid.count())): 
            widget = self.scrolledGrid.itemAt( i ).widget()
            options[str(widget.text())]= widget.isChecked()

        interfaceDesc  = copy.deepcopy( foo.getMethod(options) ) 
        # interfaceDesc ["tobeimported"] = [   'pippo'+str(count[0]), pippo   ]
        
        interfaceDesc ["tobeimported"] = [  modulename , pippo   ]
        
        view = mainwindow2( interfaceDesc, self.viewsTab  )
        self.viewsTab.addTab(view, "Go create : %s"% os.path.basename(  pippo  ))

        


        

class mainwindow2(Qt.QMainWindow):
    def __getstate__(self):
        return {}
    
    def Exit(self):
        ret = exit_confirmation()
        if ret ==  Qt.QMessageBox.Ok:
            self.paparent.removeTab( self.paparent.currentIndex() )
            self.close()
    def closetab(self,i):
        if i:
            self.viewsTab.widget(i).timer.stop()
            self.viewsTab.removeTab(i)

            
    def __init__(self, InterfaceDescription ,     parent=None):
        
        self.timerS = Qt.QTimer()
 
        self.paparent = parent

        self.InterfaceDescription = copy.deepcopy(InterfaceDescription)


        
        Qt.QMainWindow.__init__(self, parent )
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path, "mainwindow2.ui" ), self)


        
        self.actionExit.triggered.connect(self.Exit)
        self.action_load_start_parameters.triggered.connect(self.load_start_parameters)
        self.action_load_parameters.triggered.connect(self.load_parameters)
        self.action_load_parametersGuided.triggered.connect(self.load_parametersGuided)
        
        self.action_save_parameters.triggered.connect(self.save_parameters)
        
        self.viewsTab.clear()

        #### styleTabs="QTabBar::tab { width:120px; height:30px; color:#333399; font-family:Tahoma; font-size:14px; font-weight:bold; border-top:2px solid #E68B2C; margin-top:-1px; }"
        styleTabs="QTabBar::tab {  width:120px;  color:#333399; font-family:Tahoma; font-size:7px; font-weight:bold; border-top:2px solid #E68B2C; margin-top:-1px; }"
        
        self.viewsTab.setStyleSheet(styleTabs);
        view = creationtab(self)
        self.viewsTab.addTab(view, "Go create it")
        self.viewsTab.setTabsClosable(True);
        self.viewsTab.connect( self.viewsTab  , QtCore.SIGNAL("tabCloseRequested(int)"), self.closetab)

        
        dicodic = self.InterfaceDescription["CREATION_PARS"]
        # dicodic["refresher_stuff"] =  [ self.timerS ,  self.dico    ]
        self.edition_widgets={}
        self.edws = []

        map_par2w = collections.OrderedDict ()
        
        for name, par in dicodic.items():
            if not isinstance(par, Parametro):
                continue
            if 1:
                mydico = collections.OrderedDict ( [ ["mydico", dicodic]]  )
                mydico[ "refresher_stuff"] = [self.timerS ,  mydico    ]
                # print par.tipo, edition_widgets[par.tipo], par.doc
                nw = edition_widgets[par.tipo](self,  name, par , mydico  ,map_par2w= map_par2w )
                map
                self.edws.append(nw)
                
                self.edition_widgets[name]= nw
                map_par2w[par ] = nw
                nw.setToolTip( par.doc   )
                view.verticalLayout.addWidget(nw)
                
        for par,nw in map_par2w.items():
            par.pump_default()

            if par.value is not None:
                nw.take_par_value(par)

        view.Creation.setText("Create Inputs Lists accordingly")
        view.Creation .connect(view.Creation , QtCore.SIGNAL("clicked()"), self.CreateInputsLists)
        
        view.CreationGuided.setText("GUIDED MODE : Create Inputs Lists accordingly")
        view.CreationGuided .connect(view.CreationGuided , QtCore.SIGNAL("clicked()"), self.CreateInputsListsGuided)

        view.CreationInit.setText("Init from dict in the buffer (if any) ")
        view.CreationInit.connect(view.CreationInit , QtCore.SIGNAL("clicked()"), self.pick_the_dict)


        
            
        # view = creationtab(self)
        # self.viewsTab.addTab(view, "1")
        
        self.connect( self.timerS, QtCore.SIGNAL("timeout()"),  self.refresh );
        self.timerS.setSingleShot(True)

    def refresh(self):
        if DEBUG :print( " rinfresco self " ,  self)
        for edw in self.edws:
            if type(edw.par.value) != type(""):
                print( " PROBLEMA in refresh mainwindow2 ", edw.par.value)
                edw.par.value="problem"
            else:
                edw.take_par_value(edw.par )

    def pick_the_dict(self):
        global dico_copy_buffer
        if DEBUG :print(         dico_copy_buffer)
        if DEBUG :print( " era il dict ")
        if "init_initialiser" in self.InterfaceDescription:
            self.InterfaceDescription["init_initialiser"](  dico_copy_buffer )
            self.refresh()
        else:
            print( " non predisposto per init_initialiser")
        

    def CreateInputsListsGuided(self):
        self.CreateInputsLists(guided=True)

    def CreateInputsLists(self, guided=False):
        sintesi = str(self.InterfaceDescription["CREATION_PARS"]["define_new"]())
        sintesi=sintesi[-24:]
        dicodic = copy.deepcopy(self.InterfaceDescription)
        view = creationtab2(self)

        self.viewsTab.addTab(view, sintesi)
        view.populate(dicodic , guided)

    def load_start_parameters(self):
        itab = self.viewsTab.currentIndex()
        filename = str(Qt.QFileDialog.getOpenFileName(None, "select",
                 filter = "All Files (*);;Par Files (*.par);;Pick Files (*.pick)"                                 ))

        if isinstance(filename, tuple):
            filename = filename[0]

        filename=str(filename)
        
        if len(filename ):
            try:
                s=open(filename, "rb" ).read()
                dicodic = pickle.loads( s )

            except:
                s=open(filename, "r" ).read()
                dicodic = pickle.loads( bytearray(s, 'utf-8') )
                
            # s = s.replace( "XRStools.", "XRStools"+version+"." )

            for name, par in dicodic.items():
                if not isinstance(par, Parametro):
                    continue
                else:
                    if par.value is not None:
                        
                        self.edition_widgets[name].take_par_value(par )
            # f.close()

    def load_parametersGuided(self):
        self.load_parameters( guided=True)
        
    def load_parameters(self, guided=False):
        itab = self.viewsTab.currentIndex()
        filename = (Qt.QFileDialog.getOpenFileName(None, "select",
                                                      filter = "All Files (*);;Par Files (*.par);;Pick Files (*.pick)" ) )
        if isinstance(filename, tuple):
            filename = filename[0]
        filename=str(filename)
        
        if len(filename ):
            # f=open(filename, "r" )
            # dicodic = pickle.load( f )
            try:
                s=open(filename, "rb" ).read()
                dicodic = pickle.loads( s,)
                
            except:
                s=open(filename, "r" ).read()
                dicodic = pickle.loads( bytearray(s,'utf-8' ))
                
            # s =  s.replace( "XRStools.", "XRStools"+version+"." )
            
            sintesi=filename[-24:]
            view = creationtab2(self)
            self.viewsTab.addTab(view, sintesi)

            aggiorna_dics(dicodic, self.InterfaceDescription)
            view.populate(dicodic, guided)



    def _finditems(self, obj, keys, results, already_done ):
        if obj in already_done:
            return
        already_done.append(obj)
        for key in keys:
            if key in obj:
                results.append(obj[key])

        for k, v in obj.items():
            if isinstance(v,dict):
                self._finditems(v, keys, results, already_done)

            
    def save_parameters(self):
        filename = (Qt.QFileDialog.getSaveFileName(None, "select"))
        if isinstance(filename, tuple):
            filename = filename[0]
        filename=str(filename)
        
        itab = self.viewsTab.currentIndex()
        if len(filename ):
            if itab==0:
                dico2save = self.InterfaceDescription["CREATION_PARS"]
            else:
                results=[]
                dico2save = copy.deepcopy(self.viewsTab.widget(itab).dico)
                self._finditems( dico2save, ["refresher_stuff", "shadok_data",  "shadok_data_bis"], results, [] )

                for t in results:
                    for i in range(len(t)):
                        t[i]=None

                
            f=open(filename, "wb" )
            s = pickle.dumps( dico2save )
            # s =  s.replace( "XRStools"+version+".", "XRStools." )
            f.write(s)
            f.close()

            # f=open(filename, "w" )
            # pickle.dump( dico2save, f )
            # f.close()
               

            

class hdf5dialog(Qt.QDialog):
    def __getstate__(self):
        return {}
    
    def __init__(self, parent=None):
        Qt.QDialog.__init__(self, parent)
        Qt.loadUi(  os.path.join( installation_dir ,"resources" , my_relativ_path, "hdf5dialog.ui" ), self)


        
        self.hdf5File=None
    def closeEvent(self,event):
        print( " CLOSE EVENT " )
        if self.hdf5File is not None:
            self.hdf5File.close()
        event.accept()



def hdf5_filedialog(hint, fileme=1, myaction = None, modal=1):
    sphint= split_hdf5_address(hint)
    groupname="/"
    if sphint is not None and len(sphint):
        hint , groupname = sphint
    if fileme in [1]:
        if len(hint):
            
            filename =  Qt.QFileDialog.getOpenFileName(None,'Open hdf5 file and then choose groupname', hint,filter="hdf5 (*h5)\nall files ( * )"  )
        else:
            filename =  Qt.QFileDialog.getOpenFileName(None,'Open hdf5 file and then choose groupname',filter="hdf5 (*h5)\nall files ( * )"  )
    else:
        if len(hint):
            filename =  Qt.QFileDialog.getSaveFileName(None,'Open hdf5 file and then choose groupname', hint,filter="hdf5 (*h5)\nall files ( * )"  )
        else:
            filename =  Qt.QFileDialog.getSaveFileName(None,'Open hdf5 file and then choose groupname',filter="hdf5 (*h5)\nall files ( * )"  )

    if isinstance(filename, tuple):
        filename = filename[0]

            
    filename=str(filename)
    if len(filename) :
        ret =None
        if os.path.exists(filename):
            import PyMca5.PyMcaGui.HDF5Widget as HDF5Widget
            storage=[None]
            __hdf5Dialog = hdf5dialog()


            def mySlot(ddict):

                name = ddict["name"]
                if ddict["event"] =='itemDoubleClicked':
                    if myaction is not None:
                        vn = str(__hdf5Dialog.varname.text())
                        if len(vn)==0:
                            vn="grabbed"
                        myaction(filename, name, varname = vn)
                storage[0]=name


            if myaction is None:
                __hdf5Dialog.varname.hide()
            __hdf5Dialog.setWindowTitle('Select a Group  containing roi_definitions by a double click')
            __hdf5Dialog.mainLayout =  __hdf5Dialog.verticalLayout_2
            fileModel = HDF5Widget.FileModel()
            fileView = HDF5Widget.HDF5Widget(fileModel)
            hdf5File = fileModel.openFile(filename)
        
            shiftsDataset = None
            fileView.sigHDF5WidgetSignal.connect(mySlot)
            __hdf5Dialog.mainLayout.addWidget(fileView)
            __hdf5Dialog.resize(400, 700)


            indice = fileModel.index(0,0,QtCore.QModelIndex())
            fileView.expand(indice)
            if groupname !="":
                groupname.replace("//","/")
                gn_l = groupname.split("/")
                hf=h5py.File(filename,"r")
                h5 = hf
                if DEBUG :print( " gn_l " ,  gn_l)
                i=None
                for tn in gn_l:
                    if tn=="":
                        continue
                    if DEBUG :print( tn)
                    nl = list(h5.keys())
                    nl.sort()

                    if tn not in  nl:
                        break

                    i = nl.index( tn  )
                    
                    indice = fileModel.index(i,0,indice)
                    fileView.expand(indice)
                    h5=h5[tn]
                if i is not None:
                    fileView.setCurrentIndex(indice)
                hf.close()

            if modal:
                ret = __hdf5Dialog.exec_()
                hdf5File.close()
            else:
                __hdf5Dialog.hdf5file=hdf5File
                __hdf5Dialog.show()

        if ret:
            name =  storage[0]
            return filename, name
        else:
            # return filename , None
            return None , None
    else:
        
        return None,None



bacchetta32_xpm = [
"32 32 212 2",
"  	c None",
". 	c #000000",
"+ 	c #482D02",
"@ 	c #160B00",
"# 	c #FADE55",
"$ 	c #FFF3A7",
"% 	c #807149",
"& 	c #ECEDEE",
"* 	c #BAB5AE",
"= 	c #CDC3B1",
"- 	c #D7C7AC",
"; 	c #FFE720",
"> 	c #FFFFEF",
", 	c #FFFDB4",
"' 	c #FFFFFF",
") 	c #D3CEC5",
"! 	c #825E03",
"~ 	c #DEC018",
"{ 	c #FFFCA6",
"] 	c #FFFFB1",
"^ 	c #655100",
"/ 	c #0B0700",
"( 	c #372505",
"_ 	c #EFDA98",
": 	c #FFFFFA",
"< 	c #FFF87B",
"[ 	c #FFFAB6",
"} 	c #614E21",
"| 	c #634700",
"1 	c #FFF7BB",
"2 	c #FFF9B5",
"3 	c #FFF6A2",
"4 	c #FFF773",
"5 	c #694600",
"6 	c #A77710",
"7 	c #FFEE71",
"8 	c #FFF6B5",
"9 	c #FFF8A9",
"0 	c #FFF7A4",
"a 	c #FFF4AB",
"b 	c #E4B828",
"c 	c #271600",
"d 	c #FFF9AF",
"e 	c #FFF096",
"f 	c #FFF184",
"g 	c #A08B43",
"h 	c #442F00",
"i 	c #E2D5AE",
"j 	c #FFFCE8",
"k 	c #FFF090",
"l 	c #FFF5BD",
"m 	c #FFE755",
"n 	c #220E00",
"o 	c #291801",
"p 	c #FFD918",
"q 	c #FFF20D",
"r 	c #FFF523",
"s 	c #523D03",
"t 	c #77643A",
"u 	c #FFEE96",
"v 	c #BC8F05",
"w 	c #05060C",
"x 	c #090608",
"y 	c #F4D23B",
"z 	c #FFF89D",
"A 	c #E5C465",
"B 	c #1C1100",
"C 	c #8B5E05",
"D 	c #393634",
"E 	c #F7F7F7",
"F 	c #2B2A29",
"G 	c #050305",
"H 	c #DAB232",
"I 	c #FFF166",
"J 	c #FFF79D",
"K 	c #373331",
"L 	c #110F0E",
"M 	c #D5D0C8",
"N 	c #B7B7B7",
"O 	c #010101",
"P 	c #535353",
"Q 	c #F7F2D2",
"R 	c #533F0F",
"S 	c #352501",
"T 	c #7C5A07",
"U 	c #000208",
"V 	c #6B5500",
"W 	c #F5F3DD",
"X 	c #413100",
"Y 	c #281900",
"Z 	c #DCDCDC",
"` 	c #FDFAEE",
" .	c #FFF6A1",
"..	c #FEFDD4",
"+.	c #B5982A",
"@.	c #180800",
"#.	c #FFF0BD",
"$.	c #FFFFFC",
"%.	c #FFF77D",
"&.	c #FFF78B",
"*.	c #FFF874",
"=.	c #C89A0C",
"-.	c #939393",
";.	c #F3F3F3",
">.	c #CFCFCF",
",.	c #766649",
"'.	c #E5C23B",
").	c #FFEF4C",
"!.	c #1A1208",
"~.	c #090100",
"{.	c #E9CF83",
"].	c #FFF8C3",
"^.	c #FFF890",
"/.	c #FFF79E",
"(.	c #FFF3A3",
"_.	c #302100",
":.	c #2E2E2E",
"<.	c #E6E6E6",
"[.	c #FFE268",
"}.	c #FFE96B",
"|.	c #FFEA63",
"1.	c #FFDD32",
"2.	c #E1E1E1",
"3.	c #323232",
"4.	c #FFCC1C",
"5.	c #FED117",
"6.	c #92721C",
"7.	c #E0A713",
"8.	c #0B0300",
"9.	c #878787",
"0.	c #FEFEFE",
"a.	c #8D8D8D",
"b.	c #ACA59A",
"c.	c #EBE07E",
"d.	c #FFFF18",
"e.	c #625311",
"f.	c #2C2B2B",
"g.	c #F2F2F2",
"h.	c #D7D6D6",
"i.	c #C0B69C",
"j.	c #FEF8B7",
"k.	c #FFFDFB",
"l.	c #FFFFC0",
"m.	c #C7C6C6",
"n.	c #EAE9E9",
"o.	c #2D2D2D",
"p.	c #916E15",
"q.	c #FFF156",
"r.	c #FFF38C",
"s.	c #817756",
"t.	c #7B7B7B",
"u.	c #D8D7D7",
"v.	c #767575",
"w.	c #100900",
"x.	c #1A0E00",
"y.	c #755709",
"z.	c #1D1400",
"A.	c #2E2D2D",
"B.	c #D1D0D0",
"C.	c #B6B5B5",
"D.	c #040404",
"E.	c #AFAEAE",
"F.	c #C8C7C7",
"G.	c #292929",
"H.	c #695104",
"I.	c #FEEA11",
"J.	c #B6A425",
"K.	c #090704",
"L.	c #6D6D6D",
"M.	c #B7B6B5",
"N.	c #646464",
"O.	c #030303",
"P.	c #C6A62D",
"Q.	c #FFFBCA",
"R.	c #4A3810",
"S.	c #292829",
"T.	c #AAA8A8",
"U.	c #939292",
"V.	c #080808",
"W.	c #8E6D19",
"X.	c #F7D855",
"Y.	c #DCC262",
"Z.	c #020307",
"`.	c #8C8B8B",
" +	c #9D9C9D",
".+	c #252424",
"++	c #585757",
"@+	c #8E8C8C",
"#+	c #545353",
"$+	c #232323",
"%+	c #7F7D7D",
"&+	c #757373",
"*+	c #181717",
"=+	c #686768",
"-+	c #878586",
";+	c #605E5F",
">+	c #434242",
",+	c #969696",
"'+	c #B5B5B5",
")+	c #161616",
"!+	c #454545",
"~+	c #A6A5A5",
"{+	c #B1B1B1",
"]+	c #5F5F5F",
"^+	c #0E0E0E",
"/+	c #FBFBFB",
"(+	c #9B9B9B",
"_+	c #616161",
":+	c #AFAFAF",
"<+	c #949494",
"[+	c #0F0F0F",
"}+	c #313131",
". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . + @ . . . . . . . . . . . . . . . . ",
". . . . . . . . . . . . . . # $ % & * . . . = - . . . . . . . . ",
". . . . . . . . . . . . . . ; > , ' ) . ! ~ { ] ^ / . . . . . . ",
". . . . . . . . . . . . ( _ : < [ ' } . | 1 2 3 4 5 . . . . . . ",
". . . . . . . . . . . . 6 7 8 9 0 a b . c d e f g . . . . . . . ",
". . . . . . . . . . . . . h i j k l m n o p q r s . . . . . . . ",
". . . . . . . . . . . . . . t u v w x . . y z ' A . . . . . . . ",
". . . . . . . . . . . . . . B C D ' E F G H I J K L M N O . . . ",
". . . . . . . . . . . . . . . . P ' ' Q ' R S T U V ' W X Y . . ",
". . . . . . . . . . . . . . . . Z ' `  ...+.. @.#.$.%.&.*.=.. . ",
". . . . . . . . . . . . . . . -.;.>.,.'.).!.. ~.{.].^./.(._.. . ",
". . . . . . . . . . . . . . :.' <.. . . . . . . . [.}.|.1.. . . ",
". . . . . . . . . . . . . . 2.' 3.. . . . . . . . 4.5.6.7.8.. . ",
". . . . . . . . . . . . . 9.0.a.. . . . . . . . b.c.d.e.. . . . ",
". . . . . . . . . . . . f.g.h.O . . . . . . . . i.j.k.l.~.. . . ",
". . . . . . . . . . . . m.n.o.. . . . . . . . . p.q.r.s.. . . . ",
". . . . . . . . . . . t.u.v.O . . . . . . . . . w.x.y.z.. . . . ",
". . . . . . . . . . A.B.C.D.. . . . . . . . . . . . . . . . . . ",
". . . . . . . . . O E.F.G.. . . . . . . . . . . H.I.J.K.. . . . ",
". . . . . . . . . L.M.N.O.. . . . . . . . . . . P.Q.' R.. . . . ",
". . . . . . . . S.T.U.V.. . . . . . . . . . . . W.X.Y.Z.. . . . ",
". . . . . . . O.`. +.+. . . . . . . . . . . . . . . . . . . . . ",
". . . . . . . ++@+#+D.. . . . . . . . . . . . . . . . . . . . . ",
". . . . . . $+%+&+*+. . . . . . . . . . . . . . . . . . . . . . ",
". . . . . D.=+-+;+. . . . . . . . . . . . . . . . . . . . . . . ",
". . . . . >+,+'+)+. . . . . . . . . . . . . . . . . . . . . . . ",
". . . . !+~+{+]+. . . . . . . . . . . . . . . . . . . . . . . . ",
". . . ^+/+' Z . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . (+' ' _+. . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . :+' <+. . . . . . . . . . . . . . . . . . . . . . . . . . ",
". . . [+}+D.. . . . . . . . . . . . . . . . . . . . . . . . . . "]
















class Console(QtGui.QPlainTextEdit):
    def __init__(self, prompt='$> ', startup_message='', parent=None):
        QtGui.QPlainTextEdit.__init__(self, parent)
        self.prompt = prompt
        self.history = []
        self.namespace = {"QtCore":QtCore ,"QtGui":QtGui }
        self.construct = []

        self.setGeometry(50, 75, 600, 400)
        self.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)
        self.setUndoRedoEnabled(False)
        self.document().setDefaultFont(QtGui.QFont("monospace", 10, QtGui.QFont.Normal))
        self.showMessage(startup_message)

    def updateNamespace(self, namespace):
        self.namespace.update(namespace)

    def showMessage(self, message):
        self.appendPlainText(message)
        self.newPrompt()

    def newPrompt(self):
        if self.construct:
            prompt = '.' * len(self.prompt)
        else:
            prompt = self.prompt
        self.appendPlainText(prompt)
        self.moveCursor(QtGui.QTextCursor.End)

    def getCommand(self):
        doc = self.document()
        curr_line = unicode(doc.findBlockByLineNumber(doc.lineCount() - 1).text())
        curr_line = curr_line.rstrip()
        curr_line = curr_line[len(self.prompt):]
        return curr_line

    def setCommand(self, command):
        if self.getCommand() == command:
            return
        self.moveCursor(QtGui.QTextCursor.End)
        self.moveCursor(QtGui.QTextCursor.StartOfLine, QtGui.QTextCursor.KeepAnchor)
        for i in range(len(self.prompt)):
            self.moveCursor(QtGui.QTextCursor.Right, QtGui.QTextCursor.KeepAnchor)
        self.textCursor().removeSelectedText()
        self.textCursor().insertText(command)
        self.moveCursor(QtGui.QTextCursor.End)

    def getConstruct(self, command):
        if self.construct:
            prev_command = self.construct[-1]
            self.construct.append(command)
            if not prev_command and not command:
                ret_val = '\n'.join(self.construct)
                self.construct = []
                return ret_val
            else:
                return ''
        else:
            if command and command[-1] == (':'):
                self.construct.append(command)
                return ''
            else:
                return command

    def getHistory(self):
        return self.history

    def setHisory(self, history):
        self.history = history

    def addToHistory(self, command):
        if command and (not self.history or self.history[-1] != command):
            self.history.append(command)
        self.history_index = len(self.history)

    def getPrevHistoryEntry(self):
        if self.history:
            self.history_index = max(0, self.history_index - 1)
            return self.history[self.history_index]
        return ''

    def getNextHistoryEntry(self):
        if self.history:
            hist_len = len(self.history)
            self.history_index = min(hist_len, self.history_index + 1)
            if self.history_index < hist_len:
                return self.history[self.history_index]
        return ''

    def getCursorPosition(self):
        return self.textCursor().columnNumber() - len(self.prompt)

    def setCursorPosition(self, position):
        self.moveCursor(QtGui.QTextCursor.StartOfLine)
        for i in range(len(self.prompt) + position):
            self.moveCursor(QtGui.QTextCursor.Right)

    def runCommand(self):
        command = self.getCommand()
        self.addToHistory(command)

        command = self.getConstruct(command)

        if command:
            tmp_stdout = sys.stdout

            class stdoutProxy():
                def __init__(self, write_func):
                    self.write_func = write_func
                    self.skip = False
                def flush(self):
                    return None
                    
                def write(self, text):
                    if not self.skip:
                        stripped_text = text.rstrip('\n')
                        self.write_func(stripped_text)
                        QtCore.QCoreApplication.processEvents()
                    self.skip = not self.skip

            sys.stdout = stdoutProxy(self.appendPlainText)
            try:
                try:
                    result = eval(command, self.namespace, self.namespace)
                    if result != None:
                        self.appendPlainText(repr(result))
                except SyntaxError:
                    exec( command , self.namespace, self.namespace)
            except SystemExit:
                self.close()
            except:
                traceback_lines = traceback.format_exc().split('\n')
                # Remove traceback mentioning this file, and a linebreak
                for i in (3,2,1,-1):
                    traceback_lines.pop(i)
                self.appendPlainText('\n'.join(traceback_lines))
            sys.stdout = tmp_stdout
        self.newPrompt()

    def keyPressEvent(self, event):
        if event.key() in (QtCore.Qt.Key_Enter, QtCore.Qt.Key_Return):
            self.runCommand()
            return
        if event.key() == QtCore.Qt.Key_Home:
            self.setCursorPosition(0)
            return
        if event.key() == QtCore.Qt.Key_PageUp:
            return
        elif event.key() in (QtCore.Qt.Key_Left, QtCore.Qt.Key_Backspace):
            if self.getCursorPosition() == 0:
                return
        elif event.key() == QtCore.Qt.Key_Up:
            self.setCommand(self.getPrevHistoryEntry())
            return
        elif event.key() == QtCore.Qt.Key_Down:
            self.setCommand(self.getNextHistoryEntry())
            return
        elif event.key() == QtCore.Qt.Key_D and event.modifiers() == QtCore.Qt.ControlModifier:
            self.close()
        super(Console, self).keyPressEvent(event)



        
console_welcome_message = '''
   ---------------------------------------------------------------
     Welcome to a primitive Python interpreter.
      grabbed from :
         http://stackoverflow.com/questions/2758159/how-to-embed-a-python-interpreter-in-a-pyqt-widget
   ---------------------------------------------------------------
'''

def LaunchAnotherConsole(namespace = {}):
    message = console_welcome_message+"\n\nINITIAL NAMESPACE\n\n"+str(namespace)
    console = Console(startup_message=message)
    console.updateNamespace( namespace  )
    console.show();
    return console

 
def wizardMain(filename, metodi, extrawpaths, options ):
    global animation_file

    animation_file = filename
    app=Qt.QApplication([])
    w = mainwindow(metodi,extrawpaths, options )
    w.show()
    global MainWindow
    MainWindow = w
    app.exec_()       

# wizardMain("../../data/shadok_eb-003.gif", metodi )







