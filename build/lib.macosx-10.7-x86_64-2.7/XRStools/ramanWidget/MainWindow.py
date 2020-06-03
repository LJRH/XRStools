from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


try:
    from PyQt4 import Qt, QtCore, QtGui
except:
    from PyQt5 import Qt, QtCore, QtGui


from silx.gui import qt as Qt
from silx.gui import qt as QtCore
from silx.gui import qt as QtGui

# from .. import ui

from .subsetTable import subsetTable
from .scansTable import scansTable
from .acquisition import acquisition
from .edgesTable import edgesTable

from .source import source

from .. import roiSelectionWidget
from ..roiNmaSelectionGui import roiNmaSelectionWidget


import numpy as np
from XRStools import xrs_read, xrs_rois, xrs_extraction
import os
import PyMca5.PyMcaIO.specfilewrapper as SpecIO

from silx.gui.plot.PlotWindow import Plot1D  , Plot2D
from  collections import OrderedDict as odict

from scipy import optimize
import sys
     
import traceback
from six import StringIO


import yaml 
Resolver = yaml.resolver.Resolver
import re
from six import u
Resolver.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u(r"""^(?:[-+]?(?:[0-9][0-9_]*)(\.[0-9_]*)?(?:[eE][-+]?[0-9]+)?
                    |\.[0-9_]+(?:[eE][-+][0-9]+)?
                    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\.[0-9_]*
                    |[-+]?\.(?:inf|Inf|INF)
                    |\.(?:nan|NaN|NAN))$"""), re.X),
        list(u'-+0123456789.'))
import yaml
import yaml.resolver

DEBUG=0
DEBUG2=0
import pickle

from XRStools.installation_dir import installation_dir
import os

my_dir = os.path.dirname(os.path.abspath(__file__))
my_relativ_path =  my_dir [len( os.path.commonprefix([ installation_dir , my_dir  ])):]
if my_relativ_path[0]=="/":
    my_relativ_path = my_relativ_path[1:]

class MainWindow(Qt.QMainWindow) :
    def __init__(self, parent=None):
       
        super(  MainWindow, self).__init__(parent)
        print(  installation_dir  )
        print(  "resources"  )
        print(  my_relativ_path  )
        print(os.path.join(  installation_dir,"resources" , my_relativ_path ,  "MainWindow.ui" )    )
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path ,  "MainWindow.ui" ), self)
        
        sys.excepthook = excepthook

        self.tabWidget.clear()
        self.subsettable = subsetTable()
        self.scanstable  = scansTable()
        self.acquisition =  acquisition()
        self.source =  source()
        self.edges =  edgesTable()
        
        rsw = roiSelectionWidget.mainwindow()
        self.rsw = rsw
        
        self.tabWidget.addTab( rsw  ,"Spatial ROIS")

        sprsw = roiNmaSelectionWidget.mainwindow()

        rsw.user_input_signal.connect(sprsw.update_user_input)

        sprsw.user_input_signal.connect(self.source.update_user_input)

        self.sprsw = sprsw
        
        self.tabWidget.addTab( sprsw  ,"Spectral ROIS")
        
        self.tabWidget.addTab( self.source  ,"experiment")

        sca = QtGui.QScrollArea()
        sca.setWidgetResizable(True) 
        sca.setWidget(self.subsettable)
        self.tabWidget.addTab( sca  ,"analyzers subsets")
        
        sca = QtGui.QScrollArea()
        sca.setWidgetResizable(True) 
        sca.setWidget(self.scanstable)
        self.tabWidget.addTab( sca  ,"scans selection")

        self.tabWidget.addTab( self.edges  ,"edges")
        
        self.tabWidget.addTab( self.acquisition  ,"loading")
        self.acquisition.pushButton.clicked.connect(self.acquire)
        self.acquisition.pushButton_saveAnalysis.clicked.connect(self.saveAnalysis)
    
        self.lw = None
        self.lw_ex = None
        self.plots=[]
        self.plots_conf = {}
        self.emitter2plotC={}

        self.actionSave_Configuration.triggered.connect(self.saveConfiguration)
        self.actionLoad_Configuration.triggered.connect(self.loadConfiguration)

    def loadConfiguration(self):
        self.loadConfigurationOption(None)
    def loadConfigurationOption(self, option=None):

        if option is None:
            filename = QtGui.QFileDialog.getOpenFileName(None, "select", )
            if isinstance(filename, tuple):
                filename = filename[0]

            filename=str(filename)
            if len(filename)==0: return
        else:
            filename = option

        d = yaml.load(open(filename,"r"), yaml.Loader)
        
        selected_scans_d = d["selected_scans"]
        selected_scans = []
        if selected_scans_d is not None:
            for name, scan in selected_scans_d.items():
                selected_scans.append( [name]+scan)
        self.setScansSelection( selected_scans )   

        selected_subsets_d = d["selected_subsets"]
        selected_subsets = []
        if selected_subsets_d is not None:
            for name, subset in selected_subsets_d.items():
                selected_subsets.append( [subset[0],name]+subset[1:])
        self.setSubsetsSelection( selected_subsets )   

        selected_acquisition_d = d["selected_acquisition"]
        names = [ "method", "refscan", "include_elastic", "output_prefix"]
        selected_acquisition =  [     selected_acquisition_d[name] for name in names    ]
        self.setLoadingSelection(selected_acquisition)

        selected_experiment_d = d["selected_experiment"]
        names = [ "specfile_name" , "roifile_address" ]
        selected_experiment =  [     selected_experiment_d[name] for name in names    ]
        self.setExperimentSelection(selected_experiment)


        selected_edges_d = d["selected_edges"]
        formula, edges = selected_edges_d["formula"] , selected_edges_d["edges"]
        self.setEdgesSelection( formula ,  edges )


        if "plots" in d:
            self.plots_conf = d["plots"]
            self.consume_plots_definitions()

    def consume_plots_definitions(self):
        if self.plots_conf is  None:
            return

        for plotC in self.plots:
            name = plotC.name
            if name in self.plots_conf:
                defs = self.plots_conf[name]
                names = ["hfcore_shift" , "pea_center", "pea_width","pea_shape","pea_height","lin_back0", "lin_back1", "hf_factor"]
                inputs = plotC.inputs
                for i,n in zip(inputs, names):
                    if defs[n] is not None:
                        i.setText("%e"%defs[n])
                        
                D = defs
                roisDefs = odict()
                for rn in [ "range1","range2", "Output", "Norm"   ]:
                    roisDefs[rn]  = odict([["from", D[rn][0]] ,["to",D[rn][1]],["type", "energy"] ]) 
                plotC.plot.getCurvesRoiDockWidget().setRois(roisDefs)
                del self.plots_conf[name]


        
    def saveConfiguration(self):
        filename = QtGui.QFileDialog.getSaveFileName(None, "select", )
        if isinstance(filename, tuple):
            filename = filename[0]
        
        filename=str(filename)
        if len(filename)==0: return
        selected_scans       =  self.getScansSelection()
        selected_subsets     =  self.getSubsetsSelection()
        selected_acquisition =  self.getLoadingSelection()
        specfile_name , roifile_address =  self.getExperimentSelection()
        selected_edges = self.getEdgesSelection( ) # [[ "H2O",1.0]] ,  {'O':['K'] } ) 

        print( " DEVO SALVARE " )
        print  (  selected_scans  ) 
        print  (  selected_subsets  ) 
        print  (  selected_acquisition  ) 
        print  (  specfile_name , roifile_address  )
        for plotC in self.plots:
            print( " per plotC ", plotC.name)
            print  (  [ tok.text() for tok in plotC.inputs ]  ) 
            print( plotC.plot.getCurvesRoiDockWidget().getRois() )

        # @@@@@@@@@@ààà quando si ricarica se c'e' un plotC di nome corrispondente inizializzare
        # @@@@@@@@@@@@ se no si tiene da parte e si consuma quando si fa un'acquire

        ## recuperare anche output prefix
        # yaml.load( file o stringa,yaml.Loader)

        ## fare per batch senza grafica
        file = open(filename,"w")
        
        file.write("selected_scans :\n" ) 
        for scan in selected_scans:
            name, nums = scan[0], scan[1:]
            file.write( "    %s : %s\n"%( name , str(list(nums)) ) )
            
        file.write("selected_subsets :\n" ) 
        for subset in selected_subsets:
            scal, name, nums = subset[0], subset[1], subset[2:]
            file.write("    %s : %s\n"%( name , str([scal] + list(nums)) ) )

        file.write("selected_acquisition :\n" ) 
        acqui = selected_acquisition
        names = [ "method", "refscan", "include_elastic", "output_prefix"]
        vals = acqui[0], acqui[1], acqui[2], acqui[3]
        for n,v in zip(names, vals ) :
            file.write("    %s : %s\n"%( n,v) )
                
        file.write("selected_experiment :\n" ) 
        names = [ "specfile_name" , "roifile_address" ]
        vals =  specfile_name , roifile_address
        for n,v in zip(names, vals ) :
            file.write("    %s : %s\n"%( n,v) )

        file.write("selected_edges :\n" ) 
        formula, edges = selected_edges
        file.write("    %s : %s\n"%( " formula",  str(formula)) )
        file.write("    %s : %s\n"%( " edges",  str(edges)) )

            
        file.write("plots :\n" ) 
        for plotC in self.plots:
            print( " per plotC ", plotC.name)
            file.write("    %s   :\n"%plotC.name   ) 
            vals = [ tok.text() for tok in plotC.inputs ]
            names = ["hfcore_shift" , "pea_center", "pea_width","pea_shape","pea_height","lin_back0", "lin_back1", "hf_factor"] 
            for n,v in zip(names, vals ) :
                file.write("        %s : %s\n"%( n,v) )

            dizio_rois = plotC.plot.getCurvesRoiDockWidget().getRois()
            for kroi, val in dizio_rois.items():
                file.write("        %s :   [%e,%e]\n"% (kroi, val.getFrom(), val.getTo()  )  )
                
        file.close()
        
    def acquire(self):
        print(" =============== LOADING============ ")
        selected_scans       =  self.getScansSelection()
        selected_subsets     =  self.getSubsetsSelection()
        selected_acquisition =  self.getLoadingSelection()
            
        if not DEBUG:
            specfile_name , roifile_address =  self.getExperimentSelection()
            formulas, edges = self.getEdgesSelection()
            if len(edges)!=1:
                raise Exception(" So far only one edge can be processed")
            element = list(edges.keys())[0]
            if len(edges[element])!=1:
                raise Exception(" So far only one edge can be processed")
            else:
                edge = edges[element][0]
        
            forms=[]
            weights=[]
            for f,ww in formulas:
                forms.append(f)
                weights.append(float(ww) )
        
            self.lw, roinums = integrate( specfile_name,  roifile_address,  selected_scans, selected_subsets, selected_acquisition  ) 

        for p in self.plots:
            self.tabWidget.removeTab( self.tabWidget.indexOf(p))
            p.deleteLater()
            del p

                  
        toreport = {}
        for plotC in self.plots:
            vals = [ tok.text() for tok in plotC.inputs ]
            dizio_rois = plotC.plot.getCurvesRoiDockWidget().getRois()
            toreport[plotC.name]=[vals, dizio_rois ]
        self.plots=[]


        self.lw_ex_s = []
        
        for iplot,subset in enumerate(selected_subsets):

            if not DEBUG:
                scal = subset[0]
                name = subset[1]
                nums = subset[2:]

                nums = [ i for i in range(len(roinums)) if roinums[i] in nums ]

                lw_ex = xrs_extraction.edge_extraction( self.lw,forms,weights,{element:[edge]})
                lw_ex.analyzerAverage(nums, errorweighing=False)
                y = np.maximum(1.0e-10 , lw_ex.avsignals)


                print("   INITIALISATION y, ", y)
                
                plotC = plotContainer(lw_ex.eloss, y, name, lw_ex, element, edge )
                plotC.controller=self
                self.plots.append(plotC)
                self.tabWidget.addTab(plotC, name)

                if name in toreport:
                    vals, dizio = toreport[name]
                    for tok,v in zip(plotC.inputs, vals):
                        tok.setText(v)
                    plotC.plot.getCurvesRoiDockWidget().setRois(dizio)
                
                self.consume_plots_definitions()

                
            else:
                scal = subset[0]
                name = subset[1]
                nums = subset[2:]
                saved = np.load("debug%d.npy"%iplot)
                eloss, y = saved

                plotC = plotContainer(eloss, y, name)
                self.plots.append(plotC)
                self.tabWidget.addTab(plotC, name)
                

                
    def saveAnalysis(self):
        prefix = str(self.acquisition.lineEdit_outputPrefix.text())
        for C in self.plots:
            C.saveAnalysis(prefix)
                
    def setScansSelection(self,    selection):
        self.scanstable.set_selection( selection )
    def getScansSelection(self):
        res,mini,maxi = self.scanstable.get_selection()
        return res
    def setSubsetsSelection(self,    selection):
        self.subsettable.set_selection( selection )
    def getSubsetsSelection(self):
        return self.subsettable.get_selection()
       
        
    def setLoadingSelection(self,    selection):
        if selection is not None:
            self.acquisition.set_selection( selection )
    def getLoadingSelection(self):
        return self.acquisition.get_selection()
    
    def setExperimentSelection(self,    selection):
        self.source.set_selection( selection )

        
    def getExperimentSelection(self):
        return self.source.get_selection()
    
    def setEdgesSelection(self,    formulas, edges):
        self.edges.set_selection( formulas, edges )
    def getEdgesSelection(self):
        return self.edges.get_selection()
       
def integrate( specfile_name,  roifile_address,  selected_scans, selected_subsets, selected_acquisition  ) :
        print(" SETTING PATH TO LW", specfile_name)
        assert( os.path.exists(specfile_name )   )
        if not os.path.isdir(specfile_name):
            specfile_dir = os.path.dirname(specfile_name)

        print(" DEBUG   lw creo ", specfile_dir)
        lw = xrs_read.Hydra(specfile_dir)

        print("LOADING ROIS")
        myroi = xrs_rois.load_rois_fromh5_address(roifile_address)
            
        
        print("DEBUG roi ", roifile_address ) ;
        lw.set_roiObj(myroi)
        
        print(" TRIMMING KEYS")
        print  ( list(myroi.red_rois.keys()  ) ) 

        roinums = [ int(''.join(filter(str.isdigit, str(key) )))           for key in     myroi.red_rois.keys()     ]
        roinums.sort()
        #### subsets =  self.getSubsetsSelection()
        subsets =  selected_subsets
        new_subsets=[]
        scaling = np.zeros(72)
        for ss in subsets:
            
            first = ss[:2]
            last  = ss[2:]
            nl = [ i for i in range(len(roinums)) if roinums[i] in last ]

            scal,name = first
            scaling[nl] = scal
            
            if len(nl):
                nuovo = first +nl
                print(" APPENDING SET ")
                print( nuovo)
                new_subsets.append(  nuovo   )
            else:
                print( " NIENTE PER ")
                print( first)

        method, ref_scan, keep_elastic,   output_prefix  = selected_acquisition
        print(" CALCULATING COMPENSATION", ref_scan, method)
        print("DEBUG compensation factor ", ref_scan, method   ) 
        lw.get_compensation_factor(ref_scan, method=method )
        
        print(" --------------------------------------------------------")
        print(lw.cenom_dict.keys())
        print("=====================================================")

        
        print(" LOADING GROUPS OF SCANS   ")
        for scan in selected_scans:
            scan_ns = scan[1:]
            scan_name = scan[0]
            print("DEBUG  LOADING ", scan_name, scan_ns, method)
            lw.load_scan( scan_ns, method=method, direct=True, scan_type=scan_name) #  scaling = scaling)


        print(" DEBUG get spe ",  method, keep_elastic)
        lw.get_spectrum_new( method=method , include_elastic=keep_elastic)    
                        

        print(" SET detector angles")
        specfile = SpecIO.Specfile( specfile_name )
        
        scan = specfile.select(str(ref_scan))
        
        rvd  = scan.motorpos("RVD" )  
        rvu  = scan.motorpos("RVU" )  
        rvb  = scan.motorpos("RVB" )
        
        rhr  = scan.motorpos("RHR" )  
        rhl  = scan.motorpos("RHL" )  
        rhb  = scan.motorpos("RHB" )  
        
        lw.get_tths(rvd=rvd, rvu=rvu, rvb=rvb, rhr=rhr, rhl=rhl, rhb=rhb, order=[0, 1, 2, 3, 4, 5])

        return lw, roinums



class MyPlot1D(Plot1D):
    def __init__(self, parent=None):
        super(MyPlot1D, self).__init__(parent)  # , backend = "gl")
    #     self.shint = 400
        
    #     self.setSizePolicy( Qt.QSizePolicy.Fixed, Qt.QSizePolicy.Fixed  ) # ou maximum        
    # def sizeHint(self ) :
    #     return Qt.QSize( self.shint, self.shint)
    
    # def setSizeHint(self, val ) :
    #     self.shint = val
    #     self.updateGeometry()
        



# @ui.UILoadable
class plotContainer(QtGui.QWidget) :
    def __init__(self, eloss, y, name, lw_ex, element, edge, parent=None):


        super( plotContainer, self).__init__(parent)
        Qt.loadUi(  os.path.join(  installation_dir,"resources" , my_relativ_path, "plotContainer.ui" ), self)

        
        self.plot =  MyPlot1D()
        self.layout().addWidget(self.plot)
        
        plot =  self.plot
        plot.getFitAction().setVisible(True)
        plot.setVisible(True)
        plot.getFitAction().setVisible(True)
        y = np.maximum(1.0e-3 , y)
        plot.addCurve(x=eloss, y=y, legend = name, replace="True")
        plot.setYAxisLogarithmic(True)
        self.name = name
        self.lw_ex = lw_ex

        E0,E1 = eloss.min(), eloss.max()
        
        roisDefs = odict( [
            ["ICR",  odict([["from",E0-0.1*(E1-E0) ],["to", E1+0.3*(E1-E0)   ],["type","energy"]])],
            ["range1", odict([["from",E0+0.1*(E1-E0) ],["to", E0+0.3*(E1-E0)   ],["type","energy"]])],
            ["range2", odict([["from",E0+0.6*(E1-E0)],["to",E0+0.8*(E1-E0)],["type","energy"]])],
            ["Output", odict([["from",E0+0.1*(E1-E0)],["to",E0+0.9*(E1-E0)],["type","energy"]])],
            ["Norm", odict([["from",E0+0.6*(E1-E0)],["to",E0+0.7*(E1-E0)],["type","energy"]])]
        ]   )
        
        plot.getCurvesRoiDockWidget().setRois(roisDefs)
        plot.getCurvesRoiDockWidget().setVisible(True)
        ## plot.getCurvesRoiDockWidget().roiWidget.showAllMarkers(True)
        # plot.getCurvesRoiDockWidget().roiWidget._isInit = True
        plot.getCurvesRoiDockWidget().sigROISignal.connect(self.on_rois_changed)
                
        self.pushButton_guess.clicked.connect(self.do_guess)
        self.pushButton_fit.clicked.connect(self.do_fit)
        self.pushButton_plot.clicked.connect(self.plotta)
        
        self.eloss = eloss
        self.y = y 
        self.inputs = [self.lineEdit_hfshift,self.lineEdit_A0,self.lineEdit_A1,self.lineEdit_A2,self.lineEdit_A3,self.lineEdit_A4,self.lineEdit_A5,self.lineEdit_A6]
        self.element = element
        self.edge = edge
        
        self.pushButton_save.clicked.connect(self.saveAnalysis_local)
        
        
    def checkInputs(self):
        
        for tok in self.inputs:
            if len( str( tok.text()))==0:
                tok.setText("0")
        for tok in self.inputs:
            try:
                tmp = float( str(tok.text()))
            except:
                tok.setText("ERROR")
                


        
    def do_guess(self):
        roisDef = self.plot.getCurvesRoiDockWidget().getRois()
        range1 =[ roisDef["range1"].getFrom()  ,   roisDef["range1"].getTo()  ]
        range2 =[ roisDef["range2"].getFrom()  ,   roisDef["range2"].getTo()  ]
        
        print( self.plot.getCurvesRoiDockWidget().getRois())
        
        # define fitting ranges
        region1 = np.where(np.logical_and(self.eloss >= range1[0], self.eloss <= range1[1]))

        print(" REGION 1 " , region1)
        
        region2 = np.where(np.logical_and(self.eloss >= range2[0], self.eloss <= range2[1]))
        print(" REGION 2 " , region2)
        # region  = np.append(region1*weights[0],region2*weights[1])
        region  = np.append(region1,region2)

        # find indices for guessing start values from HF J_total
        
        # print(" QUI FITTO ")
        # fitfct = lambda a: np.sum( (self.y[region] - pearson7_zeroback(self.eloss,a)[region] - 
        #                             np.polyval(a[4:6],self.eloss[region]) )**2.0 )

        fact = self.y[region].max()
        # fitfct = functor2minim( self.eloss[region1], self.y[region1]/fact   ) 
        # guess1 = optimize.minimize(fitfct, [1.0,1.0,1.0,1.0  ], method='SLSQP').x
        # guess1 = optimize.minimize(fitfct, [1.0,1.0,1.0,1.0  ], method='SLSQP').x


        fitfct = functorObjectV(  self.y[region1]/fact , self.eloss[region1],  0   ) 
        bndsa = [ -np.inf]+ [ 0  for tmp in range(3)]
        bndsb =  [ np.inf for tmp in range(4)]

        guess = [1.0,1.0,1.0,1.0  ]
        guess_alt = []
        for t in self.inputs[1:5]:
            if len(str(t.text()))==0:
                break
            else:
                guess_alt.append(float(  str(t.text()) ))
        else:
            guess = guess_alt

        guess[3]/=fact
     
        soluzione    = optimize.least_squares(fitfct,guess , method='trf', bounds=[bndsa,bndsb] ) # constraints=cons)
        guess1       = soluzione.x

        guess1[3]*=fact
        print(" IL RISULTATO DEL FIT EST ", guess1)
        guess=list(guess1)+[0.0,0.0,1.0]
        
        self.lineEdit_A0.setText(str("%e"%guess[0]))
        self.lineEdit_A1.setText(str("%e"%guess[1]))
        self.lineEdit_A2.setText(str("%e"%guess[2]))
        self.lineEdit_A3.setText(str("%e"%guess[3]))
        self.lineEdit_A4.setText(str("%e"%guess[4]))
        self.lineEdit_A5.setText(str("%e"%guess[5]))
        self.lineEdit_A6.setText(str("%e"%guess[6]))
        
        # fit = fitfct.funct( guess  ,   self.eloss  )
        fit = fitfct.funct( guess  ,   self.eloss  )
        self.plot.addCurve(x=self.eloss, y=fit, legend = "guess", replace=False)
        self.guess=guess




        
    def do_fit(self):

        roisDef = self.plot.getCurvesRoiDockWidget().getRois()
        self.checkInputs()
        HFcore_shift = float(self.inputs[0].text())
        guess = [ float(tok.text())  for tok in self.inputs ][1:]

        
        result = removeCorePearson(guess,  self.eloss, self.y, roisDef,HFcore_shift, self.lw_ex, self.element, self.edge)
        guess = result["x"]

        self.lineEdit_A0.setText(str("%e"%guess[0]))
        self.lineEdit_A1.setText(str("%e"%guess[1]))
        self.lineEdit_A2.setText(str("%e"%guess[2]))
        self.lineEdit_A3.setText(str("%e"%guess[3]))
        self.lineEdit_A4.setText(str("%e"%guess[4]))
        self.lineEdit_A5.setText(str("%e"%guess[5]))
        self.lineEdit_A6.setText(str("%e"%guess[6]))
        
        self.plotta()

    def saveAnalysis_local(self):
        roisDef = self.plot.getCurvesRoiDockWidget().getRois()
        prefix =  str(self.controller.acquisition.lineEdit_outputPrefix.text())
        self.saveAnalysis_nogui( prefix, roisDef)

    def saveAnalysis(self, prefix):
        roisDef = self.plot.getCurvesRoiDockWidget().getRois()
        self.saveAnalysis_nogui( prefix, roisDef)
        
    def saveAnalysis_nogui(self, prefix, roisDef):
        y_tot=self.y
        eloss = self.eloss
        this_plot_def = roisDef

        fitfct_result = self.plotta(doplot=True)
        sqwav    =         (y_tot - fitfct_result.peapol)
        sqwaverr =         self.lw_ex.averrors
        allfit = fitfct_result.fit

        
        rangeOutput =[ this_plot_def["Output"].getFrom()  ,  this_plot_def["Output"].getTo()    ]
        rangeNorm =[ this_plot_def["Norm"].getFrom()  ,  this_plot_def["Norm"].getTo()    ]
        regionOutput = np.where(np.logical_and(eloss >= rangeOutput[0], eloss <= rangeOutput[1]))
                
        data = np.zeros((len(regionOutput[0]),3),"d")
        data[:,0] = eloss[regionOutput]
        data[:,1] = sqwav[regionOutput]
        data[:,2] = sqwaverr[regionOutput]
        
        regionNorm = np.where(np.logical_and(data[:,0] >= rangeNorm[0], data[:,0] <= rangeNorm[1]))

        if len(regionNorm):
            norm = np.trapz(data[regionNorm,1],data[regionNorm,0])
            print(" NOOOOOOOOOOOOOOOOOOOOOOOOOOOORM ", norm)
            data[:,1] /= norm
            data[:,2] /= norm
        np.savetxt(prefix+"_"+self.name+".txt",data)

        data = np.zeros((len(eloss),4),"d")
        data[:,0] = eloss
        data[:,1] = y_tot
        data[:,2] = allfit 
        data[:,3] = fitfct_result.peapol
        # np.savetxt(prefix+"_"+self.name+"_all.txt",data)

        
        return


        
    def plotta(self, doplot=True):
        
        roisDef = self.plot.getCurvesRoiDockWidget().getRois()
        self.checkInputs()        
        HFcore_shift = float(self.inputs[0].text())
        guess = [ float(tok.text())  for tok in self.inputs ][1:]        
        HF_core = np.interp(self.eloss,self.eloss+HFcore_shift,self.lw_ex.av_C[self.element][self.edge])
        fitfct_result = functorObjectV( self.y, self.eloss,  HF_core   ) 
        fitfct_result(guess)

        if doplot:
            self.plot.addCurve(x=self.eloss, y=fitfct_result.fit, legend = "Fit(Pearson+linear+CoreHF)", replace=False)
            self.plot.addCurve(x=self.eloss, y=fitfct_result.hf_fit, legend = "Core HF", replace=False)
            self.plot.addCurve(x=self.eloss, y=fitfct_result.peapol, legend = "pea+pol", replace=False)
        
        return fitfct_result
        
    def on_rois_changed(self, message):
        print(" ROIS changed ", message, self.sender())




class functor2minim:
    def __init__(self, eloss, y):
        self.eloss=eloss
        self.y = y
        print(" ELOSS ===================================================================================")
        print( eloss)

    def funct(self, a,eloss):
        pear =  pearson7_zeroback(eloss,a)
        poly =  np.polyval(a[4:6],eloss   )
        tot = pear+poly
        return tot
        
    def __call__(self, a):
        tot = self.funct( a,  self.eloss)
        diff = self.y-tot
        res = (diff*diff/self.y).sum()
        print(" RETURN ", res)
        return res


def pearson7_zeroback(x,a):
    """
    returns a pearson function but without y-offset
    a[0] = Peak position
    a[1] = FWHM
    a[2] = Shape, 1 = Lorentzian, infinite = Gaussian
    a[3] = Peak intensity
    """
    x = np.array(x)
    y = a[3] * (1.0+(2.0**(1.0/a[2])-1.0) * (2.0*(x-a[0])/a[1])**2.0)**(-a[2])
    return y


class myObject(object):
    pass

class functorObject:
    def __init__(self, y, eloss, hfcore):
        self.y      =  y
        self.eloss  =  eloss
        self.hfcore =  hfcore
    def __call__(self, x):
        pea = pearson7_zeroback(self.eloss,x[0:4])
        pol = np.polyval(x[4:6],self.eloss)
        hf  = self.hfcore*x[6]

        self.hf_fit = hf  ##
        self.fit = pea+pol+hf ##
        self.peapol=pea+pol
        
        diff = self.y-self.fit
        
        res  = (diff*diff/self.y).sum()

        # np.save( "Oeloss",self.eloss)
        # np.save( "Oy",self.y)
        # np.save( "Ohf",hf)
        # np.save( "Of" ,    +pea+pol+hf         )
        
        print(" RETURN  ",  res)
        return res
        
class functorObjectV:
    
    def __init__(self, y, eloss, hfcore):
        self.y      =  y
        self.eloss  =  eloss
        self.hfcore =  hfcore

    def funct(self, a,eloss):
        pear =  pearson7_zeroback(eloss,a)
        poly =  np.polyval(a[4:6],eloss   )
        tot = pear+poly
        return tot

        
    def __call__(self, x):
        
        pea = pearson7_zeroback(self.eloss,x[0:4])
        if len(x)==7:
            pol = np.polyval(x[4:6],self.eloss)
            hf  = self.hfcore*x[6]
        else:
            pol=0
            hf=0


        self.hf_fit = hf  ##
        self.fit = pea+pol+hf ##
        self.peapol = pea+pol
        
        diff = self.y-self.fit
        res  = diff/np.sqrt(self.y)

        # np.save( "Oeloss",self.eloss)
        # np.save( "Oy",self.y)
        # np.save( "Ohf",hf)
        # np.save( "Of" ,    +pea+pol+hf         )
        
        return res

def removeCorePearson(guess,  eloss, y_tot, roisDef,HFcore_shift, lw_ex, element, edge):
    range1 =[ roisDef["range1"].getFrom()  ,   roisDef["range1"].getTo()  ]
    range2 =[ roisDef["range2"].getFrom()  ,   roisDef["range2"].getTo()  ]
    region1 = np.where(np.logical_and(eloss >= range1[0], eloss <= range1[1]))
    region2 = np.where(np.logical_and(eloss >= range2[0], eloss <= range2[1]))
    region  = np.append(region1,region2)
    HF_core = np.interp(eloss,eloss+HFcore_shift,lw_ex.av_C[element][edge])
    
    print(" qui y_tot all inizio di removecorepearson est " , y_tot)
    
    fact = y_tot[region].max()
    guess[3] /=fact
    guess[4] /=fact
    guess[5] /=fact
    guess[6] /=fact
    y=y_tot/fact
    y_reg1 = y[region1]
    eloss_reg1 = eloss[region1]
    HF_core_reg1 = HF_core[region1]
    y_reg2 = y[region2]
    eloss_reg2 = eloss[region2]
    HF_core_reg2 = HF_core[region2]
    y_reg = y[region]
    eloss_reg = eloss[region]
    HF_core_reg = HF_core[region]
    fitfct = functorObjectV( y_reg, eloss_reg,  HF_core_reg   ) 
    bndsa = [ -np.inf]+ [ 0  for tmp in range(6)]
    bndsa[5]=-np.inf
    bndsb =  [ np.inf for tmp in range(7)]
    soluzione    = optimize.least_squares(fitfct, guess, method='trf', bounds=[bndsa,bndsb] ) # constraints=cons)
    guess = soluzione.x
    print ( "RISULTATO ", soluzione)
    guess[3] *=fact
    guess[4] *=fact
    guess[5] *=fact
    guess[6] *=fact
    fitfct_result = functorObjectV( y_tot, eloss,  HF_core   ) 
    fitfct_result(guess)
    result = {}
    result["x"] = guess
    return result


def batch( filename):
    d = yaml.load(open(filename,"r"), yaml.Loader)

    selected_scans_d = d["selected_scans"]
    selected_scans = []
    for name, scan in selected_scans_d.items():
        selected_scans.append( [name]+scan)

    selected_subsets_d = d["selected_subsets"]
    selected_subsets = []
    for name, subset in selected_subsets_d.items():
        selected_subsets.append( [subset[0],name]+subset[1:])

    selected_acquisition_d = d["selected_acquisition"]

    if selected_acquisition_d["output_prefix"] is None:
        selected_acquisition_d["output_prefix"] = ""
    
    
    names = [ "method", "refscan", "include_elastic", "output_prefix"]
    selected_acquisition =  [     selected_acquisition_d[name] for name in names    ]
        
    selected_experiment_d = d["selected_experiment"]
    names = [ "specfile_name" , "roifile_address" ]
    selected_experiment =  [     selected_experiment_d[name] for name in names    ]

    lw, roinums = integrate( selected_experiment[0], selected_experiment[1]  ,  selected_scans, selected_subsets, selected_acquisition  ) 


    selected_edges_d = d["selected_edges"]
    formula, edges = selected_edges_d["formula"] , selected_edges_d["edges"]
    if len(edges)!=1:
        raise Exception(" So far only one edge can be processed")
    element = list(edges.keys())[0]
    if len(edges[element])!=1:
        raise Exception(" So far only one edge can be processed")
    else:
        edge = edges[element][0]
        
    forms=[]
    weights=[]
    for f,ww in formula:
        forms.append(f)
        weights.append(float(ww) )


    if "plots" in d:
        plots_def = d["plots"]
    else:
        plots_def = {}
        
    for iplot,subset in enumerate(selected_subsets):
        scal = subset[0]
        name = subset[1]
        nums = subset[2:]

        if name not in plots_def:
            continue

        this_plot_def = plots_def [name]
        
        nums = [ i for i in range(len(roinums)) if roinums[i] in nums ]        
        lw_ex = xrs_extraction.edge_extraction( lw,forms,weights,{element:edge})
        lw_ex.analyzerAverage(nums, errorweighing=False)
        
        eloss = lw_ex.eloss

        
        range1 =[ this_plot_def["range1"][0]  ,  this_plot_def["range1"][1]    ]
        range2 =[ this_plot_def["range2"][0]  ,  this_plot_def["range2"][1]    ]        
        region1 = np.where(np.logical_and(eloss >= range1[0], eloss <= range1[1]))
        region2 = np.where(np.logical_and(eloss >= range2[0], eloss <= range2[1]))
        region  = np.append(region1,region2)

        HFcore_shift = this_plot_def["hfcore_shift"]
        HF_core = np.interp(eloss,eloss+HFcore_shift,lw_ex.av_C[element][edge])
        guess = [ this_plot_def[k] for k in ["pea_center","pea_width","pea_shape","pea_height","lin_back0","lin_back1","hf_factor"]    ]
        
        y_tot = np.maximum(1.0e-10 , lw_ex.avsignals)
        fact = y_tot[region].max()
        guess[3] /=fact
        guess[4] /=fact
        guess[5] /=fact
        guess[6] /=fact

        y=y_tot/fact
        
        y_reg1 = y[region1]
        eloss_reg1 = eloss[region1]
        HF_core_reg1 = HF_core[region1]
        
        y_reg2 = y[region2]
        eloss_reg2 = eloss[region2]
        HF_core_reg2 = HF_core[region2]
        
        y_reg = y[region]
        eloss_reg = eloss[region]
        HF_core_reg = HF_core[region]
        
        fitfct = functorObjectV( y_reg, eloss_reg,  HF_core_reg   ) 
        bndsa = [ -np.inf]+ [ 0  for tmp in range(6)]
        bndsa[5]=-np.inf
        bndsb =  [ np.inf for tmp in range(7)]
        soluzione    = optimize.least_squares(fitfct, guess, method='trf', bounds=[bndsa,bndsb] ) # constraints=cons)
        guess = soluzione.x
        print ( "RISULTATO ", soluzione)
        
        guess[3] *=fact
        guess[4] *=fact
        guess[5] *=fact
        guess[6] *=fact
        y=y*fact

        
        fitfct_result = functorObjectV( y, eloss,  HF_core   ) 
        fitfct_result(guess)

        all4plot      =         (y )
        fit4plot      =         (fitfct_result.fit)
        sqwav    =         (y - fitfct_result.peapol)
        sqwaverr =         lw_ex.averrors

        rangeOutput =[ this_plot_def["Output"][0]  ,  this_plot_def["Output"][1]    ]
        rangeNorm =[ this_plot_def["Norm"][0]  ,  this_plot_def["Norm"][1]    ]
        regionOutput = np.where(np.logical_and(eloss >= rangeOutput[0], eloss <= rangeOutput[1]))
        
        
        data = np.zeros((len(regionOutput[0]),3),"d")
        data[:,0] = eloss[regionOutput]
        data[:,1] = sqwav[regionOutput]
        data[:,2] = sqwaverr[regionOutput]
        
        regionNorm = np.where(np.logical_and(data[:,0] >= rangeNorm[0], data[:,0] <= rangeNorm[1]))

        if len(regionNorm):
            norm = np.trapz(data[regionNorm,1],data[regionNorm,0])
            print(" NOOOOOOOOOOOOOOOOOOOOOOOOOOOORM ", norm)
            data[:,1] /= norm
            data[:,2] /= norm
        np.savetxt(selected_acquisition_d["output_prefix"]+"_"+name+".txt",data)

        data = np.zeros((len(eloss),3),"d")
        data[:,0] = eloss
        data[:,1] = all4plot
        data[:,2] = fit4plot
        # np.savetxt(selected_acquisition_d["output_prefix"]+"_"+name+"_all",data)

   
separator = '-' * 80
def excepthook(type, value, tracebackobj):
    tbinfofile = StringIO()
    traceback.print_tb(tracebackobj, None, tbinfofile)
    # traceback.print_exc()
    # traceback.print_stack(  tracebackobj , file=tbinfofile)

    
    tbinfofile.seek(0)
    tbinfo = tbinfofile.read()
    errmsg = '%s: %s' % (str(type), str(value))
    sections = [separator, errmsg, separator, tbinfo]
    msg = '\n'.join(sections)
    msgBox = Qt.QMessageBox(None)
    msgBox.setText("An exception Occurred")
    msgBox.setInformativeText(msg)
    msgBox.setStandardButtons(  Qt. QMessageBox.Ok)
    msgBox.setDefaultButton( Qt.QMessageBox.Ok)
    ret = msgBox.exec_()
    return



def main(manageQApp = True) :
    if len(sys.argv) > 1:
        batch(sys.argv[1])
    else:
        if manageQApp:
            app=Qt.QApplication([])
        w = MainWindow()
        w.show()
        if manageQApp:
            app.exec_()
        else:
            return w


        
if __name__ =="__main__":

    if len(sys.argv) > 1:
        batch(sys.argv[1])

    else:

        lowq = list(range(24))
        lowq.extend(range(36,60))
        medq = list(range(24,36))
        highq = list(range(60,72))

        app=Qt.QApplication([])
        w = MainWindow()
        w.setScansSelection( [
            ["elastic" , 611,615,619,623 ] ,
            [ "ok0", 628],
            [ 'ok1', 612,616,620,624 ],
            [ "ok2", 613,617,621,625],
            [ "ok3", 614,618,622,626]
        ]  )

        w.setSubsetsSelection( [[3.14, "lowq" ]+lowq ,   [ 3.14, 'middleq']+ medq  , [3.14, "highq"]+ highq ]) 
        w.setLoadingSelection( ["sum",611, True, "pippo" ]  )
        w.setExperimentSelection([   "/data/id20/inhouse/data/run5_17/run7_ihr/hydra",   "Aspectral_myroi.h5:/datas/ROI"   ])
        w.setEdgesSelection( [[ "H2O",1.0]] ,  {'O':['K'] } ) 
        w.show()
        app.exec_()
        print(" SELECTION ")
        print( w.getScansSelection())
        print( w.getSubsetsSelection())
        print( w.getLoadingSelection())
        print( w.getExperimentSelection())
        print( w.getEdgesSelection())


        
        
# recordmydesktop  --v_quality 20 --s_quality 10  --fps 10 --overwrite --device plughw:0,0 -o timefinder-v4-screencast.ogv
# ffmpeg -i terza  -ss 0 -t 23 -acodec copy -vcodec copy terza_a.ogv
# ffmpeg -f concat -safe 0 -i mylist.txt -c copy del.ogv
#    alex@debian:~/src/XRStools/doc$ more ~/mylist.txt
#    file 'prima.ogv'
#    file 'seconda.ogv'
#    file 'terza_a.ogv'
#    file 'terza_b.ogv'
# avconf -i a.ogv a.webm
#                                             per prender solo il video
# avconv  -i /olddisk/home/alex/spectralRoiSelector.ogv  -map 0:0 nnmf.webm
# conversione
# ffmpe -i a.ogv -f webm a.webm
# ffmpeg -i tds2el2_1.ogv -f webm -crf 2  tds2el2_1.webm
# ffmpeg -i TUT1.ogv -vcodec copy -af "volume=enable='between(t,33.6,35.2)':volume=0" test.ogv
