from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: xrs_calctools.py

from . import xrs_utilities
import os
import warnings
from copy import deepcopy

import numpy as np
import array as arr

from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from scipy import constants

import sys

if(sys.argv[0][-12:]!="sphinx-build"):  ## otherwise the documentation cannot build : wrong docstrings imported from pylab
    from pylab import *

from scipy import signal
from scipy.ndimage import measurements
import matplotlib.pyplot as plt
from six.moves import range
from six.moves import zip

__metaclass__ = type # new style classes

A2AU_factor = constants.physical_constants['atomic unit of length'][0]*10**10

def gauss1(x,x0,fwhm):
    """
    returns a gaussian with peak value normalized to unity
    a[0] = peak position
    a[1] = Full Width at Half Maximum
    """
    y = np.exp(-np.log(2.0)*((x-x0)/fwhm*2.0)**2.0)
    return y

def gauss(x,x0,fwhm):
    # area-normalized gaussian
    sigma = fwhm/(2*np.sqrt(2*np.log(2)));
    y = np.exp(-(x-x0)**2/2/sigma**2)/sigma/np.sqrt(2*np.pi)
    return y

def gauss_areanorm(x,x0,fwhm):
    """
    area-normalized gaussian
    """
    sigma = fwhm/(2.0*np.sqrt(2.0*np.log(2.0)))
    y = np.exp(-(x-x0)**2.0/2.0/sigma**2)/sigma/np.sqrt(2.0*np.pi)
    return y

def convg(x,y,fwhm):
    """
    Convolution with Gaussian    
    """
    dx = np.min(np.absolute(np.diff(x)))
    x2 = np.arange(np.min(x)-1.5*fwhm, np.max(x)+1.5*fwhm, dx)
    xg = np.arange(-np.floor(2.0*fwhm/dx)*dx, np.floor(2.0*fwhm/dx)*dx, dx)
    yg = gauss(xg,[0,fwhm])
    yg = yg/np.sum(yg)
    y2 = spline2(x,y,x2)
    c  = np.convolve(y2,yg, mode='full')
    n  = int( np.floor(np.max(np.shape(xg))/2))
    c  = c[n:len(c)-n+1] # not sure about the +- 1 here
    f  = interpolate.interp1d(x2,c)
    return f(x)

def spline2(x,y,x2):
    """
    Extrapolates the smaller and larger valuea as a constant
    """
    xmin = np.min(x)
    xmax = np.max(x)
    imin = x == xmin
    imax = x == xmax
    f  = interpolate.interp1d(x,y, bounds_error=False, fill_value=0.0)
    y2 = f(x2)
    i     = np.where(x2<xmin)
    y2[i] = y[imin]
    i     = np.where(x2>xmax)
    y2[i] = y[imax]
    return y2

def readxas(filename):
    """
    function output = readxas(filename)%[e,p,s,px,py,pz] = readxas(filename)

    % READSTF   Load StoBe fort.11 (XAS output) data
    %
    %   [E,P,S,PX,PY,PZ] = READXAS(FILENAME)
    %
    %      E        energy transfer [eV]
    %      P        dipole transition intensity
    %      S        r^2 transition intensity 
    %      PX       dipole transition intensity along x
    %      PY       dipole transition intensity along y
    %      PZ       dipole transition intensity along z
    %
    %                   as line diagrams.
    %
    %                             T Pylkkanen @ 2011-10-17 
    """
    # Open file
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    data = []
    for line in lines[1:]:
        data.append([float(x) for x in line.replace('D', 'e').strip().split()])

    data      = np.array(data)
    data[:,0] = data[:,0]*27.211384565719481 # convert from a.u. to eV
    

    data[:,3] = 2.0/3.0*data[:,3]**2.0*e/27.211384565719481 # osc(x)
    data[:,4] = 2.0/3.0*data[:,4]**2.0*e/27.211384565719481 # osc(y)
    data[:,5] = 2.0/3.0*data[:,5]**2.0*e/27.211384565719481 # osc(z)

    return data

def broaden_diagram(e,s,params=[1.0, 1.0, 537.5, 540.0],npoints=1000):
    """
    function [e2,s2] = broaden_diagram2(e,s,params,npoints)

    % BROADEN_DIAGRAM2   Broaden a StoBe line diagram
    %
    %  [ENE2,SQW2] = BROADEN_DIAGRAM2(ENE,SQW,PARAMS,NPOINTS)
    %
    %   gives the broadened spectrum SQW2(ENE2) of the line-spectrum
    %   SWQ(ENE). Each line is substituted with a Gaussian peak,
    %   the FWHM of which is determined by PARAMS. ENE2 is a linear
    %   scale of length NPOINTS (default 1000).
    %
    %    PARAMS = [f_min f_max emin max]
    %
    %     For ENE <= e_min, FWHM = f_min.
    %     For ENE >= e_max, FWHM = f_min.
    %     FWHM increases linearly from [f_min f_max] between [e_min e_max].
    %
    %                               T Pylkkanen @ 2008-04-18 [17:37]
    """
    f_min = params[0]
    f_max = params[1]
    e_min = params[2]
    e_max = params[3]

    e2   = np.linspace(np.min(e)-10.0,np.max(e)+10.0,npoints);
    s2   = np.zeros_like(e2) 
    fwhm = np.zeros_like(e)

    # FWHM: Constant  -- Linear -- Constant
    A    = (f_max-f_min)/(e_max-e_min)
    B    = f_min - A*e_min
    fwhm = A*e + B
    inds = e <= e_min
    fwhm[inds] = f_min
    inds = e >= e_max    
    fwhm[inds] = f_max

    for i in range(len(s)):
        s2 += s[i]*gauss(e2,e[i],fwhm[i])

    return e2, s2

def broaden_linear(spec,params=[0.8, 8, 537.5, 550],npoints=1000):
    evals = spec[:,0]
    sticks= spec[:,1]
    f_min = params[0]
    f_max = params[1]
    e_min = params[2]
    e_max = params[3]
    e2    = np.linspace(np.min(evals)-10.0,np.max(evals)+10.0,npoints)
    s2    = np.zeros(len(e2))
    fwhm  = np.zeros(len(evals))
    # FWHM: Constant  -- Linear -- Constant
    A    = (f_max-f_min)/(e_max-e_min)
    B    = f_min - A*e_min
    fwhm = A*evals + B
    fwhm[evals <= e_min] = f_min
    fwhm[evals >= e_max] = f_max
    for n in range(len(sticks)):
        s2 = s2 + sticks[n]*gauss1(e2,evals[n],fwhm[n])
    spectrum = np.zeros((len(e2),2))
    spectrum[:,0] = e2
    spectrum[:,1] = s2
    return spectrum

def load_stobe_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
    """
    load a bunch of StoBe calculations, which filenames are made up of the prefix, postfix, and the counter in the between the prefix and postfix runs from 'fromnumber' to 'tonumber' in steps of 'step' (number of digits is 'stepformat')
    """
    numbers = np.linspace(fromnumber,tonumber,(tonumber-fromnumber + step)//step)
    filenames = []
    precision = '%0'+str(stepformat)+'d'
    for number in numbers:
        thenumber = precision % number
        thefilename = prefix+thenumber+postfix
        filenames.append(thefilename)
    specs = []
    for filename in filenames:
        try:
            specs.append(readxas(filename))
        except:
            print( 'found no file: ' + filename)
    return specs

def load_erkale_spec(filename):
    spec = np.loadtxt(filename)
    return spec

def load_erkale_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
    numbers = np.linspace(fromnumber,tonumber,(tonumber-fromnumber + step)//step)
    filenames = []
    precision = '%0'+str(stepformat)+'d'
    for number in numbers:
        thenumber = precision % number
        thefilename = prefix+thenumber+postfix
        filenames.append(thefilename)
    specs = []
    for filename in filenames:
        try:
            specs.append(load_erkale_spec(filename))
        except:
            print( 'found no file: ' + filename)
    return specs

def cut_spec(spec,emin=None,emax=None):
    if not emin:
        emin = spec[0,0]
    if not emax:
        emax = spec[-1,0]
    spec = spec[spec[:,0]>emin]
    spec = spec[spec[:,0]<emax]
    return spec



def readxas(filename):
    """
    function output = readxas(filename)%[e,p,s,px,py,pz] = readxas(filename)

    % READSTF   Load StoBe fort.11 (XAS output) data
    %
    %   [E,P,S,PX,PY,PZ] = READXAS(FILENAME)
    %
    %      E        energy transfer [eV]
    %      P        dipole transition intensity
    %      S        r^2 transition intensity 
    %      PX       dipole transition intensity along x
    %      PY       dipole transition intensity along y
    %      PZ       dipole transition intensity along z
    %
    %                   as line diagrams.
    %
    %                             T Pylkkanen @ 2011-10-17 
    """
    # Open file
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    data = []
    for line in lines[1:]:
        data.append([float(x) for x in line.replace('D', 'e').strip().split()])

    data      = np.array(data)
    data[:,0] = data[:,0]*27.211384565719481 # convert from a.u. to eV
    

    data[:,3] = 2.0/3.0*data[:,3]**2.0*e/27.211384565719481 # osc(x)
    data[:,4] = 2.0/3.0*data[:,4]**2.0*e/27.211384565719481 # osc(y)
    data[:,5] = 2.0/3.0*data[:,5]**2.0*e/27.211384565719481 # osc(z)

    return data

def broaden_diagram(e,s,params=[1.0, 1.0, 537.5, 540.0],npoints=1000):
    """
    function [e2,s2] = broaden_diagram2(e,s,params,npoints)

    % BROADEN_DIAGRAM2   Broaden a StoBe line diagram
    %
    %  [ENE2,SQW2] = BROADEN_DIAGRAM2(ENE,SQW,PARAMS,NPOINTS)
    %
    %   gives the broadened spectrum SQW2(ENE2) of the line-spectrum
    %   SWQ(ENE). Each line is substituted with a Gaussian peak,
    %   the FWHM of which is determined by PARAMS. ENE2 is a linear
    %   scale of length NPOINTS (default 1000).
    %
    %    PARAMS = [f_min f_max emin max]
    %
    %     For ENE <= e_min, FWHM = f_min.
    %     For ENE >= e_max, FWHM = f_min.
    %     FWHM increases linearly from [f_min f_max] between [e_min e_max].
    %
    %                               T Pylkkanen @ 2008-04-18 [17:37]
    """
    f_min = params[0]
    f_max = params[1]
    e_min = params[2]
    e_max = params[3]

    e2   = np.linspace(np.min(e)-10.0,np.max(e)+10.0,npoints);
    s2   = np.zeros_like(e2) 
    fwhm = np.zeros_like(e)

    # FWHM: Constant  -- Linear -- Constant
    A    = (f_max-f_min)/(e_max-e_min)
    B    = f_min - A*e_min
    fwhm = A*e + B
    inds = e <= e_min
    fwhm[inds] = f_min
    inds = e >= e_max    
    fwhm[inds] = f_max

    for i in range(len(s)):
        s2 += s[i]*gauss1(e2,e[i],fwhm[i])

    return e2, s2

def broaden_linear(spec,params=[0.8, 8, 537.5, 550],npoints=1000):
    """
    broadens a spectrum with a Gaussian of width params[0] below 
    params[2] and width params[1] above params[3], width increases 
    linear in between.
    returns two-column numpy array of length npoints with energy and the broadened spectrum
    """
    evals = spec[:,0]
    sticks= spec[:,1]
    f_min = params[0]
    f_max = params[1]
    e_min = params[2]
    e_max = params[3]
    e2    = np.linspace(np.min(evals)-10.0,np.max(evals)+10.0,npoints)
    s2    = np.zeros(len(e2))
    fwhm  = np.zeros(len(evals))
    # FWHM: Constant  -- Linear -- Constant
    A    = (f_max-f_min)/(e_max-e_min)
    B    = f_min - A*e_min
    fwhm = A*evals + B
    fwhm[evals <= e_min] = f_min
    fwhm[evals >= e_max] = f_max
    for n in range(len(sticks)):
        s2 = s2 + sticks[n]*gauss(e2,evals[n],fwhm[n])
    spectrum = np.zeros((len(e2),2))
    spectrum[:,0] = e2
    spectrum[:,1] = s2
    return spectrum

def load_stobe_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
    """
    load a bunch of StoBe calculations, which filenames are made up of the 
    prefix, postfix, and the counter in the between the prefix and postfix 
    runs from 'fromnumber' to 'tonumber' in steps of 'step' (number of digits 
    is 'stepformat')
    """
    numbers = np.linspace(fromnumber,tonumber,(tonumber-fromnumber + step)//step)
    filenames = []
    precision = '%0'+str(stepformat)+'d'
    for number in numbers:
        thenumber = precision % number
        thefilename = prefix+thenumber+postfix
        filenames.append(thefilename)
    specs = []
    for filename in filenames:
        try:
            specs.append(readxas(filename))
        except:
            print( 'found no file: ' + filename)
    return specs

def load_erkale_spec(filename):
    """
    returns an erkale spectrum
    """
    spec = np.loadtxt(filename)
    return spec

def load_erkale_specs(prefix,postfix,fromnumber,tonumber,step,stepformat=2):
    """
    returns a list of erkale spectra
    """
    numbers = np.linspace(fromnumber,tonumber,(tonumber-fromnumber + step)//step)
    filenames = []
    precision = '%0'+str(stepformat)+'d'
    for number in numbers:
        thenumber = precision % number
        thefilename = prefix+thenumber+postfix
        filenames.append(thefilename)
    specs = []
    for filename in filenames:
        try:
            specs.append(load_erkale_spec(filename))
        except:
            print( 'found no file: ' + filename)
    return specs

def cut_spec(spec,emin=None,emax=None):
    """
    deletes lines of matrix with first column smaller than emin and larger than emax 
    """
    if not emin:
        emin = spec[0,0]
    if not emax:
        emax = spec[-1,0]
    spec = spec[spec[:,0]>emin]
    spec = spec[spec[:,0]<emax]
    return spec


class stobe:
    """
    class to analyze StoBe results
    """
    def __init__(self,prefix,postfix,fromnumber,tonumber,step,stepformat=2):
        self.energy   = [] # array of final energy scale for all snapshots of this run
        self.signal   = [] # array of averaged and smoothed results     
        self.rawspecs = load_stobe_specs(prefix,postfix,fromnumber,tonumber,step,stepformat) # list of all raw stick spectra:rawspecs[n][energy,sticks]
        self.broadened= [] # list of broadened stick spectra 

    def cut_rawspecs(self,emin=None,emax=None):
        cutspecs = []
        for spec in self.rawspecs:
            cutspecs.append(cut_spec(spec,emin,emax))
        self.rawspecs = cutspecs

    def broaden_lin(self,params=[0.8, 8, 537.5, 550],npoints=1000):
        for spec in self.rawspecs:
            self.broadened.append(broaden_linear(spec,params,npoints))

    def sum_specs(self):
        self.energy   = self.broadened[0][:,0] # first spectrum defines energy scale
        self.signal   = np.zeros(np.shape(self.energy))
        for spec in self.broadened:
            f = interp1d(spec[:,0], spec[:,1],bounds_error=False, kind='cubic', fill_value=0.0)
            self.signal += f(self.energy)

    def norm_area(self,emin,emax):
        inds = np.where(np.logical_and(self.energy>=emin,self.energy<=emax))
        norm = np.trapz(self.signal[inds],self.energy[inds])
        self.signal = self.signal/norm
            

class erkale:
    """
    class to analyze ERKALE XRS results.
    """
    def __init__(self,prefix,postfix,fromnumber,tonumber,step,stepformat=2):
        self.energy   = [] # array of final energy scale for all snapshots of this run
        self.sqw      = [] # array of averaged and smoothed results        
        self.rawspecs = load_erkale_specs(prefix,postfix,fromnumber,tonumber,step,stepformat) # list of all raw stick spectra
        self.broadened= []
        self.norm     = [] # results of normalization

    def cut_rawspecs(self,emin=None,emax=None):
        cutspecs = []
        for spec in self.rawspecs:
            cutspecs.append(cut_spec(spec,emin,emax))
        self.rawspecs = cutspecs

    def cut_broadspecs(self,emin=None,emax=None):
        cutspecs = []
        for spec in self.broadened:
            cutspecs.append(cut_spec(spec,emin,emax))
        self.broadened = cutspecs

    def broaden_lin(self,params=[0.8, 8, 537.5, 550],npoints=1000):
        for spec in self.rawspecs:
            self.broadened.append(broaden_linear(spec,params,npoints))

    def sum_specs(self):
        self.energy   = self.broadened[0][:,0] # first spectrum defines energy scale
        self.sqw      = np.zeros(np.shape(self.energy))
        for spec in self.broadened:
            f = interp1d(spec[:,0], spec[:,1],bounds_error=False, kind='cubic', fill_value=0.0)
            self.sqw += f(self.energy)

    def norm_area(self,emin=None,emax=None):
        if not emin:
            emin = self.energy[0]
        if not emax:
            emax = self.energy[-1]

        inds = np.where(np.logical_and(self.energy>=emin,self.energy<=emax))[0]
        self.sqw = self.sqw/np.trapz(self.sqw[inds],self.energy[inds])
    
    def norm_max(self):
        pass

    def plot_spec(self):
        plt.plot(self.energy,self.sqw)
        plt.show()

################################### 
# reading function for cowan's code output

class xyzAtom:
    """ **xyzAtom**

    Class to hold information about and manipulate a single atom in xyz-style format.

    Args. :
          * name (str): Atomic symbol.
          * coordinates (np.array): Array of xyz-coordinates.
          * number (int): Integer, e.g. number of atom in a cluster.

    """
    def __init__(self,name,coordinates,number):
        self.name        = name
        self.coordinates = np.array(coordinates)
        self.x_coord     = self.coordinates[0]
        self.y_coord     = self.coordinates[1]
        self.z_coord     = self.coordinates[2]
        self.Z           = xrs_utilities.element(name)
        self.number      = number
        self.spectrum    = np.array([])

    def getNorm(self):
        return np.linalg.norm(self.coordinates)

    def getCoordinates(self):
        return self.coordinates

    def translateSelf(self, vector):
        try:
            self.coordinates += vector
            self.x_coord     += vector[0]
            self.y_coord     += vector[1]
            self.z_coord     += vector[2]
        except ValueError:
            print('Vector must be 3D np.array!') 

    def translateSelf_arb(self, lattice, lattice_inv, vector):
        rel_coords = np.dot(lattice_inv, np.array(self.coordinates))
        rel_coords += vector

        self.coordinates = np.dot(lattice,rel_coords)
        self.x_coord     = self.coordinates[0]
        self.y_coord     = self.coordinates[1]
        self.z_coord     = self.coordinates[2]

    def load_spectrum(self, file_name):
        self.spectrum = np.loadtxt(file_name)

    def normalize_spectrum(self, normrange):
        inds = np.where(np.logical_and(self.spectrum[:,0]>=normrange[0], self.spectrum[:,0]<=normrange[1]))[0]
        norm = np.trapz(self.spectrum[inds,1], self.spectrum[inds,0])
        self.spectrum[:,1] /= norm


class xyzMolecule:
    """ **xyzMolecule**

    Class to hold information about and manipulate an xyz-style molecule.

    Args.:
        * xyzAtoms (list): List of instances of the xyzAtoms class that make up the molecule.

    """

    def __init__(self,xyzAtoms,title=None):
        self.xyzAtoms = xyzAtoms
        self.title = title

    def getCoordinates(self):
        """ **getCoordinates**
        Return coordinates of all atoms in the cluster.
        """
        return [atom.getCoordinates() for atom in self.xyzAtoms]

    def getCoordinates_name(self,name):
        """ **getCoordinates_name**
        Return coordintes of all atoms with 'name'.
        """
        atoms = []
        for atom in self.xyzAtoms:
            if atom.name == name:
                atoms.append(atom.coordinates)
        return atoms

    def get_atoms_by_name(self,name):
        """ **get_atoms_by_name**
        Return a list of all xyzAtoms of a given name 'name'.
        """
        atoms = []
        for atom in self.xyzAtoms:
            if atom.name == name:
                atoms.append(atom)
            else:
                pass
        if len(atoms) == 0:
                print('Found no atoms with given name in molecule.')
                return
        return atoms

    def getGeometricCenter(self):
        """ **getGeometricCenter**
        Return the geometric center of the xyz-molecule.
        """
        for_average = np.zeros((len(self.xyzAtoms),3))
        for ii in range(len(self.xyzAtoms)):
            for_average[ii, :] = self.xyzAtoms[ii].coordinates
        return np.mean(for_average,axis = 0)

    def getGeometricCenter_arb(self, lattice, lattice_inv):
        pass

    def translateAtomsMinimumImage(self, lattice, lattice_inv):
        """ **translateAtomsMinimumImage**

        Brings back all atoms into the original box using 
        periodic boundary conditions and minimal image 
        convention.

        """
        nullatom = xyzAtom('O', np.array([0.0, 0.0, 0.0]), 0)
        for atom in self.xyzAtoms:
            new_vec = getDistVectorPBC_arb(nullatom, atom, lattice, lattice_inv)
            atom.coordinates = new_vec
            atom.x_coord     = atom.coordinates[0]
            atom.y_coord     = atom.coordinates[1]
            atom.z_coord     = atom.coordinates[2]

    def translateSelf(self,vector):
        """ **translateSelf**
        Translate all atoms of the molecule by a vector 'vector'.
        """
        for atom in self.xyzAtoms:
            atom.translateSelf(vector)

    def scatterPlot(self):
        """ **scatterPlot**
        Opens a plot window with a scatter-plot of all coordinates of the molecule.
        """
        from mpl_toolkits.mplot3d import Axes3D
        fig = figure()
        ax  = Axes3D(fig)
        x_vals = [coord[0] for coord in self.getCoordinates()]
        y_vals = [coord[1] for coord in self.getCoordinates()]
        z_vals = [coord[2] for coord in self.getCoordinates()]
        ax.scatter(x_vals, y_vals, z_vals)
        show()

    def appendAtom(self,Atom):
        """ **appendAtom**
        Add an xzyAtom to the molecule.
        """
        if isinstance(Atom,xyzAtom):
            self.xyzAtoms.append(Atom)
        elif isinstance(Atom,list):
            self.xyzAtoms.extend(Atom)

    def popAtom(self,xyzAtom):
        """ **popAtom**
        Delete an xyzAtom from the molecule.
        """
        self.xyzAtoms.remove(xyzAtom)

    def writeXYZfile(self,fname):
        """ **writeXYZfile**
        Creates an xyz-style text file with all coordinates of the molecule.
        """
        if not self.title:
            self.title = 'None'
        writeXYZfile(fname, len(self.xyzAtoms), self.title, self.xyzAtoms)


class xyzBox:

    """ **xyzBox**

    Class to hold information about and manipulate a xyz-periodic cubic box.

    Args.:
              * xyzAtoms (list): List of instances of the xyzAtoms class that make up the molecule.
              * boxLength (float): Box length.

    """

    def __init__( self, xyzAtoms, boxLength=None, title=None ):
            self.xyzMolecules = []
            self.xyzAtoms     = xyzAtoms
            self.n_atoms      = len(self.xyzAtoms)
            self.title        = title
            self.boxLength    = boxLength
            self.lattice      = None
            self.lattice_inv  = None
            self.relAtoms     = None
            self.av_spectrum  = np.array([])

    def setBoxLength( self, boxLength, angstrom=True ):
            """ **setBoxLength**
            Set the box length.
            """
            if angstrom:
                    self.boxLength = boxLength
            else:
                    self.boxLength = boxLength*constants.physical_constants['atomic unit of length'][0]*10**10

    def writeBox(self, filename):
            """ **writeBox**
            Creates an xyz-style text file with all coordinates of the box.
            """
            writeXYZfile(filename, self.n_atoms, self.title, self.xyzAtoms)

    def writeRelBox(self,filename,inclAtomNames=True):
            """ **writeRelBox**
            Writes all relative atom coordinates into a text file (useful as OCEAN input).
            """
            if not self.boxLength:
                    print('Cannot write rel. coordinates without boxLength. Need to set it first.')
                    return
            else:
                    writeRelXYZfile(filename, self.n_atoms, self.boxLength, self.title, self.xyzAtoms, inclAtomNames)

    def multiplyBoxPBC(self,numShells):
            """ **multiplyBoxPBC**
            Applies the periodic boundary conditions and multiplies the box in shells around the original.
            """
            if not self.boxLength:
                    print('Cannot multiply without boxLength. Need to set it first.')
                    return
            all_atoms = getPeriodicTestBox(self.xyzAtoms,self.boxLength,numbershells=numShells)
            self.xyzMolecules = []
            self.xyzAtoms     = all_atoms
            self.n_atoms      = self.n_atoms*(numShells*2.+1.)**3.
            self.boxLength    = self.boxLength*(numShells*2.0+1.0)

    def multiplyBoxPBC_arb( self, numShells ):
            """ **multiplyBoxPBC_arb**
            Applies the periodic boundary conditions and multiplies 
            the box in shells around the original. Works with arbitrary lattices.
            """
            if not np.any(self.lattice) and not np.any(self.lattice_inv):
                    print('Cannot multiply without lattice. Need to set it first.')
                    return

            all_atoms = getPeriodicTestBox_arb( self.xyzAtoms, self.lattice, self.lattice_inv, numbershells=numShells )
            self.xyzMolecules = []
            self.xyzAtoms     = all_atoms
            self.n_atoms      = self.n_atoms*(numShells*2.+1.)**3.
            try:
                    self.boxLength    = self.boxLength*(numShells*2.0+1.0)
                    self.lattice      = self.lattice*(numShells*2.+1.)
                    self.lattice_inv  = np.linalg.inv(self.lattice)
            except:
                    pass

    def translateAtomsMinimumImage(self, lattice, lattice_inv):
            """ **translateAtomsMinimumImage**

            Brings back all atoms into the original box using 
            periodic boundary conditions and minimal image 
            convention.

            """
            nullatom = xyzAtom('O', np.array([0.0, 0.0, 0.0]), 0)
            for atom in self.xyzAtoms:
                    new_vec = getDistVectorPBC_arb(nullatom, atom, lattice, lattice_inv)
                    atom.coordinates = new_vec
                    atom.x_coord     = atom.coordinates[0]
                    atom.y_coord     = atom.coordinates[1]
                    atom.z_coord     = atom.coordinates[2]

    def deleteTip4pCOM(self):
            """ **deleteTip4pCOM**
            Deletes the ficticious atoms used in the TIP4P water model.
            """
            for atom in self.xyzAtoms:
                    if atom.name == 'M':
                            self.xyzAtoms.remove(atom)
                            self.n_atoms -= 1

    def writeClusters(self,cenatom_name, number,cutoff,prefix,postfix='.xyz'):
            """ **writeXYZclusters**
            Write water clusters into files.
            """
            # find central atom
            cen_atom = self.get_atoms_by_name(cenatom_name)[number]
            coor1 = cen_atom.getCoordinates()
            # cut clusters and write files
            atoms = []
            atoms.append(cen_atom)
            for atom2 in self.xyzAtoms:
                    coor2 = atom2.getCoordinates()
                    if np.linalg.norm( coor1 - coor2) > 0.0 and np.linalg.norm( coor1 - coor2) <= cutoff:
                            atoms.append(atom2)
            fname = prefix + '_%03d' % number + postfix
            box = xyzBox(atoms)
            box.writeBox(fname)

    def writeClusters_arb(self,cenatom_name, number,cutoff,prefix, postfix='.xyz'):
            """ **writeXYZclusters**
            Write water clusters into files.
            """
            test_atoms = deepcopy( self.xyzAtoms )
            test_box = xyzBox(test_atoms)
            test_box.lattice = deepcopy(self.lattice)
            test_box.lattice_inv = deepcopy(self.lattice_inv)
            test_box.multiplyBoxPBC_arb( 1 )

            # find central atom
            cen_atom = self.get_atoms_by_name(cenatom_name)[number]
                        
            # cut clusters and write files
            atoms = []
            atoms.append( cen_atom )
            for atom2 in test_box.xyzAtoms:
                dist = np.linalg.norm(cen_atom.coordinates-atom2.coordinates)
                if dist > 0.0 and dist <= cutoff:
                    atoms.append(atom2)

            box2 = xyzBox(atoms)
            fname = prefix + '_%03d' % number + postfix
            box2.writeBox(fname)

    def writeH2Oclusters(self,cutoff,prefix,postfix='.xyz',o_name='O',h_name='H'):
            """ **writeXYZclusters**
            Write water clusters into files.
            """
            if not self.boxLength:
                    print('Cannot multiply without boxLength. Need to set it first.')
                    return
            # find H2O molecules
            self.get_h2o_molecules(o_name,h_name)
            # get a test box
            pbcMols = getPeriodicTestBox_molecules(self.xyzMolecules,self.boxLength,numbershells=1)
            # cut clusters and write files
            for ii, mol in enumerate(self.xyzMolecules):
                    o_atom = mol.get_atoms_by_name(o_name)[0]
                    #for o_atom, ii in zip(o_atoms,range(len(o_atoms))):
                    cluster = []
                    for molecule in pbcMols:
                            coor = molecule.getCoordinates_name(o_name)
                            if np.linalg.norm( o_atom.coordinates - coor) <= cutoff:
                                    cluster.extend(molecule.xyzAtoms)
                    fname = prefix + '_%03d' % ii + postfix
                    box = xyzBox(cluster)
                    box.writeBox(fname)

    def writeMoleculeCluster(self,molAtomList,fname,cutoff=None,numH2Omols=None,o_name='O',h_name='H',mol_center=None):
            """ **writeMoleculeCluster**
            Careful, this works only for a single molecule in water.
            """
            if not self.boxLength:
                    print('Cannot multiply without boxLength. Need to set it first.')
                    return
            # find H2O molecules
            self.get_h2o_molecules(o_name,h_name)
            # get a test box
            pbcMols = getPeriodicTestBox_molecules(self.xyzMolecules,self.boxLength,numbershells=1)
            # find the solute molecule
            cluster = findMolecule(self.xyzAtoms,molAtomList)
            # find center of mass of molecule
            if not mol_center:
                    cenom = cluster.getGeometricCenter()
            else:
                    cenom = cluster.getCoordinates_name(mol_center)
            if cutoff:
                    # use cutoff criterium
                    waters = findAllWaters(cenom,pbcMols,o_name,cutoff)
                    cluster.appendAtom(waters)
            elif numH2Omols:
                    # use number of waters to include
                    dists = getDistsFromMolecule(cenom,pbcMols,o_name=o_name)
                    inds  = np.argsort(dists)
                    for ii in range(numH2Omols):
                            cluster.appendAtom(pbcMols[inds[ii]].xyzAtoms)
            else:
                    print('Something is fishy!')
                    return
            cluster.writeXYZfile(fname)

    def writeFDMNESinput(self,fname,Filout,Range,Radius,Edge,NRIXS,Absorber):
            """ **writeFDMNESinput**
            Creates an input file to be used for q-dependent calculations with FDMNES.
            """
            writeFDMNESinput_file(self.xyzAtoms,fname,Filout,Range,Radius,Edge,NRIXS,Absorber)

    def writeOCEANinput(self,fname,headerfile,exatom,edge,subshell):
            """ **writeOCEANinput**
            Creates an OCEAN input file based on the headerfile.
            """
            if not self.boxLength:
                    print('Need box size for this function to work!')
                    return
            writeOCEANinput(fname,headerfile,self,exatom,edge,subshell)

    def getCoordinates(self):
            """ **getCoordinates**
            Return coordinates of all atoms in the cluster.
            """
            return [atom.getCoordinates() for atom in self.xyzAtoms]

    def get_atoms_by_name(self,name):
            """ **get_atoms_by_name**
            Return a list of all xyzAtoms of a given name 'name'.
            """
            # find oxygen atoms
            atoms = []
            for atom in self.xyzAtoms:
                    if atom.name == name:
                            atoms.append(atom)
                    else:
                            pass
            if len(atoms) == 0:
                    print('Found no atoms with given name in box.')
                    return
            return atoms

    def scatterPlot(self):
            """ **scatterPlot**
            Opens a plot window with a scatter-plot of all coordinates of the box.
            """
            from mpl_toolkits.mplot3d import Axes3D
            fig = figure()
            ax  = Axes3D(fig)
            x_vals = [coord[0] for coord in self.getCoordinates()]
            y_vals = [coord[1] for coord in self.getCoordinates()]
            z_vals = [coord[2] for coord in self.getCoordinates()]
            cla()
            ax.scatter(x_vals, y_vals, z_vals)
            draw()

    def get_OO_neighbors(self,Roocut=3.6):
            """ **get_OO_neighbors**
            Returns list of numbers of nearest oxygen neighbors within readius 'Roocut'.
            """
            o_atoms = self.get_atoms_by_name('O')
            if not self.boxLength:
                    return count_OO_neighbors(o_atoms,Roocut)
            else:
                    return count_OO_neighbors(o_atoms,Roocut,boxLength=self.boxLength)

    def get_OO_neighbors_pbc(self,Roocut=3.6):
            """ **get_OO_neighbors_pbc**
            Returns a list of numbers of nearest oxygen atoms, uses periodic boundary conditions.
            """
            o_atoms = self.get_atoms_by_name('O')
            return count_OO_neighbors_pbc(o_atoms,Roocut,boxLength=self.boxLength)

    def get_h2o_molecules(self,o_name='O',h_name='H'):
            """ **get_h2o_molecules**
            Finds all water molecules inside the box and collects them inside 
            the self.xyzMolecules attribute.
            """
            o_atoms  = self.get_atoms_by_name(o_name)
            h_atoms  = self.get_atoms_by_name(h_name)
            if self.boxLength:
                    self.xyzMolecules = find_H2O_molecules(o_atoms,h_atoms,boxLength=self.boxLength)
            else:
                    self.xyzMolecules = find_H2O_molecules(o_atoms,h_atoms,boxLength=None)

    def get_h2o_molecules_arb(self, o_name='O',h_name='H'):
            o_atoms  = self.get_atoms_by_name(o_name)
            h_atoms  = self.get_atoms_by_name(h_name)
            self.xyzMolecules = h2o_mols = find_H2O_molecules_PBC_arb( o_atoms, h_atoms, self.lattice, self.lattice_inv )


    def get_atoms_from_molecules(self):
            """ **get_atoms_from_molecules**
            Parses all atoms inside self.xyzMolecules into self.xyzAtoms (useful for turning
            an xyzMolecule into an xyzBox).
            """
            if not self.xyzAtoms:
                    self.xyzAtoms = []
                    for molecule in self.xyzMolecules:
                            for atom in molecule.xyzAtoms:
                                    self.xyzAtoms.append(atom)
            self.n_atoms      = len(self.xyzAtoms)

    def get_hbonds(self, Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
            """ **get_hbonds**
            Counts the hydrogen bonds inside the box, returns the number
            of H-bond donors and H-bond acceptors.
            """
            hb_accept_sum  = 0 # total number of H-bonds in box
            hb_donor_sum  = 0 # total number of donor H-bonds in box
            dbonds = []
            abonds = []
            o_atoms  = self.get_atoms_by_name('O')
            h_atoms  = self.get_atoms_by_name('H')
            if self.boxLength:
                    h2o_mols = find_H2O_molecules(o_atoms,h_atoms,boxLength=self.boxLength)
                    self.xyzMolecules = h2o_mols
            else:
                    h2o_mols = find_H2O_molecules(o_atoms,h_atoms)
                    test_h2o_mols = h2o_mols
            for mol1 in h2o_mols:
                    dbonds_1 = 0
                    abonds_1 = 0
                    for mol2 in h2o_mols:
                            donor, acceptor = countHbonds_pbc(mol1,mol2,self.boxLength,Roocut=Roocut, Rohcut=Rohcut, Aoooh=Aoooh)
                            dbonds_1 += donor
                            abonds_1 += acceptor
                            hb_donor_sum  += donor
                            hb_accept_sum += acceptor
                    dbonds.append(dbonds_1)
                    abonds.append(abonds_1)
            return dbonds, abonds #hb_donor_sum, hb_accept_sum #hbonds, dbonds, abonds, hbondspermol

    def changeOHBondlength(self,fraction, oName='O', hName='H'):
            """ **changeOHBondlength**
            Changes all OH covalent bond lengths inside the box by a fraction.
            """
            o_atoms  = self.get_atoms_by_name('O')
            h_atoms  = self.get_atoms_by_name('H')
            # find all H2O molecules
            if self.boxLength:
                    h2o_mols = find_H2O_molecules(o_atoms,h_atoms,boxLength=self.boxLength)
            else:
                    h2o_mols = find_H2O_molecules(o_atoms,h_atoms)
            # change the bond length
            new_h2o_mols = []
            for mol in h2o_mols:
                    new_h2o_mols.append(changeOHBondLength(mol, fraction, boxLength=self.boxLength, oName=oName, hName=hName))
            # redefine all molecules and atoms in box
            self.xyzMolecules = new_h2o_mols
            self.xyzAtoms = []
            for mol in self.xyzMolecules:
                    for atom in mol.xyzAtoms:
                            self.xyzAtoms.append(atom)

    def getTetraParameter(self):
            """ **getTetraParameter**
            Returns a list of tetrahedrality paprameters, according to NATURE, VOL 409, 18 JANUARY (2001).

            UNTESTED!!!
            """
            if not self.boxLength:
                    print('This only works with PBC. Need to set boxLength first.')
                    return
            else:
                    o_atoms  = self.get_atoms_by_name('O')
                    return getTetraParameter(o_atoms,self.boxLength)

    def get_angle_arb( self, atom1, atom2, atom3, degrees=True):
            """ **get_angle**
            Return angle between the three given atoms (as seen from atom2).

            """
            vec1 = getDistVectorPBC_arb(atom1, atom2, self.lattice, self.lattice_inv)
            vec2 = getDistVectorPBC_arb(atom3, atom2, self.lattice, self.lattice_inv)
            dotp = np.dot(vec1/np.linalg.norm(vec1), vec2/np.linalg.norm(vec2))
            if degrees:
                    return np.degrees( np.arccos( np.clip( dotp, -1.0, 1.0 ) ) )
            else:
                    return np.arccos( np.clip( dotp, -1.0, 1.0 ) )

    def get_angle( self, atom1, atom2, atom3, degrees=True):
            """ **get_angle**
            Return angle between the three given atoms (as seen from atom2).

            """
            vec1 = getDistVector(atom1, atom2)
            vec2 = getDistVector(atom3, atom2)
            dotp = np.dot(vec1/np.linalg.norm(vec1), vec2/np.linalg.norm(vec2))
            if degrees:
                    return np.degrees( np.arccos( np.clip( dotp, -1.0, 1.0 ) ) )
            else:
                    return np.arccos( np.clip( dotp, -1.0, 1.0 ) )
            
    def count_neighbors( self, name1, name2, cutoff_low=0.0, cutoff_high=2.0, counter_name='num_OO_shell' ):
            """ **count_neighbors**

            Counts number of neighbors (of name2) around atom of name1.

            Args:
                      * name1         (str): Name of first type of atom.
                      * name2         (str): Name of second type of atom.
                      * cutoff_low  (float): Lower cutoff (Angstrom).
                      * cutoff_high (float): Upper cutoff (Angstrom).
                      * counter_name  (str): Attribute namer under which  the result should be saved.

            """
            for atom in self.xyzAtoms:
                    if atom.name == name1:
                            cou = 0
                            atoms_2 = self.get_atoms_by_name(name2)
                            for atom_2 in atoms_2:
                                    dist = getDistancePBC_arb( atom, atom_2, self.lattice, self.lattice_inv )
                                    if dist >= cutoff_low and dist <= cutoff_high:
                                            cou += 1
                            setattr(atom, counter_name, cou)

    def count_hbonds( self, Roocut=3.6, Rohcut=2.4, Aoooh=30.0, counter_name='num_H_bonds', counter_name2='H_bond_angles'):
            """ **count_hbonds**
            Counts the number of hydrogen bonds around all oxygen atoms and sets
            that number as attribute to the accorting xyzAtom.

            """
            o_atoms  = self.get_atoms_by_name('O')
            h_atoms  = self.get_atoms_by_name('H')
            h2o_mols = find_H2O_molecules_PBC_arb( o_atoms, h_atoms, self.lattice, self.lattice_inv )

            for mol1 in h2o_mols:
                    don = 0
                    acc = 0
                    angles = []
                    for mol2 in h2o_mols:
                            d, a, ang = count_HBonds_pbc_arb( mol1, mol2, self.lattice, self.lattice_inv, Roocut=Roocut, Rohcut=Rohcut, Aoooh=Aoooh )
                            don += d
                            acc += a
                            angles.append(ang)
                    the_o_atom = mol1.get_atoms_by_name('O')[0]
                    angles = np.array(angles)
                    setattr(the_o_atom, counter_name, (don, acc))
                    setattr(the_o_atom, counter_name2, angles)

            #for atom in o_atoms:
            #	acceptor = 0
            #	donor = 0
            #	# first molecule
            #	mol_1 = []
            #	mol_1.append(atom)
            #	for h_atom in h_atoms:
            #		if getDistancePBC_arb(atom, h_atom, self.lattice, self.lattice_inv ) <= 1.2:
            #			mol_1.append(h_atom)
            #	mol1 = xyzMolecule(mol_1)
            #	for o_atom in o_atoms:
            #		# second molecule
            #		mol_2 = []
            #		mol_2.append(o_atom)
            #		for h_atom in h_atoms:
            #			if getDistancePBC_arb( o_atom, h_atom, self.lattice, self.lattice_inv) <= 1.5:
            #				mol_2.append(h_atom)
            #		mol2 = xyzMolecule(mol_2)
            #		d, a = count_HBonds_pbc_arb( mol1, mol2, self.lattice, self.lattice_inv, Roocut=3.6, Rohcut=2.4, Aoooh=30.0 )
            #		acceptor += a
            #		donor += d
            #	setattr(atom, counter_name, (donor, acceptor))

    def count_contact_pairs( self, name_1, name_2, cutoff, counter_name='contact_pair'):
            atoms_1 = self.get_atoms_by_name(name_1)
            atoms_2 = self.get_atoms_by_name(name_2)

            for atom1 in atoms_1:
                    contact_pair = []
                    contact_pair.append(atom1)
                    for atom2 in atoms_2:
                            dist = getDistancePBC_arb(atom1, atom2, self.lattice, self.lattice_inv )
                            if dist <= cutoff:
                                    contact_pair.append(atom2)
                    if len(contact_pair) == 2:
                            setattr(atom1, counter_name, 1)
                    else:
                            setattr(atom1, counter_name, 0)

    def normalize_spectrum(self, normrange):
            inds = np.where(np.logical_and(self.av_spectrum[:,0]>=normrange[0], self.av_spectrum[:,0]<=normrange[1]))[0]
            norm = np.trapz(self.av_spectrum[inds,1], self.av_spectrum[inds,0])
            self.av_spectrum[:,1] /= norm

    def normalize_arb_spectrum(self, normrange, attribute):
            spectrum = getattr(self, attribute)
            inds = np.where(np.logical_and(spectrum[:,0]>=normrange[0], spectrum[:,0]<=normrange[1]))[0]
            norm = np.trapz(spectrum[inds,1], spectrum[inds,0])
            spectrum[:,1] /= norm
            setattr(self, attribute, spectrum)

    def find_hydroniums( self, OH_cutoff=1.5 ):
            """ **find_hydroniums**
            Returns a list of hydronium molecules.

            """
            o_atoms = self.get_atoms_by_name('O')
            h_atoms = self.get_atoms_by_name('H')
            hydroniums = []
            for o_atom in o_atoms:
                    molecule = []
                    molecule.append( o_atom )
                    for h_atom in h_atoms:
                            if np.any(self.lattice) and np.any(self.lattice_inv):
                                    oh_dist = getDistancePBC_arb( o_atom, h_atom, self.lattice, self.lattice_inv )
                            elif np.any(self.boxLength):
                                    oh_dist = getDistancePbc( o_atom, h_atom, self.boxLength )
                            else:
                                    oh_dist = np.linalg.norm( o_atom.coordinates - h_atom.coordinates )
                            if oh_dist <= OH_cutoff:
                                    molecule.append( h_atom )
                    if len(molecule) == 4:
                            hydroniums.append( xyzMolecule(molecule) )
            return hydroniums

    def find_tmao_molecules_arb(self, CH_cut=1.2, CN_cut=1.6, NO_cut=1.5, CC_cut=2.5 ):
        """ **find_tmao_molecules**
        Returns a list of TMAO molecules.

        """
        o_atoms = self.get_atoms_by_name('O')
        h_atoms = self.get_atoms_by_name('H')
        c_atoms = self.get_atoms_by_name('C')
        n_atoms = self.get_atoms_by_name('N')
        tmao_mols = []
        #
        for n_atom in n_atoms:
            molecule = []
            molecule.append( n_atom )
            # find all C atoms
            for c_atom in c_atoms:
                cn_dist = getDistancePBC_arb( n_atom, c_atom, self.lattice, self.lattice_inv )
                if cn_dist <= CN_cut:
                    molecule.append(c_atom)
            # find the O atom
            for o_atom in o_atoms:
                on_dist = getDistancePBC_arb( n_atom, o_atom, self.lattice, self.lattice_inv )
                if on_dist <= NO_cut:
                    molecule.append(o_atom)
            # find the H atoms
            for c_atom in c_atoms:
                for h_atom in h_atoms:
                    ch_dist = getDistancePBC_arb( c_atom, h_atom, self.lattice, self.lattice_inv )
                    if ch_dist <= CN_cut:
                        molecule.append(h_atom)
            # check if molecule is complete
            if len(molecule) == 14:
                    tmao_mols.append(xyzMolecule(molecule))
        #
        return tmao_mols

    def find_urea_molecules_arb(self, NH_cut=1.2, CN_cut=1.6, CO_cut=1.5 ):
        """ **find_urea_molecules**
        Returns a list of Urea molecules.

        """
        o_atoms = self.get_atoms_by_name('O')
        h_atoms = self.get_atoms_by_name('H')
        c_atoms = self.get_atoms_by_name('C')
        n_atoms = self.get_atoms_by_name('N')
        urea_mols = []
        #
        for c_atom in c_atoms:
            molecule = []
            molecule.append( c_atom )
            # find the O atom
            for o_atom in o_atoms:
                oc_dist = getDistancePBC_arb( c_atom, o_atom, self.lattice, self.lattice_inv )
                if oc_dist <= CO_cut:
                    molecule.append(o_atom)
            # find the N atoms
            mol_n_atoms = []
            for n_atom in n_atoms:
                nc_dist = getDistancePBC_arb( c_atom, n_atom, self.lattice, self.lattice_inv )
                if nc_dist <= CN_cut:
                    molecule.append(n_atom)
                    mol_n_atoms.append(n_atom)
            # find the H atoms
            for n_atom in mol_n_atoms:
                for h_atom in h_atoms:
                    nh_dist = getDistancePBC_arb( n_atom, h_atom, self.lattice, self.lattice_inv )
                    if nh_dist <= NH_cut:
                        molecule.append(h_atom)
            # check if molecule is complete
            if len(molecule) == 8:
                    urea_mols.append(xyzMolecule(molecule))
        #
        return urea_mols

    def find_hydroxides( self, OH_cutoff=1.5 ):
            """ **find_hydroxides**
            Returns a list of hydroxide molecules.

            """
            o_atoms = self.get_atoms_by_name('O')
            h_atoms = self.get_atoms_by_name('H')
            hydroxides = []
            for o_atom in o_atoms:
                    molecule = []
                    molecule.append( o_atom )
                    for h_atom in h_atoms:
                            if np.any(self.lattice) and np.any(self.lattice_inv):
                                    oh_dist = getDistancePBC_arb( o_atom, h_atom, self.lattice, self.lattice_inv )
                            elif np.any(self.boxLength):
                                    oh_dist = getDistancePbc( o_atom, h_atom, self.boxLength )
                            else:
                                    oh_dist = np.linalg.norm( o_atom.coordinates - h_atom.coordinates )
                            if oh_dist <= OH_cutoff:
                                    molecule.append( h_atom )
                    if len(molecule) == 2:
                        hydroxides.append( xyzMolecule(molecule) )
            return hydroxides

    def getDistancePBC_arb(self, atom1, atom2):
            """ **getDistancePBC_arb**
            Calculates the distance of two atoms from an arbitrary 
            simulation box using the minimum image convention.

            Args:
                atom1 (obj): Instance of the xzyAtom class.
                atom2 (obj): Instance of the xzyAtom class.

            Returns:
                The distance between the two atoms.

            """
            return np.linalg.norm( getDistVectorPBC_arb(atom1, atom2, self.lattice, self.lattice_inv) )

    def getDistVectorPBC_arb(self, atom1, atom2):
        """ **getDistVectorPBC_arb**

        Calculates the distance vector between two atoms from an 
        arbitrary simulation box using the minimum image convention.

        Args:
            atom1 (obj): Instance of the xzyAtom class.
            atom2 (obj): Instance of the xzyAtom class.

        Returns:
            The distance vector between the two atoms (np.array).

        """
        dist_vec = np.array(atom2.coordinates) - np.array(atom1.coordinates)
        dist_vec -= np.dot(self.lattice, np.round(np.dot(self.lattice_inv,dist_vec)))
        return dist_vec



class xyzTrajectory:
    def __init__(self,xyzBoxes):
        self.xyzBoxes  = xyzBoxes
        try:
            self.boxLength = xyzBoxes[0].boxLength
        except:
            self.boxLength = None

    def writeRandBox(self,filename):
        ind = np.random.randint(len(self.xyzBoxes))
        self.xyzBoxes[ind].writeBox(filename)

    def loadAXSFtraj(self,filename):
        self.xyzBoxes = axsfTrajParser(filename)
        try:
            self.boxLength = self.xyzBoxes[0].boxLength
        except:
            pass

    def writeXYZtraj(self,filename):
        writeXYZtrajectory(filename,self.xyzBoxes)

    def getRDF(self,atom1='O',atom2='O',MAXBIN=1000,DELR=0.01,RHO=1.0):
        HIST = np.zeros(MAXBIN+1)
        GR   = np.zeros(MAXBIN+1)
        for box in self.xyzBoxes:
            atoms = box.get_atoms_by_name(atom1) + box.get_atoms_by_name(atom2)
            HIST += calculateRIJhist(atoms,self.boxLength,DELR=DELR,MAXBIN=MAXBIN)
        CONST = 4.0*np.pi*RHO/3.0
        for BIN in range(2,MAXBIN):
            RLOWER = (BIN-1)*DELR
            RUPPER = RLOWER + DELR
            NIDEAL = CONST * ( RUPPER**3 - RLOWER**3 )
            GR[BIN] = float(HIST[BIN]) / float(len(self.xyzBoxes)) / float(len(atoms)) / float(NIDEAL)
        return np.arange(0.0, (MAXBIN+1)*DELR, DELR),GR

    def getRDF_arb(self,atom1='O',atom2='O',MAXBIN=1000,DELR=0.01,RHO=1.0):
        HIST = np.zeros(MAXBIN+1)
        GR   = np.zeros(MAXBIN+1)
        for box in self.xyzBoxes:
            atoms1 = box.get_atoms_by_name(atom1)
            atoms2 = box.get_atoms_by_name(atom2)
            HIST += calculateRIJhist_arb(atoms1, atoms2, box.lattice, box.lattice_inv,DELR=DELR,MAXBIN=MAXBIN)
        CONST = 4.0*np.pi*RHO/3.0
        for BIN in range(2,MAXBIN):
            RLOWER = (BIN-1)*DELR
            RUPPER = RLOWER + DELR
            NIDEAL = CONST * ( RUPPER**3 - RLOWER**3 )
            GR[BIN] = float(HIST[BIN]) / float(len(self.xyzBoxes)) / float(len(atoms1)+len(atoms2)) / float(NIDEAL)
        return np.arange(0.0, (MAXBIN+1)*DELR, DELR),GR

    def getRDF2_arb( self, atom1='O', atom2='O', MAXBIN=1000, DELR=0.01, RHO=1.0 ):
        HIST = np.zeros(MAXBIN+1)
        GR   = np.zeros(MAXBIN+1)
        for box in self.xyzBoxes:
            atoms1 = box.get_atoms_by_name(atom1)
            atoms2 = box.get_atoms_by_name(atom2)
            HIST += calculateRIJhist2_arb(atoms1, atoms2, box.lattice, box.lattice_inv, DELR=DELR, MAXBIN=MAXBIN )
        V = box.lattice[:,0].dot( np.cross(box.lattice[:,1],box.lattice[:,2]))
        rhoB   = len(atoms2)/V
        Nconf  = float(len(self.xyzBoxes))
        Na     = len(atoms1)
        for BIN in range(2,MAXBIN):
            HistAB  = HIST[BIN]
            Ri2     = (BIN*DELR)**2
            GR[BIN] = HistAB/ (4.0*np.pi) / rhoB / Ri2 / DELR / Nconf / Na
        return np.arange(0.0, (MAXBIN+1)*DELR, DELR),GR


def calculateRIJhist(atoms,boxLength,DELR=0.01,MAXBIN=1000):
    HIST = np.zeros(MAXBIN+1)
    for ii in range(len(atoms)-1):
        for jj in range(ii+1,len(atoms)):
            RIJ = getDistancePbc(atoms[ii],atoms[jj],boxLength)
            BIN = int(RIJ/DELR) + 1
            if BIN <= MAXBIN:
                HIST[BIN] = HIST[BIN] + 2
    return HIST

def calculateRIJhist_arb(atoms1, atoms2, lattice, lattice_inv,DELR=0.01,MAXBIN=1000):
    HIST = np.zeros(MAXBIN+1)
    for ii in range(len(atoms1)):
        for jj in range(len(atoms2)):
            RIJ = getDistancePBC_arb( atoms1[ii], atoms2[jj], lattice, lattice_inv )
            BIN = int(RIJ/DELR) + 1
            if BIN <= MAXBIN:
                HIST[BIN] = HIST[BIN] + 2
    return HIST

def calculateRIJhist2_arb( atoms1, atoms2, lattice, lattice_inv, DELR=0.01, MAXBIN=1000 ):
    HIST = np.zeros(MAXBIN+1)
    for ii in range(len(atoms1)):
        for jj in range(len(atoms2)):
            RIJ = getDistancePBC_arb( atoms1[ii], atoms2[jj], lattice, lattice_inv )
            BIN = int(RIJ/DELR) + 1
            if BIN <= MAXBIN:
                HIST[BIN] = HIST[BIN] + 1
    return HIST



def getDistsFromMolecule(point,listOfMolecules,o_name=None):
    dists = []
    if not o_name:
        for mol in listOfMolecules:
            dists.append(np.linalg.norm(point - mol.getGeometricCenter()))
        return dists
    else:
        for mol in listOfMolecules:
            coor = mol.getCoordinates_name(o_name)[0]
            dists.append(np.linalg.norm(point - coor))
        return dists
            
def findAllWaters(point,waterMols,o_name,cutoff):
    atoms = []
    for mol in waterMols:
        coor = mol.getCoordinates_name(o_name)
        if np.linalg.norm( point - coor) <= cutoff:
            atoms.extend(mol.xyzAtoms)
    return atoms

def findMolecule(xyzAtoms,molAtomList):
    molecule = []
    for atom in xyzAtoms:
        if atom.name in molAtomList:
            molecule.append(atom)
    return xyzMolecule(molecule)

def calculateCOMlist(atomList):
    """ **calculateCOMlist**
    Calculates center of mass for a list of atoms.
    """
    r   = np.array([0.,0.,0.])
    cou = 0
    for atom in atomList:
        r   += atom.coordinates
        cou += 1
    return r/cou


def getPeriodicTestBox_molecules(Molecules,boxLength,numbershells=1):
    vectors = []
    for l in range(-numbershells,numbershells+1):
        for m in range(-numbershells,numbershells+1):
            for n in range(-numbershells,numbershells+1):
                vectors.append(np.array([l, m, n]))
    pbc_molecules = []
    for vector in vectors:
        for molecule in Molecules:
            cpAtoms = []
            for atom in molecule.xyzAtoms:
                cpAtom = deepcopy(atom)
                cpAtom.translateSelf(vector*boxLength)
                cpAtoms.append(cpAtom)
            pbc_molecules.append(xyzMolecule(cpAtoms))
    return pbc_molecules

def getPeriodicTestBox(xyzAtoms,boxLength,numbershells=1):
    vectors = []
    for l in range(-numbershells,numbershells+1):
        for m in range(-numbershells,numbershells+1):
            for n in range(-numbershells,numbershells+1):
                vectors.append(np.array([l, m, n]))
    pbc_atoms = []
    for vector in vectors:
        for atom in xyzAtoms:
            cpAtom = deepcopy(atom)
            cpAtom.translateSelf(vector*boxLength)
            pbc_atoms.append(cpAtom)
    return pbc_atoms

def getPeriodicTestBox_arb( xyzAtoms, lattice, lattice_inv, numbershells=1 ):
    vectors = []
    for l in range(-numbershells,numbershells+1):
        for m in range(-numbershells,numbershells+1):
            for n in range(-numbershells,numbershells+1):
                vectors.append(np.array([l, m, n]))

    pbc_atoms = []
    for vector in vectors:
        for atom in xyzAtoms:
            cpAtom = deepcopy(atom)
            cpAtom.translateSelf_arb(lattice, lattice_inv, vector)
            pbc_atoms.append(cpAtom)
    return pbc_atoms


def getTranslVec(atom1,atom2,boxLength):
    """ **getTranslVec**
    Returns the translation vector that brings atom2 closer to atom1 in case
    atom2 is further than boxLength away.
    """
    xcoord1 = atom1.coordinates[0]
    xcoord2 = atom2.coordinates[0]
    ycoord1 = atom1.coordinates[1]
    ycoord2 = atom2.coordinates[1]
    zcoord1 = atom1.coordinates[2]
    zcoord2 = atom2.coordinates[2]
    translVec = np.zeros(atom2.coordinates.shape)
    if xcoord1-xcoord2 > boxLength/2.0:
        translVec[0] = -boxLength/2.0
    if ycoord1-ycoord2 > boxLength/2.0:
        translVec[1] = -boxLength/2.0
    if zcoord1-zcoord2 > boxLength/2.0:
        translVec[2] = -boxLength/2.0
    return translVec

def getTranslVec_geocen(mol1COM,mol2COM,boxLength):
    """ **getTranslVec_geocen**
    """
    translVec = np.zeros((len(mol1COM),))
    if mol1COM[0]-mol2COM[0] > boxLength/2.0:
        translVec[0] = -boxLength/2.0
    if mol1COM[1]-mol2COM[0] > boxLength/2.0:
        translVec[1] = -boxLength/2.0
    if mol1COM[2]-mol2COM[0] > boxLength/2.0:
        translVec[2] = -boxLength/2.0
    return translVec

def getDistancePbc(atom1,atom2,boxLength):
    xdist  = atom1.coordinates[0] - atom2.coordinates[0]
    xdist -= boxLength*round(xdist/boxLength)
    ydist  = atom1.coordinates[1] - atom2.coordinates[1] 
    ydist -= boxLength*round(ydist/boxLength)
    zdist  = atom1.coordinates[2] - atom2.coordinates[2] 
    zdist -= boxLength*round(zdist/boxLength)
    return np.sqrt(xdist**2.0 + ydist**2.0 + zdist**2.0)

def getDistancePBC_arb(atom1, atom2, lattice, lattice_inv):
    """ **getDistancePBC_arb**

    Calculates the distance of two atoms from an arbitrary 
    simulation box using the minimum image convention.

    Args:
        atom1 (obj): Instance of the xzyAtom class.
        atom2 (obj): Instance of the xzyAtom class.
        lattice (np.array): Array with lattice vectors as columns.
        lattice_inv (np.array): Inverse of lattice.

    Returns:
        The distance between the two atoms.

    """
    return np.linalg.norm(getDistVectorPBC_arb(atom1, atom2, lattice, lattice_inv))

def getDistVectorPBC_arb(atom1, atom2, lattice, lattice_inv):
    """ **getDistVectorPBC_arb**

    Calculates the distance vector between two atoms from an 
    arbitrary simulation box using the minimum image convention.

    Args:
        atom1 (obj): Instance of the xzyAtom class.
        atom2 (obj): Instance of the xzyAtom class.
        lattice (np.array): Array with lattice vectors as columns.
        lattice_inv (np.array): Inverse of lattice.

    Returns:
        The distance vector between the two atoms (np.array).

    """
    #frac_coords1 = np.dot( lattice_inv, atom1.coordinates )
    #frac_coords2 = np.dot( lattice_inv, atom2.coordinates )

    #red_frac1 = np.array(frac_coords1) - np.floor(frac_coords1)
    #red_frac2 = np.array(frac_coords2) - np.floor(frac_coords2)

    #red_dist = np.array(red_frac2) - np.array(red_frac1)
    #return np.dot(lattice, red_dist)
    dist_vec = np.array(atom2.coordinates) - np.array(atom1.coordinates)
    dist_vec -= np.dot(lattice, np.round(np.dot(lattice_inv,dist_vec)))
    return dist_vec

def getDistance(atom1, atom2):
    return np.linalg.norm(atom2.getCoordinates()-atom2.getCoordinates())

def getDistVectorPbc(atom1,atom2,boxLength):
    xdist  = atom1.coordinates[0] - atom2.coordinates[0]
    xdist -= boxLength*round(xdist/boxLength)
    ydist  = atom1.coordinates[1] - atom2.coordinates[1] 
    ydist -= boxLength*round(ydist/boxLength)
    zdist  = atom1.coordinates[2] - atom2.coordinates[2] 
    zdist -= boxLength*round(zdist/boxLength)
    return np.array([xdist, ydist, zdist])

def getDistVector(atom1,atom2):
    xdist  = atom1.coordinates[0] - atom2.coordinates[0]
    ydist  = atom1.coordinates[1] - atom2.coordinates[1] 
    zdist  = atom1.coordinates[2] - atom2.coordinates[2] 
    return np.array([xdist, ydist, zdist])

def count_HBonds_pbc_arb( mol1, mol2, lattice, lattice_inv, Roocut=3.6, Rohcut=2.4, Aoooh=30.0 ):
    hbond_angle = 0
    hb_donor  = 0
    hb_accept = 0
    # get atoms
    mol1_o = mol1.get_atoms_by_name('O')
    mol1_h = mol1.get_atoms_by_name('H')
    mol2_o = mol2.get_atoms_by_name('O')
    mol2_h = mol2.get_atoms_by_name('H')
    # O-O dist:
    OO_dist = getDistancePBC_arb(mol1_o[0],mol2_o[0], lattice, lattice_inv)
    if OO_dist <= Roocut and OO_dist > 0.0:
        # donor bond through first hydrogen atom of mol1
        OH_dist = getDistancePBC_arb(mol1_h[0],mol2_o[0], lattice, lattice_inv)
        if OH_dist <= Rohcut:
            vec1 = getDistVectorPBC_arb(mol1_h[0], mol1_o[0], lattice, lattice_inv)
            vec2 = getDistVectorPBC_arb(mol2_o[0], mol1_o[0], lattice, lattice_inv)
            angle = np.degrees( np.arccos( np.clip( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)), -1, 1 ) ))
            if angle <= Aoooh:
                hb_donor += 1.0
                hbond_angle = angle
        # donor bond through second hydrogen atom of mol1
        try: # catch if there is a OH- involved
            OH_dist = getDistancePBC_arb(mol1_h[1],mol2_o[0], lattice, lattice_inv)
            if OH_dist <= Rohcut:
                vec1 = getDistVectorPBC_arb(mol1_h[1], mol1_o[0], lattice, lattice_inv)
                vec2 = getDistVectorPBC_arb(mol2_o[0], mol1_o[0], lattice, lattice_inv)
                angle = np.degrees( np.arccos( np.clip( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)), -1, 1 ) ))
                if angle <= Aoooh:
                    hb_donor += 1.0
                    hbond_angle = angle 
        except:
            pass
        # acceptor bond through first hydrogen atom of mol2
        OH_dist = getDistancePBC_arb(mol1_o[0],mol2_h[0], lattice, lattice_inv)
        if OH_dist <= Rohcut:
            vec1 = getDistVectorPBC_arb(mol2_h[0], mol2_o[0], lattice, lattice_inv)
            vec2 = getDistVectorPBC_arb(mol1_o[0], mol2_o[0], lattice, lattice_inv)
            angle = np.degrees( np.arccos( np.clip( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)), -1, 1 ) ))
            if angle <= Aoooh:
                hb_accept += 1.0
                hbond_angle = angle
        # acceptor bond through second hydrogen atom of mol2
        try: # catch if there is a OH- involved
            OH_dist = getDistancePBC_arb(mol1_o[0],mol2_h[1], lattice, lattice_inv)
            if OH_dist <= Rohcut:
                vec1 = getDistVectorPBC_arb(mol2_h[1], mol2_o[0], lattice, lattice_inv)
                vec2 = getDistVectorPBC_arb(mol1_o[0], mol2_o[0], lattice, lattice_inv)
                angle = np.degrees( np.arccos( np.clip( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)), -1, 1 ) ))
                if angle <= Aoooh:
                    hb_accept += 1.0
                    hbond_angle = angle
        except:
            pass
        # try/catch if there is a OH3+ involved
        # acceptor bond through THIRD hydrogen atom of mol1
        try: # catch if there is a OH- involved
            OH_dist = getDistancePBC_arb(mol1_o[0],mol2_h[2], lattice, lattice_inv)
            if OH_dist <= Rohcut:
                vec1 = getDistVectorPBC_arb(mol2_h[2], mol2_o[0], lattice, lattice_inv)
                vec2 = getDistVectorPBC_arb(mol1_o[0], mol2_o[0], lattice, lattice_inv)
                angle = np.degrees( np.arccos( np.clip( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)), -1, 1 ) ))
                if angle <= Aoooh:
                    hb_accept += 1.0
                    hbond_angle = angle
        except:
            pass
        # donor bond through THIRD hydrogen atom of mol2
        try: # catch if there is a OH- involved
            OH_dist = getDistancePBC_arb(mol1_h[2],mol2_o[0], lattice, lattice_inv)
            if OH_dist <= Rohcut:
                vec1 = getDistVectorPBC_arb(mol1_h[2], mol1_o[0], lattice, lattice_inv)
                vec2 = getDistVectorPBC_arb(mol2_o[0], mol1_o[0], lattice, lattice_inv)
                angle = np.degrees( np.arccos( np.clip( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)), -1, 1 ) ))
                if angle <= Aoooh:
                    hb_donor += 1.0
                    hbond_angle = angle
        except:
            pass

    return hb_donor, hb_accept, hbond_angle

def countHbonds_pbc(mol1,mol2,boxLength,Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
    hb_donor  = 0.0
    hb_accept = 0.0
    # get atoms
    mol1_o = mol1.get_atoms_by_name('O')
    mol1_h = mol1.get_atoms_by_name('H')
    mol2_o = mol2.get_atoms_by_name('O')
    mol2_h = mol2.get_atoms_by_name('H')
    if getDistancePbc(mol1_o[0],mol2_o[0],boxLength) <= Roocut and getDistancePbc(mol1_o[0],mol2_o[0],boxLength) > 0.0:
        # donor bond through first hydrogen atom of mol1
        if getDistancePbc(mol1_h[0],mol2_o[0],boxLength) <= Rohcut:
            vec1 = getDistVectorPbc(mol1_h[0],mol1_o[0],boxLength)
            vec2 = getDistVectorPbc(mol2_o[0],mol1_o[0],boxLength)
            if np.degrees(np.arccos( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
                hb_donor += 1.0
        # donor bond through second hydrogen atom of mol1
        if getDistancePbc(mol1_h[1],mol2_o[0],boxLength) <= Rohcut:
            vec1 = getDistVectorPbc(mol1_h[1],mol1_o[0],boxLength)
            vec2 = getDistVectorPbc(mol2_o[0],mol1_o[0],boxLength)
            if np.degrees(np.arccos( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
                hb_donor += 1.0
        # acceptor bond through first hydrogen atom of mol2
        if getDistancePbc(mol1_o[0],mol2_h[0],boxLength) <= Rohcut:
            vec1 = getDistVectorPbc(mol2_h[0],mol2_o[0],boxLength)
            vec2 = getDistVectorPbc(mol1_o[0],mol2_o[0],boxLength)
            if np.degrees(np.arccos( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
                hb_accept += 1.0
        # acceptor bond through second hydrogen atom of mol2
        if getDistancePbc(mol1_o[0],mol2_h[1],boxLength) <= Rohcut:
            vec1 = getDistVectorPbc(mol2_h[1],mol2_o[0],boxLength)
            vec2 = getDistVectorPbc(mol1_o[0],mol2_o[0],boxLength)
            if np.degrees(np.arccos( np.dot(vec2,vec1)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
                hb_accept += 1.0
    return hb_donor, hb_accept

def countHbonds(mol1,mol2, Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
    hb_donor  = 0
    hb_accept = 0
    # get the coordinates
    Aoooh = np.radians(Aoooh)
    mol1_o = mol1.getCoordinates_name('O')
    mol1_h = mol1.getCoordinates_name('H')
    mol2_o = mol2.getCoordinates_name('O')
    mol2_h = mol2.getCoordinates_name('H')
    # check O-O distance
    if np.linalg.norm(mol2_o[0] - mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0] - mol1_o[0]) <= Roocut: # check Roocut is met
        # check for donor H-bonds
        if np.linalg.norm(mol2_o[0] - mol1_h[0]) > 0.0 and np.linalg.norm(mol2_o[0] - mol1_h[0]) <= Rohcut: # check Rohcut, first H
            vec1 = mol1_h[0]-mol1_o[0]
            vec2 = mol2_o[0]-mol1_o[0]
            if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
                hb_donor += 1
        if np.linalg.norm(mol2_o[0] - mol1_h[1]) > 0.0 and np.linalg.norm(mol2_o[0] - mol1_h[1]) <= Rohcut: # check Rohcut, second H
            vec1 = mol1_h[1]-mol1_o[0]
            vec2 = mol2_o[0]-mol1_o[0]
            if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
                hb_donor += 1
        # check for acceptor H-bonds
        if np.linalg.norm(mol2_h[0] - mol1_o[0]) > 0.0 and np.linalg.norm(mol2_h[0] - mol1_o[0]) <= Rohcut: # check Rohcut, first H
            vec1 = mol1_o[0]-mol2_o[0]
            vec2 = mol2_h[0]-mol2_o[0]
            if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
                hb_accept += 1
        if np.linalg.norm(mol2_h[1] - mol1_o[0]) > 0.0 and np.linalg.norm(mol2_h[1] - mol1_o[0]) <= Rohcut: # check Rohcut, socond H
            vec1 = mol1_o[0]-mol2_o[0]
            vec2 = mol2_h[1]-mol2_o[0]
            if np.arccos((np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
                hb_accept += 1
    return hb_donor, hb_accept

def countHbonds_orig(mol1,mol2, Roocut=3.6, Rohcut=2.4, Aoooh=30.0):
    mol1_o = mol1.getCoordinates_name('O')
    mol1_h = mol1.getCoordinates_name('H')
    mol2_o = mol2.getCoordinates_name('O')
    mol2_h = mol2.getCoordinates_name('H')
    hbnoD1 = 0
    if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
        if np.linalg.norm(mol2_h[0] - mol1_o[0]) <= Rohcut:
            vec1 = mol2_h[0]-mol2_o[0]
            vec2 = mol1_o[0]-mol2_o[0]
            if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
                hbnoD1 += 1
    hbnoD2 = 0
    if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
        if np.linalg.norm(mol2_h[1]-mol1_o[0]) <= Rohcut:
            vec1 = mol2_h[1]-mol2_o[0]
            vec2 = mol1_o[0]-mol2_o[0]
            if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec2)*np.linalg.norm(vec1)))) <= Aoooh:
                hbnoD2 += 1
    hbnoA1 = 0
    if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
        if np.linalg.norm(mol2_o[0]-mol1_h[0]) <= Rohcut:
            vec1 = mol2_o[0]-mol1_o[0]
            vec2 = mol1_h[0]-mol1_o[0]
            if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
                hbnoA1 += 1
    hbnoA2 = 0
    if np.linalg.norm(mol2_o[0]-mol1_o[0]) > 0.0 and np.linalg.norm(mol2_o[0]-mol1_o[0]) <= Roocut:
        if np.linalg.norm(mol2_o[0]-mol1_h[1]) <= Rohcut:
            vec1 = mol2_o[0]-mol1_o[0]
            vec2 = mol1_h[1]-mol1_o[0]
            if np.degrees(np.arccos( np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2)))) <= Aoooh:
                hbnoA2 += 1
    dbonds = hbnoD1 + hbnoD2
    abonds = hbnoA1 + hbnoA2
    return dbonds, abonds

def getTetraParameter(o_atoms,boxLength=None):
    """
    according to NATURE, VOL 409, 18 JANUARY 2001
    """
    tetra_params = []
    for atom1 in o_atoms:
        NN_dists = []
        NN_atoms = []
        for atom2 in o_atoms:
            NN_dists.append(np.linalg.norm(atom1.coordinates - atom2.coordinates))
        order = np.argsort(NN_dists)
        for ii in order[1:5]:
            NN_atoms.append(o_atoms[ii])
        tetra_param = 0
        for j in range(0,3):
            for k in range(j+1,4):
                vec1 = getDistVectorPbc(atom1,NN_atoms[j],boxLength)
                vec2 = getDistVectorPbc(atom1,NN_atoms[k],boxLength)
                psi = np.arccos(np.dot(vec1,vec2)/np.linalg.norm(vec1)/np.linalg.norm(vec2))
                tetra_param += (np.cos(psi) + 1./3.)**2
        tetra_params.append(1-3./8.*tetra_param)
    return tetra_params

def repair_h2o_molecules_pbc(h2o_mols,boxLength):
    new_mols = []
    for mol in h2o_mols:
        o_atom  = mol.get_atoms_by_name('O')[0]
        h_atoms = mol.get_atoms_by_name('H')
        new_mol = [o_atom]
        for h_atom in h_atoms:
            cpAtom = deepcopy(h_atom)
            xdist  = o_atom.coordinates[0] - cpAtom.coordinates[0]
            ydist  = o_atom.coordinates[1] - cpAtom.coordinates[1]
            zdist  = o_atom.coordinates[2] - cpAtom.coordinates[2]
            cpAtom.translateSelf([boxLength*round(xdist/boxLength),boxLength*round(ydist/boxLength),boxLength*round(zdist/boxLength)])
            new_mol.append(cpAtom)
        new_mols.append(xyzMolecule(new_mol))
    return new_mols

def find_H2O_molecules(o_atoms,h_atoms,boxLength=None):
    h2o_molecules = []
    if not boxLength:
        warnings.warn('No box length provided, will not take PBC into account!')
        for o_atom in o_atoms:
            ho_dists = []
            for h_atom in h_atoms:
                ho_dists.append(np.linalg.norm(o_atom.coordinates - h_atom.coordinates))
            order = np.argsort(ho_dists)
            if np.linalg.norm(o_atom.coordinates - h_atoms[order[0]].getCoordinates()) <= 1.5 and np.linalg.norm(o_atom.coordinates - h_atoms[order[1]].getCoordinates()) <= 1.5:
                h2o_molecules.append(xyzMolecule([o_atom,h_atoms[order[0]],h_atoms[order[1]]]))
        return h2o_molecules
    else:
        for o_atom in o_atoms:
            ho_dists = []
            for h_atom in h_atoms:
                ho_dists.append(getDistancePbc(o_atom,h_atom,boxLength))
            order = np.argsort(ho_dists)
            h2o_molecules.append(xyzMolecule([o_atom,h_atoms[order[0]],h_atoms[order[1]]]))
        return h2o_molecules

def find_H2O_molecules_PBC_arb( o_atoms, h_atoms, lattice, lattice_inv, OH_cutoff=1.5 ):
    h2o_molecules = []
    for o_atom in o_atoms:
        ho_dists = []
        for h_atom in h_atoms:
            ho_dists.append(getDistancePBC_arb(o_atom, h_atom, lattice, lattice_inv))
        inds = np.where(np.array(ho_dists) <= OH_cutoff)[0]
        molecule = []
        molecule.append(o_atom)
        for ind in inds:
            molecule.append(h_atoms[ind])
        h2o_molecules.append(xyzMolecule(molecule))
    return h2o_molecules

def writeXYZfile(filename,numberOfAtoms, title, list_of_xyzAtoms):
        # create file
        xyz = open(filename,'w+')
        # write number of atoms
        xyz.write(str(numberOfAtoms) + ' \n')
        # write title
        if not title:
            title = 'None'
        xyz.write(title + ' \n')
        # write coordinates
        for atom in list_of_xyzAtoms:
                xyz.write('%4s %12.8f %12.8f %12.8f \n' % (atom.name, atom.x_coord, atom.y_coord, atom.z_coord))
        xyz.close()

def writeXYZtrajectory(filename,boxes):
    # create file
    xyz = open(filename,'w+')
    cou = 0
    for box in boxes:
        xyz.write(str(len(box.xyzAtoms)) + ' \n')
        xyz.write('Step: ' + str(cou) + ' , boxLength = ' + str(box.boxLength) + ' \n')
        for atom in box.xyzAtoms:
            xyz.write('%6.4s %14.8f %14.8f %14.8f \n' % (atom.name, atom.x_coord, atom.y_coord, atom.z_coord))
        cou += 1
    xyz.close()


def writeRelXYZfile(filename, n_atoms, boxLength, title, xyzAtoms, inclAtomNames=True):
    # create file
    xyz = open(filename,'w+')
    # write number of atoms
    xyz.write(str(n_atoms) + ' \n')
    # write title
    if not title:
        title = 'None'
    xyz.write(title + ' \n')
    # write coordinates
    for atom in xyzAtoms:
        if inclAtomNames:
            xyz.write('%4s %8f %8f %8f \n' % (atom.name, atom.x_coord/boxLength, atom.y_coord/boxLength, atom.z_coord/boxLength))
        else:
            xyz.write('%8f %8f %8f \n' % (atom.x_coord/boxLength, atom.y_coord/boxLength, atom.z_coord/boxLength))
    xyz.close()


def count_OO_neighbors(list_of_o_atoms,Roocut,boxLength=None):
    noo   = []
    if not boxLength:
            warnings.warn('No box length provided, will not take PBC into account!')
            for atom1 in list_of_o_atoms:
                    dists = []
                    for atom2 in list_of_o_atoms:
                            dists.append(np.linalg.norm(atom1.coordinates - atom2.coordinates))
                    noo.append(len(np.where(np.logical_and(np.sort(np.array(dists))>0.0, np.sort(np.array(dists))<=Roocut))[0]))
            return noo
    else:
        for atom1 in list_of_o_atoms:
            dists = []
            for atom2 in list_of_o_atoms:
                dists.append(getDistancePbc(atom1,atom2,boxLength))
            noo.append(len(np.where(np.logical_and(np.sort(np.array(dists))>0.0, np.sort(np.array(dists))<=Roocut))[0]))
        return noo

def count_OO_neighbors_pbc(list_of_o_atoms,Roocut,boxLength,numbershells=1):
    noo = []
    vectors = []
    for l in range(-numbershells,numbershells+1):
        for m in range(-numbershells,numbershells+1):
            for n in range(-numbershells,numbershells+1):
                vectors.append(np.array([l, m, n]))
    pbc_atoms = []
    for vector in vectors:
        for atom in list_of_o_atoms:
            cpAtom = deepcopy(atom)
            cpAtom.translateSelf(vector*boxLength)
            pbc_atoms.append(cpAtom)
    for atom1 in list_of_o_atoms:
        dists = []
        for atom2 in pbc_atoms:
            dists.append(np.linalg.norm(atom1.coordinates - atom2.coordinates))
        sdists = np.sort(np.array(dists))
        noo.append(len(np.where(np.logical_and(sdists>0.0, sdists<=Roocut))[0]))
    return noo

def boxParser(filename):
        """**parseXYZfile**
        Reads an xyz-style file.
        """
        atoms = []
        coordinates = []
        xyz = open(filename)
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        for line in xyz:
                if len(line.split())==4:
                        atom,x,y,z = line.split()[0:4]
                        atoms.append(atom)
                        coordinates.append([float(x), float(y), float(z)])
                else:
                        pass
        xyz.close()
        xyzAtoms = []
        for ii in range(n_atoms):
                xyzAtoms.append(xyzAtom(atoms[ii],coordinates[ii],ii))

        return xyzBox(xyzAtoms)

def keithBoxParser(cell_fname, coord_fname):
    """ **keithBoxParser**

    Reads structure files from Keith's SiO2 simulations.

    """

    # lattice in Angstr.
    lattice = np.loadtxt(cell_fname)*A2AU_factor

    absCoords = []
    atoms     = []
    xyz       = open(coord_fname)
    for line in xyz:
        atom,x,y,z = line.split()[0:4]
        atoms.append(atom)
        absCoords.append([float(x), float(y), float(z)])

    xyz.close()
    xyzAtoms = []
    relAtoms = []
    for ii in range(len(atoms)):
        xyzAtoms.append(xyzAtom( atoms[ii], np.array(absCoords[ii])*A2AU_factor, ii))
        relAtoms.append(xyzAtom( atoms[ii], np.dot(np.array(absCoords[ii]),np.linalg.inv(lattice)) , ii))
    box = xyzBox(xyzAtoms)
    box.lattice = lattice
    box.lattice_inv = np.linalg.inv(lattice)
    box.relAtoms = relAtoms
    return box


def axsfTrajParser(filename):
    """ **axsfTrajParser**

    """
    boxes       = []
    boxLength   = None
    alat        = [[0,0,0],[0,0,0],[0,0,0]]
    title       = None
    # open the file
    xyz     = open(filename)
    counter = 0
    line = xyz.readline()
    while line:
        if 'ANIMSTEPS' in line:
            numBoxes = int(line.split()[-1])
        if 'CRYSTAL' in line:
            line = xyz.readline()
            for ii in range(3):
                latt = xyz.readline().split()
                alat[ii] = [float(latt[0]), float(latt[1]), float(latt[2])]
        if 'PRIMCOORD' in line:
            title       = line.strip()
            line        = xyz.readline()
            n_atoms     = int(line.split()[0])
            xyzAtoms    = []
            atoms       = []
            coordinates = []
            # read in coordinates
            for ii in range(n_atoms):    
                line = xyz.readline()
                atoms.append(xrs_utilities.element(int(line.split()[0])))
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                coordinates.append([x,y,z])
                xyzAtoms.append(xyzAtom(atoms[ii],coordinates[ii],ii))
            boxes.append(xyzBox(xyzAtoms,np.amax(alat),title))
        line = xyz.readline()
    return boxes

def writeWFN1waterInput(fname,box,headerfile,exatomNo=0):
    """ **writeWFN1input**
    Writes an input for cp.x by Quantum espresso for electronic wave function minimization.
    """
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)
    inputf.write('CELL_PARAMETERS cubic \n')
    A2Bfac = constants.physical_constants['atomic unit of length'][0]*10**10
    inputf.write('%20.16f %20.16f %20.16f \n' % (box.boxLength/A2Bfac, 0.0, 0.0))
    inputf.write('%20.16f %20.16f %20.16f \n' % (0.0, box.boxLength/A2Bfac, 0.0))
    inputf.write('%20.16f %20.16f %20.16f \n' % (0.0, 0.0, box.boxLength/A2Bfac))
    inputf.write('ATOMIC_POSITIONS crystal \n')
    for atom,ii in zip(box.xyzAtoms,list(range(len(box.xyzAtoms)))):
        if ii == exatomNo:
            if atom.name == 'O':
                inputf.write('%4s %20.16f %20.16f %20.16f \n' % ('Ox', atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
            else:
                print('Not writing Ox, since chosen atom is not an oxygen atom.')
                inputf.write('%4s %20.16f %20.16f %20.16f \n' % (atom.name, atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
        else:
            inputf.write('%4s %20.16f %20.16f %20.16f \n' % (atom.name, atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
    inputf.close()
    header.close()

def writeMD1Input(fname,box,headerfile,exatomNo=0):
    """ **writeWFN1input**
    Writes an input for cp.x by Quantum espresso for electronic wave function minimization.
    """
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)
    inputf.write('CELL_PARAMETERS cubic \n')
    A2Bfac = constants.physical_constants['atomic unit of length'][0]*10**10
    inputf.write('%20.16f %20.16f %20.16f \n' % (box.boxLength/A2Bfac, 0.0, 0.0))
    inputf.write('%20.16f %20.16f %20.16f \n' % (0.0, box.boxLength/A2Bfac, 0.0))
    inputf.write('%20.16f %20.16f %20.16f \n' % (0.0, 0.0, box.boxLength/A2Bfac))
    inputf.write('ATOMIC_POSITIONS crystal \n')
    for atom,ii in zip(box.xyzAtoms,list(range(len(box.xyzAtoms)))):
        if ii == exatomNo:
            if atom.name == 'O':
                inputf.write('%4s %20.16f %20.16f %20.16f \n' % ('C', atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
            else:
                print('Not writing Ox, since chosen atom is not an oxygen atom.')
                inputf.write('%4s %20.16f %20.16f %20.16f \n' % (atom.name, atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
        else:
            inputf.write('%4s %20.16f %20.16f %20.16f \n' % (atom.name, atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
    inputf.close()
    header.close()



def writeOCEAN_XESInput(fname,box,headerfile,exatomNo=0):
    """ **writeOCEAN_XESInput**
    Writes an input for ONEAN XES calculation for 17 molecule water boxes.
    """
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)
    inputf.write('xred { \n')
    for atom,ii in zip(box.xyzAtoms,list(range(len(box.xyzAtoms)))):
        inputf.write('%20.16f %20.16f %20.16f \n' % ( atom.x_coord/box.boxLength, atom.y_coord/box.boxLength, atom.z_coord/box.boxLength))
    inputf.write('} \n')

def groBoxParser(filename,nanoMeter=True):
    """ **groBoxParser**
    Parses an gromacs GRO-style file for the xyzBox class.
    """
    atoms       = []
    coordinates = []
    boxLength   = None
    if nanoMeter:
        scale = 10.0
    else:
        scale = 1.0
    xyz     = open(filename)
    title   = xyz.readline()
    n_atoms = int(xyz.readline())
    for line in xyz:
        if len(line.split()) == 9:
            name = line.split()[1][0]
            x = float(line.split()[3])*scale
            y = float(line.split()[4])*scale
            z = float(line.split()[5])*scale
            atoms.append(name)
            coordinates.append([x, y, z])
        elif len(line.split()) == 8:
            name = line.split()[1][0]
            x = float(line.split()[2])*scale
            y = float(line.split()[3])*scale
            z = float(line.split()[4])*scale
            atoms.append(name)
            coordinates.append([x, y, z])
        elif line[0] == '#':
            pass
        elif len(line.split()) == 3:
            boxLength = float(line.split()[0])*scale
        else:
            print('Something is fishy!')
    xyz.close()
    AllxyzAtoms = []
    for ii in range(n_atoms):
        AllxyzAtoms.append(xyzAtom(atoms[ii],coordinates[ii],ii))

    return xyzBox(AllxyzAtoms,boxLength=boxLength)

def xyzTrajecParser(filename, boxLength):
    """Parses a Trajectory of xyz-files.

    Args:
        filename (str): Filename of the xyz Trajectory file.

    Returns:
        A list of xzyBoxes.
    """
    boxes       = []
    # read the file
    xyz     = open(filename)
    n_atoms = int(xyz.readline())
    title   = xyz.readline()
    # headerlines
    headerlines = 2
    lines_per_box = n_atoms + headerlines
    # reopen file and start from scratch
    xyz.close()
    xyz = open(filename)
    counter = 0
    for line in xyz:
        if counter%lines_per_box == 0:
            # new box 
            n_atoms = int(line)
            atoms = []
            coordinates = []
            cou = 0
        if counter%lines_per_box == 1:
            title = line
        if counter%lines_per_box in list(range(lines_per_box))[2::]:
            atom,x,y,z = line.split()[0:4]
            #print(atom,x,y,z)
            atoms.append(xyzAtom(atom,[float(x), float(y), float(z)],cou))
            #coordinates.append([float(x), float(y), float(z)])
            cou += 1
        if counter%lines_per_box == lines_per_box-1:
            boxes.append(xyzBox(atoms,boxLength=boxLength))
        counter += 1
    return boxes

def vaspBoxParser(filename):
    """ **groTrajecParser**
    Parses an gromacs GRO-style file for the xyzBox class.
    """
    all_lines   = open(filename).readlines()
    atom_kinds  = [ii for ii in all_lines[5].split()]
    num_atom_kind = len(atom_kinds)
    num_atoms = np.array([int(ii) for ii in all_lines[6].split()])
    headerlines = 2
    box_latt_lines = 3
    headerlines2 = 3
    linesPerBox = headerlines + box_latt_lines + headerlines2 + np.sum(num_atoms)
    cell1 = np.array([float(ii) for ii in all_lines[2].split()])
    cell2 = np.array([float(ii) for ii in all_lines[3].split()])
    cell3 = np.array([float(ii) for ii in all_lines[4].split()]) 
    cell = np.array([cell1, cell2, cell3])
    coordinates = []
    atoms = []
    coor_start = headerlines + box_latt_lines + headerlines2
    coor_stop  = headerlines + box_latt_lines + headerlines2+np.sum(num_atoms)
    cou = 0
    for zz,kind in enumerate(atom_kinds):
        for ind in range(num_atoms[zz]):
            # print(ind, kind, cou)
            rel_coords = np.array([float(ii) for ii in all_lines[coor_start+cou].split()])
            abs_coords = cell.dot(rel_coords)
            atoms.append(xyzAtom(kind, abs_coords, cou ))
            cou += 1
    box = xyzBox(atoms,boxLength=0.0)
    box.lattice = np.array([cell1, cell2, cell3])
    box.lattice_inv = np.linalg.inv(box.lattice)
    return box

def vaspTrajecParser(filename, min_boxes=0, max_boxes=1000):
    """ **groTrajecParser**
    Parses an gromacs GRO-style file for the xyzBox class.
    """
    boxes       = []
    atoms       = []
    coordinates = []
    all_lines   = open(filename).readlines()
    atom_kinds  = [ii for ii in all_lines[5].split()]
    num_atom_kind = len(atom_kinds)
    num_atoms = np.array([int(ii) for ii in all_lines[6].split()])
    headerlines = 2
    box_latt_lines = 3
    headerlines2 = 3
    linesPerBox = headerlines + box_latt_lines + headerlines2 + np.sum(num_atoms)
    del(all_lines)
    xyz = open(filename)
    counter = 0
    for line in xyz:
        if counter//linesPerBox>= min_boxes:
            if counter%linesPerBox == 0:
                #print('hehe1')
                # write a new box if there is something to write
                if len(coordinates) > 0:
                    coorcou = 0
                    for ii_kind in range(num_atom_kind):
                        for jj_num in range(num_atoms[ii_kind]):
                            atoms.append(xyzAtom( ([atom_kinds[ii_kind]]*num_atoms[ii_kind])[jj_num] ,coordinates[coorcou],counter))
                            coorcou += 1
                    box = xyzBox(atoms,boxLength=0.0)
                    box.lattice = np.array([cell1, cell2, cell3])
                    box.lattice_inv = np.linalg.inv(box.lattice)
                    boxes.append(box)
                # if first box, start here
                else:
                    pass
                # new box 
                title = line.strip()
                atoms = []
                coordinates = []
                atom_type_counter = np.zeros_like(num_atoms)
                if counter/linesPerBox >= max_boxes:
                    return boxes
            if counter%linesPerBox == 1:
                title2 = line.strip()
            if counter%linesPerBox == 2:
                cell1 = np.array([float(ii) for ii in line.split()])
            if counter%linesPerBox == 3:
                cell2 = np.array([float(ii) for ii in line.split()])
            if counter%linesPerBox == 4:
                cell3 = np.array([float(ii) for ii in line.split()]) 
            if counter%linesPerBox == 5:
                atom_kinds  = [ii for ii in line.split()]
            if counter%linesPerBox == 6:
                num_atoms = np.array([int(ii) for ii in line.split()])
            if counter%linesPerBox == 7:
                title3 = line.strip()
            if not counter%linesPerBox in [0,1,2,3,4,5,6,7]:
                #print(counter, line.split())
                cell = np.array([cell1, cell2, cell3])
                rel_coords = np.array([float(ii) for ii in line.split()])
                abs_coords = cell.dot(rel_coords)
                coordinates.append(abs_coords)
        counter += 1
    return boxes
    
def groTrajecParser(filename,nanoMeter=True):
    """ **groTrajecParser**
    Parses an gromacs GRO-style file for the xyzBox class.
    """
    boxes       = []
    atoms       = []
    coordinates = []
    boxLength   = None
    if nanoMeter:
        scale = 10.0
    else:
        scale = 1.0
    # read the header
    xyz     = open(filename)
    title   = xyz.readline()
    n_atoms = int(xyz.readline())
    xyz.close()
    # how many lines per box
    headerlines = 2
    taillines   = 1
    linesPerBox = headerlines + n_atoms + taillines
    # read the rest
    xyz     = open(filename)
    counter = 0
    for line in xyz:
        if counter%linesPerBox == 0:
            # new box 
            title = line
            atoms = []
            coordinates = []
        if counter%linesPerBox == 1:
            n_atoms = int(line)
        if counter in list(range(linesPerBox))[2::]:
            if len(line.split()) == 9:
                atoms.append(line.split()[1][0])
                x = float(line.split()[3])*scale
                y = float(line.split()[4])*scale
                z = float(line.split()[5])*scale
                coordinates.append([x, y, z])
            elif len(line.split()) == 8:
                atoms.append(line.split()[1][0])
                x = float(line.split()[2])*scale
                y = float(line.split()[3])*scale
                z = float(line.split()[4])*scale
                coordinates.append([x, y, z])
            elif len(line.split()) == 6:
                atoms.append(line.split()[1][0])
                x = float(line.split()[3])*scale
                y = float(line.split()[4])*scale
                z = float(line.split()[5])*scale
                coordinates.append([x, y, z])
        if counter > 1 and counter%(linesPerBox-1) == 0:
            boxLength = float(line.split()[0])*scale
            AllxyzAtoms = []
            for ii in range(n_atoms):
                AllxyzAtoms.append(xyzAtom(atoms[ii],coordinates[ii],ii))
            boxes.append(xyzBox(AllxyzAtoms,boxLength=boxLength))
            counter = -1
        counter += 1
    return boxes

def changeOHBondLength(h2oMol, fraction, boxLength=None, oName='O', hName='H'):
    o_atom = h2oMol.get_atoms_by_name(oName)
    h_atom = h2oMol.get_atoms_by_name(hName)
    # get OH bond-vectors
    ohVectors = []
    if not boxLength:
        warnings.warn('No box length provided, will not take PBC into account!')
        for atom in h_atom:
            ohVectors.append(getDistVector(atom,o_atom[0]))
    else:
        for atom in h_atom:
            ohVectors.append(getDistVectorPbc(atom,o_atom[0],boxLength))
    # get new OH bond vectors
    new_ohVectors = []
    for vector in ohVectors:
        new_ohVectors.append( vector*fraction)
    # make new molecule
    newmol =[deepcopy(o_atom[0])]
    for vector in new_ohVectors:
        newmol.append(xyzAtom(hName,o_atom[0].getCoordinates()+vector,1))
    return xyzMolecule(newmol)


def parseXYZfile(filename):
        """**parseXYZfile**
        Reads an xyz-style file.
        """
        atoms = []
        coordinates = []
        xyz = open(filename)
        n_atoms = int(xyz.readline())
        title = xyz.readline()
        for line in xyz:
                atom,x,y,z = line.split()
                atoms.append(atom)
                coordinates.append([float(x), float(y), float(z)])
        xyz.close()
        return n_atoms, title, atoms, coordinates

def alterGROatomNames(filename,oldName, newName):
    import re
    with open(filename, "r") as sources:
        lines = sources.readlines()
    with open(filename, "w") as sources:
        for line in lines:
            sources.write(re.sub(r'%s'%oldName, '%s'%newName, line))

def writeFDMNESinput_file(xyzAtoms,fname,Filout,Range,Radius,Edge,NRIXS,Absorber):
    """ **writeFDMNESinput_file**
    Writes an input file to be used for FDMNES.
    """
    # create file
    inp = open(fname,'w+')
    # write some header line
    inp.write('! Fdmnes indata file \n')
    inp.write('\n')
    # Filout
    inp.write(' Filout \n')
    inp.write('  ' + Filout + ' \n')
    inp.write('\n')
    # Range
    inp.write(' Range \n')
    inp.write('  ')
    for rr in [str(ii) for ii in Range]:
        inp.write(rr + ' ')
    inp.write(' \n')
    inp.write('\n')
    # Radius
    inp.write(' Radius \n')
    inp.write('  ' + str(Radius) + ' \n')
    inp.write('\n')
    # Edge
    inp.write(' Edge \n')
    inp.write('  ' + Edge + ' \n')
    inp.write('\n')
    # NRIXS
    if NRIXS:
        inp.write(' NRIXS \n')
        inp.write('  ')
        for nn in [str(ii) for ii in NRIXS]:
            inp.write(nn + ' ')
        inp.write(' \n')
        inp.write('\n')
    # Absorber
    inp.write(' Absorber \n')
    inp.write('  ' + str(Absorber) + ' \n')
    inp.write('\n')
    # molecule
    inp.write(' molecule \n')
    inp.write('  1.0 1.0 1.0 90.0 90.0 90.0 \n')
    # Atoms
    for atom in xyzAtoms:
        inp.write('%4d %10f %10f %10f \n' % (xrs_utilities.element(atom.name), atom.x_coord, atom.y_coord, atom.z_coord))
    inp.write('\n')
    # Convolution
    inp.write(' Convolution \n')
    inp.write('\n')
    # End
    inp.write(' End \n')
    inp.write('\n')

def writeOCEANinput(fname,headerfile,xyzBox,exatom,edge,subshell):
    """ **writeOCEANinput**
    """
    # write everything that is in the header
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)

    inputf.write('\n')

    # write the cell
    const = constants.physical_constants['atomic unit of length'][0]*10**10
    inputf.write('acell { %10.8f %10.8f %10.8f } \n'%(xyzBox.boxLength/const, xyzBox.boxLength/const, xyzBox.boxLength/const))
    inputf.write('\n')
    # write natoms
    inputf.write('natom %d \n' %len(xyzBox.xyzAtoms))
    inputf.write('\n')
    # write typat
    inputf.write('typat { \n')
    for atom in xyzBox.xyzAtoms:
        if atom.name == 'H':
            inputf.write('%d ' %1)
        if atom.name == 'O':
            inputf.write('%d ' %2)
        if atom.name == 'N':
            inputf.write('%d ' %3)
        if atom.name == 'C':
            inputf.write('%d ' %4)
        if atom.name == 'Cl':
            inputf.write('%d ' %3)
        if atom.name == 'Na':
            inputf.write('%d ' %4)
    inputf.write('\n')
    inputf.write('} \n')
    inputf.write('\n')
    # write nedges
    ind = 0
    for atom in xyzBox.xyzAtoms:
        if atom.name == exatom:
            ind += 1
    inputf.write('nedges %d \n' %ind)
    inputf.write('\n')
    # write edges
    inputf.write('edges { \n')
    ind = 1
    for atom in xyzBox.xyzAtoms:
        if atom.name == exatom:
            inputf.write('%d %d %d \n' %(ind,edge,subshell))
        ind += 1
    inputf.write('} \n')
    inputf.write('\n')
    # write coordinates
    inputf.write('xred { \n')
    for atom in xyzBox.xyzAtoms:
        inputf.write('%16f %16f %16f \n' % (atom.x_coord/xyzBox.boxLength, atom.y_coord/xyzBox.boxLength, atom.z_coord/xyzBox.boxLength))
    inputf.write('} \n')
    inputf.write('\n')
    inputf.close()

def writeOCEANinput_new(fname,headerfile,xyzBox,exatom,edge,subshell):
    """ **writeOCEANinput**
    """
    # write everything that is in the header
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)

    inputf.write('\n')

    # write the cell
    const = constants.physical_constants['atomic unit of length'][0]*10**10
    inputf.write('acell { %10.8f %10.8f %10.8f } \n'%(xyzBox.boxLength/const, xyzBox.boxLength/const, xyzBox.boxLength/const))
    inputf.write('\n')
    # write natoms
    inputf.write('natom %d \n' %len(xyzBox.xyzAtoms))
    inputf.write('\n')
    # write typat
    inputf.write('typat { \n')
    all_atom_kinds = []
    for atom in xyzBox.xyzAtoms:
        if not atom.name in all_atom_kinds:
            all_atom_kinds.append(atom.name)
        inputf.write('%d '%(all_atom_kinds.index(atom.name)+1))
        # if atom.name == 'H':
        #     inputf.write('%d ' %1)
        # if atom.name == 'O':
        #     inputf.write('%d ' %2)
        # if atom.name == 'N':
        #     inputf.write('%d ' %6)
        # if atom.name == 'C':
        #     inputf.write('%d ' %5)
        # if atom.name == 'Cl':
        #     inputf.write('%d ' %3)
        # if atom.name == 'Na':
        #     inputf.write('%d ' %4)
    inputf.write('\n')
    inputf.write('} \n')
    inputf.write('\n')
    # write nedges
    #ind = 0
    #for atom in xyzBox.xyzAtoms:
    #    if atom.name == exatom:
    #        ind += 1
    #inputf.write('nedges %d \n' %ind)
    #inputf.write('\n')
    # write edges
    inputf.write('edges { \n')
    #ind = 1
    #for atom in xyzBox.xyzAtoms:
    #    if atom.name == exatom:
    inputf.write('%d %d %d \n' %(-xrs_utilities.element(exatom),edge,subshell))
    #    ind += 1
    inputf.write('} \n')
    inputf.write('\n')
    # write coordinates
    inputf.write('xred { \n')
    for atom in xyzBox.xyzAtoms:
        inputf.write('%16f %16f %16f \n' % (atom.x_coord/xyzBox.boxLength, atom.y_coord/xyzBox.boxLength, atom.z_coord/xyzBox.boxLength))
    inputf.write('} \n')
    inputf.write('\n')
    inputf.close()

def writeOCEANinput_arb(fname, headerfile, xyzBox, exatom, edge,subshell):
    """ **writeOCEANinput**
    """
    # write everything that is in the header
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)

    inputf.write('\n')

    # write the cell
    const = constants.physical_constants['atomic unit of length'][0]*10**10
    inputf.write('acell { 1.0 1.0 1.0 }' )
    inputf.write('\n')
    # write rprim
    inputf.write('rprim { \n')
    inputf.write('%10.8f %10.8f %10.8f \n'%(xyzBox.lattice[0,0]/const, xyzBox.lattice[0,1]/const, xyzBox.lattice[0,2]/const))
    inputf.write('%10.8f %10.8f %10.8f \n'%(xyzBox.lattice[1,0]/const, xyzBox.lattice[1,1]/const, xyzBox.lattice[1,2]/const))
    inputf.write('%10.8f %10.8f %10.8f \n'%(xyzBox.lattice[2,0]/const, xyzBox.lattice[2,1]/const, xyzBox.lattice[2,2]/const))
    inputf.write('}')
    inputf.write('\n')
    # write natoms
    inputf.write('natom %d \n' %len(xyzBox.xyzAtoms))
    inputf.write('\n')
    # write typat
    inputf.write('typat { \n')
    all_atom_kinds = []
    for atom in xyzBox.xyzAtoms:
        if not atom.name in all_atom_kinds:
            all_atom_kinds.append(atom.name)
        inputf.write('%d '%(all_atom_kinds.index(atom.name)+1))
    #inputf.write('typat { \n')
    #for atom in xyzBox.xyzAtoms:
    #    if atom.name == 'H':
    #        inputf.write('%d ' %4)
    #    if atom.name == 'O':
    #        inputf.write('%d ' %1)
    #    if atom.name == 'N':
    #        inputf.write('%d ' %2)
    #    if atom.name == 'C':
    #        inputf.write('%d ' %3)
    #    if atom.name == 'Cl':
    #        inputf.write('%d ' %5)
    #    if atom.name == 'Na':
    #        inputf.write('%d ' %3)
    inputf.write('\n')
    inputf.write('} \n')
    inputf.write('\n')
    # write ntypat
    inputf.write('ntypat %d \n' %len(all_atom_kinds))
    inputf.write('\n')
    #
    inputf.write('znucl { ' )
    for ii in all_atom_kinds:
        inputf.write(' %d '%xrs_utilities.element(ii))
    inputf.write('} \n')
    inputf.write('\n')
    # write nedges
    #ind = 0
    #for atom in xyzBox.xyzAtoms:
    #    if atom.name == exatom:
    #        ind += 1
    #inputf.write('nedges %d \n' %ind)
    #inputf.write('\n')
    # write edges
    inputf.write('edges { \n')
    #ind = 1
    #for atom in xyzBox.xyzAtoms:
    #    if atom.name == exatom:
    inputf.write('%d %d %d \n' %(-xrs_utilities.element(exatom),edge,subshell))
    #    ind += 1
    inputf.write('} \n')
    inputf.write('\n')
    # write coordinates
    inputf.write('xred { \n')
    for atom in xyzBox.xyzAtoms:
        coords = np.dot(xyzBox.lattice_inv,atom.coordinates)
        inputf.write('%16f %16f %16f \n' % (coords[0], coords[1], coords[2]))
    inputf.write('} \n')
    inputf.write('\n')
    inputf.close()

def writeOCEANinput_full(fname,xyzBox,exatom,edge,subshell):
    """ Writes a complete OCEAN input file.

    Args:
          * fname     (str): Filename for the input file to be written.
          * xyzBox (xyzBox): Instance of the xyzBox class to be converted into an OCEAN input file.
          * exatom    (str): Atomic symbol for the excited atom.
          * edge      (int): Integer defining which shell to excite (e.g. 0 for K-shell, 1 for L, etc.).
          * subshell  (int): Integer defining which sub-shell to excite ( e.g. 0 for s, 1 for p, etc.).

    """
    # some pre-defined input parameters
    std_input = {}
    std_input['control']       = 'control 0 \n'
    std_input['core']          = 'core 28 \n'
    std_input['para_prefix']   = 'para_prefix{ mpirun -machinefile ../nodelist -n 28} \n'
    std_input['dft']           = 'dft { obf } \n'
    std_input['trace_tol']     = 'trace_tol{ 1.0d-10 } \n'
    std_input['core_off']      = 'core_offset 130 \n'
    std_input['ngkpt']         = 'ngkpt{ 1 1 1 } \n'
    std_input['paw.nkpt']      = 'paw.nkpt{ 1 1 1 } \n'
    std_input['nkpt']          = 'nkpt{ 2 2 2 } \n'
    std_input['obkpt']         = 'obkpt{ 1 1 1 } \n'
    std_input['ham_kpoints']   = 'ham_kpoints{ 2 2 2 } \n'
    std_input['nbands']        = 'nbands 1600 \n'
    std_input['paw.nbands']    = 'paw.nbands 1600 \n'
    std_input['obf.nbands']    = 'obf.nbands 1000 \n'
    std_input['rprim']         = 'rprim { \n 1 0 0 \n 0 1 0 \n 0 0 1 \n } \n'
    std_input['ppdir']         = 'ppdir { \n \'/users/sahle/pseudos\' \n } \n'
    std_input['ecut']          = 'ecut 70 \n'
    std_input['toldfe']        = 'toldfe 1.0d-8 \n'
    std_input['toldfr']        = 'tolwfr 1.0d-16 \n'
    std_input['nstep']         = 'nstep 600 \n'
    std_input['mixing']        = 'mixing 0.1 \n'
    std_input['diemac']        = 'diemac 2.0 \n'
    std_input['atomic_paw']    = 'paw.fill {8 o.fill } \n' +'paw.opts {8 o.opts } \n'
    std_input['paw.shells']    = 'paw.shells { 4.5 } \n'
    std_input['cnbse.rad']     = 'cnbse.rad { 4.5 } \n'
    std_input['cnbse.broaden'] = 'cnbse.broaden { 0.05 } \n'
    std_input['scfac']         = 'scfac 1.0 \n'
    std_input['cks.normal']    = 'cks.normal { .true. } \n'
    std_input['xas']           = 'cnbse.mode { xas } \n'
    std_input['cnbse.ways']    = 'cnbse.ways { 3 } \n'
    std_input['cnbse.xmesh']   = 'cnbse.xmesh{ 24 24 24 } \n'

    # open the filename
    f = open(fname,'w')

    # start with printing all standard input parameters
    for key in std_input:
        f.write(std_input[key])

    # write the cell parameters in atomic units
    bl = xyzBox.boxLength/0.52917721092
    f.write( 'acell { ' + str(bl) + ' ' + str(bl) + ' ' + str(bl) + '}\n' )

    # write number of types of atoms
    unique_atoms = []
    for atom in xyzBox.xyzAtoms:
        if atom.name not in unique_atoms:
            unique_atoms.append(atom.name)
    f.write( 'ntypat ' + str(len(unique_atoms)) + ' \n' )

    # write znucl
    znucl_str = 'znucl { ' 
    for ii in range(len(unique_atoms)):
        znucl_str = znucl_str + str(xrs_utilities.element(unique_atoms[ii])) + ' '
    f.write( znucl_str + '}\n' )

    # write pp_list
    f.write( 'pp_list { \n  \n } \n' )

    # write number of atoms in calculation
    f.write( 'natom ' + str(len(xyzBox.xyzAtoms)) + '\n' )

    # write typat
    typat_str = 'typat { '
    list_of_unique_atomic_numbers = [xrs_utilities.element(ii) for ii in unique_atoms]
    for atom in xyzBox.xyzAtoms:
        ind = int(np.where(np.array(list_of_unique_atomic_numbers) == xrs_utilities.element(atom.name))[0] +1)
        typat_str = typat_str + str(ind) + ' '
    f.write( typat_str + '}\n' )

    # write nedges (number of edges to be calculated)
    excited_inds = list(np.where(np.array([atom.name for atom in xyzBox.xyzAtoms]) == exatom)[0])
    f.write( 'nedges ' + str(len(excited_inds)) + ' \n' )

    # write edges
    excited_inds = list(np.where(np.array([atom.name for atom in xyzBox.xyzAtoms]) == exatom)[0])
    edges_str = 'edges { '
    for ind in excited_inds:
        edges_str = edges_str + str(ind+1) +' '+ str(int(edge))+' '+ str(int(subshell))+' '+ '\n'+' '
    f.write( edges_str + '}\n' )

    # write reduced atomic coordinates
    xred_str = 'xred { \n'
    for atom in xyzBox.xyzAtoms:
        xred_str = xred_str + str(atom.x_coord) +' '+ str(atom.y_coord) +' '+ str(atom.z_coord) + ' \n'
    f.write( xred_str + '}\n' )

def parsePwscfFile(fname):
    """ **parsePwscfFile**

    Parses a PWSCF file and returns a xyzBox object.

    Args:
        fname (str): Absolute filename of OCEAN input file.

    Returns:
        xyzBox object


    """
    nat    = 0
    ntyp   = 0
    celldm = 0
    infile = open(fname)
    while True:
        line = infile.readline()
        if not line:
            break
        if '&SYSTEM' in line:
            system_read = True
            while system_read:
                line = infile.readline()
                if 'nat' in line:
                    nat = int(line.split()[-1])
                if 'ntyp' in line:
                    ntyp = int(line.split()[-1])
                if 'celldm' in line:
                    celldm = float(line.split()[-1])
                if nat>0 and ntyp>0 and celldm>0:
                    system_read = False
        if 'CELL_PARAMETERS' in line:
            rprim1 = infile.readline()
            rprim2 = infile.readline()
            rprim3 = infile.readline()
            #lattice = np.array([[rprim1.split()[0], rprim1.split()[1], rprim1.split()[2]], 
            #                    [rprim2.split()[0], rprim2.split()[1], rprim2.split()[2]], 
            #                    [rprim3.split()[0], rprim3.split()[1], rprim3.split()[2]] ] ,dtype=float)
            lattice = np.array([[rprim1.split()[0], rprim2.split()[0], rprim3.split()[0]], 
                                [rprim1.split()[1], rprim2.split()[1], rprim3.split()[1]], 
                                [rprim1.split()[2], rprim2.split()[2], rprim3.split()[2]] ] ,dtype=float)

        if 'ATOMIC_POSITIONS' in line:
            xred     = []
            atomType = [] 
            for ii in range(nat):
                line = infile.readline()
                atomType.append(str(line.split()[0]))
                xred.append(np.array( [line.split()[1], line.split()[2], line.split()[3]  ],dtype=float ))

    lattice[:,0] *= celldm * A2AU_factor
    lattice[:,1] *= celldm * A2AU_factor
    lattice[:,2] *= celldm * A2AU_factor
    xabs  = [np.dot(lattice, xred[ii]) for ii in range(len(xred))]

    xyzAtoms = []
    for ii in range(len(xabs)):
            xyzAtoms.append( xyzAtom(atomType[ii],[xabs[ii][0], xabs[ii][1], xabs[ii][2]], ii)  )

    box = xyzBox( xyzAtoms )

    box.lattice = lattice
    box.lattice_inv = np.linalg.inv(lattice)
    box.relAtoms     = xred

    return box


def parseOCEANinputFile(fname):
    """ **parseOCEANinputFile**

    Parses an OCEAN input file and returns lattice vectors, 
    atom names, and relative atom positions.

    Args:
          * fname (str): Absolute filename of OCEAN input file.
          * atoms (list): List of elemental symbols in the same
            order as they appear in the input file.

    Returns:
        * lattice (np.array): Array of lattice vectors.
        * rel_coords (np.array): Array of relative atomic coordinates.
        * oceaatoms (list): List of atomic names.

    """
    infile = open(fname)
    while True:
        line = infile.readline()
        if not line:
            break
        if 'acell' in line:
            acell = np.array( [line.split()[2], line.split()[3], line.split()[4]] ,dtype=float )
        if 'rprim' in line:
            rprim1 = infile.readline()
            rprim2 = infile.readline()
            rprim3 = infile.readline()
            lattice = np.array([[rprim1.split()[0], rprim1.split()[1], rprim1.split()[2]], 
                                [rprim2.split()[0], rprim2.split()[1], rprim2.split()[2]], 
                                [rprim3.split()[0], rprim3.split()[1], rprim3.split()[2]] ] ,dtype=float)
        if 'xred' in line:
            xred = []
            if len(line.split())==1 or len(line.split())==2:
                line = infile.readline()
            elif len(line.split())== 5:
                xred.append(np.array( [line.split()[-3], line.split()[-2], line.split()[-1]  ],dtype=float ))
                line = infile.readline()
            while True:
                if not '}' in line:
                    xred.append(np.array( [line.split()[0], line.split()[1], line.split()[2]  ],dtype=float ))
                    line =infile.readline()
                elif not line.split()[0]=='}' and line.split()[-1] == '}':
                    xred.append(np.array( [line.split()[0], line.split()[1], line.split()[2]  ],dtype=float ))
                    break
                elif line.split()[0]=='}' and line.split()[-1] == '}':
                    break
                else:
                    print('>>>>>> check your OCEAN input file!')
        if 'ntypat' in line:
            ntypat = line.split()[-1]
        if 'znucl' in line and not line.strip().startswith('#'):
            t1 = [ii for ii in line.split()[1:] if ii!='{']
            znucl = np.array([ii for ii in t1 if ii!='}'], dtype=int)
        if 'typat' in line and 'ntypat' not in line:
            typat = []
            if len(line.split())==1 or len(line.split())==2:
                line = infile.readline()
            elif len(line.split())> 2:
                print(line)
                # line_list = line.split()
                # try:
                #     line_list.remove('typat')
                # except:
                #     pass
                # try:
                #     line_list.remove('{')
                # except:
                #     pass
                # try:
                #     line_list.remove('}')
                # except:
                #     pass
                typat.append(np.array( [line.split()[ii] for ii in range(len(line.split()))  ],dtype=int ))
                #typat = np.asarray(line_list, dtype=int)
                line = infile.readline()
                break
            while True:
                if not '}' in line:
                    typat.append(np.array( [line.split()[ii] for ii in range(len(line.split()))  ],dtype=int ))
                    line =infile.readline()
                elif not line.split()[0]=='}' and line.split()[-1] == '}':
                    typat.append(np.array( [line.split()[ii] for ii in range(len(line.split())-1)  ],dtype=int ))
                    break
                elif line.split()[0]=='}' and line.split()[-1] == '}':
                    break
                else:
                    print('>>>>>> check your OCEAN input file!')
            typat=typat[0]

    lattice[:,0] *= acell[0] * A2AU_factor
    lattice[:,1] *= acell[1] * A2AU_factor
    lattice[:,2] *= acell[2] * A2AU_factor

    atoms_name = [xrs_utilities.element(int(znucl[ii])) for ii in range(len(znucl))]
    atoms = []
    #print('>>>>>>>> ' , typat )
    #print('========', atoms_name)
    for at in typat:
        atoms.append(atoms_name[at-1])

    xabs  = [np.dot(lattice, xred[ii]) for ii in range(len(xred))]

    xyzAtoms = []
    for ii in range(len(xabs)):
            xyzAtoms.append( xyzAtom(atoms[ii],[xabs[ii][0], xabs[ii][1], xabs[ii][2]], ii)  )

    box = xyzBox( xyzAtoms )

    box.lattice = lattice
    box.lattice_inv = np.linalg.inv(lattice)
    box.relAtoms     = xred

    return box







    

def translateOcean2FDMNES_p1(ocean_in, fdmnes_out, header_file):
	# parse the OCEAN input
	box = parseOCEANinputFile( ocean_in )

	# read the header
	header = open(header_file)
	inputf = open(fdmnes_out,'w+')
	for line in header:
		inputf.write(line)

	inputf.write('\n')

	# unit cell lengths
	a = np.linalg.norm(box.lattice[0,:])#/0.52917721092
	b = np.linalg.norm(box.lattice[1,:])#/0.52917721092
	c = np.linalg.norm(box.lattice[2,:])#/0.52917721092

	# unit cell angles
	alpha = xrs_utilities.vangle(box.lattice[1,:], box.lattice[2,:])
	beta  = xrs_utilities.vangle(box.lattice[0,:], box.lattice[2,:])
	gamma = xrs_utilities.vangle(box.lattice[0,:], box.lattice[1,:])

	# write the cell size and angles

	# write the atom numbers and relative coordinates
	inputf.write(' %10.8f %10.8f %10.8f %8.6f %8.6f %8.6f \n'%(a,b,c,alpha,beta,gamma))
	for ii in range(len(box.xyzAtoms)):	
		name = xrs_utilities.element(box.xyzAtoms[ii].name)
		xyz = np.dot(box.lattice_inv,box.xyzAtoms[ii].getCoordinates())
		inputf.write( '%4d %10.8f %10.8f %10.8f \n'%(name, xyz[0], xyz[1], xyz[2] ) ) 

	inputf.write('\n')
	inputf.write('  Convolution \n')
	inputf.write('\n')	
	inputf.write('  End \n')
	inputf.write('\n')

def sorter(elem):
    return elem.dist_to_center

def writeFEFFinput_arb(fname, headerfile, xyzBox, exatom, edge):
    """ **writeFEFFinput_arb**
    """
    # write everything that is in the header
    header = open(headerfile)
    inputf = open(fname,'w+')
    for line in header:
        inputf.write(line)

    inputf.write('\n')

    # find some numbers
    numat = len([atom.name for atom in xyzBox.xyzAtoms])
    numex = len(xyzBox.get_atoms_by_name(exatom))
    tags = list( set([atom.name for atom in xyzBox.xyzAtoms]))
    pots = np.arange(len(tags)) + 1
    exatom_ind = tags.index(exatom) + 1

    # write ATOMS
    inputf.write('ATOMS \n')
    inputf.write(' * x y z ipot tag distance \n')

    cen_atom = xyzBox.xyzAtoms[0]
    dists = []
    for atom in xyzBox.xyzAtoms:
        atom.dist_to_center = np.linalg.norm(atom.coordinates - cen_atom.coordinates)
        dists.append(atom.dist_to_center)

    inds = np.argsort(dists)
    for ind in inds:
        atom = xyzBox.xyzAtoms[ind]
        name   = atom.name
        tag    = tags.index(name)+1
        coords = atom.coordinates
        dist   = atom.dist_to_center
        if ind == 0:
            tag = 0
        inputf.write('%16f %16f %16f %d %s %5f\n' % (coords[0], coords[1], coords[2], tag, name, dist))

    inputf.write('END \n')
    inputf.write('\n')
    inputf.close()
