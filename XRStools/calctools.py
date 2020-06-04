from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: stobe_analyze.py

import os

import numpy as np
import array as arr

from itertools import groupby
from scipy.integrate import trapz
from scipy.interpolate import interp1d

from pylab import *
from scipy import signal
from scipy.ndimage import measurements
import matplotlib.pyplot as plt
from six.moves import range

__metaclass__ = type # new style classes

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
        s2 = s2 + sticks[n]*gauss(e2,evals[n],fwhm[n])
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
    class to analyze ERKALE results
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


#### calculations
# cut clusters from xyz files
# cut clusters from gro files

################################### 
# reading function for cowan's code output

def read_rcg_born(fname,q=9.2*0.529,egauss=0.2,elore=0.1):
    #% function [e,spectr]=read_rcg_born(fname,q);
    #% fname = filename (outg11 with plane wave born collision calculation)
    #% q     = momentum transfer in at.units
    #% egauss = gaussian component of the resolution
    #% elore  = lorentzian component of the resolution
    #%   S. Huotari 2012


    return e,spectr

# function [e,spectr]=read_rcg_born(fname,q,egauss,elore);


# fid=fopen(fname);
# ei=[]; y=[];

# conti=1;
# while conti &   ~feof(fid),
#   s=fgetl(fid);
#   conti=isempty(findstr('k-values',s));
# end
# conti=1;
# while conti & ~feof(fid),
#   s=fgetl(fid);
#   conti=isempty(findstr('k-values',s));
# end
#       s=fgetl(fid); % empty string
#       conti=1;k=[];
#       while conti,
#           s=fgetl(fid);
#           conti=length(str2num(s)); 
#           k=[k str2num(s)];
#       end

# while ~feof(fid)
#   s=fgetl(fid);
#   if ~isempty(findstr('delta E',s));
#       s=fgetl(fid); % empty line
#       s=fgetl(fid); disp(s);
#       sn=str2num(s);
#       ei=[ei sn(8)];
#       s=fgetl(fid); % the text "0generalised.."
#       conti=1;ytmp=[];
#       while conti,
#           s=fgetl(fid);
#           conti=length(str2num(s)); 
#           ytmp=[ytmp str2num(s)];
#       end      
#       y=[y;ytmp];
#   end
# end
# fclose(fid);
# ei=ei/8.065-1.45;
# %q=9.2*0.529; 
# [q,qi]=min((k-q).^2)
# e=[min(ei)-2:0.01:max(ei)+2]';
# spectr=zeros(size(e));
# for ii=1:length(ei),
#   spectr=spectr+convg(e,lorentz2([y(ii,qi) elore ei(ii)],e),egauss);
# end









