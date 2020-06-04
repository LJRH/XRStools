from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from six.moves import range
#!/usr/bin/python
# Filename: xrs_utilities.py

#/*##########################################################################
#
# The XRStools software package for XRS spectroscopy
#
# Copyright (c) 2013-2014 European Synchrotron Radiation Facility
#
# This file is part of the XRStools XRS spectroscopy package developed at
# the ESRF by the DEC and Software group and contains practical functions, 
# most of which are translated from Matlab functions from the University of
# Helsinki Electronic Structure Laboratory.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
__author__ = "Christoph J. Sahle - ESRF"
__contact__ = "christoph.sahle@esrf.fr"
__license__ = "MIT"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"

import os
import numpy as np
import matplotlib.pyplot as plt
from . import xrs_utilities

data_installation_dir = os.path.join( os.path.dirname(os.path.abspath(__file__)), 'resources',  'data')

# tabulated data / form factors
Si_bragg_data = np.loadtxt(os.path.join(data_installation_dir, 'si_bragg_params.dat'  ))
Si_f2         = np.loadtxt(os.path.join(data_installation_dir, 'si_f2.dat'  ))
Ge_bragg_data = np.loadtxt(os.path.join(data_installation_dir, 'ge_bragg_params.dat'  ))
Ge_f2         = np.loadtxt(os.path.join(data_installation_dir, 'ge_f2.dat'  ))

# constants
c = 299792458.0             # m/sec
h = 4.13566727e-15          # eV/s
e_mass = 9.10938291e-31     # kg
e_charge = 1.60217657e-19   # C
epsilon_0 = 8.854187817e-12 #F/m
r_electron = 1.0/(4.0*np.pi*epsilon_0)*e_charge**2/e_mass/c**2
Si_latt = 5.43095           # Angstr.
Ge_latt = 5.65734992136     # Angstr.

def get_crystal_data(crystal):
    """ **get_crystal_data**
    Returns mass density, atomic mass number, atomic density, and
    poisson ratio for given crystal ('Si' and 'Ge' so far only).

    Args:
    -----
        crystal (str): Keyword for crystal used (so far only 'Si' and 'Ge').

    Returns:
    --------
        mass_density, mass_number, atom_density, poisson_ratio (floats):
            mass density, atomic mass number, atomic density, and
            poisson ratio
    """
    # Si params
    mass_density_si  = 2329000 # g/m**3
    mass_number_si   = 28
    atom_density_si  = mass_density_si*6.022e23/mass_number_si
    poisson_ratio_si = 0.22

    # Ge params
    mass_density_ge  = 5323000 # g/m**3
    mass_number_ge   = 73
    atom_density_ge  = mass_density_ge*6.022e23/mass_number_ge
    poisson_ratio_ge = 0.27
    if crystal=='Si':
        return mass_density_si, mass_number_si, atom_density_si, poisson_ratio_si
    elif crystal=='Ge':
        return mass_density_ge, mass_number_ge, atom_density_ge, poisson_ratio_ge
    else:
        print('No parameters for %s crystals, yet.'%crystal)
        return

def get_bragg_data(crystal):
    """ **get_bragg_data**
    Returns an array of hkl-indices and delE/E_p data.

    Args:
    -----
        crystal (str): Keyword for crystal used (so far only 'Si' and 'Ge').

    Returns:
    --------
        array (nd.array): hkl-indices and delE/E_p data.

    """
    if crystal == 'Si':
        return Si_bragg_data
    elif crystal == 'Ge':
        return Ge_bragg_data
    else:
        print('No parameters for %s crystals, yet.'%crystal)
        return

def get_lattice_constant(crystal):
    """ **get_lattice_constant**
    Returns an array of hkl-indices and delE/E_p data.

    Args:
    -----
        crystal (str): Keyword for crystal used (so far only 'Si' and 'Ge').

    Returns:
    --------
        latt (float): Lattice constant of given crystal.

    """
    latt = {'Si': 5.43095,
            'Ge': 5.65734992136,
            'SIXOP': 5.430919,
            'SIKOH': 5.430707,
            'LIF': 4.027,
            'INSB': 6.4784,
            'C': 6.708,
            'DIA': 3.57,
            'LI': 3.41,
            'TCNE': 9.736,
            'CU': 3.61,
            'PB': 4.95,
            'NA': 4.2906,
            'AL': 4.0495 }
    try:
        return latt[crystal]
    except:
        print('No parameters for %s crystals, yet.'%crystal)
        return

def get_f2(crystal, energy):
    """ **get_f2**
    Returns the anomalous form factor f" for given element and energy.

    Args:
    -----
        crystal (str): String for element (so far only 'Ge', 'Si').
        energy (float,array): Energy in [keV].

    Returns:
    --------
        f2 (float, array): form factor f2.
    """
    if crystal == 'Si':
        return np.interp(energy*1e3,Si_f2[:,0], Si_f2[:,1])
    elif crystal == 'Ge':
        return np.interp(energy*1e3,Ge_f2[:,0], Ge_f2[:,1])
    else:
        print('No parameters for %s crystals, yet.'%crystal)
        return

def find_column( crystal, hkl ):
    """ **find_column**
    Returns the row-index for a given crystal and hkl for the
    Bragg data file.

    Args:
    -----
        crystal (str): Keyword describing the crystal (so far only 'Ge' and 'Si').
        hkl (array): HKL-indices for wanted lattice plane.

    Returns:
    --------
        ind (int): Row-index.

    """
    data = get_bragg_data(crystal)
    return np.where( np.logical_and(np.logical_and( data[:,0]==hkl[0], data[:,1]==hkl[1]), data[:,2]==hkl[2] ) )[0]
        
def get_delE_over_Ep( crystal, hkl ):
    data = get_bragg_data(crystal)
    ind  = find_column( crystal, hkl )
    return data[ind,4]

def get_thetaB( energy, crystal, hkl ):
    latt = get_lattice_constant(crystal)
    lam = h*c/energy*1e7
    sqrt_hkl = np.sqrt( hkl[0]**2 + hkl[1]**2 + hkl[2]**2 )
    thetaB = 180./np.pi * np.arcsin((lam * sqrt_hkl)/ (2.0* latt))
    return thetaB
    
def get_ax(rowlandR, thetaB, asymmetry):
    return 2.0 * rowlandR * np.sin((thetaB-asymmetry)*np.pi/180.0 )

def get_delTh_darwin_microrad( energy, crystal, hkl ):
    delEOverEp = get_delE_over_Ep( crystal, hkl )
    thetaB = get_thetaB( energy, crystal, hkl )
    return delEOverEp*np.tan(thetaB*np.pi/180.0)*1.0e6

def get_delTh_darwin_meV( energy, crystal, hkl ):
    delEOverEp = get_delE_over_Ep( crystal, hkl )
    return delEOverEp*energy*1e6

def get_source_contrib_x(energy, source_hx, crystal, hkl, rowlandR, asymmetry):
    thetaB = get_thetaB(energy, crystal, hkl)
    ax = get_ax(rowlandR, thetaB, asymmetry)
    return energy*source_hx/1000.0/ax/np.tan(thetaB*np.pi/180.0)**2*1000000 # in meV

def get_source_contrib_y(energy, source_hy, crystal, hkl, rowlandR, asymmetry):
    thetaB = get_thetaB(energy, crystal, hkl)
    ax = get_ax(rowlandR, thetaB, asymmetry)
    return energy/8.0*(source_hy/1000/ax)**2*1000000

def get_source_contrib_z(energy, source_v, crystal, hkl, rowlandR, asymmetry):
    thetaB = get_thetaB(energy, crystal, hkl)
    ax = get_ax(rowlandR, thetaB, asymmetry)
    return source_v/1000/ax/np.tan(thetaB*np.pi/180.0)*energy*1000000

def get_pixel_contrib(energy, pixels, crystal, hkl, rowlandR, asymmetry):
    thetaB = get_thetaB(energy, crystal, hkl)
    ax = get_ax(rowlandR, thetaB, asymmetry)
    return pixels/1000/(2.0*rowlandR*np.sin((thetaB-asymmetry)*np.pi/180)+2.0*rowlandR*np.sin((thetaB+asymmetry)*np.pi/180))/np.tan(thetaB*np.pi/180)*energy*1000000

def get_Johann_contrib(energy, maskD, crystal, hkl, rowlandR, asymmetry):
    thetaB = get_thetaB(energy, crystal, hkl)
    ax = get_ax(rowlandR, thetaB, asymmetry)
    return energy*0.5*(maskD/2.0/ax)**2/np.tan(thetaB*np.pi/180.0)/np.tan((thetaB-asymmetry)*np.pi/180.0)*1000000

def get_offRowland_contrib(energy, z ,maskD, crystal, hkl, rowlandR):
    thetaB = get_thetaB(energy, crystal, hkl)
    return energy*z*maskD/(2*rowlandR*np.sin(np.radians(thetaB))+z)/(2*rowlandR*np.sin(np.radians(thetaB)))/2.0/np.tan(np.radians(thetaB))*1000000

def get_absorption_length(crystal, energy):
    mass_density, mass_number, atom_density, poisson_ratio = get_crystal_data(crystal)    
    lam = h*c/energy*1e7
    f2  = get_f2_data(crystal, energy)
    return (2.0*r_electron * atom_density*lam*0.0000000001*f2)**(-1)*1000000

def get_stress_contrib(energy, rowlandR, crystal, hkl):
    mass_density, mass_number, atom_density, poisson_ratio = get_crystal_data(crystal)
    a_length = get_absorption_length(crystal, energy)
    thetaB = get_thetaB(energy, crystal, hkl)
    f2  = get_f2_data(crystal, energy)
    return energy*a_length*np.sin(np.radians(thetaB))/1000.0/(2.0*rowlandR)*np.abs(1.0/np.tan(np.radians(thetaB))**2-2*poisson_ratio)*1000000

def get_diced_resolution(energy, crystal, hkl, source_hx, source_hy, source_z, pixel_s, rowlandR, maskD, asymmetry):
    delE_darwin = get_delTh_darwin_meV( energy, crystal, hkl )
    delE_sx     = get_source_contrib_x(energy, source_hx, crystal, hkl, rowlandR, asymmetry)
    delE_sy     = get_source_contrib_y(energy, source_hy, crystal, hkl, rowlandR, asymmetry)
    delE_sz     = get_source_contrib_z(energy, source_z, crystal, hkl, rowlandR, asymmetry)
    delE_pixel  = get_pixel_contrib(energy, pixel_s, crystal, hkl, rowlandR, asymmetry)
    delE_Johann = get_Johann_contrib(energy, maskD, crystal, hkl, rowlandR, asymmetry)
    return np.sqrt( delE_darwin**2 + delE_sx**2 + delE_sy**2 + delE_sz**2 + delE_pixel**2 + delE_Johann**2  )

def get_asymmetry_factor(asymmetry, energy, crystal, hkl):
    thetaB = get_thetaB(energy, crystal, hkl)
    return np.sin(np.radians(thetaB + asymmetry))/( np.sin( np.radians(thetaB - asymmetry)) )

def get_resolution_4bounce(energy, crystal, hkl, asymmetry):
    delTh = get_delTh_darwin_microrad( energy, crystal, hkl )
    asym  = get_asymmetry_factor(asymmetry, energy, crystal, hkl)
    thetaB = get_thetaB(energy, crystal, hkl)
    return delTh/np.sqrt(asym)/asym*energy/np.tan(np.radians(thetaB))

def get_sigma_v(energy, crystal, hkl, asymmetry, source_dist):
    # additional contribution to the divergence of the X-ray beam
    # at distance souce_dist from the mono, due to the asymmetry
    delTh = get_delTh_darwin_microrad( energy, crystal, hkl )
    asym  = get_asymmetry_factor(asymmetry, energy, crystal, hkl)
    return delTh*(asym**2 -1.0)/np.sqrt(asym**3)*source_dist 

def get_deltaTh_max(energy, crystal, hkl, asymmetry):
    # 4-bounce
    delTh = get_delTh_darwin_microrad( energy, crystal, hkl )
    asym  = get_asymmetry_factor(asymmetry, energy, crystal, hkl)
    return delTh*(asym**2-1.0)/np.sqrt(asym)

def find_angle(hkl1, hkl2, verbose=True):
    """ **find_angle**
    Find angle between two reflections HKL1 and HKL2.

    Args:
    -----
    hkl1 (list, array): First reflection.
    hkl2 (list, array): Second reflection.

    Returns:
    angle (float): Angle between two reflection in degrees.

    """
    angle = xrs_utilities.vangle(hkl1, hkl2)
    if verbose:
        print('The angle between ',hkl1,' and ',hkl2,' is: %f'%angle, 'degrees.')
    return angle
        


class analyzer:
    """ Class for estimating analyzer resolutions.
    """
    def __init__(self, hkl, energy, crystal='Si',asymmetry=0.0, maskD=100.0, rowlandR=500 ):

        # constants
        self.c = 299792458.0         # m/sec
        self.h = 4.13566727e-15      # eV/s
        self.Si_latt = 5.43095       # Angstr.
        self.Ge_latt = 5.65734992136 # Angstr.

        # params

    def find_deltaE_over_Ep(crystal, hkl):
        pass
        
    def find_angle(self, hkl1, hkl2, verbose=True):
        """ **find_angle**
        Find angle between two reflections HKL1 and HKL2.

        Args:
        -----
        hkl1 (list, array): First reflection.
        hkl2 (list, array): Second reflection.

        Returns:
        angle (float): Angle between two reflection in degrees.

        """
        angle = find_angle(hkl1, hkl2, verbose=verbose)
        return angle
        


