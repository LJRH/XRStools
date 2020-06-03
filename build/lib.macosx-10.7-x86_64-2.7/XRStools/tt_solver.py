from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import numpy as np
from . import xrs_utilities
from scipy.integrate import odeint
from six.moves import range

def TTsolver1D(energy, hkl=[6,6,0], crystal='Si', R=1.0, dev=np.arange(-250.0,450.0,1.0), alpha=0.0, chitable_prefix='/home/christoph/sources/XRStools/data/chitables/chitable_'):
    """ **TTsolver**
    Solves the Takagi-Taupin equation for a bent crystal.
    
    This function is based on a Matlab implementation by S. Huotari of M. Krisch's
    Fortran programs.

    Args:
    -----
    energy (float): Fixed nominal (working) energy in keV.
    hkl (array): Reflection order vector, e.g. [6, 6, 0]
    crystal (str): Crystal used (can be silicon 'Si' or 'Ge')
    R (float): Crystal bending radius in m.
    dev (np.array): Deviation parameter (in arc. seconds) for 
                    which the reflectivity curve should be calculated.
    alpha (float): Crystal assymetry angle.

    Returns:
    --------
    refl (np.array): Reflectivity curve.
    e (np.array): Deviation from Bragg angle in meV.
    dev (np.array): Deviation from Bragg angle in microrad.

    """
    # load dielectric susceptibility data (tabulated)
    chi = np.loadtxt(chitable_prefix + crystal.lower() + str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2])) + '.dat')
    if len(chi[:,0]) == 1:
        print( 'Will only  calculate for the following energy: ' + '%.4f' % chi[0,0] + ' keV!!!')
    else:
        if energy < np.min(chi[:,0]) or energy > np.max(chi[:,0]):
            print( 'Energy outside of values defined in Chi-table.')
            return
    # interpolate 
    chi0 = complex(np.interp(energy,chi[:,0],chi[:,1]),np.interp(energy,chi[:,0],chi[:,2]))
    chih = complex(np.interp(energy,chi[:,0],chi[:,3]),np.interp(energy,chi[:,0],chi[:,4]))
    # set the stress tensor values depending on crystal used
    if crystal.upper() == 'SI':
        s13 = -0.278
    elif crystal.upper() == 'GE':
        s13 = -0.273
    else:
        print( 'Poisson ratio for this crystal not defined')
        return
    s15 = -0.0 # s15/s11
    # scattering angle in degree
    th  = xrs_utilities.braggd(hkl,energy,crystal)
    # dspace in m
    dsp = xrs_utilities.dspace(hkl,crystal)/10.0*1e-9
    # wavelength in m
    lam = 12.3984191/energy/10.0*1e-9
    # debye-waller factor
    dwf = 1.0 # dwf = 0.899577 
    # meridional bending radius
    radius = R 
    # sagittal bending radius
    rsag   = R*np.sin(np.radians(th))**2.0 
    # thickness in m
    thick  = 500.0*1e-6 
    # asymmetry in radians
    alpha    = np.radians(alpha)
    # deviation parameter in arcsec
    dev      = dev/3600.0/180.0*np.pi
    # gamma0,gammah = cos<n,K_0,h > , n = inward normal of crystal surface, K_0,h = wave vector
    gammah = -np.sin(np.arcsin(lam/(2.0*dsp)) + alpha) # Krisch et al. convention
    gamma0 =  np.sin(np.arcsin(lam/(2.0*dsp)) - alpha) # Krisch et al. convention
    gamma  = gammah/gamma0
    a0 = np.sqrt(1-gamma0**2.0)
    ah = np.sqrt(1-gammah**2.0)
    beta  = gamma0/np.abs(gammah)
    # polarization factor
    cpol   = 1.0 
    # penetration depth
    mu     = -2.0*np.pi/lam*chi0.imag
    tdepth = 1.0/mu/(1.0/np.abs(gamma0)+1.0/np.abs(gammah))
    lex    = lam*np.sqrt(gamma0*np.abs(gammah))/(np.pi*chih.real)
    y0     = chi0.imag*(1.0+beta)/(2.0*np.sqrt(beta)*chih.real)
    c1     = cpol*dwf* complex(1.0,-chih.imag/chih.real)
    #abbreviation concerning the deviation parameter y
    abb0 = -np.sqrt(beta)/2.0/chih.real
    abb1 = chi0.real*(1.0+beta)/(2.0*np.sqrt(beta)*chih.real)
    #abbreviations concerning the deformation field
    abb2 = gamma0*gammah*(gamma0-gammah)
    abb3 = 1.0 + 1.0/(gamma0*gammah)
    abb4 = s13*(1.0 + radius/rsag)
    abb5 = (ah - a0)/(gamma0 - gammah)*s15 
    abb6 = 1.0/(np.abs(cpol)*chih.real*np.cos(np.arcsin(lam/(2.0*dsp)))*radius)
    abb7 = 2.0*np.abs(cpol)*chih.real*np.cos(np.arcsin(lam/(2.0*dsp)))/gamma0
    # spherical diced crystal, 1-m bending radius, nearly backscattering conditions, strain gradient
    sgbeta = abb6*(abb2*(abb3 - abb4 + abb5))
    # number of steps along reflectivity curve
    nstep=len(dev)
    # reflectivity curve
    refl  = np.zeros_like(dev)

    OLDMETHOD = 0
    if OLDMETHOD:
        # loop over all steps along the reflectivity curve
        for l in range(nstep):
            # deviation parameter
            abb8   = -2.0*np.sin(2.0*np.arcsin(lam/(2.0*dsp)))*dev[l]
            T = np.arange(np.max([-10.0*tdepth, -thick]),0.0,1e-8)
            Y = odeint(odefctn,np.array([0.0, 0.0]),T,args=(abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1)) 
            # normalized reflectivity at this point
            refl[l] = np.sum(Y[-1,:]**2.0)
    else:
        # deviation parameter
        abb8   = -2.0*np.sin(2.0*np.arcsin(lam/(2.0*dsp)))*dev
        # dev axis (complex)
        YY = np.zeros([nstep],"D")
        # small step size
        ministep = tdepth/10000.0
        ssrk     = ministep/2
        xpoints  = np.arange(np.max([-10.0*tdepth, -thick]),0, ministep)
        for xpos in xpoints[:-1]:
            Yp0  = odefctn_CN( YY,              xpos+0*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            Yp1  = odefctn_CN( YY+1.0*ssrk*Yp0, xpos+1*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            Yp2  = odefctn_CN( YY+1.0*ssrk*Yp1, xpos+1*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            Ypb  = odefctn_CN( YY+2.0*ssrk*Yp2, xpos+2*ssrk, abb0, abb1, abb7, abb8, lex, sgbeta, y0, c1)
            YY   =  YY  + ministep*(  Yp0 + 2.0*(Yp1+Yp2) + Ypb  )/6.0
            refl = YY.real*YY.real+YY.imag*YY.imag

    # deviation in degree
    dev = dev/4.848136811e-06/3600.0 
    # deviation in meV
    e   = -(xrs_utilities.energy(xrs_utilities.dspace(hkl,crystal),th+dev)-energy)*1.0e6
    # deviation in microrad
    dev = dev*1.e3
    return refl, e, dev


def odefctn_CN(yCN,t,abb0,abb1,abb7,abb8N,lex,sgbeta,y0,c1):
    fcomp = 1.0/(complex(0,-lex)) * (-2.0*((abb0*(abb8N + abb7*sgbeta*t) + abb1) + complex(0,y0))*(yCN) + c1*(1.0 + yCN* yCN) )
    return fcomp

def odefctn(y,t,abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1):
    """
    #%    [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
    #%    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
    #%    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
    #%    no options are set.   
    """
    #print 'shape of y is ' , np.shape(y), np.shape(t)
    fcomp = 1.0/(complex(0,-lex)) * (-2.0*((abb0*(abb8 + abb7*sgbeta*t) + abb1) + complex(0,y0))*(y[0] + complex(0,y[1])) + c1*(1.0 + (y[0] + complex(0,y[1]))**2.0))
    return fcomp.real,fcomp.imag


def taupgen(e, hkl = [6,6,0], crystals = 'Si', R = 1.0, dev = np.arange(-50.0,150.0,1.0), alpha = 0.0):
    """
    % TAUPGEN          Calculates the reflectivity curves of bent crystals
    %
    % function [refl,e,dev]=taupgen_new(e,hkl,crystals,R,dev,alpha);
    %
    %              e = fixed nominal energy in keV
    %            hkl = reflection order vector, e.g. [1 1 1]
    %       crystals = crystal string, e.g. 'si' or 'ge'
    %              R = bending radius in meters
    %            dev = deviation parameter for which the 
    %                  curve will be calculated (vector) (optional)
    %          alpha = asymmetry angle 
    % based on a FORTRAN program of Michael Krisch
    % Translitterated to Matlab by Simo Huotari 2006, 2007
    % Is far away from being good matlab writing - mostly copy&paste from
    % the fortran routines. Frankly, my dear, I don't give a damn. 
    % Complaints -> /dev/null
    """
    prefix = data_installation_dir+'/'
    path = '/home/christoph/sources/XRStools/data/chitables/chitable_' # prefix + 'data/chitables/chitable_' # path to chitables
    # load the according chitable (tabulated)
    hkl_string = str(int(hkl[0])) + str(int(hkl[1])) + str(int(hkl[2]))
    filestring = path + crystals.lower() + hkl_string + '.dat'
    chi = np.loadtxt(filestring)

    # good for 1 m bent crystals in backscattering
    ystart = -50.0 # start value of angular range in arcsecs
    yend   = 150.0 # end value of angular range in arcsecs
    ystep  = 1.0   # step width in arcsecs

    if len(chi[:,0]) == 1:
        print( ' I will only  calculate for the following energy: ' + '%.4f' % chi[0,0] + ' keV!!!')
    else:
        if e < np.min(chi[:,0]) or e > np.max(chi[:,0]):
            print( 'Energy outside of the range in ' + filestring)
            return

        chi0r = np.interp(e,chi[:,0],chi[:,1])
        chi0i = np.interp(e,chi[:,0],chi[:,2])
        chihr = np.interp(e,chi[:,0],chi[:,3])
        chihi = np.interp(e,chi[:,0],chi[:,4])

    th = braggd(hkl,e,crystals)
    lam = 12.3984191/e/10.0 # wavelength in nm

    reflcorr = 0.0
    chi0 = complex(chi0r,chi0i)
    chih = complex(chihr,chihi)

    if crystals.upper() == 'SI':
        s13 = -0.278
    elif crystals.upper() == 'GE':
        s13 = -0.273
    else:
        print( 'Poisson ratio for this crystal not defined')
        return

    s15 = -0.0 # s15/s11
    dsp = dspace(hkl,crystals)/10.0 # dspace

    dwf    = 1.0 # dwf = 0.899577 # debye-waller factor
    radius = R # meridional bending radius
    rsag   = R*np.sin(np.radians(th))**2.0 # sagittal bending radius
    thick  = 500.0 # thickness in micrometers #rsag = R

    lam      = lam*1e-9
    dsp      = dsp*1e-9
    alpha    = np.radians(alpha) # alpha in rad
    thick    = thick*1e-6
    ystart   = ystart/3600.0/180.0*np.pi
    yend     = yend/3600.0/180.0*np.pi
    ystep    = ystep/3600.0/180*np.pi
    dev      = dev/3600.0/180.0*np.pi
    reflcorr = reflcorr/3600.0/180.0*np.pi

    thetab = np.arcsin(lam/(2.0*dsp))
    cpol   = 1.0 # cpol=0.5*(1+cos(2*thetab).^2) # cpol=cos(2*thetab).^2

    # gamma0 = sin(thetab+alpha) # normal convention
    # gammah = -sin(thetab-alpha) # normal convention
    gammah = -np.sin(thetab + alpha) # Krisch et al. convention (really!)
    gamma0 = np.sin(thetab - alpha) # Krisch et al. convention (I'm not kidding!!)

    beta  = gamma0/np.abs(gammah)
    gamma = gammah/gamma0

    a0 = np.sqrt(1-gamma0**2.0)
    ah = np.sqrt(1-gammah**2.0)

    mu = -2.0*np.pi/lam*chi0i

    tdepth = 1.0/mu/(1.0/np.abs(gamma0)+1.0/np.abs(gammah))

    lex = lam*np.sqrt(gamma0*np.abs(gammah))/(np.pi*chihr)

    y0 = chi0i*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

    pfried = -chihi/chihr

    c1 = cpol*dwf* complex(1.0,pfried)

    #abbreviation concerning the deviation parameter y
    abb0 = -np.sqrt(beta)/2.0/chihr
    abb1 = chi0r*(1.0+beta)/(2.0*np.sqrt(beta)*chihr)

    #abbreviations concerning the deformation field

    abb2 = gamma0*gammah*(gamma0-gammah)
    abb3 = 1.0 + 1.0/(gamma0*gammah)
    abb4 = s13*(1.0 + radius/rsag)
    abb5 = (ah - a0)/(gamma0 - gammah)*s15 
    abb6 = 1.0/(np.abs(cpol)*chihr*np.cos(thetab)*radius)
    abb7 = 2.0*np.abs(cpol)*chihr*np.cos(thetab)/gamma0

    #   a spectrometer based on a spherical diced analyzer crystal with a 1-m bending radius in nearly backscattering conditions utilizing a strain gradient beta
    sgbeta = abb6*(abb2*(abb3 - abb4 + abb5))

    nstep=len(dev)
    eta  = np.zeros_like(dev)
    abb8z = np.zeros_like(dev)
    refl  = np.zeros_like(dev)
    refl1 = np.zeros_like(dev)
    refl2 = np.zeros_like(dev)
    for l in range(nstep):
        # actual value of the deviation angle
        # dev[l] = ystart + (l - 1)*ystep

        # deviation parameter
        abb8   = -2.0*np.sin(2.0*thetab)*dev[l]
        eta[l] = (dev[l]*np.sin(2.0*thetab)+np.abs(chi0.real)/2.0*(1.0-gamma))/(np.abs(cpol)*np.sqrt(np.abs(gamma))*np.sqrt(chih*chih))
        eta[l] = eta[l].real

        ndiff = 2
        xend = 0
        x = np.max([-10.0*tdepth, -thick])
        y = np.array([0.0, 0.0])
        h = xend
        abb8z[l] = abb8

        # in this point call the subroutine
        #     [T,Y] = ODE23(ODEFUN,TSPAN,Y0,OPTIONS,P1,P2,...) passes the additional
        #    parameters P1,P2,... to the ODE function as ODEFUN(T,Y,P1,P2...), and to
        #    all functions specified in OPTIONS. Use OPTIONS = [] as a place holder if
        #    no options are set.   
        #print 'the fucking shape of y is ', np.shape(y)
        T = np.arange(x,xend,1e-8)
        Y = odeint(odefctn,y,T,args=(abb0,abb1,abb7,abb8,lex,sgbeta,y0,c1)) 

        # normalized reflectivity at this point
        refl[l] = np.sum(Y[-1,:]**2.0)
        refl1[l] = Y[-1,0]
        refl2[l] = Y[-1,1]

    de = dev * e * 1.0e6 /np.tan(thetab)

    lam    = lam *1.0e+09        
    dsp    = dsp*1.0e+09        
    alpha  = alpha/np.pi*180.0        
    ystart = ystart/4.848136811e-06
    yend   = yend/4.848136811e-06   
    ystep  = ystep/4.848136811e-06
    dev    = dev/4.848136811e-06 # dev in arcsecs
    
    dev = dev/3600.0 # in degrees
    thb = th
    th  = thb + dev
    e0  = e
    e   = energy(dspace(hkl,crystals),th)-e0
    e = e*1e6

    dev = dev*3600.0 # back to arcsecs

    return refl,e,dev,e0








