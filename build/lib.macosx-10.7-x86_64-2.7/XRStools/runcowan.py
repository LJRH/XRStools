from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#!/usr/bin/python
# Filename: runcowan.py

import os
import numpy as np
from scipy import interpolate, signal, integrate, constants, optimize, ndimage
import pylab
from six.moves import range

# the scheme:
# take rcn-file
# run RCN2.sh 
# edit rcf to rcg
# run RCG2.sh
# create rac-file
# run RAC2.sh
# create/take plo-file
# run PLO2.sh
# read .dat with readracah
# square difference, minimize that

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

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)    
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)               
    #print size, sizey    
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc) 

def gauss(x,a):
    """
    returns a gaussian with peak value normalized to unity
    a[0] = peak position
    a[1] = Full Width at Half Maximum
    """
    y = np.exp(-np.log(2.0)*((x-a[0])/a[1]*2.0)**2.0)
    return y

def lorentz(x,x0,fwhm):
    """
    % OUTPUT = LORENTZ(X,X0,FWHM)
    % X      = x-scale (row or column vector)
    % X0     = peak position
    % FWHM   = Full Width at Half Maximum of the Lorentzian

    """
    y = (((x-x0)/(fwhm/2.0))**2.0+1.0)**(-1.0)
    return y

def lorentz2(x,x0,w):
    y = 0.5*w/((x-x0)**2.0+(0.5*w)**2.0)/np.pi
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
    n  = int( np.floor(np.max(np.shape(xg))//2))
    c  = c[n:len(c)-n+1] # not sure about the +- 1 here
    f  = interpolate.interp1d(x2,c)
    return f(x)

def readracah(fname,degauss=0.5,delorentz=0.4):
    lines = open(fname,'r').readlines()
    foundsticks = []
    A = []
    B = []
    counter = 0
    for line in lines:
        counter += 1
        if 'Sticks' in line:
            break
    for line in lines[:counter-2]:
        A.append([float(line.strip().split()[0]), float(line.strip().split()[1])])    
    for line in lines[counter+1:]:
        B.append([float(line.strip().split()[0]), float(line.strip().split()[1])])
    espectr = np.array([col[0] for col in A])
    yspectr = np.array([col[1] for col in A])
    estick  = np.array([col[0] for col in B])
    ystick  = np.array([col[1] for col in B])
    # e       = espectr
    e = np.arange( np.min(estick)-10.0, np.max(estick)+10,np.mean(np.diff(espectr)))
    y       = np.zeros_like(e)
    for ii in range(len(estick)):
        y += lorentz(e,estick[ii],delorentz)*ystick[ii]
    y = convg(e,y,degauss)
    return e,y,estick,ystick,espectr,yspectr

def dqtox400(tendq,dt=0,ds=0):
    dq   = tendq/10
    x400 = np.sqrt(30.0)*(6.0*dq-7.0/2.0*dt)
    x420 = -5.0/2.0*np.sqrt(42.0)*dt
    x220 = -np.sqrt(70.0)*ds
    return x400,x420,x220

def x400todq(x400,x420=0,x220=0):
    dq    = x400/np.sqrt(30.0)/6.0-7.0/30.0*x420/np.sqrt(42.0)
    tendq = 10.0*dq
    ds    = -x220/np.sqrt(70.0)
    dt    = -2.0/5.0*x420/np.sqrt(42.0)
    return tendq,dt,ds

def createracinput(fname1,fname2,tendq):
    x400 = dqtox400(tendq)[0]
    fid1 = open(fname1,'r').readlines()
    fid2 = open(fname2,'w+')
    for line in fid1:
        if not 'BRANCH 4+' in line: 
            fid2.write(line)
        else:
            string = '    BRANCH 4+ > 0 0+ %.2f\n' % x400
            fid2.write(string)
    fid2.close()

def writercginput_r1(fnamef,fnameg,scscale):
    fidf = open(fnamef,'r').readlines()
    fidg = open(fnameg,'w+')
    string = '   10    1    0   14    2    4    1    1 SHELL03000000 SPIN03000000 INTER8      \n'
    fidg.write(string)
    string = '    0                         %2d99%2d%2d            8065.47800     0000000        \n' % (scscale,scscale,scscale)
    fidg.write(string)
    for line in fidf[2:]:
        fidg.write(line)
    fidg.close()

def writercginput_r2(fnamef,fnameg,scscale):
    fidf = open(fnamef,'r').readlines()
    fidg = open(fnameg,'w+')
    string = '   10    1    0   14    2    6    2    2 SHELL03000000 SPIN03000000 INTER8      \n'
    fidg.write(string)
    string = '    0                         %2d99%2d%2d            8065.47800     0000000        \n' % (scscale,scscale,scscale)
    fidg.write(string)
    string = '    1     2 1 12 1 10         00      9 00000000 3 8065.4790 .00       1\n'
    fidg.write(string)
    for line in fidf[3:]:
        fidg.write(line)
    fidg.close()

def writercginput_r3(fnamef,fnameg,scscale):
    fidf = open(fnamef,'r').readlines()
    fidg = open(fnameg,'w+')
    string = '   10    1    0   14    2    4    3    3 SHELL03000000 SPIN03000000 INTER8      \n'
    fidg.write(string)
    string = '    0                         %2d99%2d%2d            8065.47800     0000000        \n' % (scscale,scscale,scscale)
    fidg.write(string)
    for line in fidf[2:]:
        if '//R1//' in line:
            s = 'Fe2+ 3P06 3D06      Fe2+ 3p05 3D07         1.21660( 3P//R3// 3D)-0.991HR -91 -96\n'
            fidg.write(s)
        else:
            fidg.write(line)
    fidg.close()

def for_fitfe2_r1(a,e,y):
    scscale = a[1]#int(np.around(80)) # a[1]
    tendq   = a[0]
    eshift  = 1.5
    if scscale<0.0 or scscale>99.0 or tendq<-10.0 or tendq>7.0:
        d = 1000000
    os.system('../batch/RCN2.sh fe3_r1')
    writercginput_r1('fe3_r1.rcf','fe3_r1.rcg',scscale)
    os.system('../batch/RCG2.sh fe3_r1')    
    createracinput('rac_Oh_template.rac','fe3_r1.rac',tendq)
    os.system('../batch/RAC2.sh fe3_r1')
    os.system('../batch/PLO2.sh fe3_r1')
    et,yt,estick,ystick,espectr,yspectr = readracah('fe3_r1.dat',degauss=1.5,delorentz=0.4)
    espectr = espectr-eshift
    f  = interpolate.interp1d(espectr, yspectr,bounds_error=False,fill_value=0.0)
    yspectr = f(e)
    yspectr = yspectr/np.trapz(yspectr,e)*np.trapz(y,e)
    pylab.plot(e,yspectr,e,y)
    pylab.show(block=False)
    #return np.sum((yt-y)**2.0)
    return yspectr


def make_10dq_series(tendqs,scscale,degauss=1.0,delorentz=0.4):
    """
    calculates spectra for a series of 10dq values
    writes a file with columns: energy, 10dq1, 10dq2, ... 
    """
    sticks   = {}
    spectra  = {}
    spectra2 = {}
    # run the first tendq to get the output shape etc.
    for tendq in tendqs:
        os.system('../batch/RCN2.sh fe3_r2')
        writercginput_r2('fe3_r2.rcf','fe3_r2.rcg',scscale)
        os.system('../batch/RCG2.sh fe3_r2')    
        createracinput('rac_Oh_template_r2.rac','fe3_r2.rac',tendq)
        os.system('../batch/RAC2.sh fe3_r2')
        os.system('../batch/PLO2.sh fe3_r2')
        et,yt,estick,ystick,espectr,yspectr = readracah('fe3_r2.dat',degauss=degauss,delorentz=delorentz)
        sticks[str(tendq)]   = estick,ystick
        spectra[str(tendq)]  = et,yt
        spectra2[str(tendq)] = espectr,yspectr
    return sticks, spectra, spectra2

def write_dict_to_txt(dictionary,prefix,postfix):
    for key in dictionary:
        fname = prefix + key + postfix
        data  = np.zeros((len(dictionary[key][0]),len(dictionary[key])))
        for ii in range(len(dictionary[key])):        
            data[:,ii] = dictionary[key][ii]
        np.savetxt(fname,data)

# make calculations for Christopher Weiss
tendqs  = np.arange(1.0,5.1,0.1)
scscale = 80.0 
sticks, spectra, spectra2 = make_10dq_series(tendqs,scscale,degauss=1.0,delorentz=0.4)

from pylab import *
ion()
thelegend = []
for key in spectra:
    thelegend.append(key)
    plot(spectra[key][0],spectra[key][1])

legend(thelegend)


#exp = np.loadtxt('/home/christoph/data/fe_data_alex/2/fe3oct_lq.dat')
#e = exp[:,0]
#y = exp[:,1]

#fitfunc = lambda p, x, y: for_fitfe2_r1(p,x,y) # Target function
#errfunc = lambda p, x, y: fitfunc(p, x, y) - y # Distance to the target function
#p0 = [2.1,70] # Initial guess for the parameters
#p1, success = optimize.leastsq(errfunc, p0[:], args=(e, y))

#a = for_fitfe2_r1(p1,e,y)

#p1 = [2.5,70]

#pylab.plot(e,y,e,for_fitfe2_r1(p1,e,y))
#pylab.show()

#print p1, success

#optimize.leastsq(for_fitfe2_r1, [1.3, 2.1], args=(e,y))


# work in xrstools/cowans/WORK directory

# take rcn-file
# run RCN2.sh 
#os.system('../batch/RCN2.sh fe3_r1')
# edit rcf to rcg
#writercginput_r1('fe3_r1.rcf','fe3_r1.rcg',80)
# run RCG2.sh
#os.system('../batch/RCG2.sh fe3_r1')
# create rac-file
#createracinput('rac_Oh_template.rac','fe3_r1.rac',1.3)
# run RAC2.sh
#os.system('../batch/RAC2.sh fe3_r1')
# create/take plo-file
# run PLO2.sh
#os.system('../batch/PLO2.sh fe3_r1')
# read .dat with readracah
#e,y,estick,ystick,espectr,yspectr = readracah('fe3_r1.dat',degauss=1.5,delorentz=0.4)
# square difference, minimize that
#pylab.plot(e,y,espectr,yspectr)
#pylab.show()

#function [e,y,estick,ystick,espectr,yspectr]=readracah(fname,degauss,delorentz,delorentz_split);
#%function [e,y,estick,ystick,espectr,yspectr]=readracah(fname,degauss,delorentz,delonrentz_split);
#% degauss = fwhm of gaussian broadening (eV); delorentz=fwhm of lorentzian broadening
#
#try, degauss; catch, degauss=0.5; end
#try, delorentz; catch, delorentz=0.4; end
#try, delorentz_split; catch; delorentz_split=[]; end
#if length(delorentz_split)~=length(delorentz)-1, error('Need also lorentz conv split'); end
#
#fid=fopen(fname,'r');
#foundsticks=[];A=[];B=[];
#while ~feof(fid) & length(foundsticks)==0,
#  s=fgetl(fid);
#  foundsticks=strfind(s,'Sticks');
#  if length(foundsticks)==0,
#      A=[A;str2num(s)];
#  end
#end
#while ~feof(fid);
#    s=fgetl(fid);
#    B=[B;str2num(s)];
#end
#fclose(fid);
#espectr=A(:,1);yspectr=A(:,2);
#estick=B(:,1); ystick=B(:,2);
#%e=espectr;
#e=[min(estick)-10:mean(diff(espectr)):max(estick)+10]';
#delor=ones(size(estick))*delorentz(1);
#y=zeros(size(e));
#for ii=1:length(delorentz_split)
#    delor(find(estick>delorentz_split(ii)))=delorentz(ii+1);
#%    ystick(find(estick>delorentz_split(ii)))=ystick(find(estick>delorentz_split(ii)))*1.3;
#end
#
#for ii=1:length(estick);
#   y=y+lorentz2(e,estick(ii),delor(ii))*ystick(ii);
#end
#y=convg(e,y,degauss);


#function d=fitti4(a,e,y);
#%scscale=a(1);
#scscale=70;
#tendq=a(1);
#eshift=a(2);
#scscale=round(scscale);
#if (scscale<1 | scscale>99 | tendq<1 | tendq>3), d=1000000;
#else
#!c:\cowan\batch\rcn2 ti4
#writercginp;
#!c:\cowan\batch\rcg2 ti4
#createracinput('ti4orig.rac','ti4.rac',tendq);
#!c:\cowan\batch\rac2 ti4
#!c:\cowan\batch\plo2 ti4
#[et,yt]=readracah('ti4_r1.dat',1.0,[0.2 0.9]*2,[466.7]);
#et=et-eshift;
#yt=cnan(interp1(et,yt,e));
#y=y/isum(e,y)*isum(e,yt);
#plot(e,y,e,yt);drawnow
#d=sum((yt-y).^2);
#end








