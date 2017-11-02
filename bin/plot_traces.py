#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
import sys
import argparse

def u(wave,wavemin,wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, nargs="*", default = None, required = True,
                    help = 'path to psf files')

args = parser.parse_args()

waveref=None
xref=None
yref=None
nw=50
u=np.linspace(-1,1,nw)
for filename in args.psf :
    print(filename)
    psf=pyfits.open(filename)
    wmin=psf[0].header["WAVEMIN"]
    wmax=psf[0].header["WAVEMAX"]
    wave=np.linspace(wmin,wmax,nw)
    xcoef=psf[0].data
    ycoef=psf[1].data
    if waveref is None :
        waveref=wave
    
    fibers=np.arange(xcoef.shape[0])
    
    x=np.zeros((fibers.size,nw))
    y=np.zeros((fibers.size,nw))
    
    wave_is_diff=np.max(np.abs(wave-waveref))>0.05
    
    for fiber in fibers :
        #if fiber%10==0: print(fiber)
        if wave_is_diff :
            x[fiber] = np.interp(waveref,wave,legval(u, xcoef[fiber]))
            y[fiber] = np.interp(waveref,wave,legval(u, ycoef[fiber]))
        else :
            x[fiber] = legval(u, xcoef[fiber])
            y[fiber] = legval(u, ycoef[fiber])
            
    if xref is None :
        xref=x
        yref=y
        plt.figure("traces")
        for fiber in range(yref.shape[0]) :
            plt.plot(xref[fiber],yref[fiber])
    else :
        plt.figure("dx")
        for fiber in range(xref.shape[0]) :
            plt.plot(waveref,x[fiber]-xref[fiber])
        plt.figure("dy")
        for fiber in range(yref.shape[0]) :
            plt.plot(waveref,y[fiber]-yref[fiber])


plt.show()
