#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
from teststand.graph_tools         import parse_fibers
import sys
import argparse
import fitsio

def u(wave,wavemin,wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, nargs="*", default = None, required = True,
                    help = 'path to psf files')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('--image', type=str, default = None, required = False,
                    help = 'overplot traces on image')
parser.add_argument('--lines', type=str, default = None, required = False,
                    help = 'coma separated list of lines')
parser.add_argument('--vmin', type=float, default = None, required = False,
                    help = 'min value for image display')
parser.add_argument('--vmax', type=float, default = None, required = False,
                    help = 'max value for image display')

args = parser.parse_args()

fibers = parse_fibers(args.fibers)


lines=None
if args.lines :
    lines=list()
    for tmp in args.lines.split(",") :
        lines.append(float(tmp))
    print("lines=",lines)

waveref=None
xref=None
yref=None
nw=50
u=np.linspace(-1,1,nw)
nfibers=0
for filename in args.psf :
    print(filename)
    psf=pyfits.open(filename)
    kx = "XTRACE"
    ky = "YTRACE"
    if not kx in psf :
        kx = "XCOEFF"
        ky = "YCOEFF"
        
    wmin=psf[kx].header["WAVEMIN"]
    wmax=psf[kx].header["WAVEMAX"]
    wave=np.linspace(wmin,wmax,nw)
    xcoef=psf[kx].data
    ycoef=psf[ky].data
    if waveref is None :
        waveref=wave
    
    
    if fibers is None :
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
        #print(fiber,"x=",x[fiber])
        #print(fiber,"y=",y[fiber])
        
        
    if xref is None :
        xref=x
        yref=y
        if(len(args.psf)==1) :
            plt.figure("traces")
            
            if args.image is not None :
                img=fitsio.read(args.image)
                vmax=1000
                for l in range(5) :
                    vmax=np.median(img[img>vmax])
                vmin=0
                if args.vmin is not None :
                    vmin=args.vmin
                if args.vmax is not None :
                    vmax=args.vmax
                plt.imshow(img,origin=0,vmin=vmin,vmax=vmax,aspect="auto")

            for fiber in range(yref.shape[0]) :
                
                if lines is not None :
                    for line in lines :
                        xl=np.interp(line,wave,legval(u, xcoef[fiber]))
                        yl=np.interp(line,wave,legval(u, ycoef[fiber]))
                        plt.plot(xl,yl,"x",color="white")

                color=None
                if args.image is not None: color="white"
                plt.plot(xref[fiber],yref[fiber],lw=1,color=color)
    else :
        plt.figure("delta")
        plt.subplot(2,1,1)
        for fiber in range(xref.shape[0]) :
            plt.plot(waveref,x[fiber]-xref[fiber],color=plt.cm.rainbow(fiber/float(len(fibers))))
        plt.xlabel("wavelength")
        plt.ylabel("dx")
        plt.grid()
        plt.subplot(2,1,2)        
        for fiber in range(yref.shape[0]) :
            plt.plot(waveref,y[fiber]-yref[fiber],color=plt.cm.rainbow(fiber/float(len(fibers))))
        plt.xlabel("wavelength")
        plt.ylabel("dy")
        plt.grid()
        
plt.show()
