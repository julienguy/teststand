#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
import desimodel.io
import desispec.io
import sys
import argparse
import string
from desispec.xytraceset import XYTraceSet

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True, nargs="*",
                    help = 'path of psf boot file')
args        = parser.parse_args()

plt.figure("wavelength_range")

for filename in args.psf :
    tset = desispec.io.read_xytraceset(filename)
    wavemin=tset.y_vs_wave_traceset._xmin
    wavemax=tset.y_vs_wave_traceset._xmax
    
    fibers=np.arange(tset.nspec)
    
    ymax = tset.npix_y-1
    wb=[]
    we=[]
    wave=np.linspace(wavemin-10,wavemax+10,int(wavemax-wavemin))
    for fiber in fibers :
        y = tset.y_vs_wave(fiber,wave)
        wb.append(np.interp(0,y,wave))
        we.append(np.interp(ymax,y,wave))

    if wavemin<3700 : arm="b"
    elif wavemax>9000 : arm="z"
    else : arm = "r"
    
    if arm=="b" :
        camid=0
        color="b"
    elif arm=="r" :
        camid=1
        color="g"
    else :
        camid=2
        color="r"

    tfibers = fibers+1+1*(fibers>1)
    
    plt.subplot(2,3,4+camid)
    plt.plot(tfibers,wb,"o",color=color)
    if arm=="b" :
        plt.ylabel("Min. wavelength (A)")
    plt.xlabel("Fiber")
    plt.grid()
    
    plt.subplot(2,3,1+camid,title="%s cam"%arm)
    plt.plot(tfibers,we,"o",color=color)
    if arm=="b" :
        plt.ylabel("Max. wavelength (A)")
    plt.grid()
        
plt.show()

