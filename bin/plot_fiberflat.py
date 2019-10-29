#!/usr/bin/env python

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--infile', type = str, default = None, required = True,
                    help = 'path to fiber flat file')

args        = parser.parse_args()
h=fits.open(args.infile)


h.info()
wave=h["WAVELENGTH"].data

plt.figure()
for fiber in range(h[0].data.shape[0]) :
    ok=np.where(h["MASK"].data[fiber]==0)[0]
    plt.plot(wave,h[0].data[fiber],color="gray",alpha=0.1)
    if ok.size>1 :
      plt.plot(wave[ok],h[0].data[fiber,ok])

if "FIBERMAP" in h :
    plt.figure()
    mflat=np.median(h[0].data,axis=1)
    x=h["FIBERMAP"].data["FIBERASSIGN_X"]
    y=h["FIBERMAP"].data["FIBERASSIGN_Y"]
    plt.scatter(x,y,c=mflat)
    plt.colorbar()

if False :
    plt.figure()

    #meanspec=np.mean(h[0].data,axis=0)
    meanspec=np.ones((h[0].data.shape[1]))
    ratio=h[0].data/meanspec
    ivar=h[1].data*meanspec**2

    wrebin=100
    rwave=wave[:(wave.size//wrebin)*wrebin].reshape((wave.size//wrebin),wrebin).mean(-1)

    fbin=50
    nfbin=500//fbin
    for f in range(nfbin) :
        color=plt.cm.rainbow(f/float(nfbin))
        mratio=np.mean(ratio[fbin*f:fbin*(f+1)],axis=0)
        mivar=np.sum(ivar[fbin*f:fbin*(f+1)],axis=0)
        rmratio=mratio[:(wave.size//wrebin)*wrebin].reshape((wave.size//wrebin),wrebin).mean(-1)
        rmivar=mivar[:(wave.size//wrebin)*wrebin].reshape((wave.size//wrebin),wrebin).sum(-1)
        err=1./np.sqrt(rmivar)
        plt.fill_between(rwave,rmratio-1-err,rmratio-1+err,color=color,alpha=0.2)
        plt.plot(rwave,rmratio-1,c=color,lw=2)
    plt.xlabel("wavelength (A)")
    plt.ylabel("fiber flat syst. error")
plt.grid()    
plt.show()
