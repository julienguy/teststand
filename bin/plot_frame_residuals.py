#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
from desispec.log import get_logger
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--frame', type = str, default = None, required = True,
                    help = 'path of frame fits file')

log         = get_logger()
args        = parser.parse_args()
frame_file  = pyfits.open(args.frame)
flux        = frame_file[0].data
ivar        = frame_file[1].data
wave        = frame_file["WAVELENGTH"].data


for amp in range(2) :
    if amp==0 :
        fibers=np.arange(0,10)
    else :
        fibers=np.arange(10,20)

    sw          = np.sum(ivar[fibers],axis=0)
    mflux       = np.sum(ivar[fibers]*flux[fibers],axis=0)/(sw+(sw==0))
    chi2        = np.sum(ivar[fibers]*(flux[fibers]-mflux)**2,axis=0)
    rebin       = 8
    chi2        = chi2[:(chi2.size/rebin)*rebin].reshape(chi2.size/rebin,rebin).sum(-1)
    chi2pdf     = chi2/(fibers.size-1)/rebin
    rwave       = wave[:(wave.size/rebin)*rebin].reshape(wave.size/rebin,rebin).mean(-1)
    plt.plot(rwave,chi2pdf)

plt.xlabel("wavelength")
plt.ylabel("chi2/ndf")
plt.grid()
plt.show()

