#!/usr/bin/env python

import sys,os
import numpy as np

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt



for filename in sys.argv[1:] :
    h=pyfits.open(filename)
    wave=h["WAVELENGTH"].data
    cal=h["FLUXCALIB"].data[100] # electrons/A  /  (1e-17  ergs/s/cm2/A)
    cal *= 1e17 # electrons/A  /  ( ergs/s/cm2/A)
    exptime=h[0].header["EXPTIME"]
    cal /= exptime # electrons  /  (ergs/cm2)
    # reference effective collection area
    area = 8.678709421*1e4 # cm2
    cal /= area # electrons  /  ergs
    hplanck = 6.62606957e-34 #J.s
    hplanck *= 1e7 # erg.s
    cspeed = 2.99792458e8 # m/s
    cspeed *= 1e10 # A/s
    energy = hplanck*cspeed/wave # erg ( erg.s * A/s / A)
    cal *= energy # electrons /photons
    plt.plot(wave,cal,label=filename)

#plt.legend()

plt.show()


