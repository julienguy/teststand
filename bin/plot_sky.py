#!/usr/bin/env python

import sys,os
import numpy as np

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt



for filename in sys.argv[1:] :
    h=pyfits.open(filename)
    wave=h["WAVELENGTH"].data
    sky=h["SKY"].data
    ivar=h["IVAR"].data

    fiber=100
    plt.figure("sky")
    #for fiber in range(sky.shape[0]) :
    plt.plot(wave,sky[fiber])
    plt.figure("ivar")
    #for fiber in range(sky.shape[0]) :
    plt.plot(wave,ivar[fiber])
#plt.legend()

plt.show()


