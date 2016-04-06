#!/usr/bin/env python

import numpy as np
import pyfits
import pylab
from numpy.polynomial.legendre import legval
#from specter.psf.gausshermite import GaussHermitePSF
import sys

def u(wave,wavemin,wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

if len(sys.argv)<2 :
    print sys.argv[0],"psf-?-boot.fits"
    sys.exit(12)

psf=pyfits.open(sys.argv[1])
has_specex_psf=(len(sys.argv)>2)

print psf[0].header
print psf.info()
wavemin=psf[0].header["WAVEMIN"]
wavemax=psf[0].header["WAVEMAX"]
xcoef=psf[0].data
ycoef=psf[1].data
sigma=psf[2].data
print xcoef.shape
nfibers=xcoef.shape[0]

wave=np.linspace(wavemin,wavemax,100)

pylab.figure()
a0=pylab.subplot(1,3,1)
a1=pylab.subplot(1,3,2)
a2=pylab.subplot(1,3,3)
for fiber in range(nfibers) :
    x = legval(u(wave,wavemin,wavemax), xcoef[fiber])
    y = legval(u(wave,wavemin,wavemax), ycoef[fiber])
    a0.plot(x,y)
    a1.plot(y,wave)
a0.set_xlabel("X CCD")
a0.set_ylabel("Y CCD")
a1.set_xlabel("Y CCD")
a1.set_ylabel("Wavelength [A]")
a2.plot(sigma)
a2.set_xlabel("fiber #")
a2.set_ylabel("PSF sigma")


pylab.show()
