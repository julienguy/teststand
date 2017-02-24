#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
import sys
import argparse
import string

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf file')

args        = parser.parse_args()

psf=pyfits.open(args.psf)
spots=psf["SPOTS"].data
ok=np.where((spots["STATUS"]==1))[0]

plt.figure("spots")
plt.subplot(1,2,1)
plt.plot(spots["X"],spots["Y"],"o",c="b")
plt.plot(spots["X"][ok],spots["Y"][ok],"o",c="r")
plt.xlabel("XCCD")
plt.ylabel("YCCD")
plt.subplot(1,2,2)
ok1=np.where((spots["EFLUX"]>0))[0]
ok2=np.where((spots["STATUS"]==1)&(spots["EFLUX"]>0))[0]
plt.plot(spots["WAVE"][ok1],spots["FLUX"][ok1]/spots["EFLUX"][ok1],"o",c="b")
plt.plot(spots["WAVE"][ok2],spots["FLUX"][ok2]/spots["EFLUX"][ok2],"o",c="r")
plt.xlabel("wavelength")
plt.ylabel("S/N")
plt.grid()
plt.show()

