#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
from teststand.graph_tools import parse_fibers
import sys
import argparse
import string

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, nargs="*", default = None, required = True,
                    help = 'path of psf file')
parser.add_argument('--fibers', type = str, required=False, default=None, help= "defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)")

args        = parser.parse_args()

plt.figure("spots")
a1=plt.subplot(1,2,1)
a2=plt.subplot(1,2,2)

fibers=parse_fibers(args.fibers)

for filename in args.psf :
    psf=pyfits.open(filename)
    spots=psf["SPOTS"].data
    if fibers is None :
        fibers = np.unique(spots["FIBER"])

    for fiber in fibers :
        ok1=np.where((spots["FIBER"]==fiber))[0]
        ok2=np.where((spots["STATUS"]==1)&(spots["FIBER"]==fiber))[0]
        a1.plot(spots["X"][ok1],spots["Y"][ok1],"x",c="gray")
        a1.plot(spots["X"][ok2],spots["Y"][ok2],"o")


        ok1=np.where((spots["EFLUX"]>0)&(spots["FIBER"]==fiber))[0]
        ok2=np.where((spots["STATUS"]==1)&(spots["EFLUX"]>0)&(spots["FIBER"]==fiber))[0]
        #a2.plot(spots["WAVE"][ok1],spots["FLUX"][ok1]/spots["EFLUX"][ok1],"x",c="gray")
        #a2.plot(spots["WAVE"][ok2],spots["FLUX"][ok2]/spots["EFLUX"][ok2],"o")
        a2.plot(spots["WAVE"][ok1],spots["FLUX"][ok1],"x",c="gray")
        a2.plot(spots["WAVE"][ok2],spots["FLUX"][ok2],"o")
        
a1.set_xlabel("XCCD")
a1.set_ylabel("YCCD")
a2.set_xlabel("wavelength")
#a2.set_ylabel("S/N")
a2.set_ylabel("FLUX")
a2.grid()        
plt.show()

