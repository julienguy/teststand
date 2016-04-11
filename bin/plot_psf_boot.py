#!/usr/bin/env python

import numpy as np
import pyfits
import pylab
from numpy.polynomial.legendre import legval
import desimodel.io
import desispec.io
#from specter.psf.gausshermite import GaussHermitePSF
import sys
import argparse

def u(wave,wavemin,wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf boot file')
parser.add_argument('-f','--fibermap',type = str, default = None, required = False,
                    help = 'path to fibermap file to compare with truth for simulations')
parser.add_argument('-a','--arm',type = str, default = None, required = False,
                    help = 'camera arm : b, r or z; to compare with truth for simulations')
args        = parser.parse_args()

if args.fibermap is not None :    
    if args.arm is None :
        print "need to know which camera arm b, r or z"
        sys.exit(12)
if args.arm is not None :
    if not args.arm in ['b','r','z'] :
        print "camera arm must be b, r or z"
        sys.exit(12)

psf=pyfits.open(args.psf)

print psf[0].header
print psf.info()
wavemin=psf[0].header["WAVEMIN"]
wavemax=psf[0].header["WAVEMAX"]
xcoef=psf[0].data
ycoef=psf[1].data
sigma=psf[2].data
print xcoef.shape
nspecs=xcoef.shape[0]

wave=np.linspace(wavemin,wavemax,100)

pylab.figure()
a0=pylab.subplot(1,3,1)
a1=pylab.subplot(1,3,2)
a2=pylab.subplot(1,3,3)
for spec in range(nspecs) :
    x = legval(u(wave,wavemin,wavemax), xcoef[spec])
    y = legval(u(wave,wavemin,wavemax), ycoef[spec])
    a0.plot(x,y)
    a1.plot(y,wave)
a0.set_xlabel("X CCD")
a0.set_ylabel("Y CCD")
a1.set_xlabel("Y CCD")
a1.set_ylabel("Wavelength [A]")
a2.plot(sigma)
a2.set_xlabel("spec #")
a2.set_ylabel("PSF sigma")

if args.fibermap is not None :
    
    if args.arm is None :
        print "need to know which camera arm b, r or z"
        sys.exit(12)
    
    
    fm, fmhdr = desispec.io.read_fibermap(args.fibermap, header=True)
    fibers = fm["FIBER"][:nspecs]
    psf = desimodel.io.load_psf(args.arm)
    


    """
    print "loading simpix traces",args.simpix
    simpix=pyfits.open(args.simpix)
    simpix.info()
    xcoef_truth=simpix[1].data
    xwavemin_truth=simpix[1].header["WAVEMIN"]
    xwavemax_truth=simpix[1].header["WAVEMAX"]
    ycoef_truth=simpix[2].data
    ywavemin_truth=simpix[2].header["WAVEMIN"]
    ywavemax_truth=simpix[2].header["WAVEMAX"]
    """
    fig = pylab.figure()
    a0=pylab.subplot(1,2,1)
    a1=pylab.subplot(1,2,2)
    for spec,fiber in enumerate(fibers) : 
        x = legval(u(wave,wavemin,wavemax), xcoef[spec])
        y = legval(u(wave,wavemin,wavemax), ycoef[spec])
        x_truth = psf.x(int(fiber),wave)
        y_truth = psf.y(int(fiber),wave)
        a0.plot(wave,x-x_truth)
        a1.plot(wave,y-y_truth)
    a0.set_xlabel("Wavelength [A]")
    a0.set_ylabel("delta X CCD")
    a1.set_xlabel("Wavelength [A]")
    a1.set_ylabel("delta Y CCD")
    fig.savefig(args.psf.replace(".fits",".png"))
pylab.show()
