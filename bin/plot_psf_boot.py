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
parser.add_argument('--fig',type = str, default = None, required = False,
                    help = 'figure filename')
parser.add_argument('--batch', action = 'store_true',help="do not display result")
parser.add_argument('--sim', action = 'store_true',help="compare with simulation truth")


args = parser.parse_args()


psf=pyfits.open(args.psf)
cam=psf[0].header["CAMERA"].strip()
arm=cam[0].lower()
if not arm in ['b','r','z'] :
    print "camera arm must be b, r or z, and read '%s' in psf header"%arm
    sys.exit(12)

#print psf[0].header
#print psf.info()
wavemin=psf[0].header["WAVEMIN"]
wavemax=psf[0].header["WAVEMAX"]
print "wavemin,wavemax=",wavemin,wavemax
xcoef=psf[0].data
ycoef=psf[1].data
sigma=psf[2].data
print "xcoef.shape=",xcoef.shape
print "ycoef.shape=",ycoef.shape
nspec=xcoef.shape[0]
print "nspec=",nspec

wave=np.linspace(wavemin,wavemax,100)
mwave=np.mean(wave)


fig = pylab.figure()
nx=3
ny=2

pcount=1
a0=pylab.subplot(ny,nx,pcount) ; pcount +=1
a1=pylab.subplot(ny,nx,pcount) ; pcount +=1
a2=pylab.subplot(ny,nx,pcount) ; pcount +=1
a3=pylab.subplot(ny,nx,pcount) ; pcount +=1


mx=[]
for spec in range(nspec) :

    x = legval(u(wave,wavemin,wavemax), xcoef[spec])
    y = legval(u(wave,wavemin,wavemax), ycoef[spec])
    a0.plot(x,y)
    a1.plot(y,wave)
    mx.append(legval(u(mwave,wavemin,wavemax), xcoef[spec]))


a0.set_xlabel("X CCD")
a0.set_ylabel("Y CCD")
a1.set_xlabel("Y CCD")
a1.set_ylabel("Wavelength [A]")
a2.plot(sigma)
a2.set_xlabel("spec #")
a2.set_ylabel("PSF sigma")
mx=np.array(mx)
dx=(mx-np.roll(mx,1))[1:]
print dx[:12]
a3.plot(dx,"o-")
a3.set_xlabel("spec #")
a3.set_ylabel("Delta X CCD @%dA"%int(mwave))

if args.sim :
    
    if args.fibermap is not None :
        fm, fmhdr = desispec.io.read_fibermap(args.fibermap, header=True)
        fibers = fm["FIBER"][:nspec]
    else :
        fibers = np.arange(nspec)
        print "assuming it's the first %d fibers in the sims (if wrong, rerun with --fibermap option)"%nspec
    psf = desimodel.io.load_psf(arm)
    


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
    
    a0=pylab.subplot(ny,nx,pcount) ; pcount +=1
    a1=pylab.subplot(ny,nx,pcount) ; pcount +=1
    for spec,fiber in enumerate(fibers) : 
        print "spec #%d fiber #%d"%(spec,fiber)
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

if args.fig is not None :
    fig.savefig(args.fig)

if not args.batch :
    pylab.show()

