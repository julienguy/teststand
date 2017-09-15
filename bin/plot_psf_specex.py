#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import pylab
from numpy.polynomial.legendre import legval
import desimodel.io
import desispec.io
import sys
import argparse
import string

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
args        = parser.parse_args()

psf=pyfits.open(args.psf)

xhdu=psf["XTRACE"]
yhdu=psf["YTRACE"]
wavemin=xhdu.header["WAVEMIN"]
wavemax=xhdu.header["WAVEMAX"]
xcoef=xhdu.data
ycoef=yhdu.data

psfhdu=psf["PSF"]
cam=psfhdu.header["CAMERA"].strip().replace("'","").strip()
arm=cam[0]
print("CAMERA=",cam,"ARM=",arm)


params=psfhdu.data["PARAM"]
for i in range(params.size) :
    params[i]=params[i].strip()
ghsigx_index= np.where(params=="GHSIGX")[0][0]
ghsigy_index= np.where(params=="GHSIGY")[0][0]
fibermin=int(psfhdu.header["FIBERMIN"])
fibermax=int(psfhdu.header["FIBERMAX"])
legdeg=psfhdu.header["LEGDEG"]
table=psfhdu.data
ghsigx_coef=table["COEFF"][ghsigx_index]
ghsigy_coef=table["COEFF"][ghsigy_index]


print("wavemin,wavemax=",wavemin,wavemax)
nspec=xcoef.shape[0]

wave=np.linspace(wavemin,wavemax,100)
refwave=int(np.mean(wave))

fig=pylab.figure()

nx=2
if args.sim :
    ny=3
else :
    ny=3

pcount=1
a0=pylab.subplot(ny,nx,pcount) ; pcount+=1
a5=pylab.subplot(ny,nx,pcount) ; pcount+=1
a1=pylab.subplot(ny,nx,pcount) ; pcount+=1
a3=pylab.subplot(ny,nx,pcount) ; pcount+=1
a2=pylab.subplot(ny,nx,pcount) ; pcount+=1
a4=pylab.subplot(ny,nx,pcount) ; pcount+=1

for spec in range(nspec) :
    x = legval(u(wave,wavemin,wavemax), xcoef[spec])
    y = legval(u(wave,wavemin,wavemax), ycoef[spec])
    a0.plot(x,y)
    ghsigx = legval(u(wave,wavemin,wavemax), ghsigx_coef[spec])
    a1.plot(wave,ghsigx)
    ghsigy = legval(u(wave,wavemin,wavemax), ghsigy_coef[spec])
    a2.plot(wave,ghsigy)
    ghsigx = legval(u(refwave,wavemin,wavemax), ghsigx_coef[spec])
    a3.plot(spec,ghsigx,"o")
    ghsigy = legval(u(refwave,wavemin,wavemax), ghsigy_coef[spec])
    a4.plot(spec,ghsigy,"o")
    a5.plot(y,wave)

a0.set_xlabel("X CCD")
a0.set_ylabel("Y CCD")
a5.set_xlabel("Y CCD")
a5.set_ylabel("Wavelength [A]")

#a1.set_xlabel("Wavelength [A]")
a1.set_ylabel("Cross-dispersion sigma (pixels)")

#a3.set_xlabel("Fiber #, @%dA"%refwave)
#a3.set_ylabel("Cross-dispersion sigma (pixels)")

a2.set_xlabel("Wavelength [A]")
a2.set_ylabel("Dispersion/resolution sigma (pixels)")

a4.set_xlabel("Fiber #, @%dA"%refwave)
#a4.set_ylabel("Dispersion/resolution sigma (pixels)")






if args.sim :
    
    
    if args.fibermap is not None :
        fm, fmhdr = desispec.io.read_fibermap(args.fibermap, header=True)
        fibers = fm["FIBER"][:nspec]
    else :
        fibers = np.arange(nspec)
        print("assuming it's the first %d fibers in the sims (if wrong, rerun with --fibermap option)"%nspec)
    
    psf = desimodel.io.load_psf(arm)
    

    
    a0=pylab.subplot(ny,nx,pcount) ; pcount+=1
    a1=pylab.subplot(ny,nx,pcount) ; pcount+=1
    
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
    

if args.fig is not None :
    fig.savefig(args.fig)

if not args.batch :
    pylab.show()

