#!/usr/bin/env python

import numpy as np
import pyfits
import specter.psf
import sys
import argparse
import string
import os.path
def readpsf(filename) :
    try :
        psftype=pyfits.open(filename)[0].header["PSFTYPE"]
    except KeyError :
        psftype=""
    print("PSF Type=",psftype)
    if psftype=="GAUSS-HERMITE" :
        return specter.psf.GaussHermitePSF(filename)
    elif psftype=="SPOTGRID" :
        return specter.psf.SpotGridPSF(filename)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf', type = str, default = None, required = True,
                    help = 'path of psf file')
parser.add_argument('--fiber', type = int, default = None, required = True,
                    help = 'fiber')
parser.add_argument('--wave', type = float, default = 6000., required = True,
                    help = 'wavelength')
parser.add_argument('-o','--output', type = str, default = None, required = True,
                    help = 'path to output fits image')


args        = parser.parse_args()

psf=readpsf(args.psf)
xy=psf.xy(args.fiber,args.wave)
hw=4.
n1d=2*hw*8+1
x=np.tile(np.linspace(-hw,hw,n1d),(n1d,1))
y=x.T
fpix=psf._value(x+xy[0],y+xy[1],args.fiber,args.wave)
ds=(x[0,1]-x[0,0])*(y[1,0]-x[0,0])
print(ds)
fpix *= 1./(np.sum(fpix*ds))

pyfits.writeto(args.output,fpix,clobber=True)
