#!/usr/bin/env python

import numpy as np
import pyfits
import pylab
import specter.psf
import sys
import argparse
import string

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf boot file')
parser.add_argument('-o','--output', type = str, default = None, required = True,
                    help = 'path to output fits file')

args        = parser.parse_args()

psf=specter.psf.GaussHermitePSF(args.psf)
wave=np.linspace(psf.wmin+200,psf.wmax-200,15)

zoom=16
image=np.zeros((psf.npix_y/zoom,psf.npix_x/zoom))
for fiber in range(psf.nspec) :
    for wavelength in wave :
        print fiber,wavelength
        x,y = psf.xy(fiber,wavelength)
        xslice,yslice,pix = psf.xypix(fiber,wavelength)
        
        xslice = slice(xslice.start+int(x/zoom-x),xslice.stop+int(x/zoom-x))
        yslice = slice(yslice.start+int(y/zoom-y),yslice.stop+int(y/zoom-y))    
        try :
            image[yslice,xslice] += pix
        except ValueError :
            print "failed for fiber wave=",fiber,wavelength
            pass
      
pyfits.writeto(args.output,image,clobber=True)
