#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import pylab
import specter.psf
import sys
import argparse
import string
import scipy.interpolate

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf boot file')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output image (png) file')


args        = parser.parse_args()

psf=specter.psf.GaussHermitePSF(args.psf)
wmin=max(3650,psf.wmin+500)
wmax=min(10200,psf.wmax-500)


wave=np.linspace(wmin,wmax,15)

zoom_out=32

image=np.zeros((psf.npix_y//zoom_out,psf.npix_x//zoom_out))
fibers=10+20*np.arange(25)
for fiber in fibers :
    for wavelength in wave :
        print(fiber,wavelength)
        x,y = psf.xy(fiber,wavelength)
        xslice,yslice,pix = psf.xypix(fiber,wavelength)
        #dx=xslice.stop-xslice.start
        #dy=yslice.stop-yslice.start
        dx=(x/zoom_out-x)
        dy=(y/zoom_out-y)
        idx=int(dx+0.5)
        idy=int(dy+0.5)
        
        ipix=np.zeros(pix.shape)
        x=np.arange(xslice.start,xslice.stop)+idx
        y=np.arange(yslice.start,yslice.stop)+idy
        f=scipy.interpolate.interp2d(x,y,pix) 
        ipix = f(x+(idx-dx), y+(idy-dy))
        xslice = slice(xslice.start+idx,xslice.stop+idx)
        yslice = slice(yslice.start+idy,yslice.stop+idy)   
        
        try :
            image[yslice,xslice] += ipix
        except ValueError :
            print("failed for fiber wave=",fiber,wavelength)
            pass

fig=pylab.figure()
pylab.imshow(image,origin=0,interpolation="nearest",extent=(0,psf.npix_x,0,psf.npix_y),vmin=0,vmax=0.15)
if args.output is not None :
    fig.savefig(args.output)
pylab.show()

#pyfits.writeto(args.output,image,clobber=True)
