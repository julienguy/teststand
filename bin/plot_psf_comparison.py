#!/usr/bin/env python

import numpy as np
import pyfits
import pylab
import sys
import argparse
import string

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf1', type = str, default = None, required = True,
                    help = 'path of psf file')
parser.add_argument('--psf2', type = str, default = None, required = True,
                    help = 'path of psf second file')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output image (png) file')


args        = parser.parse_args()

psf1=pyfits.open(args.psf1)[0].data
psf2=pyfits.open(args.psf2)[0].data


pylab.figure()
pylab.subplot(1,2,1)
pylab.imshow(psf1,origin=0,interpolation="nearest",vmin=0)
pylab.subplot(1,2,2)
pylab.imshow(psf2,origin=0,interpolation="nearest",vmin=0)

pylab.figure()
x1=np.arange(psf1.shape[1])-psf1.shape[1]/2.
x2=np.arange(psf2.shape[1])-psf2.shape[1]/2.
pylab.plot(x1,psf1[psf1.shape[0]/2,:])
pylab.plot(x2,psf2[psf2.shape[0]/2,:])
pylab.xlabel("X CCD (cross-dispersion axis)")
pylab.figure()
y1=np.arange(psf1.shape[0])-psf1.shape[0]/2.
y2=np.arange(psf2.shape[0])-psf2.shape[0]/2.
pylab.plot(y1,psf1[:,psf1.shape[1]/2])
pylab.plot(y2,psf2[:,psf2.shape[1]/2])
pylab.xlabel("Y CCD (dispersion axis)")

pylab.show()

#pyfits.writeto(args.output,image,clobber=True)
