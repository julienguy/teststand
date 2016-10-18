#!/usr/bin/env python


import sys,string
import astropy.io.fits as pyfits
import argparse
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--image', type = str, default = None, required = True, nargs="*",
                    help = 'path of image fits files')
parser.add_argument('-o','--outfile', type = str, default = None, required = True, 
                    help = 'output median image filename')
parser.add_argument('--hdu',type = str, default = 0, required = False, 
                    help = 'header HDU (int or string)')

args        = parser.parse_args()

try :
    hdu = int(args.hdu)
except ValueError:
    hdu = args.hdu

print "read images ..."
images=[]
shape=None
for filename in args.image :
    print filename
    image=pyfits.open(filename)[hdu].data.astype("float32")
    if shape is None :
        shape=image.shape
    images.append(image.ravel())

print "compute median image ..."
medimage=np.median(images,axis=0).reshape(shape)

print "write ..."
hdulist=pyfits.HDUList([pyfits.PrimaryHDU(medimage)])
i=0
for filename in args.image :
    hdulist[0].header["INPUT%03d"%i]=filename
    i+=1
hdulist.writeto(args.outfile,clobber="True")


