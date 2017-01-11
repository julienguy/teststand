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
parser.add_argument('--with-primary-header',action="store_true")

args        = parser.parse_args()

try :
    hdu = int(args.hdu)
except ValueError:
    hdu = args.hdu

print "read images ..."
images=[]
shape=None
primary_header=None
image_header=None
for filename in args.image :
    print filename
    fitsfile=pyfits.open(filename)
    
    image=fitsfile[hdu].data.astype("float32")
    if shape is None :
        shape=image.shape
    images.append(image.ravel())
    
    if primary_header is None and args.with_primary_header :
        primary_header=fitsfile[0].header
    if image_header is None :
        image_header=fitsfile[hdu].header

print "compute median image ..."
medimage=np.median(images,axis=0).reshape(shape)

print "write ..."
if not args.with_primary_header :
    hdulist=pyfits.HDUList([pyfits.PrimaryHDU(medimage)])
    hdu=0
else :
    hdulist=pyfits.HDUList([pyfits.PrimaryHDU(),pyfits.ImageHDU(medimage,name=hdu)])
    hdulist[0].header=primary_header
    hdu=1

hdulist[hdu].header=image_header
i=0
for filename in args.image :
    hdulist[hdu].header["INPUT%03d"%i]=filename
    i+=1
hdulist.writeto(args.outfile,clobber="True")


