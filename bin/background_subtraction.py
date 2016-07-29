#!/usr/bin/env python

import astropy.io.fits as pyfits
import argparse
import sys
import numpy as np
import re
import scipy.interpolate
#import pylab

def parse_sec_keyword(value):
    m = re.search('\[(\d+):(\d+)\,(\d+):(\d+)\]', value)
    if m is None:
        raise ValueError, 'unable to parse {} as [a:b, c:d]'.format(value)
    return  map(int, m.groups())

parser=argparse.ArgumentParser(description="Reformat DESI raw data")
parser.add_argument('-i','--infile', type = str, default = None, required=True,
                    help = 'path to input DESI raw data file')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,
                    help = 'path to output DESI raw data file')
parser.add_argument('-s','--step', type = int, default = 200, required=False,
                    help = 'step in pixels of the median filter')

args = parser.parse_args()

ifile=pyfits.open(args.infile)
image = ifile[0].data
header = ifile[0].header

bkg=0.*image

# do it independently per amplifier

for amp in range(1,5) :
    xmin,xmax,ymin,ymax = parse_sec_keyword(header["CCDSEC%d"%amp])
    print amp,xmin,xmax,ymin,ymax
    xmin-=1
    ymin-=1
    nx=xmax-xmin
    ny=ymax-xmin
    xbins=np.linspace(xmin,xmax,int((xmax-xmin)/args.step)).astype(int)
    ybins=np.linspace(ymin,ymax,int((ymax-ymin)/args.step)).astype(int)

    bkg_grid=np.zeros((ybins.size-1,xbins.size-1))
    for i in range(xbins.size-1) :
        for j in range(ybins.size-1) :
            bkg_grid[j,i]=np.median(image[ybins[j]:ybins[j+1],xbins[i]:xbins[i+1]])
            print i,j,bkg_grid[j,i]
    
    x=xbins[:-1]+(xbins[1]-xbins[0])/2.
    y=ybins[:-1]+(ybins[1]-ybins[0])/2.
    f = scipy.interpolate.interp2d(x, y, bkg_grid, kind='linear')
    bkg[ymin:ymax,xmin:xmax] = f(np.arange(xmin,xmax),np.arange(ymin,ymax))
    print "done amp",amp

"""
xmin,xmax,ymin,ymax = parse_sec_keyword(header["CCDSEC4"])
xsplit=xmin
ysplit=ymin

x=xbins[:-1]+(xbins[1]-xbins[0])/2.
y=ybins[:-1]+(ybins[1]-ybins[0])/2.

f = scipy.interpolate.interp2d(x, y, bkg, kind='cubic')
bkg = f(np.arange(image.shape[1]),np.arange(image.shape[0]))
"""
ifile[0].data -= bkg
ifile.append(pyfits.ImageHDU(bkg,name="BKG"))
ifile.writeto(args.outfile,clobber=True)
