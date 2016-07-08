#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import argparse
import sys
import os.path
import desimodel.io

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required=True,
                        help = 'path to main-peaks-vacuum.data')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,
                        help = 'path to output file ')

args = parser.parse_args()

wave=[]
ion=[]
intensity=[]
ifile=open(args.infile)
for line in ifile.readlines() :
    if line[0]=="#" :
        continue
    vals=line.strip().split()
    if len(vals)==3 :
        wave.append(float(vals[0]))
        ion.append((vals[1]))
        intensity.append((float(vals[2])))
    else :
        print "WARNING IGNORE LINE '%s'"%line.strip()
    ifile.close()
wave=np.array(wave)
ion=np.array(ion)
intensity=np.array(intensity)

ok=np.where(intensity>0)[0]
wave=wave[ok]
ion=ion[ok]
intensity=intensity[ok]
intensity *= 100000./np.max(intensity)

# write this
cols=[]
cols.append(fits.Column(name='WAVE', format='D', array=wave))
cols.append(fits.Column(name='PHOTONS', format='D', array=intensity))
cols.append(fits.Column(name='ION', format='8A', array=ion))
cols = fits.ColDefs(cols)
hdulist=fits.HDUList([fits.PrimaryHDU()])
hdulist.append(fits.BinTableHDU.from_columns(cols))
hdulist.writeto(args.outfile,clobber=True)

print "wrote %d lines in %s"%(wave.size,args.outfile)
