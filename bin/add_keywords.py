#!/usr/bin/env python


import sys,string
import argparse
import numpy as np
import astropy.io.fits as pyfits

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True,
                    help = 'path to file to modify')
                  


args        = parser.parse_args()

h=pyfits.open(args.infile)
h[0].header["FLAVOR"]="none"
h.writeto(args.infile,clobber=True)

