#!/usr/bin/env python

import sys
import argparse
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from desispec.qproc.util import parse_fibers
from desiutil.log import get_logger
import os.path
			
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--frame', type = str, default = None, required = True, nargs="*",help = 'path to one or several frame fits files')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('--wmin',type=float,default=3700,required=False,help="min wavelength")
parser.add_argument('--wmax',type=float,default=9700,required=False,help="max wavelength")

log         = get_logger()
args        = parser.parse_args()
fig         = plt.figure()
fibers      = parse_fibers(args.fibers)

def median_flux_in_bins(wave,flux,ivar,bins) :
    x=[]
    y=[]
    yivar=[]
    for i in range(bins.size-1) :
        j=np.where((wave>=bins[i])&(wave<bins[i+1]))[0]
        if j.size<2 : continue
        x.append(np.median(wave[j]))
        y.append(np.median(flux[j]))
        yivar.append(2./np.pi*np.sum(ivar[j]))
    return np.array(x),np.array(y),np.array(yivar)

wavestep=10
bins=np.linspace(args.wmin,args.wmax,int((args.wmax-args.wmin)/wavestep))
bwave=bins[:-1]+(bins[1]-bins[0])/2.
print("# expnum exptime expreq nd fiber flux ivar")
for filename in args.frame :
    
    expnum=int(os.path.basename(filename).split("-")[2].split(".")[0])

    h  = pyfits.open(filename)
    header = h[0].header
    
    if fibers is None :
        fibers = np.arange(h[0].data.shape[0])

    for fiber in fibers :
        flux=h[0].data[fiber]
        ivar=h[1].data[fiber]
        fwave=h["WAVELENGTH"].data[fiber]
        if 1 : # first a median per bin
            bwave,bflux,bivar=median_flux_in_bins(fwave[ivar!=0],flux[ivar!=0],ivar[ivar!=0],bins)            
            mean_flux      = np.sum(bflux)/bflux.size
            mean_flux_var  = np.sum(1./bivar)/(bflux.size)**2
            mean_flux_ivar  = 1./mean_flux_var
        
        
        print("%d %f %f %d %02d %g %g"%(header["EXPNUM"],header["EXPTIME"],header["EXPREQ"],header["NDNUM"],fiber,mean_flux,mean_flux_ivar))
        
