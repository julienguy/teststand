#!/usr/bin/env python

import sys
import argparse
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from teststand.graph_tools         import plot_graph,parse_fibers
from desiutil.log                  import get_logger
import os.path
			
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--frame', type = str, default = None, required = True, nargs="*",help = 'path to one or several frame fits files')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('--wmin',type=float,default=3700,required=False,help="min wavelength")
parser.add_argument('--wmax',type=float,default=9700,required=False,help="max wavelength")
parser.add_argument('--sb',action="store_true",help="remove flux of side bands")

log         = get_logger()
args        = parser.parse_args()
fig         = plt.figure()
fibers      = parse_fibers(args.fibers)






def median_flux(wave,flux,rwave) :
    x=[]
    y=[]
    for i in range(bins.size-1) :
        j=np.where((wave>=bins[i])&(wave<bins[i+1]))[0]
        if j.size<2 : continue
        x.append(np.median(wave[j]))
        y.append(np.median(flux[j]))        
    return np.array(x),np.array(y)
#BUG ???? frame-b1-00002917.fits 18 n= 17.1304347826

wavestep=10
bins=np.linspace(args.wmin,args.wmax,int((args.wmax-args.wmin)/wavestep))
bwave=bins[:-1]+(bins[1]-bins[0])/2.
rwave=None
print("# expnum exptime expreq nd fiber flux")
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

        
        
        # subtract side region
        #zero=np.median(flux[(ivar>0)&(fwave>5500)&(fwave<6000)])
        #flux -= zero

        if True :
            x,y=median_flux(fwave[ivar!=0],flux[ivar!=0],bins)            
            if False :
                print(filename,fiber)
                plt.figure()
                plt.plot(fwave[ivar!=0],flux[ivar!=0],'.')
                plt.plot(x,y,"o")
                plt.show()
                
        else :
            if rwave is None :
                rwave=fwave[(fwave>args.wmin)&(fwave<args.wmax)]
            y=np.interp(rwave,fwave[ivar>0],flux[ivar>0],left=0,right=0)
        
        # apply non-linearity correction !!!
        # this is for b1-B
        # corr = (1-6.84e-6*y+1.346e-10*y**2)
        # y *= corr
        
        sflux=np.mean(y)
        sbflux=0
        if args.sb :
            sbflux = np.median(flux[(ivar>0)&(((fwave>args.wmin-200)&(fwave<args.wmin))|((fwave>args.wmax)&(fwave<args.wmax+200)))])
            sflux -= sbflux
        
        print("%d %f %f %d %02d %g %g"%(header["EXPNUM"],header["EXPTIME"],header["EXPREQ"],header["NDNUM"],fiber,sflux,sbflux))
        
