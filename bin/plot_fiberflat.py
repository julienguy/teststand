#!/usr/bin/env python

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import argparse
import os.path

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--infile', type = str, default = None, required = True, nargs="*",
                    help = 'path to fiber flat file')
parser.add_argument('--spectra', action = 'store_true')

args        = parser.parse_args()

x=[]
y=[]
z=[]

for filename in args.infile :
    
    h=fits.open(filename)
    h.info()
    wave=h["WAVELENGTH"].data

    if args.spectra :
        figname=os.path.basename(filename).split(".")[0]
        plt.figure(figname)
        for fiber in range(h[0].data.shape[0]) :
            ok=np.where(h["MASK"].data[fiber]==0)[0]
            plt.plot(wave,h[0].data[fiber],color="gray",alpha=0.1)
            if ok.size>1 :
              plt.plot(wave[ok],h[0].data[fiber,ok])
            nbad=sum(h[0].data[fiber]>1.5)
            if nbad>20:
                print("fiber {} as {} flux points with flat>1.5".format(fiber,nbad))
            nbad=sum(h[0].data[fiber]<0.6)
            if nbad>20:
                print("fiber {} as {} flux points with flat<0.6".format(fiber,nbad))
    if "FIBERMAP" in h :
        h[0].data *= (h["MASK"].data == 0 )
        z.append( np.median(h[0].data,axis=1) )
        x.append( h["FIBERMAP"].data["FIBERASSIGN_X"] )
        y.append( h["FIBERMAP"].data["FIBERASSIGN_Y"] )


if len(x)>0 :
    x=np.hstack(x)
    y=np.hstack(y)
    z=np.hstack(z)
    plt.figure("focalplane")
    plt.scatter(x,y,c=z,vmin=0.7,vmax=1.2)
    plt.colorbar()


  
plt.show()
