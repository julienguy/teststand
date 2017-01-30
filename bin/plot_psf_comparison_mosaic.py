#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import pylab
import specter.psf
import sys
import argparse
import string
import os.path
def readpsf(filename) :
    try :
        psftype=pyfits.open(filename)[0].header["PSFTYPE"]
    except KeyError :
        psftype=""
    print("PSF Type=",psftype)
    if psftype=="GAUSS-HERMITE" :
        return specter.psf.GaussHermitePSF(filename)
    elif psftype=="SPOTGRID" :
        return specter.psf.SpotGridPSF(filename)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf1', type = str, default = None, required = True,
                    help = 'path of psf file')
parser.add_argument('--psf2', type = str, default = None, required = True,
                    help = 'path of psf second file')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output image (png) file')


fibers1=[0,10,19]
fibers2=[0,250,490]
fibers2=[0,10,19]
waves=[5500,6700,7600]



args        = parser.parse_args()

psf1=readpsf(args.psf1)
psf2=readpsf(args.psf2)
name1=os.path.basename(args.psf1)
name2=os.path.basename(args.psf2)



f,a = pylab.subplots(len(waves),len(fibers1),sharex=True, sharey=True)

for i in range(len(fibers1)) :

    fiber1=fibers1[i]
    fiber2=fibers2[i]

    print("fiber1 %d fiber2 %d"%(fiber1,fiber2))

    for j in range(len(waves)) :
        wave=waves[len(waves)-1-j]
        
        xy1=psf1.xy(fiber1,wave)
        xy2=psf2.xy(fiber2,wave)
        
        print("for psf1, xy=",xy1)
        print("for psf2, xy=",xy2)

        hw=5.
        n1d=51
        x=np.tile(np.linspace(-hw,hw,51),(n1d,1))
        y=x.T
        fpix1=psf1._value(x+xy1[0],y+xy1[1],fiber1,wave)
        fpix2=psf2._value(x+xy2[0],y+xy2[1],fiber2,wave)
        fpix1 /= np.sum(fpix1)
        fpix2 /= np.sum(fpix2)

        aa=a[j,i]
        #aa.imshow(fpix1,origin=0,interpolation="nearest",extent=(-hw,hw,-hw,hw),aspect="auto")
        levels = np.array([0.1,0.5,0.9])*np.max(fpix1)
        aa.contour(x,y,fpix1,colors='b',levels=levels)
        aa.contour(x,y,fpix2,colors='r',levels=levels)
        
        if False :
            #aa.set_title("#%d / #%d"%(fiber1,fiber2))
            aa.text(-hw+0.3,hw-2,"%s Fiber #%d"%(name1,fiber1),fontsize=10,color="b")
            aa.text(-hw+0.3,hw-1.2,"%s Fiber #%d"%(name2,fiber2),fontsize=10,color="r")
        if True :
            aa.text(-hw+0.3,-hw+1.7,"x y psf:",fontsize=10,color="k")
            aa.text(-hw+0.3,-hw+0.3,"%4.1f %4.1f %s"%(xy1[0],xy1[1],name1),fontsize=10,color="b")
            aa.text(-hw+0.3,-hw+1.,"%4.1f %4.1f %s"%(xy2[0],xy2[1],name2),fontsize=10,color="r")
        #aa.legend(loc="upper left",fontsize="small")
        
        if i==0 :
            aa.set_ylabel("%dA"%wave)
        if j==0 :
            aa.set_title("fibers %d&%d"%(fiber1,fiber2))
        #aa.set_xlim([-hw,hw])
        #aa.set_ylim([-hw,hw])

f.subplots_adjust(hspace=0,wspace=0)    

if args.output is not None :
   fig.savefig(args.output)
 
pylab.show()
