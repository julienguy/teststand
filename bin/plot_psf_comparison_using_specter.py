#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import specter.psf
import sys
import argparse
import string
import os.path
from   scipy.signal import fftconvolve

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
parser.add_argument('--fiber', type = int, default = None, required = True,
                    help = 'fiber for psf1')
parser.add_argument('--fiber2', type = int, default = None, required = False,
                    help = 'fiber for psf2 (default=fiber1)')
parser.add_argument('--wavelength', type = float, default = 6000., required = False,
                    help = 'wavelength')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output image (png) file')
parser.add_argument('--no-pixel-convolution', action = "store_true",
                    help = 'do not convolve PSFs with pixel size')
parser.add_argument('--transpose', action = "store_true",
                    help = 'flip x and y for psf2')


args        = parser.parse_args()

if args.fiber2 is None :
    fiber2=args.fiber
else :
    fiber2=args.fiber2

psf1=readpsf(args.psf1)
psf2=readpsf(args.psf2)
xy1=psf1.xy(args.fiber,args.wavelength)
xy2=psf2.xy(fiber2,args.wavelength)
print("for psf1, xy=",xy1)
print("for psf2, xy=",xy2)



hw=5.
n1d=51
x1d=np.linspace(-hw,hw,51)
x=np.tile(x1d,(n1d,1))
y=x.T
fpix1=psf1._value(x+xy1[0],y+xy1[1],args.fiber,args.wavelength)
fpix2=psf2._value(x+xy2[0],y+xy2[1],fiber2,args.wavelength)
if args.transpose :
    fpix2=fpix2.T
fpix1 /= np.sum(fpix1)
fpix2 /= np.sum(fpix2)

if not args.no_pixel_convolution :
    
    print("convolve PSF with pixel size")
    kernel=(np.abs(x1d)<0.5).astype(float)
    kernel/=np.sum(kernel)
    for i in range(n1d) :
        fpix1[i]=fftconvolve(fpix1[i],kernel, mode='same')
    for j in range(n1d) :
        fpix1[:,j]=fftconvolve(fpix1[:,j],kernel, mode='same')
    for i in range(n1d) :
        fpix2[i]=fftconvolve(fpix2[i],kernel, mode='same')
    for j in range(n1d) :
        fpix2[:,j]=fftconvolve(fpix2[:,j],kernel, mode='same')

fpix1 /= np.sum(fpix1)
fpix2 /= np.sum(fpix2)



mx1=np.sum(fpix1*x)
my1=np.sum(fpix1*y)
sigx1=np.sqrt(np.sum(fpix1*(x-mx1)**2))
sigy1=np.sqrt(np.sum(fpix1*(y-my1)**2))
mx2=np.sum(fpix2*x)
my2=np.sum(fpix2*y)
sigx2=np.sqrt(np.sum(fpix2*(x-mx2)**2))
sigy2=np.sqrt(np.sum(fpix2*(y-my2)**2))

print("psf1 sigx=%f sigy=%f"%(sigx1,sigy1))
print("psf2 sigx=%f sigy=%f"%(sigx2,sigy2))
print("sigx1/sigx2=%f sigy1/sigy2=%f"%(sigx1/sigx2,sigy1/sigy2))







plt.figure()
plt.figtext(0.25, 0.96,os.path.basename(args.psf1), fontsize='medium', color='b', ha ='center')
plt.figtext(0.75, 0.96,os.path.basename(args.psf2), fontsize='medium', color='r', ha ='center')

a=plt.subplot(2,2,1)
plt.imshow(fpix1,origin=0,interpolation="nearest",extent=(-hw,hw,-hw,hw))
plt.text(-hw+0.3,-hw+0.8,"fiber #%d lambda=%dA"%(args.fiber,args.wavelength),fontsize=10,color="white")
plt.text(-hw+0.3,-hw+0.1,"(x,y)=(%4.1f,%4.1f)"%(xy1[0],xy1[1]),fontsize=10,color="white")
plt.subplot(2,2,2)
plt.imshow(fpix2,origin=0,interpolation="nearest",extent=(-hw,hw,-hw,hw))
plt.text(-hw+0.3,-hw+0.8,"fiber #%d lambda=%dA"%(fiber2,args.wavelength),fontsize=10,color="white")
plt.text(-hw+0.3,-hw+0.1,"(x,y)=(%4.1f,%4.1f)"%(xy2[0],xy2[1]),fontsize=10,color="white")
a=plt.subplot(2,2,3,title="x prof.")

plt.plot(x[n1d//2,:],fpix1[n1d//2,:],c="b",lw=2)
plt.plot(x[n1d//2,:],fpix2[n1d//2,:],"--",c="r",lw=2)
plt.xlabel("x ccd")
plt.gca().axes.get_yaxis().set_ticks([])
a=plt.subplot(2,2,4,title="y prof.")
plt.plot(y[:,n1d//2],fpix1[:,n1d//2],c="b",lw=2)
plt.plot(y[:,n1d//2],fpix2[:,n1d//2],"--",c="r",lw=2)
plt.xlabel("y ccd")
plt.gca().axes.get_yaxis().set_ticks([])
# do the ratio of psf2 integrals

ratio = np.sum(fpix1*fpix2)/np.sum(fpix2**2)
print("single line photometric error=%f"%np.abs(ratio-1.))



plt.show()

#pyfits.writeto(args.output,image,clobber=True)
