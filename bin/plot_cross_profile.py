#!/usr/bin/env python

import numpy as np
import pyfits
import pylab
from numpy.polynomial.legendre import legval
import desimodel.io
import desispec.io
#from specter.psf.gausshermite import GaussHermitePSF
import sys
import argparse
from numpy.polynomial.legendre import legval, legfit

def u(wave,wavemin,wavemax) :
    return 2.*(wave-wavemin)/(wavemax-wavemin)-1.

def invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, sigmacoef, fiber, npix_y, wave_of_y, width=7) :
 
    #   Wavelength array used in 'invert_legendre_polynomial'
    wave                = np.linspace(wavemin, wavemax, 100)
    #   Determines value of Y, so we can know its coeficient and then its position
    y_of_wave           = legval(u(wave, wavemin, wavemax), ycoef[fiber])
    coef                = legfit(u(y_of_wave, 0, npix_y), wave, deg=ycoef[fiber].size)
    wave_of_y[fiber]    = legval(u(np.arange(npix_y).astype(float), 0, npix_y), coef)
    #   Determines wavelength intensity (x) based on Y
    x_of_y              = legval(u(wave_of_y[fiber], wavemin, wavemax), xcoef[fiber])
    sigma_of_y          = legval(u(wave_of_y[fiber], wavemin, wavemax), sigmacoef[fiber])
    
    #   Ascertain X by using low and high uncertainty
    x1_of_y             = np.floor(x_of_y).astype(int) - width/2
    x2_of_y             = np.floor(x_of_y).astype(int) + width/2 + 2
    return (x1_of_y, x_of_y, x2_of_y, sigma_of_y)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf boot file')
parser.add_argument('-i','--image',type = str, default = None, required = True,
                    help = 'path to preprocessed image')
parser.add_argument('--nbins',type = int, default = 4, required = False,
                    help = 'number of wavelength bins')

args=parser.parse_args()

psf=pyfits.open(args.psf)
cam=psf[0].header["CAMERA"].strip()
arm=cam[0].lower()
if not arm in ['b','r','z'] :
    print("camera arm must be b, r or z, and read '%s' in psf header"%arm)
    sys.exit(12)

wavemin=psf[0].header["WAVEMIN"]
wavemax=psf[0].header["WAVEMAX"]
print("wavemin,wavemax=",wavemin,wavemax)
xcoef=psf[0].data
ycoef=psf[1].data
sigmacoef=psf[2].data
print("xcoef.shape=",xcoef.shape)
print("ycoef.shape=",ycoef.shape)
nspec=xcoef.shape[0]
print("nspec=",nspec)

image_file=pyfits.open(args.image)
flux        = image_file[0].data
#   Inverse variance of the image's value
flux_ivar   = image_file[1].data
#   Variance based on inverse variance's size
flux_var    = np.zeros(flux_ivar.shape)
mask        = (flux_ivar > 0)
flux_var[mask] = 1./flux_ivar[mask]
#   Number of pixels in an image 
#   We are going to extract one flux per fiber per Y pixel (total = nfibers x npix_y)
npix_y  = flux.shape[0]
npix_x  = flux.shape[1]
width   = 5

wavebins=np.linspace(wavemin,wavemax,args.nbins+1)
wave_of_y=np.zeros((nspec, npix_y))

oversampling=4
x=np.linspace(-width,width,2*width*oversampling+1)

#nspec=3

f,a = pylab.subplots(args.nbins,nspec, sharex=True, sharey=True)
for spec in range(nspec) :
    
    

    x1_of_y, xc_of_y, x2_of_y, sigma_of_y = invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, sigmacoef, spec, npix_y, wave_of_y, width)
    for ibin in range(args.nbins) :

        
            
        wmin=wavebins[ibin]
        wmax=wavebins[ibin+1]
        ymin=max(0,int(legval(u(wmin,wavemin,wavemax), ycoef[spec])))
        ymax=min(npix_y,int(legval(u(wmax,wavemin,wavemax), ycoef[spec])))

        print(spec,wmin,wmax,ymin,ymax)

        
        sw=np.zeros(x.size)
        swf=np.zeros(x.size)
        swm=np.zeros(x.size)
        


        for y in range(ymin,ymax) :
            prof=flux[y, x1_of_y[y]:x2_of_y[y]]
            ivar=flux_ivar[y, x1_of_y[y]:x2_of_y[y]]
            dx=np.arange(x1_of_y[y],x2_of_y[y])-xc_of_y[y]
            model=np.exp(-dx**2/2./sigma_of_y[y]**2)
            model *= np.sum(prof)/np.sum(model)
            
            tmp_w=np.interp(x,dx,ivar)
            tmp_wf=np.interp(x,dx,ivar*prof)
            tmp_wm=np.interp(x,dx,ivar*model)
            sum_tmp_w=np.sum(tmp_w)
            if sum_tmp_w==0 :
                continue
            mflux=np.sum(tmp_wf)/sum_tmp_w
            if mflux>10 :
                sw+=tmp_w
                swf+=tmp_wf
                swm+=tmp_wm
        
        prof=swf/(sw+(sw==0))
        mprof=swm/(sw+(sw==0))
        
        sprof=np.sum(prof)
        print(spec,ibin,"sprof=",sprof)
        if sprof>0 :
            a[ibin,spec].plot(x,prof/sprof,c="b")
            a[ibin,spec].plot(x,mprof/sprof,c="r")
            
        #a[ibin,spec].grid()
f.subplots_adjust(hspace=0)    
f.subplots_adjust(vspace=0)    
#pylab.setp([aa.get_xticklabels() for aa in a[0, :]], visible=False)
#pylab.setp([aa.get_yticklabels() for aa in a[:, 1]], visible=False)
pylab.show()

