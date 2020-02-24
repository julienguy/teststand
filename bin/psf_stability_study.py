#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import specter.psf
import sys
import argparse
import string
import os.path
from desispec.qproc.util import parse_fibers
from desiutil.log                  import get_logger

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
parser.add_argument('--psf', type = str, nargs = "*", default = None, required = True,
                    help = 'path of psf files')
parser.add_argument('--fibers', type = str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output ascii file')
parser.add_argument('--plot', action='store_true',help="plot result")


args        = parser.parse_args()
log = get_logger()
psfs=[]
for filename in args.psf :
    log.info("reading %s"%filename)
    psfs.append(readpsf(filename))

wmin=psfs[0]._wmin_all
wmax=psfs[0]._wmax_all
nw=20
waves=np.linspace(wmin+(wmax-wmin)/nw/2.,wmax-(wmax-wmin)/nw/2.,nw)
fibers=parse_fibers(args.fibers)
if fibers is None :
        fibers = np.arange(psfs[0].nspec)


res_y=[]
res_x=[]
res_emission_line_rms=[]
res_continuum_rms=[]
res_fiber=[]
res_wave=[]
res_x_rms=[]
res_y_rms=[]

use_trace = True

for fiber in fibers :
    for wave in waves :
        images = []
        i0 = []
        i1 = []
        xc = []
        yc = []
        for psf in psfs :
            xx, yy, ccdpix = psf.xypix(fiber,wave)
            xc.append(psf.x(fiber,wave))
            yc.append(psf.y(fiber,wave))
            #print("fiber,wave",fiber,wave)
            #print("xx",xx)
            #print("yy",yy)
            
            #plt.imshow(ccdpix)
            #plt.show()


            
            images.append(ccdpix)
            i1.append(xx.start)
            i0.append(yy.start)
        
        mi0 = int(np.min(i0))
        mi1 = int(np.min(i1))
        n=len(images)
        # add a margin and apply offset if necessary
        nimages=np.zeros((n,images[0].shape[0]+3,images[0].shape[1]+3))
        for j in range(n) :
            nimages[j,i0[j]-mi0:i0[j]-mi0+images[j].shape[0],i1[j]-mi1:i1[j]-mi1+images[j].shape[1]] = images[j]
        
        images=nimages
        #log.info(images.shape)
        n=images.shape[0]
        mimage=np.mean(images,axis=0)
        x=np.tile(np.arange(mimage.shape[0]),(mimage.shape[1],1))    # check with visual inspection of ccd image 
        y=np.tile(np.arange(mimage.shape[1]),(mimage.shape[0],1)).T  # y is wavelength axis
        
        delta_ratio_emission_line = np.zeros(n)
        delta_ratio_continuum = np.zeros(n)
        if not use_trace :
            xc = np.zeros(n)
            yc = np.zeros(n)
        for j in range(n) :
            delta_ratio_emission_line[j] = np.sum(images[j]*mimage)/np.sum(images[j]**2)-1
            pmimage=np.sum(mimage,axis=0) # projection to get 1D PSF along cross-dispersion for continuum fit normalization
            pimage=np.sum(images[j],axis=0) 
            delta_ratio_continuum[j] = np.sum(pimage*pmimage)/np.sum(pimage**2)-1
            use_trace = False # checked has no impact
            if not use_trace :
                xc[j] = np.sum(x*images[j])/np.sum(images[j])
                yc[j] = np.sum(y*images[j])/np.sum(images[j])
            
            
            
            if False and np.abs(delta_ratio_emission_line[j])>0.008 :
                plt.figure()
                a=plt.subplot(1,2,1)
                plt.imshow(mimage,origin=0,interpolation="nearest")
                a=plt.subplot(1,2,2)
                plt.imshow(images[j],origin=0,interpolation="nearest")
                log.warning("delta_ratio_emission_line fiber=%d wave=%d img=%d = %f"%(fiber,wave,j,delta_ratio_emission_line[j]))
                plt.show()
        
        rms2d=np.sqrt(np.sum(delta_ratio_emission_line**2))*np.sqrt(n/(n-1))
        rms1d=np.sqrt(np.sum(delta_ratio_continuum**2))*np.sqrt(n/(n-1))
        log.info("fiber=%d wave=%d rms2d=%f rms1d=%f"%(fiber,wave,rms2d,rms1d))
        res_y.append(mi0)
        res_x.append(mi1)
        res_emission_line_rms.append(rms2d)
        res_continuum_rms.append(rms1d)
        res_x_rms.append(np.std(xc))
        res_y_rms.append(np.std(yc))
        res_fiber.append(fiber)
        res_wave.append(wave)
        
res_x=np.array(res_x)
res_y=np.array(res_y)
res_emission_line_rms=np.array(res_emission_line_rms)
res_continuum_rms=np.array(res_continuum_rms)
res_fiber=np.array(res_fiber)
res_wave=np.array(res_wave)
res_x_rms=np.array(res_x_rms)
res_y_rms=np.array(res_y_rms)

if args.output :
    file=open(args.output,"w")
    file.write("# fiber wavelength rms_emission_line_flux rms_continuum_flux xccd yccd\n")
    for i in range(res_x.size) :
        file.write("%d %f %f %f %f %f\n"%(res_fiber[i],res_wave[i],res_emission_line_rms[i],res_continuum_rms[i],res_x[i],res_y[i]))
    file.close()
                

if args.plot :
    plt.figure()
    plt.subplot(2,2,1)
    for fiber in fibers :
        plt.plot(res_wave[res_fiber==fiber],res_emission_line_rms[res_fiber==fiber],"-")
    plt.xlabel("Wavelength")
    plt.ylabel("flux error for emission lines")

    plt.subplot(2,2,2)
    for fiber in fibers :
        plt.plot(res_wave[res_fiber==fiber],res_continuum_rms[res_fiber==fiber],"-")
    plt.xlabel("Wavelength")
    plt.ylabel("flux error for continuum")

    plt.subplot(2,2,3)
    for fiber in fibers :
        plt.plot(res_wave[res_fiber==fiber],res_x_rms[res_fiber==fiber],"-")
    plt.xlabel("Wavelength")
    plt.ylabel("x error (rms, Pixels)")

    plt.subplot(2,2,4)
    for fiber in fibers :
        plt.plot(res_wave[res_fiber==fiber],res_y_rms[res_fiber==fiber],"-")
    plt.xlabel("Wavelength")
    plt.ylabel("y error (rms, Pixels)")

    plt.show()



