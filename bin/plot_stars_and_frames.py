#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import astropy.io.fits as pyfits
import astropy.units  as units
import matplotlib.pyplot as plt

from desispec.io import read_frame
from desispec.interpolation import resample_flux
from desispec.io.filters import load_legacy_survey_filter


band="R"
fluxunits = 1e-17 * units.erg / units.s / units.cm**2 / units.Angstrom

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="Display spectra, looping over targets if targetid not set, and optionally show best fit from redrock"
)
parser.add_argument('--cframes', type = str, default = None, required = True, nargs="*", 
                    help = 'path to cframe fits files')

parser.add_argument('--stdstars', type = str, default = None, required = True, 
                    help = 'path to stdstars fits files')

args        = parser.parse_args()

stars_filename = args.stdstars #"stdstars-0-00051001.fits"
frame_filenames = args.cframes #["cframe-r0-00051001.fits","cframe-b0-00051001.fits","cframe-z0-00051001.fits"]

h=pyfits.open(stars_filename)
h.info()
fibers=h["FIBERS"].data
table=h["METADATA"].data
print("std stars fibers=",fibers)
model_wave = h["WAVELENGTH"].data
model_flux = h["FLUX"].data

frames=[]
for frame_filename in frame_filenames :
    frame = read_frame(frame_filename)
    selection=np.intersect1d(frame.fibermap["FIBER"],fibers)
    frame = frame[selection]
    frames.append(frame)
    

for i,fiber in enumerate(fibers) :
    j=np.where(frame.fibermap['FIBER']==fiber)[0][0]
    print("fiber={}, i={}, j={}".format(fiber,i,j))
    photsys = frame.fibermap['PHOTSYS'][j]
    filter_response=load_legacy_survey_filter("R",photsys)
    model_mag=filter_response.get_ab_magnitude(model_flux[i]*fluxunits,model_wave)
    fiber_mag=-2.5*np.log10(frame.fibermap['FLUX_R'][j])+22.5
    print("model mag={:4.2f} fiber mag={:4.2f}".format(model_mag,fiber_mag))
    a=0
    for frame in frames :
        mflux = resample_flux(frame.wave,model_wave,model_flux[i])
        rflux = frame.R[j].dot(mflux)
        plt.plot(frame.wave,frame.flux[j])
        plt.plot(frame.wave,rflux,c="k",alpha=0.6)
        if a==0 :
            a=np.sum(frame.flux[j]*rflux)/np.sum(rflux**2)
            print("scale a={}".format(a))
        plt.plot(frame.wave,rflux*a,c="gray",alpha=0.6)
        


        
        
    plt.show()

