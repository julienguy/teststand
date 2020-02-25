#!/usr/bin/env python

import sys,os
import numpy as np

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import desimodel.io
from specter.throughput import ObjType

# example: plot_calib.py fluxcalib-*-00044565.fits

cmap = plt.get_cmap("tab10")
colors = cmap(np.linspace(0,1,10))
labels=[]

count=0

for filename in sys.argv[1:] :

    
    h=pyfits.open(filename)
    head=h[0].header
    if "CAMERA" in head :
        cam = head["CAMERA"]
        band=cam[0]
        if len(cam)>1 :
            s=int(cam[1])
            label="SP{}".format(int(cam[1]))
            count=s
        else :
            label=None
            count = (count+1)%10
    else :
        label=None
        count = (count+1)%10
        
    color = colors[count]
        
    
    sim = "DOSVER" in head and head["DOSVER"].strip()=="SIM"
    
    
    
    wave=h["WAVELENGTH"].data
    cal=h["FLUXCALIB"].data
    if len(cal.shape)>1 :
        cal=cal[cal.shape[0]//2] # electrons/A  /  (1e-17  ergs/s/cm2/A)
    cal *= 1e17 # electrons/A  /  ( ergs/s/cm2/A)
    if "EXPTIME" in head :
        exptime=head["EXPTIME"]
    else :
        exptime=1. # average calib vector
    cal /= exptime # electrons  /  (ergs/cm2)
    # reference effective collection area
    area = 8.678709421*1e4 # cm2
    cal /= area # electrons  /  ergs
    hplanck = 6.62606957e-34 #J.s
    hplanck *= 1e7 # erg.s
    cspeed = 2.99792458e8 # m/s
    cspeed *= 1e10 # A/s
    energy = hplanck*cspeed/wave # erg ( erg.s * A/s / A)
    cal *= energy # electrons /photons

    
    if sim: label="SIM"
    
    if label in labels :
        label=None
    else :
        labels.append(label)

    plt.plot(wave,cal,label=label,color=color)


if True :
    wave=np.linspace(3500,10000.,(10000-3500+1))
    for camera in ['b', 'r', 'z']:
        thru = desimodel.io.load_throughput(camera)
        t=thru(wave,objtype=ObjType.STAR, airmass=1.1)
        ii=np.where(t>0)[0]
        b=ii[0]
        e=ii[-1]
        label=None
        if camera=="b" : label="DESIMODEL"
        plt.plot(wave[b:e],t[b:e],c="gray",alpha=0.5,label=label)

        
plt.legend(loc="upper left")
plt.xlabel("Wavelength (A)")
plt.ylabel("Throughput (inc. fiber aperture loss)")

plt.show()


