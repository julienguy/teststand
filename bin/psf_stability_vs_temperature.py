#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import string
import os.path
from desispec.log import get_logger

def read_temperatures(filename) :
    x=np.loadtxt(filename).T
    file=open(filename)
    line=file.readlines()[0]
    file.close()
    keys=line.strip().replace("# ","").split(" ")
    vals={}
    for i,k in enumerate(keys) :
        vals[k]=x[i]
    if "DAY" in vals :
        vals["DAY"]=vals["DAY"].astype(int)
    if "EXPNUM" in vals :
        vals["EXPNUM"]=vals["EXPNUM"].astype(int)
    return keys,vals # to keep ordering we return the original key list

def read_psf_properties(filename) :
    x=np.loadtxt(filename).T
    file=open(filename)
    line=file.readlines()[0]
    file.close()
    keys=line.strip().replace("# ","").split(" ")
    vals={}
    for i,k in enumerate(keys) :
        vals[k]=x[i]
    if "EXPNUM" in vals :
        vals["EXPNUM"]=vals["EXPNUM"].astype(int)
    return keys,vals # to keep ordering we return the original key list

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psfprop', type = str, nargs = "*", default = None, required = True,
                    help = 'path of psf properties ASCII files')
parser.add_argument('--temp', type = str, default = None, required = True,
                    help = 'path to temperature file')

parser.add_argument('-o','--output', type = str, default = None, required = False, help = 'path to output ascii file')


args        = parser.parse_args()
log = get_logger()

temp_keys , temps = read_temperatures(args.temp)

ofile = open(args.output,"w")

text = '''## J. Guy 2017/03/27
## Version 1
## Code : psf_stability_vs_temperature.py --psfprop propertie-psf-?1-* --temp modified-temperatures.txt -o psf_stability_vs_temperature_version_1.txt
##
## The file gives a list of PSF properties per exposure and several environmental parameters.
## The list of keys per line of this ASCII file is given at the end of this description.
##
## Details on some parameters:
## EXPNUM is exposure number
## CAMID = 0,1,2 for b1(BLUE), r1(RED), z1(NIR)
## FIBER goes from 0 to 19 (skipping the broken fiber)
## orientation : X is the cross-dispersion direction (increasing X with increasing fiber number)
##               Y is the wavelength dispersion direction (increasing Y with increasing wavelength)
## DX,DY are trace coordinate offsets along X,Y with respect to an arbitrary average for fiber FIBER and wavelength WAVE
## CX,CY are just for debugging (centroid of PSF in the PSF stamp used internally for computation of PSF shape)
## SX,SY are sigma values along X,Y defined by : SX = sum((X-meanX)**2*PSF)/sum(PSF) , with meanX=sum(X*PSF)/sum(PSF)
## EBIAS is the bias on emission line flux measurement EBIAS = sum(meanPSF*PSF)/sum(PSF**2)-1 (function of SX and SY for a gaussian, but not a gaussian)
## DAY is YYYYMMDD (convenient to find back raw data directory)
## TIME in second is elapsed time with arbitrary date
## HOUR is hour during the day at WINLIGHT (i.e. < 24)
## EXPREQ is requested exp. time (a few seconds)
## BLUTEMP REDTEMP NIRTEMP temperatures in deg. C measured on spectrograph bench (from table in HDU 'SPECTCON1')
## EXPSHUT HARTL HARTR BLUHUMID REDHUMID NIRHUMID (from table in HDU 'SPECTCON1')
## PLCTEMP1 PLCTEMP2 temperatures in deg. C measured on injection bench (from header in HDU 'PLC')
## RCCDTEMP1 RCCDTEMP2 temperatures in deg. C measured on R CCD (keywords CCDTEMP1 CCDTEMP2 from header in HDU 'R1')
## ZCCDTEMP1 ZCCDTEMP2 temperatures in deg. C measured on Z CCD (keywords CCDTEMP1 CCDTEMP2 from header in HDU 'Z1')
## MXXXTEMP is a spline fit of temperature XXXTEMP with time to reduce statistical noise and glitches
## DHHXXXTEMP is a "delayed" temperature based on MXXXTEMP, with a characteristic time of HH hours (from 1 to 5 hours)
##
## list of keys: 
'''

ofile.write(text)

first = True
for filename in args.psfprop :
    print("reading %s"%filename)
    psf_keys , props = read_psf_properties(filename)

    expnum=props["EXPNUM"][0]
    temp_index=np.where(temps["EXPNUM"]==expnum)[0]
    if temp_index.size==0 :
        log.warning("didn't find temperature info for EXPNUM=%d"%props["EXPNUM"])
        continue
    temp_index=temp_index[0]
    
    if first :
        line="#"
        for k in psf_keys :
            line=line+" %s"%k
        for k in temp_keys :
            line=line+" %s"%k
        ofile.write("%s\n"%line)
    first=False

    line=""
    for k in psf_keys :
        line += " %s"%np.mean(props[k])
    for k in temp_keys :
        line += " %s"%temps[k][temp_index]
    ofile.write("%s\n"%line)
ofile.close()
