#!/usr/bin/env python

import sys
import os
import fitsio
import numpy as np
import argparse
from desispec.preproc import _parse_sec_keyword
from desispec.calibfinder import CalibFinder
from astropy.time import Time

def mean_and_rms(vals,gain=1.) :
    vals=vals.astype(float)*gain
    mean=np.median(vals)
    for i in range(2) :
        ok=np.abs(vals-mean)<3.5*3. # we fix the cutoff at 3.5 sigma, assuming a fixed rdnoise value of 3 electrons for stability purpose (we know that the noise is between 2.5 and 4.5)
        mean=np.mean(vals[ok])
    rms=np.std(vals[ok])
    return mean,rms

def mean(vals,gain=1.) :
    m,r = mean_and_rms(vals,gain)
    return m

def rms(vals,gain=1.) :
    m,r = mean_and_rms(vals,gain)
    return r


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ifile2', metavar='ifile', type=str, nargs='*',
                    help='input fits file(s)')
parser.add_argument('-i','--ifile', type = str, default = None, required = False, nargs="*",
                    help = 'input fits file(s)')
parser.add_argument('--camera',type = str, default = 0, required = True, 
                    help = 'camera b0,r7 ...')
#parser.add_argument('-k','--keys',type = str, required = False, nargs="*" , help = 'dump keys from headers',default=[])
parser.add_argument('-o','--outfile',type = str, default = None, required = False, help = 'output filename')
parser.add_argument('--nobias', action='store_true', help="do not do a bias correction")
parser.add_argument('--bias', type=str, required=False,default=None,help="use this bias instead of the default one")
parser.add_argument('--gradient', action='store_true', help="use difference of adjacent pixels for rms")
parser.add_argument('--flavor', type=str,required=False,default="zero", help="flavor")



args        = parser.parse_args()

filenames = args.ifile2
if args.ifile is not None :
   filenames += args.ifile 

line="EXPID MJD-OBS DATE-OBS "
for a,amp in enumerate(['A','B','C','D']) :
    line+="MEAN_COL_OVERSCAN_{} RMS_COL_OVERSCAN_{} MEAN_ROW_OVERSCAN_{} RMS_ROW_OVERSCAN_{} MEAN_CCD_{} RMS_CCD_{} ".format(amp,amp,amp,amp,amp,amp) 



#x = np.zeros((len(filenames),4*6)).astype(float)

xx=[]
for f,filename in enumerate(filenames) :
    print('reading {} hdu={}'.format(filename,args.camera))
    try: 
        pheader=fitsio.read_header(filename)
        if "flavor" in pheader :
            if pheader["flavor"].lower().strip() != args.flavor :
                print("ignore image with flavor='{}'".format(pheader["flavor"]))
                continue
        else :
            print("warning: no flavor keyword in header")
        
        img,header=fitsio.read(filename,args.camera,header=True)
        try :
            cfinder=CalibFinder([header,pheader])
        except :
            print("warning: could not find calib")
            cfinder=None
    
    except OSError :
        print("# failed to read",filename)
        continue
    if "EXPID" in pheader :
        expid = pheader["EXPID"]
    else :
        expid = 0.
    if "MJD-OBS" in pheader :
        mjdobs = pheader["MJD-OBS"]
    else :
        mjdobs = 0.
    if "DATE-OBS" in pheader :
        try : 
            dateobs = Time(pheader["DATE-OBS"]).mjd
        except :
            dateobs = 0.
    else :
        dateobs = 0.
    img=img.astype(float)
    sub = None
    if not args.nobias :
        filename=args.bias
        if filename is None and cfinder is not None and cfinder.haskey("BIAS") :
            filename=cfinder.findfile("BIAS")
        if filename is not None :
            print("subtracting bias",filename)
            bias=fitsio.read(filename)
            sub = img - bias
    if sub is None :
        sub = img # we don't do bias subtraction
           
    rms_scale = 1.
    if args.gradient :
        tmp = sub[:,1:] - sub[:,:-1]
        sub[:,1:] = tmp
        sub[:,0]  = 0
        rms_scale = 1./np.sqrt(2.)
    i=0
    x=np.zeros(3+4*6).astype(float)
    
    x[i] = expid ; i += 1
    x[i] = mjdobs ; i += 1
    x[i] = dateobs ; i += 1

    for a,amp in enumerate(['A','B','C','D']) :
        gain = 1.
        if cfinder is not None and cfinder.haskey("GAIN"+amp) :
            gain=cfinder.value("GAIN"+amp)
        x[i] = mean(img[_parse_sec_keyword(header["BIASSEC"+amp])],gain=gain); i+=1
        x[i] = rms_scale*rms(sub[_parse_sec_keyword(header["BIASSEC"+amp])],gain=gain); i+=1
        if "ORSEC"+amp in header :
            x[i] = mean(img[_parse_sec_keyword(header["ORSEC"+amp])],gain=gain); i+=1
            x[i] = rms_scale*rms(sub[_parse_sec_keyword(header["ORSEC"+amp])],gain=gain); i+=1
        else :
            i += 2
        x[i] = mean(img[_parse_sec_keyword(header["DATASEC"+amp])],gain=gain); i+=1
        x[i] = rms_scale*rms(sub[_parse_sec_keyword(header["DATASEC"+amp])],gain=gain); i+=1
        print("assuming gain= {:3.2f} for amp {}, ccd rms= {:3.2f}".format(gain,amp,x[i-1]))
        sys.stdout.flush()
    xx.append(x)

if args.outfile is not None :
    xx=np.vstack(xx)
    np.savetxt(args.outfile,xx,header=line)
    print("wrote",args.outfile)


