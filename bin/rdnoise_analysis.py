#!/usr/bin/env python

import sys
import os
import fitsio
import numpy as np
import argparse
from desispec.preproc import _parse_sec_keyword
from desispec.calibfinder import CalibFinder

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
parser.add_argument('-o','--outfile',type = str, required = True, help = 'output filename')
parser.add_argument('--nobias', action='store_true', help="do not do a bias correction")
parser.add_argument('--gradient', action='store_true', help="use difference of adjacent pixels for rms")
parser.add_argument('--flavor', type=str,required=False,default="zero", help="flavor")



args        = parser.parse_args()

filenames = args.ifile2
if args.ifile is not None :
   filenames += args.ifile 

line="EXPID "
for a,amp in enumerate(['A','B','C','D']) :
    line+="MEAN_COL_OVERSCAN_{} RMS_COL_OVERSCAN_{} MEAN_ROW_OVERSCAN_{} RMS_ROW_OVERSCAN_{} MEAN_CCD_{} RMS_CCD_{} ".format(amp,amp,amp,amp,amp,amp) 



#x = np.zeros((len(filenames),4*6)).astype(float)

xx=[]
for f,filename in enumerate(filenames) :
    print('reading {} hdu={}'.format(filename,args.camera))
    try: 
        pheader=fitsio.read_header(filename)
        if pheader["flavor"].lower().strip() != args.flavor :
            print("ignore image with flavor='{}'".format(pheader["flavor"]))
            continue
                  
        img,header=fitsio.read(filename,args.camera,header=True)
        cfinder=CalibFinder([header,pheader])
        
    except OSError :
        print("# failed to read",filename)
        continue
        
    expid=pheader["EXPID"]
    
    img=img.astype(float)
    sub = None
    if not args.nobias :
        if cfinder.haskey("BIAS") :
            filename=cfinder.findfile("BIAS")
            print("subtractiong bias",filename)
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
    x=np.zeros(1+4*6).astype(float)
    
    x[i] = expid ; i += 1


    for a,amp in enumerate(['A','B','C','D']) :
        gain = 1.
        if cfinder.haskey("GAIN"+amp) :
            gain=cfinder.value("GAIN"+amp)
            print("assuming gain for amp {} = {}".format(amp,gain))
        x[i] = mean(img[_parse_sec_keyword(header["BIASSEC"+amp])],gain=gain); i+=1
        x[i] = rms_scale*rms(sub[_parse_sec_keyword(header["BIASSEC"+amp])],gain=gain); i+=1
        x[i] = mean(img[_parse_sec_keyword(header["ORSEC"+amp])],gain=gain); i+=1
        x[i] = rms_scale*rms(sub[_parse_sec_keyword(header["ORSEC"+amp])],gain=gain); i+=1
        x[i] = mean(img[_parse_sec_keyword(header["CCDSEC"+amp])],gain=gain); i+=1
        x[i] = rms_scale*rms(sub[_parse_sec_keyword(header["CCDSEC"+amp])],gain=gain); i+=1
        print("ccd rms {} = {:3.2f}".format(amp,x[1+a*4+5]))
        sys.stdout.flush()
    xx.append(x)
xx=np.vstack(xx)

np.savetxt(args.outfile,xx,header=line)
print("wrote",args.outfile)

"""

    
        #ii=_parse_sec_keyword(header['CCDSEC%d'%amp0])
    sys.exit(12)
    
    n0=img.shape[0]
    n1=img.shape[1]
    m=200
    #print("stamp size={}x{}".format(n0//2-2*m,n1//2-2*m))
    if bias is not None :
        img = img.astype(float)
        img -= bias
    for amp in 'ABCD' :
        if amp=='A' :
            quad=img[m:n0//2-m,m:n1//2-m]
        elif amp=='B' :
            quad=img[m:n0//2-m,n1//2+m:-m]
        elif amp=='C':
            quad=img[n0//2+m:-m,m:n1//2-m]
        elif amp=='D' :
            quad=img[n0//2+m:-m,n1//2+m:-m]
        shape=quad.shape
        quad=quad.astype(float)
        med=np.median(quad)
        mean=np.mean(quad[np.abs(quad-med)<20])
        vals.append(mean)
        overscan_rdnoise.append(header["OBSRDN"+amp])
        rms=np.std(quad[np.abs(quad-mean)<20])
        central_rdnoise.append(rms)


    line=""
    for k in args.keys :
        if k in pheader :
            line += " "+str(pheader[k])
        elif k in header :
            line += " "+str(header[k])
        else :
            line += " 0"
    for i in range(4) :
        line+=" {}".format(vals[i])
    for i in range(4) :
        line+=" {}".format(central_rdnoise[i])
    for i in range(4) :
        line+=" {}".format(overscan_rdnoise[i])
    print(line)
    sys.stdout.flush()

"""
