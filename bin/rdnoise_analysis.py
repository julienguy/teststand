#!/usr/bin/env python

import sys
import os
import fitsio
import numpy as np
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('ifile2', metavar='ifile', type=str, nargs='*',
                    help='input fits file(s)')
parser.add_argument('-i','--ifile', type = str, default = None, required = False, nargs="*",
                    help = 'input fits file(s)')
parser.add_argument('--hdu',type = str, default = 0, required = False, 
                    help = 'header HDU (int or string)')
parser.add_argument('-k','--keys',type = str, required = False, nargs="*" , help = 'dump keys from headers',default=[])
parser.add_argument('-b','--bias',type = str, required = False, default=None , help = 'bias')


args        = parser.parse_args()

try :
    hdu = int(args.hdu)
except ValueError:
    hdu = args.hdu
if hdu is None :
    hdu = 0

filenames = args.ifile2
if args.ifile is not None :
   filenames += args.ifile 

if len(sys.argv)<2 :

    print(sys.argv[0],"image.fits")

line="#"
for k in args.keys :
    line+=" "+k
line+=" MEANA MEANB MEANC MEAND"
line+=" RMSA RMSB RMSC RMSD"
line+=" OSRMSA OSRMSB OSRMSC OSRMSD"
line+=" (electrons)"
print(line)

bias=None
if args.bias :
    bias = fitsio.read(args.bias)

for filename in filenames :
    #print('reading {} hdu={}'.format(filename,hdu))
    try: 
        pheader=fitsio.read_header(filename)
        img,header=fitsio.read(filename,hdu,header=True)
    except OSError :
        print("# failed to read",filename)
        continue
    overscan_rdnoise=[]
    central_rdnoise=[]
    vals=[]
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

