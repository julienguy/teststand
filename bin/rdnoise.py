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
line+=" RMS1 RMS2 RMS3 RMS4 (ADUs)"
print(line)

for filename in filenames :
    #print('reading {} hdu={}'.format(filename,hdu))
    try: 
        pheader=fitsio.read_header(filename)
        img,header=fitsio.read(filename,hdu,header=True)
    except OSError :
        print("# failed to read",filename)
        continue
    vals=[]
    n0=img.shape[0]
    n1=img.shape[1]
    m=600
    #print("stamp size={}x{}".format(n0//2-2*m,n1//2-2*m))
    
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
        
        # remove smooth gradients
        x=np.arange(quad.shape[1])
        prof=np.median(quad,axis=0)
        prof=np.poly1d(np.polyfit(x,prof,1))(x)
        for i in range(prof.size) :
            quad[:,i] -= prof[i]
            
        x=np.arange(quad.shape[0])
        prof=np.median(quad,axis=1)
        prof=np.poly1d(np.polyfit(x,prof,1))(x)
        for i in range(prof.size) :
            quad[i,:] -= prof[i]
        
        quad=quad.ravel()
        med=np.median(quad)
        quad -= med
        rms=1.4826*np.median(np.abs(quad))
        nsig=2.7
        rms=np.std(quad[np.abs(quad)<nsig*rms])
        #print(amp,rms)
        vals.append(rms)

        if False :
            import matplotlib.pyplot as plt
            print("rms=",rms)
            

            tmp=quad.copy()
            tmp[tmp>nsig*rms]=0
            tmp[tmp<-nsig*rms]=0
            tmp=tmp.reshape(shape)
            
            plt.figure()
            plt.imshow(tmp,origin=0)
            plt.colorbar()
            
            plt.figure()
            bins=np.linspace(-nsig*rms,nsig*rms,200)
            h,bins=np.histogram(tmp[tmp!=0],bins=bins)
            x=bins[:-1]+(bins[1]-bins[0])/2.
            i=(h>0)
            x=x[i]
            h=h[i]
            plt.plot(x,h,"o-")
            y=np.exp(-x**2/2/rms**2)
            y *= np.sum(h)/np.sum(y)
            plt.plot(x,y,"-")
            plt.show()



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
    print(line)
    sys.stdout.flush()

