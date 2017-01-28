#!/usr/bin/env python

import astropy.io.fits as pyfits
import argparse
import re
import sys
import numpy as np

def parse_sec_keyword(value):
    m = re.search('\[(\d+):(\d+)\,(\d+):(\d+)\]', value)
    if m is None:
        raise ValueError, 'unable to parse {} as [a:b, c:d]'.format(value)
    return  map(int, m.groups())

def parse_date_obs(value):
    m = re.search('(\d+)-(\d+)-(\d+)T', value)
    if m is None:
        raise ValueError, 'unable to parse {}'.format(value)
    vals=map(int, m.groups())
    yy=vals[0]
    mm=vals[1]
    dd=vals[2]
    return yy*10000+mm*100+dd

parser=argparse.ArgumentParser(description="Reformat DESI raw data")
parser.add_argument('-i','--infile', type = str, default = None, required=True,
                    help = 'path to input DESI raw data file')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,
                    help = 'path to output DESI raw data file')

args = parser.parse_args()

print("Will read %s, reformat it and write it to %s"%(args.infile,args.outfile))
print("This is a minimal reformatting of ICS raw data to adapt to DESI soft")

ffile=pyfits.open(args.infile)
ffile.info()


header=ffile[0].header
try :
    date_obs=header["DATE-OBS"]
except KeyError :
    print("no DATE-OBS, cannot guess ICS version, quit")
    sys.exit(1)
date=parse_date_obs(date_obs)
print("Observation date : %d"%date)
ICS_VERSION=0
if date <= 20161027    : ICS_VERSION=0
elif  date <= 20161107 : ICS_VERSION=1
elif  date <= 20161110 : ICS_VERSION=2
elif  date <= 20161123 : ICS_VERSION=3
else                   : ICS_VERSION=4
print("Guessed ICS_VERSION : %d"%ICS_VERSION)

if ICS_VERSION==3 :
    print("Do not do anything for ICS_VERSION>=3")
    sys.exit(0)
if ICS_VERSION==1 :
    print("ERROR : I haven't implemented the formatting for ICS_VERSION=1")
    sys.exit(1)

extnames=[]
cameras=[]
if ICS_VERSION<=2 :
    extnames=["CCDS1R"]
    cameras=["R1"]
elif ICS_VERSION==3 :
    # but code won't do anything
    extnames=["R1"]
    cameras=extnames
elif ICS_VERSION==4 :
    extnames=["B1","R1","Z1"]
    cameras=["b1","r1","z1"]
else :
    print("ERROR : don't know about ICS_VERSION=%d"%ICS_VERSION)
    sys.exit(0)

for extname,camera in zip(extnames,cameras) :
    print("Looking at extname %s"% extname)
    hdu=ffile[extname]
    header=hdu.header
    if ICS_VERSION<=2 :
        header["CAMERA"]=camera
        header["EXTNAME"]=camera
    
    if ICS_VERSION==0 and camera=="R1" :
        print("WARNING : changing DATASEC,CCDSEC,BIASEC for",extname)
        header["DATASEC1"]='[10:2064,4:2065]'
        d1xmin,d1xmax,d1ymin,d1ymax = parse_sec_keyword(header["DATASEC1"])
        header["PRESEC1"]='[1:%d,%d:%d]'%(d1xmin-1,d1ymin,d1ymax)
        header["CCDSEC1"]='[1:%d,1:%d]'%(d1xmax-d1xmin+1,d1ymax-d1ymin+1)
        header["BIASSEC1"]='[%d:2121,%d:%d]'%(d1xmax+1,d1ymin,d1ymax)
        c1xmin,c1xmax,c1ymin,c1ymax = parse_sec_keyword(header["CCDSEC1"])
        p1xmin,p1xmax,p1ymin,p1ymax = parse_sec_keyword(header["PRESEC1"])
        b1xmin,b1xmax,b1ymin,b1ymax = parse_sec_keyword(header["BIASSEC1"])

        header["DATASEC2"]='[2179:4230,%d:%d]'%(d1ymin,d1ymax)
        d2xmin,d2xmax,d2ymin,d2ymax = parse_sec_keyword(header["DATASEC2"])
        header["PRESEC2"]='[%d:%d,%d:%d]'%(d2xmax,header["NAXIS1"],d2ymin,d2ymax)    
        header["BIASSEC2"]='[%d:%d,%d:%d]'%(b1xmax+1,d2xmin-1,d2ymin,d2ymax)    
        header["CCDSEC2"]='[%d:%d,%d:%d]'%(c1xmax+1,c1xmax+1+d2xmax-d2xmin,1,d2ymax-d2ymin+1)    
        c2xmin,c2xmax,c2ymin,c2ymax = parse_sec_keyword(header["CCDSEC2"])
        p2xmin,p2xmax,p2ymin,p2ymax = parse_sec_keyword(header["PRESEC2"])
        b2xmin,b2xmax,b2ymin,b2ymax = parse_sec_keyword(header["BIASSEC2"])

        header["DATASEC3"]='[%d:%d,2090:4152]'%(d1xmin,d1xmax)
        d3xmin,d3xmax,d3ymin,d3ymax = parse_sec_keyword(header["DATASEC3"])    
        header["PRESEC3"]='[%d:%d,%d:%d]'%(p1xmin,p1xmax,d3ymin,d3ymax)
        header["BIASSEC3"]='[%d:%d,%d:%d]'%(b1xmin,b1xmax,d3ymin,d3ymax)
        header["CCDSEC3"]='[%d:%d,%d:%d]'%(c1xmin,c1xmax,c1ymax+1,c1ymax+1+d3ymax-d3ymin)
        c3xmin,c3xmax,c3ymin,c3ymax = parse_sec_keyword(header["CCDSEC3"])

        header["DATASEC4"]='[%d:%d,%d:%d]'%(d2xmin,d2xmax,d3ymin,d3ymax)
        header["PRESEC4"]='[%d:%d,%d:%d]'%(p2xmin,p2xmax,d3ymin,d3ymax)
        header["BIASSEC4"]='[%d:%d,%d:%d]'%(b2xmin,b2xmax,d3ymin,d3ymax)
        header["CCDSEC4"]='[%d:%d,%d:%d]'%(c2xmin,c2xmax,c3ymin,c3ymax)
    elif ICS_VERSION==1 :
        print("ERROR : not implemented")
    else :
        print("do not change data sec")


    print("WARNING : changing gains")
    reference_gain=1.
    header["GAIN1"]=reference_gain
    header["GAIN2"]=reference_gain
    header["GAIN3"]=reference_gain
    header["GAIN4"]=reference_gain
    

    if ICS_VERSION<=1 and camera=="R1" :

        print("WARNING : Flip image along Y for camera %s and ICS VERSION=%d"%(camera,ICS_VERSION))
        hdu.data = hdu.data[::-1,:]
        print("now modify all SEC info")
        xmin,xmax,ymin,ymax = parse_sec_keyword(header["CCDSEC4"])        
        ny_ccd=ymax
        ny_input=header["NAXIS2"]

        header_copy=header.copy()
        for amp in range(1,5) :
            for sec in ["PRESEC","DATASEC","BIASSEC","CCDSEC"] :
                key="%s%d"%(sec,amp)
                xmin,xmax,ymin,ymax = parse_sec_keyword(header_copy[key])
                if sec == "CCDSEC" :
                    flipped_ymax=ny_ccd-ymin+1
                    flipped_ymin=ny_ccd-ymax+1
                else :
                    flipped_ymax=ny_input-ymin+1
                    flipped_ymin=ny_input-ymax+1

                key="%s%d"%(sec,amp) # DO NOT CHANGE AMP NAMES, CONFUSING
                old=header[key]
                header[key]='[%d:%d,%d:%d]'%(xmin,xmax,flipped_ymin,flipped_ymax)
                print("%s %s -> %s"%(key,old,header[key]))

    


    if ICS_VERSION==4 and camera.lower() == "b1" :
         print("WARNING : rearrange amplifiers  for camera=%s and ICS VERSION=%d"%(camera,ICS_VERSION))
         pixels=hdu.data
         nx=pixels.shape[1]
         ny=pixels.shape[0]
         # get amps
         # defined as 
         # A B
         # C D
         C=pixels[:ny/2,:nx/2]
         D=pixels[:ny/2,nx/2:]
         A=pixels[ny/2:,:nx/2]
         B=pixels[ny/2:,nx/2:]
         
         # flip B and C, both axes
         B=B[::-1,::-1]
         C=C[::-1,::-1]
         
         # reassemble to get
         # A C
         # B D
         AC=np.concatenate((A,C),axis=1)
         BD=np.concatenate((B,D),axis=1)
         hdu.data=np.concatenate((BD,AC),axis=0)
         
         print("WARNING, still need to figure out wavelength direction, fiber direction for amps C and D for camera=%s and ICS VERSION=%d"%(camera,ICS_VERSION))
         print("WARNING, WOULD NEED TO CHANGE NAMES OF AMPLIFIERS")
         
         
         
            
         
         

    if ICS_VERSION==4 and camera.lower() == "z1" :
         print("WARNING : flip along X (fibers) for camera=%s and ICS VERSION=%d"%(camera,ICS_VERSION))
         hdu.data=hdu.data[:,::-1]
         print("WARNING, WOULD NEED TO CHANGE NAMES OF AMPLIFIERS")

ffile.writeto(args.outfile,clobber=True)
print("wrote",args.outfile)
print("done")
