#!/usr/bin/env python

import astropy.io.fits as pyfits
import argparse
import re
import sys

def parse_sec_keyword(value):
    m = re.search('\[(\d+):(\d+)\,(\d+):(\d+)\]', value)
    if m is None:
        raise ValueError, 'unable to parse {} as [a:b, c:d]'.format(value)
    return  map(int, m.groups())


parser=argparse.ArgumentParser(description="Reformat DESI raw data")
parser.add_argument('-i','--infile', type = str, default = None, required=True,
                    help = 'path to input DESI raw data file')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,
                    help = 'path to output DESI raw data file')
parser.add_argument('-c','--camera', type = str, default = None, required=True,
                    help = 'camera')
parser.add_argument('--hdu', type = int, default = 0, required=False,
                    help = 'HDU to consider')

args = parser.parse_args()

ffile=pyfits.open(args.infile)
hdu=ffile[args.hdu]
header=hdu.header
header["EXTNAME"]=args.camera.upper()

if args.camera.lower().find("b")==0 :
    print "warning, from sims ..."
    header["GAIN1"]=1.0
    header["GAIN2"]=1.0
    header["GAIN3"]=1.0
    header["GAIN4"]=1.0
    header["PRESEC1"]='[1:4,1:2048]'
    header["DATASEC1"]='[5:2052,1:2048]'
    header["BIASSEC1"]='[2053:2102,1:2048]'
    header["CCDSEC1"]='[1:2048,1:2048]'
    header["PRESEC2"]='[4201:4204,1:2048]'
    header["DATASEC2"]='[2153:4200,1:2048]'
    header["BIASSEC2"]='[2103:2152,1:2048]'
    header["CCDSEC2"]='[2049:4096,1:2048]'
    header["PRESEC3"]='[1:4,2049:4096]'
    header["DATASEC3"]='[5:2052,2049:4096]'
    header["BIASSEC3"]='[2053:2102,2049:4096]'
    header["CCDSEC3"]='[1:2048,2049:4096]'
    header["PRESEC4"]='[4201:4204,2049:4096]'
    header["DATASEC4"]='[2153:4200,2049:4096]'
    header["BIASSEC4"]='[2103:2152,2049:4096]'
    header["CCDSEC4"]='[2049:4096,2049:4096]'
elif args.camera.lower().find("r")==0: 
    print "tuned to first teststand image"
    
    '''
    parse keywords like BIASSECB='[7:56,51:4146]' into python slices

    python and FITS have almost opposite conventions,
      * FITS 1-indexed vs. python 0-indexed
      * FITS upperlimit-inclusive vs. python upperlimit-exclusive
      * FITS[x,y] vs. python[y,x]

    i.e. BIASSEC2='[7:56,51:4146]' -> (slice(50,4146), slice(6,56))
    '''
    
    # [x:x,y:y]

    header["GAIN1"]=1.0
    header["GAIN2"]=1.0
    header["GAIN3"]=1.0
    header["GAIN4"]=1.0
    
    header["DATASEC1"]='[10:2064,4:2065]'
    d1xmin,d1xmax,d1ymin,d1ymax = parse_sec_keyword(header["DATASEC1"])
    header["PRESEC1"]='[1:%d,%d:%d]'%(d1xmin-1,d1ymin,d1ymax)
    header["CCDSEC1"]='[1:%d,1:%d]'%(d1xmax-d1xmin+1,d1ymax-d1ymin+1)
    header["BIASSEC1"]='[%d:2100,%d:%d]'%(d1xmax+1,d1ymin,d1ymax)
    c1xmin,c1xmax,c1ymin,c1ymax = parse_sec_keyword(header["CCDSEC1"])
    p1xmin,p1xmax,p1ymin,p1ymax = parse_sec_keyword(header["PRESEC1"])
    b1xmin,b1xmax,b1ymin,b1ymax = parse_sec_keyword(header["BIASSEC1"])

    header["DATASEC2"]='[2137:4187,%d:%d]'%(d1ymin,d1ymax)
    d2xmin,d2xmax,d2ymin,d2ymax = parse_sec_keyword(header["DATASEC2"])
    header["PRESEC2"]='[%d:%d,%d:%d]'%(d2xmax,header["NAXIS1"],d2ymin,d2ymax)    
    header["BIASSEC2"]='[%d:%d,%d:%d]'%(b1xmax+1,d2xmin-1,d2ymin,d2ymax)    
    header["CCDSEC2"]='[%d:%d,%d:%d]'%(c1xmax+1,c1xmax+1+d2xmax-d2xmin,1,d2ymax-d2ymin+1)    
    c2xmin,c2xmax,c2ymin,c2ymax = parse_sec_keyword(header["CCDSEC2"])
    p2xmin,p2xmax,p2ymin,p2ymax = parse_sec_keyword(header["PRESEC2"])
    b2xmin,b2xmax,b2ymin,b2ymax = parse_sec_keyword(header["BIASSEC2"])
    
    header["DATASEC3"]='[%d:%d,2136:4197]'%(d1xmin,d1xmax)
    d3xmin,d3xmax,d3ymin,d3ymax = parse_sec_keyword(header["DATASEC3"])    
    header["PRESEC3"]='[%d:%d,%d:%d]'%(p1xmin,p1xmax,d3ymin,d3ymax)
    header["BIASSEC3"]='[%d:%d,%d:%d]'%(b1xmin,b1xmax,d3ymin,d3ymax)
    header["CCDSEC3"]='[%d:%d,%d:%d]'%(c1xmin,c1xmax,c1ymax+1,c1ymax+1+d3ymax-d3ymin)
    c3xmin,c3xmax,c3ymin,c3ymax = parse_sec_keyword(header["CCDSEC3"])
    
    header["DATASEC4"]='[%d:%d,%d:%d]'%(d2xmin,d2xmax,d3ymin,d3ymax)
    header["PRESEC4"]='[%d:%d,%d:%d]'%(p2xmin,p2xmax,d3ymin,d3ymax)
    header["BIASSEC4"]='[%d:%d,%d:%d]'%(b2xmin,b2xmax,d3ymin,d3ymax)
    header["CCDSEC4"]='[%d:%d,%d:%d]'%(c2xmin,c2xmax,c3ymin,c3ymax)
    

    for k in ["DATASEC1","PRESEC1","CCDSEC1","BIASSEC1","DATASEC2","PRESEC2","CCDSEC2","BIASSEC2","DATASEC3","PRESEC3","CCDSEC3","BIASSEC3","DATASEC4","PRESEC4","CCDSEC4","BIASSEC4"]:
        print k,"=",header[k]
    
    if 0 :
        # from sims 
        header["GAIN1"]=1.0
        header["GAIN2"]=1.0
        header["GAIN3"]=1.0
        header["GAIN4"]=1.0
        header["PRESEC1"]='[1:7,1:2064]'
        header["DATASEC1"]='[8:2064,1:2064]'
        header["BIASSEC1"]='[2065:2114,1:2064]'
        header["CCDSEC1"]='[1:2057,1:2064]'
        header["PRESEC2"]='[4222:4228,1:2064]'
        header["DATASEC2"]='[2165:4221,1:2064]'
        header["BIASSEC2"]='[2115:2164,1:2064]'
        header["CCDSEC2"]='[2058:4114,1:2064]'
        header["PRESEC3"]='[1:7,2065:4128]'
        header["DATASEC3"]='[8:2064,2065:4128]'
        header["BIASSEC3"]='[2065:2114,2065:4128]'
        header["CCDSEC3"]='[1:2057,2065:4128]'
        header["PRESEC4"]='[4222:4228,2065:4128]'
        header["DATASEC4"]='[2165:4221,2065:4128]'
        header["BIASSEC4"]='[2115:2164,2065:4128]'
        header["CCDSEC4"]='[2058:4114,2065:4128]'


elif args.camera.lower().find("z")==0:
    print "warning, from sims ..."
    header["GAIN1"]=1.0
    header["GAIN2"]=1.0
    header["GAIN3"]=1.0
    header["GAIN4"]=1.0
    header["PRESEC1"]='[1:7,1:2064]'
    header["DATASEC1"]='[8:2064,1:2064]'
    header["BIASSEC1"]='[2065:2114,1:2064]'
    header["CCDSEC1"]='[1:2057,1:2064]'
    header["PRESEC2"]='[4222:4228,1:2064]'
    header["DATASEC2"]='[2165:4221,1:2064]'
    header["BIASSEC2"]='[2115:2164,1:2064]'
    header["CCDSEC2"]='[2058:4114,1:2064]'
    header["PRESEC3"]='[1:7,2065:4128]'
    header["DATASEC3"]='[8:2064,2065:4128]'
    header["BIASSEC3"]='[2065:2114,2065:4128]'
    header["CCDSEC3"]='[1:2057,2065:4128]'
    header["PRESEC4"]='[4222:4228,2065:4128]'
    header["DATASEC4"]='[2165:4221,2065:4128]'
    header["BIASSEC4"]='[2115:2164,2065:4128]'
    header["CCDSEC4"]='[2058:4114,2065:4128]'

ffile.writeto(args.outfile,clobber=True)
print "wrote",args.outfile
print "done"
