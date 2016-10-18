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

args = parser.parse_args()

print "will read %s, reformat it and write it to %s"%(args.infile,args.outfile)
print "this is a minimal reformatting of ICS raw data to adapt to DESI soft"

ffile=pyfits.open(args.infile)
ffile.info()

print "loop on image hdus and look for those with EXTNAME starting with 'CCD'" 
for hdu in ffile :
    header=hdu.header
    try :
        extname=header["EXTNAME"]
    except KeyError :
        continue
    if extname.find("CCD")<0 :
        continue
    print "reformat extention",extname
    
    camera="UNKNOWN"
    if extname == "CCDS1B" :
        camera="B1"
    elif extname == "CCDS1R" :
        camera="R1"
    elif extname == "CCDS1Z" :
        camera="Z1"
    
    
    # I need to do this because desi code use first character to find camera arm
    header["CAMERA"]=camera
    
    print "add CAMERA=%s"%header["CAMERA"]
    
    if extname == "CCDS1R" :
        print "changing DATASEC,CCDSEC,BIASEC for",extname
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
    else :
        print "no change of DATASEC,CCDSEC,BIASEC for ",extname
    
    for k in ["GAIN1","GAIN2","GAIN3","GAIN4"] :
        gain=0
        try :
            print "read %s=%s in header"%(k,header[k])
            gain=header[k]
        except KeyError:
            print "WARNING : missing",k,"set value=1"
            header[k]=1.            
        if gain>3. :
            print "WARNING : this gain value seems too large"
            
    if extname == "CCDS1R"  :
        print "WARNING : changing gains"
        reference_gain=1.
        header["GAIN1"]=reference_gain
        header["GAIN2"]=reference_gain
        header["GAIN3"]=reference_gain
        header["GAIN4"]=reference_gain
    
    for k in ["GAIN1","GAIN2","GAIN3","GAIN4"] :
        print "%s=%f"%(k,header[k])
    
    if extname == "CCDS1R"  :
        print "flip %s image !"%extname
        hdu.data = hdu.data[::-1,:]
        print "now modify all SEC info"
        xmin,xmax,ymin,ymax = parse_sec_keyword(header["CCDSEC4"])        
        ny_ccd=ymax
        ny_input=header["NAXIS2"]
        
        header_copy=header.copy()
        #newamps={1:3,2:4,3:1,4:2}
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
                #key="%s%d"%(sec,newamps[amp])
                key="%s%d"%(sec,amp) # DO NOT CHANGE AMP NAMES, CONFUSING
                old=header[key]
                header[key]='[%d:%d,%d:%d]'%(xmin,xmax,flipped_ymin,flipped_ymax)
                print "%s %s -> %s"%(key,old,header[key])
    for a in range(1,5) :
        key="GAIN%d"%a
        if not key in header :
            header[key]=1.
            print "WARNING !! MADE UP",key,header[key]
    
    
      
ffile.writeto(args.outfile,clobber=True)
print "wrote",args.outfile
print "done"
