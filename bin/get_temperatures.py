#!/usr/bin/env python
import astropy.io.fits as pyfits  
import numpy as np 
import sys
import datetime
import dateutil.parser

file=open("temperatures.txt","w")

skeys=["DAY","EXPNUM","TIME","HOUR","EXPREQ","BLUTEMP","REDTEMP","NIRTEMP","PLCTEMP1","PLCTEMP2","EXPSHUT","HARTL","HARTR","BLUHUMID","REDHUMID","NIRHUMID","RCCDTEMP1","RCCDTEMP2","ZCCDTEMP1","ZCCDTEMP2"]

line="#"
for k in skeys :
    line=line+" %s"%k
file.write("%s\n"%line)

dtref=dateutil.parser.parse("2017-01-01T00:00:00.000000-05:00")
for filename in sys.argv[1:] :
    h=pyfits.open(filename)
    speckeys=h["SPECTCON1"].columns.names
    if not "REDTEMP" in speckeys :
        continue
    d={}
    d["EXPNUM"]=h[0].header["EXPNUM"]
    d["EXPREQ"]=h[0].header["EXPREQ"]
    dt = dateutil.parser.parse(h[0].header["DATE-OBS"])
    d["DAY"]=int(dt.strftime('%Y%m%d'))
    d["TIME"]=(dt.timestamp()-dtref.timestamp())
    offset=0 # I think to the time zone is wrong here
    d["HOUR"]=int(dt.strftime('%H'))+float(dt.strftime('%M'))/60.+float(dt.strftime('%S'))/3600.+offset

    if "R1" in h :
        d["RCCDTEMP1"]=h["R1"].header["CCDTEMP1"]
        d["RCCDTEMP2"]=h["R1"].header["CCDTEMP2"]
    else :
        d["RCCDTEMP1"]=0.
        d["RCCDTEMP2"]=0.
    if "Z1" in h :
        d["ZCCDTEMP1"]=h["Z1"].header["CCDTEMP1"]
        d["ZCCDTEMP2"]=h["Z1"].header["CCDTEMP2"]
    else :
        d["ZCCDTEMP1"]=0.
        d["ZCCDTEMP2"]=0.
    if "PLC" in h :
        temps=h["PLC"].header["TEMPS"].split(",")
        d["PLCTEMP1"]=float(temps[0])
        d["PLCTEMP2"]=float(temps[1])
    else :
        d["PLCTEMP1"]=0.
        d["PLCTEMP2"]=0.
    for k in ["BLUTEMP","REDTEMP","NIRTEMP","BLUHUMID","REDHUMID","NIRHUMID"] :
        d[k]=h["SPECTCON1"].data[k][0]
    d["EXPSHUT"]=int(h["SPECTCON1"].data["EXPSHUT"][0]=="CLOSED")
    d["HARTL"]=int(h["SPECTCON1"].data["HARTL"][0]=="OPEN")
    d["HARTR"]=int(h["SPECTCON1"].data["HARTR"][0]=="OPEN")
    
    
    line=""
    for k in skeys :
        line=line+" %s"%str(d[k])
    print(line)
    file.write("%s\n"%line)
file.close()

    
