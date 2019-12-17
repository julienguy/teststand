#!/usr/bin/env python


import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table



timeindex=dict()
xpix=dict()
ypix=dict()

for number,filename in enumerate(sys.argv[1:]) :
    
    table=Table.read(filename,format="csv")
    ii = np.where(table["FID_ID"]>0)[0]
    plt.plot(table["XPIX"][ii],table["YPIX"][ii],".")

   
    for i in ii :
        fid=table["FID_ID"][i]
        if int(fid) not in timeindex.keys() :
            timeindex[fid] = [ number ,]
            xpix[fid] = [ table["XPIX"][i], ]
            ypix[fid] = [ table["YPIX"][i], ]
        else :
            timeindex[fid].append(number)
            xpix[fid].append(table["XPIX"][i])
            ypix[fid].append(table["YPIX"][i])
    
    
if "PIN_ID" in table.dtype.names and number == 0 :
    ii=table["PIN_ID"]==1027
    if np.sum(ii)>0 :
        plt.plot(table["XPIX"][ii],table["YPIX"][ii],"X")

if len(sys.argv[1:])>0 :
    plt.figure("XPIX")
    for fid in xpix.keys() :
        plt.plot(timeindex[fid],xpix[fid]-np.mean(xpix[fid]),"o-")
    plt.figure("YPIX")
    for fid in ypix.keys() :
        plt.plot(timeindex[fid],ypix[fid]-np.mean(ypix[fid]),"o-")
    

plt.show()

