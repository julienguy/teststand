#!/usr/bin/env python


import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table



tmp=np.loadtxt("/global/homes/j/jguy/project/focalplane/fvc/new/expid-adc.txt").T

adc1=dict()
adc2=dict()

for i in range(tmp.shape[1]) :
    adc1[tmp[0,i]] = tmp[1,i]
    adc2[tmp[0,i]] = tmp[2,i]
    



timeindex=dict()
xpix=dict()
ypix=dict()
xfp=dict()
yfp=dict()
adcangle1=dict()
adcangle2=dict()




for number,filename in enumerate(sys.argv[1:]) :

    expid = int(filename.split("-")[1].split(".")[0])
    a1 = adc1[expid]
    a2 = adc2[expid]
    print(filename,a1,a2)
    
    
    table=Table.read(filename,format="csv")
    ii = np.where((table["FID_ID"]>0)&(table["PIN_ID"]==1027))[0]
    plt.figure("fp")
    plt.plot(table["XPIX"][ii],table["YPIX"][ii],".")

   
    for i in ii :
        fid=table["FID_ID"][i]
        if int(fid) not in timeindex.keys() :
            timeindex[fid] = [ number ,]
            xpix[fid] = [ table["XPIX"][i], ]
            ypix[fid] = [ table["YPIX"][i], ]
            if "XFP" in table.dtype.names :
                xfp[fid] = [ table["XFP"][i], ]
                yfp[fid] = [ table["YFP"][i], ]
            adcangle1[fid] = [ a1 ]
            adcangle2[fid] = [ a2 ]
            
        else :
            timeindex[fid].append(number)
            xpix[fid].append(table["XPIX"][i])
            ypix[fid].append(table["YPIX"][i])
            if "XFP" in table.dtype.names :
                xfp[fid].append(table["XFP"][i])
                yfp[fid].append(table["YFP"][i])
            adcangle1[fid].append(a1)
            adcangle2[fid].append(a2)


    if number == 0 :
        if "XFP" in table.dtype.names :
            x=table["XMETRO"][ii]
            y=table["YMETRO"][ii]
            dx=(table["XFP"][ii]-table["XMETRO"][ii])
            dy=(table["YFP"][ii]-table["YMETRO"][ii])
            plt.figure("quiver")
            plt.plot(x,y,"o")
            plt.plot(table["XFP"][ii],table["YFP"][ii],"o")

            jj = np.where(table["FID_ID"]>0)[0]
            plt.plot(table["XFP"][jj],table["YFP"][jj],".")
            plt.quiver(x,y,dx,dy)
            plt.show()
            sys.exit(12)
            
        
            
#adcangle1=np.array(adcangle1)
#adcangle2=np.array(adcangle2)


if "PIN_ID" in table.dtype.names and number == 0 :
    ii=table["PIN_ID"]==1027
    if np.sum(ii)>0 :
        plt.figure("fp")
        plt.plot(table["XPIX"][ii],table["YPIX"][ii],"X")

if len(sys.argv[1:])>0 :
    plt.figure("XPIX")
    for fid in xpix.keys() :
        #plt.plot(timeindex[fid],xpix[fid]-np.mean(xpix[fid]),"o-")
        plt.plot(np.array(adcangle1[fid])-np.array(adcangle2[fid]),xpix[fid]-np.mean(xpix[fid]),"o")
    plt.xlabel("adc1-adc2")
    plt.ylabel("pixels")
    
    plt.figure("YPIX")
    for fid in ypix.keys() :
        #plt.plot(timeindex[fid],ypix[fid]-np.mean(ypix[fid]),"o-")
        plt.plot(np.array(adcangle1[fid])-np.array(adcangle2[fid]),ypix[fid]-np.mean(ypix[fid]),"o")
    plt.xlabel("adc1-adc2")
    plt.ylabel("pixels")

    plt.figure("YPIX2")
    for fid in ypix.keys() :
        #plt.plot(timeindex[fid],ypix[fid]-np.mean(ypix[fid]),"o-")
        plt.plot((np.array(adcangle1[fid])+np.array(adcangle2[fid]))%360,ypix[fid]-np.mean(ypix[fid]),"o")
    plt.xlabel("adc1+adc2")
    plt.ylabel("pixels")
     
    if len(xfp.keys())>0 :
        plt.figure("XFP")
        for fid in xfp.keys() :
            plt.plot(timeindex[fid],xfp[fid]-np.mean(xfp[fid]),"o-")
            break
        plt.figure("YFP")
        for fid in yfp.keys() :
            plt.plot(timeindex[fid],yfp[fid]-np.mean(yfp[fid]),"o-") 
            break

plt.show()

