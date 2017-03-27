#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt

def read(filename) :
    x=np.loadtxt(filename).T
    file=open(filename)
    for line in file.readlines() :
        if line[0]=="#" and line[1]!="#" :
            keys=line.strip().replace("# ","").split(" ")
            vals={}
            for i,k in enumerate(keys) :
                vals[k]=x[i]
            break
    return vals

d=read("psf_stability_vs_temperature_b1.txt")

#for cam,tempkey in zip(["b1","r1","z1"],["BLUTEMP","REDTEMP","NIRTEMP"]) :
#for cam,tempkey in zip(["b1","r1"],["BLUTEMP","REDTEMP"]) :
for cam,tempkey in zip(["b1",],["BLUTEMP",]) :

    if cam == "b1" :
        camid=0
    elif cam == "r1" :
        camid=1
    elif cam == "z1" :
        camid=2
    ok=(d["CAMID"]==camid)
    print("n selected=",np.sum(ok))
    istime=False
    isexpnum=False
    x=d["EXPNUM"][ok] ; xlabel="EXPOSURE" ; isexpnum=True
    #x=d["TIME"][ok]/3600. ; x -= np.min(x) ;  xlabel="Time (hours, integrated, arbitrary offset)" ; istime=True
    #x=d["HOUR"][ok]; xlabel="Hour (in the day)"
    #x=d["DAY"]-20170100; xlabel="Day"

    print("days = ",np.unique(d["DAY"]).astype(int))
    
    breaks = []
    if istime :
        moffset=20
        print ("remove large offsets")
        for i in range(1,x.size) :
            if x[i]-x[i-1]>moffset :
                offset=moffset-(x[i]-x[i-1])
                x[i:] += offset
                breaks.append(x[i]-moffset/2.)
    
    ii=np.argsort(x)
    x=x[ii]
    
    ebias=d["EBIAS"][ok][ii]
    ebias-=np.mean(ebias)
    sx=d["SX"][ok][ii]
    sy=d["SY"][ok][ii]
    dx=d["DX"][ok][ii]
    dy=d["DY"][ok][ii]
    cx=d["CX"][ok][ii]
    cy=d["CY"][ok][ii]
    temp=d[tempkey][ok][ii]
    blutemp=d["BLUTEMP"][ok][ii]
    redtemp=d["REDTEMP"][ok][ii]
    nirtemp=d["NIRTEMP"][ok][ii]
    temp2key="PLCTEMP1"
    temp3key="ZCCDTEMP2"
    temp2=d[temp2key][ok][ii]
    temp3=d[temp3key][ok][ii]
    expnum=d["EXPNUM"][ok][ii]
    humid=(d["BLUHUMID"]+d["REDHUMID"]+d["NIRHUMID"])[ok][ii]/3.
    msx=np.mean(sx)
    msy=np.mean(sy)
    ebias2=1/( 2./np.sqrt((1+(sx/msx)**2)*(1+(sy/msy)**2)) )-1.
    
    ######################################################
    plt.figure("psf-shape-%s"%cam,figsize=(19,12))
    ######################################################
    a=plt.subplot(3,1,1)
    plt.plot(x,sx,"o",label=r"$\sigma_X$")
    plt.plot(x,sy,"o",label=r"$\sigma_Y$")
    ymin=min(np.min(sx),np.min(sy))
    ymax=max(np.max(sx),np.max(sy))
    deltay=(ymax-ymin)
    plt.ylim([ymin-0.1*deltay,ymax+0.3*deltay])
    
    plt.xlabel(xlabel)
    plt.ylabel("PSF shape")
    plt.grid()    
    for b in breaks : plt.axvline(x=b, linewidth=20, color = 'gray', alpha=0.5)
    plt.legend(loc="upper left")
    ######################################################
    a=plt.subplot(3,1,2)
    #plt.plot(x,ebias2,"o",c="g",alpha=0.5,label=r"$(\sigma_X \sigma_Y)/< \sigma_X \sigma_Y > -1$")
    plt.plot(x,ebias,"o",c="b",label="bias on emission line flux")
    ymin=min(np.min(ebias),np.min(ebias2))
    ymax=max(np.max(ebias),np.max(ebias2))
    deltay=(ymax-ymin)
    plt.ylim([ymin-0.1*deltay,ymax+0.1*deltay])
    plt.ylim([-0.03,0.03])
    
    plt.xlabel(xlabel)
    plt.ylabel("PSF shape")
    plt.grid()
    for b in breaks : plt.axvline(x=b, linewidth=20, color = 'gray', alpha=0.5)
    plt.legend(loc="upper left")
    ######################################################
    a=plt.subplot(3,1,3)
    plt.plot(x,temp,"o",color="red",label=tempkey)
    plt.plot(x,temp2,"o",color="orange",label=temp2key)
    #plt.plot(x,humid,"o",color="purple",label="HUMIDITY")
    plt.plot(x,temp3,"o",color="purple",label=temp3key)
    plt.xlabel(xlabel)
    plt.ylabel("temperature (deg. Celsius)")
    plt.legend(loc="upper left")
    #plt.ylim([17.,20.5])
    plt.grid()
    for b in breaks : plt.axvline(x=b, linewidth=20, color = 'gray', alpha=0.5)
    ######################################################
    if False :
        plt.figure("trace-shift-%s"%cam,figsize=(19,12))
        a=plt.subplot(2,1,1)
        plt.plot(x,dx,"o",label=r"$\Delta x$")
        plt.plot(x,dy,"o",label=r"$\Delta y$")
        plt.xlabel(xlabel)
        plt.ylabel("trace offset (pix.)")
        plt.legend(loc="lower left")
        for b in breaks : plt.axvline(x=b, linewidth=20, color = 'gray', alpha=0.5)
        plt.grid()
        a=plt.subplot(2,1,2)
        plt.plot(x,temp,"o",color="red",label=tempkey)
        plt.plot(x,temp2,"o",color="orange",label=temp2key)
        
        #plt.plot(x,blutemp,"o",color="b",label="BLUTEMP")
        #plt.plot(x,redtemp,"o",color="g",label="REDTEMP")
        #plt.plot(x,nirtemp,"o",color="r",label="NIRTEMP")
        plt.xlabel(xlabel)
        plt.ylabel("temperature (deg. Celsius)")
        plt.legend(loc="upper left")
        #plt.ylim([17.,20.5])
        for b in breaks : plt.axvline(x=b, linewidth=20, color = 'gray', alpha=0.5)
        plt.grid()
    ######################################################
    
    ######################################################
    if False :
        plt.figure("relation-with-temperature-%s"%cam)
        ok=(expnum<3500)|(expnum>3530)
        
        plt.plot(temp2[ok],ebias[ok],"o")
        plt.xlabel("temperature (deg. Celsius)")
        plt.ylabel("bias on emission line flux")
        coef=np.polyfit(temp2[ok],ebias[ok],deg=1)
        pol=np.poly1d(coef)
        plt.plot(temp2,pol(temp2),c="r")
        print("slope = %f per degree Celcius"%coef[0])
        plt.grid()
    ######################################################
    if False :
        plt.figure("relation-with-centers-%s"%cam)
        #ok=(expnum<3500)|(expnum>3530)
        ok=(expnum>2000)&(expnum<2700)        
        plt.plot(cx[ok],ebias[ok],"o")
        plt.xlabel("barycenter of PSF in stamp")
        plt.ylabel("bias on emission line flux")
        plt.grid()
    ######################################################
    
plt.show()

