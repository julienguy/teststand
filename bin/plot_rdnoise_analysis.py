#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os

# BLUE CAMERA (before ICS transfo) :
# (complex tranfo, see etc/ics_transfo.py)
# A C
# B D
# 
# RED CAMERA (before ICS transfo data[::-1,::-1])
# B A
# D C
# 
# NIR CAMERA (before ICS transfo data[,::-1]) :
# D C
# B A

gains={}
detectors={}
saclay_rdnoise={}

# SM2
spectro="SM02"
detectors["b1"]="sn22802"
detectors["r1"]="M1-46"
detectors["z1"]="M1-42"
gains["b1"]=[1.28,1.29,1.32,1.31]
gains["r1"]=[1.66,1.50,1.60,1.60]
gains["z1"]=[1.48,1.49,1.67,1.67]
saclay_rdnoise["b1"]=[4.55,3.74,3.69,3.40]
saclay_rdnoise["r1"]=[2.65,2.33,2.56,2.62]
saclay_rdnoise["z1"]=[3.64,2.70,3.72,2.76]
# SM3
spectro="SM03"
detectors["b2"]="sn22794"
detectors["r2"]="M1-12"
detectors["z2"]="M1-22"
gains["b2"]=[1.28,1.27,1.27,1.29]
gains["r2"]=[1.77,1.67,1.50,1.53]
gains["z2"]=[1.45,1.50,1.63,1.52]
saclay_rdnoise["b2"]=[4.90,3.86,3.90,3.77]
saclay_rdnoise["r2"]=[3.24,2.88,3.37,2.68]
saclay_rdnoise["z2"]=[3.72,4.52,3.91,3.23]
# SM4
title="SM04 (20190130)"
spectro="sm4"
detectors["b3"]="sn22797"
detectors["r3"]="M1-49"
detectors["z3"]="M1-53" # WILL CHANGE!
gains["b3"]= [ 1.32, 1.29 ,1.30 , 1.30] # [1.32,1.30,1.29,1.30]
gains["r3"]= [ 1.68, 1.51, 1.52, 1.60 ] # 1.51 1.68 1.52 1.60
gains["z3"]= [ 1.44, 1.47, 1.61, 1.53] # SN5 NIR CRYO HERE !!!!!!!!!!
saclay_rdnoise["b3"]=[ 5.01,4.14 ,4.39 , 3.67] # [5.01,4.39,4.14,3.67]
saclay_rdnoise["r3"]=[ 3.45, 2.87, 4.22, 3.23] # 2.87 3.45 3.23 4.22
saclay_rdnoise["z3"]=[ 3.06, 2.91, 3.24, 2.74] # SN5 NIR CRYO HERE !!!!!!!!!!


labels="ABCD"


#for cam in ['b1','r1','z1'] :
for cam in ['r3','z3'] :
 
    print(cam)
    #filename="rdnoise-pairs-{}.txt".format(cam)
    filename="levels-rms-{}-{}.txt".format(spectro,cam)
    t=np.loadtxt(filename).T
    
    expnum=t[0]
    mjd=t[1]
    time=(mjd-mjd[0])*24.*60. # min
    offset=t[2:6]
    rdnoise=t[6:10]
    overscan_rdnoise=t[10:14]
    colors="bgrk"
    #print(rdnoise)
    #print(overscan_rdnoise)
    
    if 1 :
        name="rdnoise-{}-{}".format(cam,spectro)
        fig=plt.figure(name)
        ax0=plt.subplot(2,1,1)
        ax1=plt.subplot(2,1,2)
        
        for a in range(4) :
            color=colors[a]
            gain=gains[cam][a]
            noise=rdnoise[a]
            osnoise=overscan_rdnoise[a]
            off=offset[a]
            sn=saclay_rdnoise[cam][a]
            x=time
            ax0.plot(x,noise,"-",c=color,label=labels[a])
            ax1.plot(x,off,"-",c=color)
            label=None
            if a==0 :
                label=labels[a]+" overscan"
            ax0.plot(x,osnoise,"--",c=color,label=label)
            #ax0.plot(x,np.sqrt(osnoise**2+off),"--",c=color)
            
#,label="{} noise={:3.2f} e (with G={:3.2f} e/ADU)".format(labels[a-2],noise,gain))
            print("{} {} Read noise={:3.2f} elec (Saclay: {:3.2f}) Overscan noise={:3.2f}; assuming G={:3.2f} elec/ADU".format(detectors[cam],labels[a],noise[-1],sn,osnoise[-1],gain))
        #ax1.set_xlabel("exp. num")
        ax1.set_xlabel("minutes")
        ax0.set_ylim([2.,6.])
        
        ax0.set_ylabel("read noise (electrons)")
        ax0.grid()
        ax0.legend(title="{} {} {}".format(title,cam,detectors[cam]),loc="upper right")
        ax1.set_ylabel("mean level (electrons)")
        ax1.set_ylim([0.,10.])
        ax1.grid()
        fig.savefig(name+".png")
   
plt.show()

