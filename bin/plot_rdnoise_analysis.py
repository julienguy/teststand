#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, nargs="*", default = None, required = True,
                    help = 'path to ASCII files produced with the script rdnoise_analysis.py')
parser.add_argument('-o','--outfile', type = str, default = None, required = False,
                    help = 'save figure in this file')

args = parser.parse_args()

for filename in args.infile :
    
    file=open(filename)
    lines=file.readlines()
    file.close()
    keys=lines[0].replace("#","").split()
    print(keys)
    x=np.loadtxt(filename).T
    print(x.shape)
    t=dict()
    for i in range(x.shape[0]) :
        t[keys[i]]=x[i]
    
    nfiles=x.shape[1]
    print("nfiles=",nfiles)
    
    colors="bgrk"
    
    #print("overscan_rdnoise.shape=",overscan_rdnoise.shape)
    
    if 1 :
        name=os.path.basename(filename).split(".")[0]
        fig=plt.figure(name,figsize=(8,12))
        ax0=plt.subplot(4,2,1)
        ax1=plt.subplot(2,1,2)
        
        i=1
        for a,amp in enumerate(["A","B","C","D"]) :
            color=colors[a]
            x=np.arange(nfiles)
            ax=plt.subplot(4,2,i); i+=1
            ax.plot(x,t["RMS_CCD_"+amp],"-",c=color,label="CCD")
            ax.plot(x,t["RMS_COL_OVERSCAN_"+amp],"--",c=color,label="Overscan cols")
            ax.plot(x,t["RMS_ROW_OVERSCAN_"+amp],":",c=color,label="Overscan rows")
            ax.grid()
            ax.legend(title="Amplifier {} rms".format(amp))
            if amp=="D" : ax.set_xlabel("expnum - offset")
            ax.set_ylabel("electrons/pixel")
            
            ax=plt.subplot(4,2,i); i+=1
            #offset=np.min(t["MEAN_COL_OVERSCAN_"+amp])
            offset=t["MEAN_COL_OVERSCAN_"+amp][-1]
            ax.plot(x,t["MEAN_CCD_"+amp]-offset,"-",c=color,label="CCD")
            ax.plot(x,t["MEAN_COL_OVERSCAN_"+amp]-offset,"--",c=color,label="Overscan cols")
            ax.plot(x,t["MEAN_ROW_OVERSCAN_"+amp]-offset,":",c=color,label="Overscan rows")
            ax.legend(title="Amplifier {} mean".format(amp))
            ax.grid()
            if amp=="D" : ax.set_xlabel("expnum - offset")

            print("Amplifier {} rdnoise = {:3.2f} electrons/pixel".format(amp,np.mean(t["RMS_CCD_"+amp][-3:])))
            
        
        fig.savefig(name+".png")
   
plt.show()
