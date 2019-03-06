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
parser.add_argument('-t','--title', type = str, default = None, required = False,
                    help = 'Figure title')
parser.add_argument('--x-axis',type=str, default='EXPID', required=False, 
    help = 'x-axis, can be EXPID, MJD-OBS, DATE-OBS')

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
        if not args.title is None:
            plt.suptitle(r'$\mathrm{'+args.title+'}$',fontsize=20)
        plt.subplots_adjust(top=0.95)
        
        i=1
        for a,amp in enumerate(["A","B","C","D"]) :
            color=colors[a]
            x = t[args.x_axis]
            print(x)
            ax=plt.subplot(4,2,i); i+=1
            ax.plot(x,t["RMS_CCD_"+amp],"-",c=color,label="CCD")
            ax.plot(x,t["RMS_COL_OVERSCAN_"+amp],"--",c=color,label="Overscan cols")
            ax.plot(x,t["RMS_ROW_OVERSCAN_"+amp],":",c=color,label="Overscan rows")
            ax.grid()
            ax.legend(title="Amplifier {} rms".format(amp))
            if amp=="D" : ax.set_xlabel(args.x_axis)
            ax.set_ylabel("electrons/pixel")
            
            ax=plt.subplot(4,2,i); i+=1
            #offset=np.min(t["MEAN_COL_OVERSCAN_"+amp])
            offset=t["MEAN_COL_OVERSCAN_"+amp][-1]
            ax.plot(x,t["MEAN_CCD_"+amp]-offset,"-",c=color,label="CCD")
            ax.plot(x,t["MEAN_COL_OVERSCAN_"+amp]-offset,"--",c=color,label="Overscan cols")
            ax.plot(x,t["MEAN_ROW_OVERSCAN_"+amp]-offset,":",c=color,label="Overscan rows")
            ax.legend(title="Amplifier {} mean".format(amp))
            ax.grid()
            if amp=="D" : ax.set_xlabel(args.x_axis)

            print("|| Amplifier {} rdnoise || {:3.2f} electrons/pixel||".format(amp,np.mean(t["RMS_CCD_"+amp][-3:])))
            

        if args.outfile is None :
            args.outfile = name+".png"
        print("saving image as",args.outfile)
        fig.savefig(args.outfile)
   
plt.show()

