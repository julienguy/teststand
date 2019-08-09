#!/usr/bin/env python


import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
from desispec.io import read_spectra
               
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*", 
                    help = 'path to spectra file(s)')
parser.add_argument('-t','--targetid', type = int, default = None, required = False, 
                    help = 'plot specific targetid')

#log         = get_logger()
args        = parser.parse_args()

spectra = []
for filename in args.infile :
    spectra.append(read_spectra(filename))

if args.targetid is not None :
    targetids=[args.targetid,]
else :
    targetids=np.unique(spectra[0].fibermap["TARGETID"])
    
for tid in targetids :
    print("TARGETID {}",tid)
    for spec in spectra :
        jj=np.where(spec.fibermap["TARGETID"]==tid)[0]
        for k in ['FIRST_EXPID','LAST_EXPID','NUM_EXPID','EXPID'] :
            if k in spec.fibermap.colnames :
                for j in jj :
                    print("{} {}".format(k,spec.fibermap[k][j]))

        for j in jj :
            for b in spec._bands :
                
                i=np.where(spec.ivar[b][j]*(spec.mask[b][j]==0)>0)[0]
                plt.plot(spec.wave[b][i],spec.flux[b][j,i])
                c=np.polyfit(spec.wave[b][i],spec.flux[b][j,i],3)
                pol=np.poly1d(c)(spec.wave[b][i])
                print("mean chi2=",np.mean(spec.ivar[b][j,i]*(spec.flux[b][j,i]-pol)**2))
    plt.show()


plt.show()


