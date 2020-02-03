#!/usr/bin/env python

import os
import sys
import numpy as np

import fitsio
import matplotlib.pyplot as plt

snr=[]
chi2=[]
x=[]
y=[]
xmin=12.
xmax=-12.
for filename in sys.argv[1:] :

    t=fitsio.read(filename,"METADATA")
    print(t.dtype.names)
    plt.plot(t["DATA_G-R"],t["MODEL_G-R"],"o",label='{}'.format(os.path.basename(filename)))
    x.append(t["DATA_G-R"])
    y.append(t["MODEL_G-R"])
    if "BLUE_SNR" in t.dtype.names :
        snr.append(t["BLUE_SNR"])
    chi2.append(t["CHI2DOF"])

x=np.hstack(x)
y=np.hstack(y)
xmin=np.min(x)
xmax=np.max(x)
snr=np.hstack(snr)

print("Std of color difference                   = {:4.3f}".format(np.sqrt(np.mean((x[y!=0]-y[y!=0])**2))))

selection=(x>0.25)&(y!=0)&(snr>0)&(np.abs(x-y)<0.3)
print("Std of color difference for G-R(data)>0.3 = {:4.3f}".format(np.sqrt(np.mean((x[selection]-y[selection])**2))))
#plt.plot(x[selection],y[selection],".",c="k")
c=np.polyfit(x[selection],y[selection],1)
pol=np.poly1d(c)




plt.plot([xmin,xmax],[xmin,xmax],c="gray")
ii=np.argsort(x[selection])
plt.plot(x[selection][ii],pol(x[selection])[ii],"--",c="gray")

rms=np.sqrt(np.mean((pol(x[selection])-y[selection])**2))
print("RMS of linear fit = {}".format(rms))

plt.xlabel("DATA G-R")
plt.ylabel("MODEL G-R")
plt.legend()


if len(snr)>0 :
    plt.figure()
    
    plt.plot(snr[x>0.25],y[x>0.25]-x[x>0.25],"o")
    plt.xlabel("BLUE SNR")
    plt.ylabel("MODEL-DATA G-R")

if True :
    plt.figure()
    chi2=np.hstack(chi2)
    plt.plot(chi2[x>0.25],y[x>0.25]-x[x>0.25],"o")
    plt.xlabel("CHI2/NDF")
    plt.ylabel("MODEL-DATA G-R")
plt.show()


