#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import specter.psf
import sys
import argparse
import string
import os.path
from teststand.graph_tools import parse_fibers
from desispec.log                  import get_logger

def readpsf(filename) :
    try :
        psftype=pyfits.open(filename)[0].header["PSFTYPE"]
    except KeyError :
        psftype=""
    print("PSF Type=",psftype)
    if psftype=="GAUSS-HERMITE" :
        return specter.psf.GaussHermitePSF(filename)
    elif psftype=="SPOTGRID" :
        return specter.psf.SpotGridPSF(filename)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf', type = str, nargs = "*", default = None, required = True,
                    help = 'path of psf files')
parser.add_argument('--temp', type = str, default = None, required = True,
                    help = 'path to temperature file')

parser.add_argument('-o','--output', type = str, default = None, required = False, help = 'path to output ascii file')
parser.add_argument('--plot', action='store_true',help="plot result")
parser.add_argument('--wave', type = float, required=True, help= "wavelength")
parser.add_argument('--fiber', type = int, required=True, help= "fiber")


args        = parser.parse_args()
log = get_logger()

x=np.loadtxt(args.temp).T
file=open(args.temp)
line=file.readlines()[0]
file.close()
keys=line.strip().replace("# ","").split(" ")
vals={}
for i,k in enumerate(keys) :
    vals[k]=x[i]
vals["DAY"]=vals["DAY"].astype(int)
vals["EXPNUM"]=vals["EXPNUM"].astype(int)

tkeys=["BLUTEMP","REDTEMP","NIRTEMP","PLCTEMP1","PLCTEMP2"]
psfs=[]
exp1=[]
exp2=[]
temps={}
for k in tkeys :
  temps[k]=[]
day = []
for filename in args.psf :
    log.info("reading %s"%filename)
    psfs.append(readpsf(filename))
    x=os.path.basename(filename).replace(".fits","").split("-")
    e1=int(x[3])
    e2=int(x[4])
    print(e1,e2)
    ok=np.where((vals["EXPNUM"]>=e1)&(vals["EXPNUM"]<=e2))[0]
    for k in tkeys :
        temps[k].append(np.mean(vals[k][ok]))
    exp1.append(e1)
    exp2.append(e2)
    day.append(np.mean(vals["DAY"][ok]))

day=np.array(day)
for k in temps.keys() :
    temps[k]=np.array(temps[k])

print(np.unique(day))

images = []
i0 = []
i1 = []


for psf in psfs :
    xx, yy, ccdpix = psf.xypix(args.fiber,args.wave)
    images.append(ccdpix)
    i1.append(xx.start)
    i0.append(yy.start)

mi0 = int(np.min(i0))
mi1 = int(np.min(i1))
n=len(images)
# add a margin and apply offset if necessary
nimages=np.zeros((n,images[0].shape[0]+3,images[0].shape[1]+3))
for j in range(n) :
    nimages[j,i0[j]-mi0:i0[j]-mi0+images[j].shape[0],i1[j]-mi1:i1[j]-mi1+images[j].shape[1]] = images[j]        
images=nimages
n=images.shape[0]
mimage=np.mean(images,axis=0)
delta_ratio_emission_line = np.zeros(n)
delta_ratio_continuum = np.zeros(n)
delta_x = np.zeros(n)
delta_y = np.zeros(n)
sigma_x = np.zeros(n)
sigma_y = np.zeros(n)

y=np.tile(np.arange(mimage.shape[0]),(mimage.shape[1],1))    # which is which ???
x=np.tile(np.arange(mimage.shape[1]),(mimage.shape[0],1)).T
if 0 :
    plt.figure()
    plt.subplot(2,1,1)
    plt.imshow(x,origin=0)
    plt.subplot(2,1,2)
    plt.imshow(y)
    plt.show()



mx=np.sum(x*mimage)/np.sum(mimage)
my=np.sum(y*mimage)/np.sum(mimage)
#msx=np.sqrt(np.sum(x**2*mimage)/np.sum(mimage)-mx**2)
#msy=np.sqrt(np.sum(y**2*mimage)/np.sum(mimage)-my**2)


for j in range(n) :
    delta_ratio_emission_line[j] = np.sum(images[j]*mimage)/np.sum(images[j]**2)-1
    pmimage=np.sum(mimage,axis=0) # projection to get 1D PSF along cross-dispersion for continuum fit normalization
    pimage=np.sum(images[j],axis=0) 
    delta_ratio_continuum[j] = np.sum(pimage*pmimage)/np.sum(pimage**2)-1
    xj = np.sum(x*images[j])/np.sum(images[j])
    yj = np.sum(y*images[j])/np.sum(images[j])
    delta_x[j] = xj - mx
    delta_y[j] = yj - my
    sigma_x[j] = np.sqrt(np.sum(x**2*images[j])/np.sum(images[j]) - xj**2)
    sigma_y[j] = np.sqrt(np.sum(y**2*images[j])/np.sum(images[j]) - yj**2)
    

if args.output is not None :
    file = open(args.output,"w")
    file.write("# wave=%d fiber=%d\n"%(args.wave,args.fiber))
    line="# day first_expnum last_expnum"
    for k in tkeys :
        line += " %s"%k
    line+=" delta_ratio_emission_line delta_ratio_continuum delta_x delta_y sigma_x sigma_y"
    file.write("%s\n"%line)
    for j in range(n) :
        line="%d %d %d"%(day[j],exp1[j],exp2[j])
        for k in tkeys :
            line=line+" %s"%temps[k][j]
        line=line+" %f %f %f %f %f %f"%(delta_ratio_emission_line[j],delta_ratio_continuum[j],delta_x[j],delta_y[j],sigma_x[j],sigma_y[j])
        file.write("%s\n"%line)
    file.close()
    print("wrote",args.output)

if args.plot :    
    temp=temps["PLCTEMP1"]
    plt.figure()
    ny=3
    nx=2
    a=1
    plt.subplot(ny,nx,a) ; a+=1
    plt.plot(temp,delta_x,"o")
    plt.grid()
    plt.ylabel("delta x (pixels)")
    plt.subplot(ny,nx,a) ; a+=1
    plt.plot(temp,delta_y,"o")
    plt.grid()
    plt.ylabel("delta y (pixels)")
    plt.subplot(ny,nx,a) ; a+=1
    plt.plot(temp,sigma_x,"o")
    plt.grid()
    plt.ylabel("sigma x (pixels)")
    plt.subplot(ny,nx,a) ; a+=1
    plt.plot(temp,sigma_y,"o")
    plt.grid()
    plt.ylabel("sigma y (pixels)")
    plt.subplot(ny,nx,a) ; a+=1
    plt.plot(temp,delta_ratio_emission_line,"o")
    plt.grid()
    plt.ylabel("emission line flux ratio")
    plt.subplot(ny,nx,a) ; a+=1
    plt.plot(temp,delta_ratio_continuum,"o")
    plt.grid()
    plt.ylabel("continuum flux ratio")
    plt.show()
