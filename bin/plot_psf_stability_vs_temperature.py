#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input', type = str, default = None, required = True, help = 'path to ascii file')
parser.add_argument('--title', type = str, default = None, required = False, help = 'title')
parser.add_argument('--temp', type = str, default = None, required = True, help = 'BLUTEMP, REDTEMP, NIRTEMP, PLCTEMP1 or PLCTEMP2')
args        = parser.parse_args()


x=np.loadtxt(args.input).T
i=0
day=x[i] ; i+=1
e1=x[i] ; i+=1
e2=x[i] ; i+=1

temps={}
for k in ["BLUTEMP","REDTEMP","NIRTEMP","PLCTEMP1","PLCTEMP2"] :
    temps[k]=x[i] ; i+=1

# hack
temps["EXPNUM"]=e1

delta_ratio_emission_line=x[i] ; i+=1
delta_ratio_continuum=x[i] ; i+=1
delta_x=x[i] ; i+=1
delta_y=x[i] ; i+=1
sigma_x=x[i] ; i+=1
sigma_y=x[i] ; i+=1

temp=temps[args.temp]
label=args.temp

uday=np.unique(day)
    
plt.figure(args.title,figsize=(14.,9.))
#plt.text(args.title,7.,8.)

ny=3
nx=2
a=1
plt.subplot(ny,nx,a) ; a+=1
for d in uday :
    plt.plot(temp[day==d],delta_x[day==d],"o",label="%d"%d)
plt.grid()
plt.ylabel("delta x (pixels)")
plt.legend(loc="upper right",fontsize="small")
plt.subplot(ny,nx,a) ; a+=1
for d in uday :
    plt.plot(temp[day==d],delta_y[day==d],"o")
plt.grid()
plt.ylabel("delta y (pixels)")
plt.subplot(ny,nx,a) ; a+=1
for d in uday :
    plt.plot(temp[day==d],sigma_x[day==d],"o")
plt.grid()
plt.ylabel("sigma x (pixels)")
plt.subplot(ny,nx,a) ; a+=1
for d in uday :
    plt.plot(temp[day==d],sigma_y[day==d],"o")
plt.grid()
plt.ylabel("sigma y (pixels)")
plt.subplot(ny,nx,a) ; a+=1
for d in uday :
    plt.plot(temp[day==d],delta_ratio_emission_line[day==d],"o")
plt.grid()
plt.xlabel(label)
plt.ylabel("emission line flux ratio")
plt.subplot(ny,nx,a) ; a+=1
for d in uday :
    plt.plot(temp[day==d],delta_ratio_continuum[day==d],"o")
plt.grid()
plt.xlabel(label)
plt.ylabel("continuum flux ratio")
plt.show()
