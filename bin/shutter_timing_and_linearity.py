#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from teststand.graph_tools         import plot_graph,parse_fibers
from desispec.log                  import get_logger
import os.path

def mypolfit(x,y,w,deg,force_zero_offset=False) :
    
    n=deg+1
    
    if len(x)<n or np.sum(w>0)<n :
        log.error("not enough data to fit a polynomial of degree %s"%deg)
        log.error("x=%s"%str(x))
        log.error("y=%s"%str(y))
        log.error("w=%s"%str(w))
        sys.exit(12)
        
    
    swxn=np.zeros(2*n)
    for i in range(2*n) :
        swxn[i]=np.sum(w*x**i)
    
    A=np.zeros((n,n))
    for i in range(n) :
        for j in range(n) :
            A[i,j]=swxn[i+j]
    
    B=np.zeros((n))  
    for i in range(n) :
        B[i]=np.sum(w*x**i*y)

    
    if force_zero_offset :
        A[0,0] += 1e8

    try :
        Ai=np.linalg.inv(A)
    except np.linalg.linalg.LinAlgError :
        log.error("cannot fit the polynomial, there must be an error in the data set")
        log.error(sys.exc_info())
        sys.exit(12)
        
    p=Ai.dot(B)
    err=np.sqrt(np.diag(Ai).copy())
    return p[::-1],err[::-1] # to agree with np.polyfit ordering




def readfile(filename) :
    print("reading",filename)
    x=np.loadtxt(filename).T
    file=open(filename)
    keys=file.readlines()[0].strip().replace("# ","").split(" ")
    #print(keys)
    file.close()
    res={}
    for i,k in enumerate(keys) :
        res[k]=x[i].copy()
    return res

def non_linearity_break(args) :
    if args.camera == "b1" and args.amp == "A" :
        return 1.5e4
    elif args.camera == "b1" and args.amp == "B" :
        return 1.e4
    elif args.camera == "r1" and args.amp == "B" :
        return 2.4e4
    elif args.camera == "z1" and args.amp == "B" :
        return 8.e3
    elif args.camera == "z1" and args.amp == "D" :
        return 10.e3
    else :
       print("Non linearity correction for %s %s not implemented !!"%(args.camera,args.amp))
       return 1.e4
   
def non_linearity_correction_b1_A(meas_flux) :
    nl1=+3.7e-06
    nl2=+2.e-06
    thres=1.5e4
    non_linear_correction = (1. + (meas_flux<thres)*nl1*meas_flux + (meas_flux>thres)*(nl2*(meas_flux-thres)+nl1*thres))
    return non_linear_correction

def non_linearity_correction_b1_B(meas_flux) :
    nl1=-3.9e-06
    nl2=+1.8e-06
    thres=1e4
    non_linear_correction = (1. + (meas_flux<thres)*nl1*meas_flux + (meas_flux>thres)*(nl2*(meas_flux-thres)+nl1*thres))
    return non_linear_correction

def non_linearity_correction_r1_B(meas_flux) :
    nl1=-0.e-06
    nl2=+0e-06
    thres=2.4e4
    non_linear_correction = (1. + (meas_flux<thres)*nl1*meas_flux + (meas_flux>thres)*(nl2*(meas_flux-thres)+nl1*thres))
    return non_linear_correction

def non_linearity_correction_z1_B(meas_flux) :
    nl1=-6.85e-06
    nl2=0.
    thres=7.e3
    non_linear_correction = (1. + (meas_flux<thres)*nl1*meas_flux + (meas_flux>thres)*(nl2*(meas_flux-thres)+nl1*thres))
    return non_linear_correction

def non_linearity_correction_z1_D(meas_flux) :
    nl1=-3.5e-06
    nl2=+0e-06
    thres=10.e3
    non_linear_correction = (1. + (meas_flux<thres)*nl1*meas_flux + (meas_flux>thres)*(nl2*(meas_flux-thres)+nl1*thres))
    return non_linear_correction


def non_linearity_correction(meas_flux,args) :
    if args.camera == "b1" and args.amp == "A" :
        return non_linearity_correction_b1_A(meas_flux)
    elif args.camera == "b1" and args.amp == "B" :
        return non_linearity_correction_b1_B(meas_flux)
    elif args.camera == "r1" and args.amp == "B" :
        return non_linearity_correction_r1_B(meas_flux)
    elif args.camera == "z1" and args.amp == "B" :
        return non_linearity_correction_z1_B(meas_flux)
    elif args.camera == "z1" and args.amp == "D" :
        return non_linearity_correction_z1_D(meas_flux)
    print("Non linearity correction for %s %s not implemented !!"%(args.camera,args.amp))
    return np.ones(meas_flux.shape)
##############################################################################################


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input', type = str, default = None, required = True, help = 'ASCII file from meanflux_for_shutter_timing_and_linearity.py')
parser.add_argument('--camera', type = str, default = None, required = True, help = 'camera (b1, r1 or z1)')
parser.add_argument('--amp', type = str, default = None, required = True, help = 'amplifier (A,B,C or D)')
parser.add_argument('--fit-non-lin', action = 'store_true', help = 'try and fit non linear correction')


args = parser.parse_args()

x=readfile(args.input)
x["nd"]=x["nd"].astype(int)
x["fiber"]=x["fiber"].astype(int)

if False :
    plt.figure()
    plt.plot(x["expreq"],x["exptime"]-x["expreq"],"o") # scary correlation
    plt.show()


if False :
    ok=np.where(x["expreq"]<100)[0]
    for k in x.keys() :
        x[k]=x[k][ok]


# hack, use exptime
# x["expreq"]=x["exptime"]

if False : # hack, use only some fibers to avoid potential gain mismatch issue 
    ok=np.where(x["fiber"]==12)[0]
    for k in x.keys() :
        x[k]=x[k][ok]
if False : # hack, 
    ok=np.where((x["nd"]==1))[0]
    for k in x.keys() :
        x[k]=x[k][ok]

xall={}
for k in x.keys() :
    xall[k]=np.zeros(x[k].size)
    for i in range(x[k].size) :
        xall[k][i]=x[k][i]+0



if True :
    print("apply a non linear correction")
    x["flux"] *= non_linearity_correction(x["flux"],args)

scale=1.e-3 # only for plots
threshold=non_linearity_break(args)
print("non linearity break (used for threshold)=",threshold)
if True : # remove brightest data possibly non-linear    
    print("FLUX threshold = %g"%threshold)
    plt.plot(x["expreq"],x["flux"]*scale,"o",c="r")
    plt.plot(x["expreq"],x["flux"]*scale*(x["flux"]<threshold),"o",c="b")
    plt.xlabel(r"EXPREQ (sec)")
    plt.ylabel(r"FLUX $\times 10^{%d}$"%np.log10(scale))

    non_linearity_mask=np.where(x["flux"]<threshold)[0]
    for k in x.keys() :
        x[k]=x[k][non_linearity_mask]



nds=np.unique(x["nd"])
fibers=np.unique(x["fiber"])
expreqs=np.unique(x["expreq"])

print("NDs=",nds)
print("fibers=",fibers)
print("EXPREQs=",expreqs)

force_same_time = True

deltat=0

#ndtrans=None
ndtrans={1:1.,2:10**-0.5,3:0.1,4:0.01}
for nd in nds :
    ok=np.where((x["nd"]==nd))[0]
    if ok.size==0 : continue
    x["flux"][ok] /= ndtrans[nd]


weight = 1./(0.01*x["flux"])**2

fibertrans=None
old_deltat=0.
for loop in range(100) :

    
    if len(fibers)>1 : # do a fiber flat per ND !
        fibertrans_corr={}
        st = 0.
        s1 = 0.
        for fiber in fibers :
            ok=np.where((x["fiber"]==fiber))[0]
            if ok.size==0 : continue
            if loop==0 or not force_same_time :
                coef=np.polyfit(x["expreq"][ok],x["flux"][ok],w=weight[ok],deg=1)
                fibertrans_corr[fiber]=np.poly1d(coef)(np.mean(x["expreq"][ok]))  # the value at the mean exptime
            else :
                # force same deltat
                truetime=x["expreq"][ok]-deltat
                fibertrans_corr[fiber]=np.sum(weight[ok]*x["flux"][ok]*truetime)/np.sum(weight[ok]*truetime**2)
            st += fibertrans_corr[fiber]
            s1 += 1

        # force a mean fiber transmission of one
        mean_fibertrans = st/s1
        for fiber in fibers :    
            fibertrans_corr[fiber] /= mean_fibertrans

        # correct for the fiber transmission
        for fiber in fibers :
            ok=np.where((x["fiber"]==fiber))[0]
            if ok.size==0  : continue
            x["flux"][ok] /= fibertrans_corr[fiber]
            weight[ok] *= fibertrans_corr[fiber]**2
        
        if fibertrans is None :
            fibertrans = fibertrans_corr
        else :
            for fiber in fibers :
                fibertrans[fiber] *= fibertrans_corr[fiber]
    

            


    ############################


    
    if len(nds)>1 : # fit neutral densities
        
        # for each fiber and neutral density, fit polynomial
        ndtrans_corr={}
        for nd in nds :
            ndtrans_corr[nd]={}
            ok=np.where((x["nd"]==nd))[0]
            if ok.size==0 : continue
            if loop==0 or not force_same_time :
                coef=np.polyfit(x["expreq"][ok],x["flux"][ok],w=weight[ok],deg=1)  
                ndtrans_corr[nd]=np.poly1d(coef)(np.mean(x["expreq"][ok]))  # the value at the mean exptime
            else :
                truetime=x["expreq"][ok]-deltat
                ndtrans_corr[nd]=np.sum(weight[ok]*x["flux"][ok]*truetime)/np.sum(weight[ok]*truetime**2)
        if 1 in ndtrans_corr :
            ref=ndtrans_corr[1]
        else :
            ref=ndtrans_corr[2]/0.3
        for nd in nds :
            ndtrans_corr[nd] /= ref
                
        if ndtrans is None :
            ndtrans = ndtrans_corr
        else :
            for nd in nds :
                ndtrans[nd] *= ndtrans_corr[nd]
     
        # correct for the nd transmission
        for nd in nds :
            ok=np.where((x["nd"]==nd))[0]
            if ok.size==0 : continue
            x["flux"][ok] /= ndtrans_corr[nd]
            weight[ok] *= ndtrans_corr[nd]**2

    ######################################
    
    expreq=x["expreq"]
    flux=x["flux"]
    ii=np.argsort(expreq)
    expreq=expreq[ii]
    flux=flux[ii]
    coef=np.polyfit(expreq,flux,deg=1)

    # flux = a*t+b
    # flux=0 for t= -b/a
    deltat=-coef[1]/coef[0]
    #print("iter #%d delta exptime = %f"%(loop,deltat))
    #print("iter #%d ndtrans_corr = %s"%(loop,str(ndtrans_corr)))
    if np.abs(deltat-old_deltat)<0.000001 :
        break
    old_deltat=deltat

print("####### RESULTS #######")
for nd in nds :
    line="ND#%d trans="%nd
    line+=" %5.4f"%ndtrans[nd]
    print(line)
for fiber in fibers :
    line="Fiber #%d trans="%fiber
    line+=" %5.4f"%fibertrans[fiber]
    print(line)
print("delta exptime = %f"%(deltat))
print("nloop = %d"%(loop))
print("#######################")

#ndtrans[3]*=(1-0.04)
#for nd in nds :
#    print("HACKED ND#%d trans= %f"%(nd,ndtrans[nd]))

plt.figure("exposure-time-%s-%s"%(args.camera,args.amp))
plt.subplot(3,1,1)
for nd in nds :
    ok=np.where(x["nd"]==nd)[0]
    plt.plot(x["expreq"][ok],x["flux"][ok]*scale,"o",label="ND #%d"%nd)
#plt.xlabel("EXPREQ")
plt.ylabel(r"FLUX $\times 10^{%d}$"%np.log10(scale))
plt.legend(loc="upper left")


pol=np.poly1d(coef)
tt=np.linspace(deltat,np.max(expreq),3)
plt.plot(tt,pol(tt)*scale,"-",color="k")
xlim=[deltat-0.5,np.max(expreq)*1.1]
plt.xlim(xlim)
plt.grid()

plt.subplot(3,1,2)
#plt.plot(expreq,(flux-pol(expreq))*scale,"o")
plt.plot(tt,0*tt,"-",color="k")
for nd in nds :
    ok=np.where(x["nd"]==nd)[0]
    plt.plot(x["expreq"][ok],(x["flux"][ok]-pol(x["expreq"][ok]))*scale,"o",label="ND #%d"%nd) 
plt.xlim(xlim)
#plt.xlabel("EXPREQ")
plt.ylabel(r"(FLUX-FIT) $\times 10^{%d}$"%np.log10(scale))
plt.grid()


plt.subplot(3,1,3)
#plt.plot(expreq,flux/pol(expreq)-1,"o")
plt.plot(tt,0*tt,"-",color="k")
for nd in nds :
    ok=np.where(x["nd"]==nd)[0]
    offset=0.05*(int(nd)-3.)
    plt.plot(x["expreq"][ok]+offset,x["flux"][ok]/pol(x["expreq"][ok])-1,"o",label="ND #%d"%nd) 

plt.xlim(xlim)
plt.xlabel("EXPREQ (sec)")
plt.ylabel("FLUX/FIT-1")
plt.grid()

####################################################################
# now look at non-linearity
####################################################################

# restart with the whole dataset
flux=xall["flux"]
modelflux=xall["expreq"]-deltat
# apply ND and fiber trans
for fiber in fibers :
    ok=np.where((xall["fiber"]==fiber))[0]
    corr = fibertrans[fiber]
    #print("CORR",fiber,corr)
    modelflux[ok] *= corr
for nd in nds :
    ok=np.where((xall["nd"]==nd))[0]
    corr = ndtrans[nd]
    #print("CORR",nd,corr)
    modelflux[ok] *= corr

# fit the flux of the LED (for this, we apply the non-linearity correction)
w=np.ones(flux.size)
coef,err=mypolfit(modelflux[flux<threshold],flux[flux<threshold]*non_linearity_correction(flux[flux<threshold],args),w=w[flux<threshold],deg=1,force_zero_offset=True)
modelflux *= coef[-2] # apply slope which is defined by unknown illumination

# sort
ii=np.argsort(modelflux)
modelflux=modelflux[ii]
flux=flux[ii]
w=w[ii]

if args.fit_non_lin :
    ok1=np.where(flux<threshold)[0]
    coef1,err1=mypolfit(modelflux[ok1],flux[ok1],w=w[ok1],deg=2,force_zero_offset=True)
    pol1=np.poly1d(coef1)
    non_linear_coef1=coef1[0]
    print("non_linear_coef1=",non_linear_coef1)
    ok2=np.where(flux>=threshold)[0]
    coef2,err2=mypolfit(modelflux[ok2],flux[ok2],w=w[ok2],deg=2,force_zero_offset=False)
    pol2=np.poly1d(coef2)
    coef2[-1] += (pol1(threshold)-pol2(threshold)) # continuity
    pol2=np.poly1d(coef2)
    non_linear_coef2=coef2[0]
    print("non_linear_coef2=",non_linear_coef2)


npy=4
npx=1
plt.figure("linearity-%s-%s"%(args.camera,args.amp))
plt.subplot(npy,npx,1)
#plt.plot(modelflux*scale,flux*scale,"o")
plt.plot(modelflux*scale,modelflux*scale,"-",color="k")
if args.fit_non_lin :
    plt.plot(modelflux[ok1]*scale,pol1(modelflux[ok1])*scale,"--",color="r")
    plt.plot(modelflux[ok2]*scale,pol2(modelflux[ok2])*scale,"--",color="r")
if True :
    for nd in nds :
        ok=np.where(xall["nd"][ii]==nd)[0]
        plt.plot(modelflux[ok]*scale,flux[ok]*scale,"o",label="ND #%d"%nd)
        
plt.ylabel(r"FLUX$\times 10^{%d}$"%np.log10(scale))
plt.xlabel(r"LINEAR FLUX MODEL$\times 10^{%d}$"%np.log10(scale))
plt.legend(loc="upper left")
plt.grid()

plt.subplot(npy,npx,2)
if args.fit_non_lin :
    plt.plot(modelflux[ok1]*scale,(pol1(modelflux[ok1])-modelflux[ok1])*scale,"--",color="r")
    plt.plot(modelflux[ok2]*scale,(pol2(modelflux[ok2])-modelflux[ok2])*scale,"--",color="r")
plt.plot(modelflux*scale,0*modelflux,"-",color="k")
#plt.plot(modelflux*scale,(flux-modelflux)*scale,"o")

if True :
    for nd in nds :
        ok=np.where(xall["nd"][ii]==nd)[0]
        #ok=ok[ii]
        plt.plot(modelflux[ok]*scale,(flux[ok]-modelflux[ok])*scale,"o",label="ND #%d"%nd)
        if False :
            for fiber in fibers :
                ok=np.where((xall["nd"][ii]==nd)&(xall["fiber"][ii]==fiber))[0]
                plt.plot(modelflux[ok]*scale,(flux[ok]-modelflux[ok])*scale,"-",c="gray")

plt.xlabel(r"LINEAR FLUX MODEL$\times 10^{%d}$"%np.log10(scale))   
plt.ylabel(r"(FLUX-MODEL) $\times 10^{%d}$"%np.log10(scale))
plt.grid()

plt.subplot(npy,npx,3)
#plt.plot(modelflux*scale,flux/modelflux-1,"o")
if args.fit_non_lin :
    plt.plot(modelflux[ok1]*scale,(pol1(modelflux[ok1])/modelflux[ok1]-1),"--",color="r")
    plt.plot(modelflux[ok2]*scale,(pol2(modelflux[ok2])/modelflux[ok2]-1),"--",color="r")
plt.plot(modelflux*scale,0*modelflux,"-",color="k")
for nd in nds :
    ok=np.where(xall["nd"][ii]==nd)[0]
    #ok=ok[ii]
    plt.plot(modelflux[ok]*scale,flux[ok]/modelflux[ok]-1,"o",label="ND #%d"%nd)
plt.ylabel("FLUX/MODEL -1")
plt.xlabel(r"LINEAR FLUX MODEL$\times 10^{%d}$"%np.log10(scale))
plt.grid()



plt.subplot(npy,npx,4)

if args.fit_non_lin :
    coef1,err1=mypolfit(flux[ok1],modelflux[ok1],w=w[ok1],deg=2,force_zero_offset=True)
    pol1b=np.poly1d(coef1)
    non_linear_coef1=coef1[0]
    coef2,err2=mypolfit(flux[ok2],modelflux[ok2],w=w[ok2],deg=2,force_zero_offset=False)
    pol2b=np.poly1d(coef2)
    coef2[-1] += (pol1b(threshold)-pol2b(threshold)) # continuity
    pol2b=np.poly1d(coef2)
    non_linear_coef2=coef2[0]
    print("non_linear_coef1 (flux->model)=",non_linear_coef1)
    print("non_linear_coef2 (flux->model)=",non_linear_coef2)
if args.fit_non_lin :
    plt.plot(flux[ok1]*scale,pol1b(flux[ok1])/flux[ok1]-1,"--",color="r")
    plt.plot(flux[ok2]*scale,pol2b(flux[ok2])/flux[ok2]-1,"--",color="r")
plt.plot(flux*scale,non_linearity_correction(flux,args)-1,"-",color="g")
plt.plot(flux*scale,0*modelflux,"-",color="k")
for nd in nds :
    ok=np.where(xall["nd"][ii]==nd)[0]
    #ok=ok[ii]
    plt.plot(modelflux[ok]*scale,modelflux[ok]/flux[ok]-1,"o",label="ND #%d"%nd)

plt.ylabel("MODEL/FLUX -1")
plt.xlabel(r"MEAS. FLUX$\times 10^{%d}$"%np.log10(scale))
plt.grid()

plt.show()
