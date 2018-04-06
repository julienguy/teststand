#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from teststand.graph_tools         import plot_graph,parse_fibers
from desiutil.log                  import get_logger
import os.path


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

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input', type = str, default = None, required = True, help = 'ASCII file from meanflux_for_shutter_timing_and_linearity.py')
parser.add_argument('--title', type = str, default = None, required = False)
parser.add_argument('--fit-non-lin', action = 'store_true', help = 'try and fit non linear correction')
parser.add_argument('--tmax', type = float , default = None)


args = parser.parse_args()

x=readfile(args.input)
x["err"]=1./np.sqrt(x["ivar"])
x["nd"]=x["nd"].astype(int)
x["fiber"]=x["fiber"].astype(int)


xall={}
for k in x.keys() :
    xall[k]=np.zeros(x[k].size)
    for i in range(x[k].size) :
        xall[k][i]=x[k][i]+0

nds=np.unique(x["nd"])
fibers=np.unique(x["fiber"])
expreqs=np.unique(x["expreq"])
print("NDs=",nds)
print("fibers=",fibers)
print("EXPREQs=",expreqs)

# init
###########################################
deltat=0
ndtrans={1:1.,2:10**-0.5,3:0.1,4:0.01}
fibertrans={}
for fiber in fibers :
    fibertrans[fiber]=1.
###########################################

def compute_corrflux() :
    x["corrflux"]=np.zeros(x["flux"].shape)
    x["corrfluxivar"]=np.zeros(x["flux"].shape)
    for fiber in fibers :
        for nd in nds :
            ok=np.where((x["fiber"]==fiber)&(x["nd"]==nd))[0]
            if ok.size==0 : continue
            x["transcorr"][ok]=1/(fibertrans[fiber]*ndtrans[nd])
            x["corrflux"][ok]=x["flux"][ok]*(x["nonlin"][ok]*x["transcorr"][ok])
            x["corrfluxivar"][ok]=x["ivar"][ok]/(x["nonlin"][ok]*x["transcorr"][ok])**2
            
    nbad=np.sum(x["corrflux"]==0)
    if nbad>0 :
        print("error nbad=",nbad)
        sys.exit(12)

def u(flux) :
    return np.arctan(flux/10000.)

nonlindeg=2
def fit_nonlin_corr(deltat,x) :
    
    mcoef=np.zeros(nonlindeg+1)
    mcoef[-1]=1.
    mpol=np.poly1d(mcoef)
    nloop=5
    for loop in range(nloop) :
        coefs=[]
        chi2=0.
        ndata=0.
        for fiber in np.unique(x["fiber"]) :
            for nd in  np.unique(x["nd"]) :
                ok=np.where((x["fiber"]==fiber)&(x["nd"]==nd))[0]
                flux=x["flux"][ok]
                ivar=x["ivar"][ok]            
                truetime=x["expreq"][ok]-deltat
                trueflux=flux*mpol(u(flux))
                coef=np.sum(ivar*truetime*trueflux)/np.sum(ivar*truetime**2)            
                mflux=coef*truetime
                            
                #plt.plot(flux,flux-mflux,"o-")
                coef = np.polyfit(u(flux),mflux/flux,nonlindeg)
                err = mflux/flux**2/np.sqrt(ivar)
                coefs.append(coef)
                pol=np.poly1d(coef)
                
                #if loop == nloop-1 : plt.plot(flux,mflux/flux,"o")
                #if loop == nloop-1 : plt.errorbar(flux,mflux/flux,err,fmt="o")
                
        mcoef=np.median(coefs,axis=0)
        mcoef /= mcoef[-1]
        print(loop,mcoef)
        mpol=np.poly1d(mcoef)
        #f=np.linspace(1,np.max(x["flux"]),100)        
        #plt.plot(f,mpol(u(f)),color="gray",alpha=0.2)
    #plt.plot(f,mpol(u(f)),color="gray")
    x["nonlin"]=mpol(u(x["flux"]))
    
def fit_nonlin_corr2(deltat,x,plot=False) :
     truetime=x["expreq"]-deltat
     flux=x["flux"]
     err=1./np.sqrt(x["ivar"])
     corrflux=x["corrflux"]
     corrfluxivar=x["corrfluxivar"]
     scale=np.sum(corrfluxivar*truetime*corrflux)/np.sum(corrfluxivar*truetime**2)     
     y=scale*truetime/x["transcorr"]/flux
     ye=scale*truetime/x["transcorr"]*(err/flux**2)
     coef = np.polyfit(u(flux),y,nonlindeg)
     #print(scale,coef[-1])
     correction = coef[-1]
     coef /= correction
     y /= correction
     ye /= correction
     
     print("fit_nonlin_corr2 coef=",coef)
     pol=np.poly1d(coef)
     x["nonlin"]=pol(u(flux))
     chi2=np.sum((y-pol(u(flux)))**2/ye**2)
     ndata=y.size
     print("chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata))

     if plot :
         title="nonlin"
         if args.title is not None : title += "-%s"%args.title
         fig=plt.figure(title)
         plt.errorbar(flux,y,ye,fmt="o")
         f=np.linspace(1,np.max(x["flux"]),100)        
         plt.plot(f,pol(u(f)),color="gray")
         plt.xlabel("mean flux in spectrum")
         plt.ylabel("non-linearity correction")
         plt.grid()
         fig.savefig(title+".png")
         
def fit_fibertrans(fibertrans,fibers,deltat,x) :
    if len(fibers)<2 : return
    for fiber in fibers :
        ok=np.where((x["fiber"]==fiber))[0]
        if ok.size==0 : continue
        truetime=x["expreq"][ok]-deltat
        weight=x["corrfluxivar"][ok]
        correction =np.sum(weight*x["corrflux"][ok]*truetime)/np.sum(weight*truetime**2)
        fibertrans[fiber] *= correction
    
    st=0.
    s=0.
    for fiber in fibers :
        st += fibertrans[fiber]
        s  += 1.
    mean_fibertrans = st/s
    for fiber in fibers :    
        fibertrans[fiber] /= mean_fibertrans

def fit_ndtrans(ndtrans,nds,deltat,x) :
    if len(nds)<1 : return
        
    # for each fiber and neutral density, fit polynomial
    for nd in nds :
        ok=np.where((x["nd"]==nd))[0]
        if ok.size==0 : continue
        truetime=x["expreq"][ok]-deltat
        weight=x["corrfluxivar"][ok]
        correction=np.sum(weight*x["corrflux"][ok]*truetime)/np.sum(weight*truetime**2)
        ndtrans[nd] *= correction
    
    # normalization
    if 1 in ndtrans :
        ref=ndtrans[1]
    else :
        ref=ndtrans[2]/0.3
    for nd in nds :
        ndtrans[nd] /= ref
            
def fitall(previous_deltat=0) :
    compute_corrflux()
    deltat=previous_deltat
    old_deltat=0.
    for loop in range(100) :
        
        fit_fibertrans(fibertrans,fibers,deltat,x)
        compute_corrflux()
        fit_ndtrans(ndtrans,nds,deltat,x)
        compute_corrflux()

        if loop<5 : continue
    
        if args.fit_non_lin :
            fit_nonlin_corr2(deltat,x)
            compute_corrflux()
        
        # fit deltat , need to add weights !!!  
        coef=np.polyfit(x["expreq"],x["corrflux"],w=x["corrfluxivar"],deg=1)
        deltat=-coef[1]/coef[0]
        compute_corrflux()
        
        print("loop #%d, deltat="%loop,deltat)
    
        if np.abs(deltat-old_deltat)<0.000001 :
            break
        old_deltat=deltat
    return deltat

def fittrans() :
    compute_corrflux()    
    for loop in range(5) :        
        fit_fibertrans(fibertrans,fibers,deltat,x)
        compute_corrflux()
        fit_ndtrans(ndtrans,nds,deltat,x)
        compute_corrflux()
   

x["nonlin"]=np.ones(x["flux"].size)
x["transcorr"]=np.ones(x["flux"].size)
x["ivar"]= 1./( 1./x["ivar"] + (0.002*x["flux"])**2 )

#mask=(x["flux"]<35000) # for b1
mask=(x["flux"]<40000) # for r1

#mask=(x["expreq"]<100)
for k in x.keys() : x[k]=x[k][mask]

deltat=-0.6
deltats=[deltat]

compute_corrflux()
fit_nonlin_corr(deltat,x)

for sloop in range(4) :
    fit_nonlin_corr2(deltat,x)
    fittrans()
for sloop in range(20) :
    fit_nonlin_corr2(deltat,x)
    deltat=fitall(deltat)
    deltats.append(deltat)
    
fit_nonlin_corr2(deltat,x,plot=True)

plt.figure("converge")
plt.plot(deltats,"o")
    
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
print("#######################")



if 1 :
    title="dt"
    if args.title is not None : title += "-%s"%args.title
    fig=plt.figure(title)
    
    a0=plt.subplot(1,1,1)
    
    c=np.polyfit(x["expreq"],x["corrflux"],1,w=x["corrfluxivar"])
    a=c[0]
    for nd in nds :
        ok=np.where(x["nd"]==nd)[0]
        flux=x["nonlin"][ok]*x["flux"][ok]
        oflux=x["flux"][ok]
        err=x["err"][ok]
        corr = x["corrflux"][ok]/flux
        exptime = x["expreq"][ok]
        dt  = corr/a*flux - exptime
        dte = corr/a*err
        a0.errorbar(flux,dt,dte,fmt="o")
        
        
    a0.set_xlabel("mean flux in spectra")
    a0.set_ylabel("exposure time correction (sec)")
    #a1.set_ylabel("flux/model")
    a0.grid()
    a0.set_ylim([0,1])
    
    fig.savefig(title+".png")
    
    plt.show()
    exit(12)
