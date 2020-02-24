#!/usr/bin/env python


import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import fitsio
from desispec.io import read_spectra
from desispec.interpolation import resample_flux
from prospect import plotframes
from astropy.table import Table
import redrock.templates

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*", 
                    help = 'path to spectra file(s)')
parser.add_argument('-t','--targetid', type = int, default = None, required = False, 
                    help = 'plot specific targetid')
parser.add_argument('--rebin',type = int, default = None, required = False, 
                    help = 'rebin')
parser.add_argument('--zbest',type = str, default = None, required = False, 
                    help = 'zbest file')
parser.add_argument('--spectype',type = str, default = None, required = False, 
                    help = 'spectype to select')


args        = parser.parse_args()

if  args.zbest is not None :
    #- Load redrock templates
    templates = dict()
    for filename in redrock.templates.find_templates():
        tx = redrock.templates.Template(filename)
        templates[(tx.template_type, tx.sub_type)] = tx

spectra = []
for filename in args.infile :
    spectra.append(read_spectra(filename))

targetids=None
if ( targetids is None ) and ( args.targetid is not None ) :
    targetids=[args.targetid,]

zbest=None
if ( targetids is None ) and ( args.zbest is not None ) and ( args.spectype is not None ):
    zbest=Table.read(args.zbest,"ZBEST")
    print(zbest.dtype.names)
    print(np.unique(zbest["SPECTYPE"]))
    #selection = np.where((zbest["SPECTYPE"]==bytes(args.spectype, 'utf-8'))&(zbest["ZWARN"]==0))[0]
    selection = np.where((zbest["SPECTYPE"]==args.spectype)&(zbest["ZWARN"]==0))[0]
    targetids=np.unique(spectra[0].fibermap["TARGETID"][selection])
    
if targetids is None :        
    targetids=np.unique(spectra[0].fibermap["TARGETID"])


lines = {
    'Halpha'      : 6562.8,
    'Hbeta'       : 4862.68,
    'MgII(2804)'  : 2803.5324,
    'CIII(1909)'  : 1909.,
    'CIV(eff)'    : 1549.06,
    'SiIV(1394)'  : 1393.76018,
    'LYA'         : 1215.67,
    'LYB'         : 1025.72
}


           
for tid in targetids :
    line="TARGETID={}".format(tid)

    model_flux=dict()
    if zbest is not None :
        j=np.where(zbest["TARGETID"]==tid)[0][0]
        line += " ZBEST={} SPECTYPE={} ZWARN={}".format(zbest["Z"][j],zbest["SPECTYPE"][j],zbest["ZWARN"][j])

        tx = templates[(zbest['SPECTYPE'][j], zbest['SUBTYPE'][j])]
        for band in spectra[0].bands:
            model_flux[band] = np.zeros(spectra[0].wave[band].shape)
            coeff = zbest['COEFF'][j][0:tx.nbasis]
            model = tx.flux.T.dot(coeff).T
            print(spectra[0].wave[band].shape,tx.wave.shape,model.shape)
            mx = resample_flux(spectra[0].wave[band], tx.wave*(1+zbest['Z'][j]), model)
            k=np.where(spectra[0].fibermap["TARGETID"]==tid)[0][0]
            model_flux[band] = spectra[0].R[band][k].dot(mx)
            
    fig=plt.figure(figsize=[10,6])    
    ax = fig.add_subplot(111)   
    print(line)
    for spec in spectra :
        jj=np.where(spec.fibermap["TARGETID"]==tid)[0]
        for k in ['FIRST_EXPID','LAST_EXPID','NUM_EXPID','EXPID'] :
            if k in spec.fibermap.colnames :
                for j in jj :
                    print("{} {}".format(k,spec.fibermap[k][j]))

        for j in jj :
            for b in spec._bands :
                
                i=np.where(spec.ivar[b][j]*(spec.mask[b][j]==0)>0)[0]
                
                if args.rebin is not None and args.rebin>0:
                    rwave=np.linspace(spec.wave[b][0],spec.wave[b][-1],spec.wave[b].size//args.rebin)
                    rflux,rivar = resample_flux(rwave,spec.wave[b],spec.flux[b][j],ivar=spec.ivar[b][j]*(spec.mask[b][j]==0))
                    plt.plot(rwave,rflux)
                else :
                    plt.plot(spec.wave[b][i],spec.flux[b][j,i])
                
                c=np.polyfit(spec.wave[b][i],spec.flux[b][j,i],3)
                pol=np.poly1d(c)(spec.wave[b][i])
                print("mean chi2=",np.mean(spec.ivar[b][j,i]*(spec.flux[b][j,i]-pol)**2))

    if zbest is not None :
        for band in spectra[0].bands:
            plt.plot(spectra[0].wave[band],model_flux[band],"-",alpha=0.6)
            for elem in lines :
                line=(1+zbest['Z'][j])*lines[elem]
                if line>spectra[0].wave[band][0] and line<spectra[0].wave[band][-1] :
                    plt.axvline(line,color="red",linestyle="--",alpha=0.4)
                    y=np.interp(line,spectra[0].wave[band],model_flux[band])
                    #plt.plot([line,line],[0,y],"--",alpha=0.8,color="red")
                    plt.text(line+100,y*1.1,elem,color="red")
    plt.xlabel("wavelength [A]")
    props = dict(boxstyle='round', facecolor='yellow', alpha=0.2)
    bla="TID = {}\n".format(zbest['TARGETID'][j])
    bla+="Z  = {:4.3f}".format(zbest['Z'][j])
    plt.text(0.9,0.9,bla,fontsize=12, bbox=props,transform = ax.transAxes,verticalalignment='top', horizontalalignment='right')
    plt.tight_layout()
    plt.show()


plt.show()


