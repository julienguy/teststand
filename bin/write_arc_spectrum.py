#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import argparse
import sys

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required=True,
                        help = 'path to file "DESI_OGSE_SpecSources.fits"')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,
                        help = 'path to output file ')
parser.add_argument('--lamps', type = str, default = None, required=True,
                        help = 'Ne,Kr,Cd,HgAr (one or several, separated by ",")')
parser.add_argument('--scale', type = float, default = 1e11, required=False,
                        help = 'scale factor')
parser.add_argument('--threshold', type = float, default = 7., required=False,
                        help = 'threshold for finding peaks')
parser.add_argument('--plot', action="store_true",
                        help = 'plot spectrum')

args = parser.parse_args()

lamps=args.lamps.split(",")
if len(lamps)<1 :
    print "error when parsing '%s'"%args.lamps
    sys.exit(12)


hdulist=fits.open(args.infile)
table=hdulist[1].data
wave=table["wavelength"]*10. # nm -> A

line_wave=None
line_flux=None


integral_to_maxval=[]

for lamp in lamps :
    try :
        flux=table[lamp]
    except KeyError :
        print "cannot find this lamp '%s' in %s"%(lamp,hdulist[1].columns.names)
        sys.exit(12)
        
    # need to find peaks, to avoid issue of resolution
    # do it very very crudly    
    diff_plus=flux-np.roll(flux,1)
    diff_minus=flux-np.roll(flux,-1)
    ii=np.arange(flux.size)
    margin=10 # avoid doubtful edges
    peaks=np.where((diff_plus>0)&(diff_minus>0)&(ii>margin)&(ii<flux.size-margin))[0]
    # apply a threshold 
    med_all=np.median(flux)
    med_peak=np.median(flux[peaks])
    threshold=args.threshold*(med_peak-med_all)+med_all
    peaks=np.where((diff_plus>0)&(diff_minus>0)&(ii>margin)&(ii<flux.size-margin)&(flux>threshold))[0]
    
    if args.plot :
        import pylab
        pylab.plot(wave,flux,c="b",alpha=0.3)
        #pylab.plot(wave,med_all*np.ones(wave.shape))
        #pylab.plot(wave,threshold*np.ones(wave.shape))
        #pylab.plot(wave[peaks],flux[peaks],"o")
        #pylab.show()

    if line_wave is None :
        line_wave=wave[peaks]
        line_flux=flux[peaks]
    else :
        line_wave=np.append(line_wave,wave[peaks])
        line_flux=np.append(line_flux,flux[peaks])
    
    # evaluate the integral based on the brightest line
    i=np.argmax(flux[peaks])
    ii=np.where(flux==flux[peaks[i]])[0]
    hw=8
    ib=max(0,ii-hw)
    ie=min(flux.size,ii+hw)
    if args.plot :
        pylab.plot(wave[ii],flux[ii],"o",c="k",ms=12)
        pylab.plot(wave[ib:ie],flux[ib:ie],"-",c="k")
    totflux=np.sum(flux[ib:ie])
    integral_to_maxval.append(totflux/flux[peaks[i]])

if args.plot :
    pylab.plot(line_wave,line_flux,"o",c="r")
integral_to_maxval = np.median(integral_to_maxval)

print "integral_to_maxval = ",integral_to_maxval
print "scale              = ",args.scale

line_flux *= integral_to_maxval # now estimate of the integrated flux in lines
line_flux *= args.scale # arbitrary scale to get something decent

# sort according to wave
ii=np.argsort(line_wave)
line_wave=line_wave[ii]
line_flux=line_flux[ii]


# write this
cols=[]
cols.append(fits.Column(name='WAVE', format='D', array=line_wave))
cols.append(fits.Column(name='ELECTRONS', format='D', array=line_flux))
cols = fits.ColDefs(cols)
hdulist=fits.HDUList([fits.PrimaryHDU()])
hdulist.append(fits.BinTableHDU.from_columns(cols))
hdulist.writeto(args.outfile,clobber=True)

print "wrote",args.outfile

if args.plot :
    pylab.xlabel("wavelength")
    pylab.ylabel("flux [??]")
    pylab.show()   
