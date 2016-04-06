#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import argparse
import sys
import os.path

def read_known_lines() :
    # find data file in this repo
    current_file_path=os.path.realpath(__file__)
    known_lines_filename="%s/../data/main-peaks.data"%os.path.dirname(current_file_path)
    if not os.path.isfile(known_lines_filename) :
        print "failed to find file main-peaks.data from program path '%s'"%current_file_path
        print "was trying to open '%s'"%known_lines_filename
        sys.exit(12)
    file=open(known_lines_filename)
    wave=[]
    element=[]
    for line in file.readlines() :
        if line[0]=="#" :
            continue
        vals=line.strip().split(" ")
        if len(vals)!=2 :
            print "ignore line '%s'"%line.strip()
        try :
            w=float(vals[0])
            e=vals[1]
            wave.append(w)
            element.append(e)
        except :
            print "failed to read properly line '%s'"%line.strip()
    known_lines={}
    for ee in np.unique(element) :
        waves=[]
        for w,e in zip(wave,element) :
            if e==ee :
                waves.append(w)
        known_lines[ee]=waves

    return known_lines

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
parser.add_argument('--tolerance', type = float, default = 30., required=False,
                        help = 'tolerance for matching bright peaks in A')

args = parser.parse_args()

lamps=args.lamps.split(",")
if len(lamps)<1 :
    print "error when parsing '%s'"%args.lamps
    sys.exit(12)

known_lines = read_known_lines()


hdulist=fits.open(args.infile)
table=hdulist[1].data
wave=table["wavelength"]*10. # nm -> A

line_wave=None
line_flux=None
matched_line_wave=None
matched_line_flux=None


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

    lamp_line_wave=wave[peaks]
    lamp_line_flux=flux[peaks]
    
    # now try to match this to known lines
    matched_lamp_line_wave=[]
    matched_lamp_line_flux=[]
    
    # get elements we can find in this lamp
    elements=[]
    for e in known_lines :
        if lamp.find(e)>=0 :
            elements.append(e)
    print "for lamp '%s' : elements = '%s'"%(lamp,str(elements))
    
    for w,f in zip(lamp_line_wave,lamp_line_flux) :
        found=False
        for e in elements :
            i=np.argmin(np.abs(known_lines[e]-w))
            dist=abs(known_lines[e][i]-w)
            if dist<args.tolerance :
                print "found match %f <-> %f for %s"%(w,known_lines[e][i],e)
                matched_lamp_line_wave.append(known_lines[e][i])
                matched_lamp_line_flux.append(f)
                found=True
                break
        if not found :
            print "line at %f in lamp %s not found in list of known bright lines"%(w,lamp)
                
    


    if line_wave is None :
        line_wave=lamp_line_wave
        line_flux=lamp_line_flux
        matched_line_wave=matched_lamp_line_wave
        matched_line_flux=matched_lamp_line_flux
    else :
        line_wave=np.append(line_wave,lamp_line_wave)
        line_flux=np.append(line_flux,lamp_line_flux)
        matched_line_wave=np.append(matched_line_wave,matched_lamp_line_wave)
        matched_line_flux=np.append(matched_line_flux,matched_lamp_line_flux)
    
    # evaluate the integral based on the brightest line
    i=np.argmax(lamp_line_flux)
    ii=np.where(flux==flux[peaks[i]])[0]
    hw=8
    ib=max(0,ii-hw)
    ie=min(flux.size,ii+hw)
    if args.plot :
        pylab.plot(wave[ii],flux[ii],"o",c="gray",ms=12)
        pylab.plot(wave[ib:ie],flux[ib:ie],"-",c="gray")
    totflux=np.sum(flux[ib:ie])
    integral_to_maxval.append(totflux/flux[peaks[i]])

if args.plot :
    pylab.plot(line_wave,line_flux,"o",c="r",alpha=0.4)
    pylab.plot(matched_line_wave,matched_line_flux,"o",c="g")

integral_to_maxval = np.median(integral_to_maxval)

print "integral_to_maxval = ",integral_to_maxval
print "scale              = ",args.scale

matched_line_flux *= integral_to_maxval # now estimate of the integrated flux in lines
matched_line_flux *= args.scale # arbitrary scale to get something decent

# sort according to wave
ii=np.argsort(matched_line_wave)
matched_line_wave=matched_line_wave[ii]
matched_line_flux=matched_line_flux[ii]


# write this
cols=[]
cols.append(fits.Column(name='WAVE', format='D', array=matched_line_wave))
cols.append(fits.Column(name='ELECTRONS', format='D', array=matched_line_flux))
cols = fits.ColDefs(cols)
hdulist=fits.HDUList([fits.PrimaryHDU()])
hdulist.append(fits.BinTableHDU.from_columns(cols))
hdulist.writeto(args.outfile,clobber=True)

print "wrote",args.outfile

if args.plot :
    pylab.xlabel("wavelength")
    pylab.ylabel("flux [??]")
    pylab.show()   
