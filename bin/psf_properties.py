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
    if psftype=="GAUSS-HERMITE" :
        return specter.psf.GaussHermitePSF(filename)
    elif psftype=="SPOTGRID" :
        return specter.psf.SpotGridPSF(filename)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf', type = str, nargs = "*", default = None, required = True,
                    help = 'path of psf files')

parser.add_argument('--refpsf', type = str, default = None, required = True,
                    help = 'reference psf')

parser.add_argument('--fibers', type = str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('--wave', type = int, default = None, required = True,
                    help = 'wavelength')


args        = parser.parse_args()
log = get_logger()

refpsf=readpsf(args.refpsf)
        
#wmin=refpsf[0]._wmin_all
#wmax=refpsf[0]._wmax_all
#nw=20
#waves=np.linspace(wmin+(wmax-wmin)/nw/2.,wmax-(wmax-wmin)/nw/2.,nw).astype(int)

fibers=parse_fibers(args.fibers)
if fibers is None :
        fibers = np.arange(refpsf.nspec)

use_trace = True

reftx={} # x from trace
refty={} # y from trace
refcx={} # barycenter of stamp
refcy={} # barycenter of stamp
refsx={} # sigma x
refsy={} # sigma y
refimage={} # stamp
xpix=None
ypix=None

for fiber in fibers :
    xx, yy, img = refpsf.xypix(fiber,args.wave)
    if xpix is None :
        # check with visual inspection of ccd image (?)
        xpix=np.tile(np.arange(img.shape[0]),(img.shape[1],1))    
        ypix=np.tile(np.arange(img.shape[1]),(img.shape[0],1)).T

    refimage[fiber]=img
    refcx[fiber]= np.sum(xpix*img)/np.sum(img)
    refcy[fiber]= np.sum(ypix*img)/np.sum(img)
    refsx[fiber]= np.sqrt(np.sum((xpix-refcx[fiber])**2*img)/np.sum(img))
    refsy[fiber]= np.sqrt(np.sum((ypix-refcy[fiber])**2*img)/np.sum(img))
    reftx[fiber]= refpsf.x(fiber,args.wave)
    refty[fiber]= refpsf.y(fiber,args.wave)
    psfs=[]

for filename in args.psf :

    ofilename="propertie-%s"%(os.path.basename(filename).replace(".fits",".txt"))
    if os.path.isfile(ofilename) :
        continue
    
    log.info("reading %s"%filename)
    psf = readpsf(filename)
    
    ofile = open(ofilename,"w")
    ofile.write("# EXPNUM CAMID FIBER WAVE DX DY CX CY SX SY EBIAS\n")
    
    cam    = os.path.basename(filename).split("-")[1][0]
    if cam=="b" : camid=0
    elif cam=="r" : camid=1
    elif cam=="z" : camid=2
    
    expnum = int(os.path.basename(filename).split("-")[2].replace(".fits",""))
    
    for fiber in fibers :

                
        tx=psf.x(fiber,args.wave)
        ty=psf.y(fiber,args.wave)
        dtx=tx-reftx[fiber]
        dty=ty-refty[fiber]
        
        psf._cache={} # reset cache !!
        psf.coeff['X']._coeff[fiber][0] -= dtx
        psf.coeff['Y']._coeff[fiber][0] -= dty
        
        xx, yy, img = psf.xypix(fiber,args.wave)
        cx = np.sum(xpix*img)/np.sum(img)
        cy = np.sum(ypix*img)/np.sum(img)
        sx = np.sqrt(np.sum((xpix-cx)**2*img)/np.sum(img))
        sy = np.sqrt(np.sum((ypix-cy)**2*img)/np.sum(img))
        ebias = np.sum(refimage[fiber]*img)/np.sum(img**2)-1.
        
        line = "%d %d %d %d"%(expnum,camid,fiber,args.wave)
        line += " %f %f %f %f %f %f %f"%(dtx,dty,cx,cy,sx,sy,ebias)

        #print(line)
        sys.stdout.flush()
        ofile.write("%s\n"%line)

    ofile.close()
