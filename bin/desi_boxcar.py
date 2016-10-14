#!/usr/bin/env python

import argparse
import astropy.io.fits as pyfits

from teststand.boxcar_extraction   import boxcar
from teststand.resample            import resample_to_same_wavelength_grid
from teststand.graph_tools         import plot_graph,show_graph
from desispec.log                       import get_logger

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'path of psf fits file to get wavelength from')
parser.add_argument('-i','--image', type = str, default = None, required = True,
                    help = 'path of image fits file')
parser.add_argument('-o','--outframe', type = str, default = None, required = True,
                    help = 'path of output frame file')
parser.add_argument('-n','--nfibers', type = int, default = None, required = False,
                    help = 'number of fibers (default=all)')
parser.add_argument('--show', action='store_true',
                    help = 'plot result')
parser.add_argument('-r','--resample', action='store_true',
                    help = 'resample to save wavelength grid')
parser.add_argument('--plot', type=str,
                    help = 'defines from_to which fiber to work on. (ex: --plot=50,60 means that only fibers 50 to 60 included will be plotted')

log         = get_logger()
args        = parser.parse_args()
psf         = pyfits.open(args.psf)
image_file  = pyfits.open(args.image)

spectra, ivar, wave = boxcar(psf, image_file, args.nfibers)

if args.resample :
    log.info("Starting resampling...")
    spectra, ivar, wave = resample_to_same_wavelength_grid(spectra, ivar, wave)
    log.info("Data resampled.")

hdulist = pyfits.HDUList([pyfits.PrimaryHDU(spectra),
                        pyfits.ImageHDU(ivar,name="IVAR"),
                        pyfits.ImageHDU(wave,name="WAVELENGTH")])
                        #pyfits.ImageHDU(rdata, name="RESOLUTION")])
hdulist.writeto(args.outframe,clobber=True)

if args.plot or args.show :
    frame = pyfits.open(args.outframe)
    if not (args.plot) :
        plot_graph(frame)
    else :
        nb = args.plot.split(',')
        if ((len(nb) is 2) or (nb[0].isdigit == True and nb[1].isdigit == True)) :
            plot_graph(frame, args.nfibers, nb[0], nb[1])
        else :
            log.error("--plot parsing error. correct format is : --plot=nb_from,nb_to")
    log.info("Showing.......")
    show_graph()

log.info("Script done")
