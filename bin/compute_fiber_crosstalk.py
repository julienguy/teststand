#!/usr/bin/env python



"""
Run DESI qproc on a given exposure
"""


import argparse
import sys,os
import time
import numpy as np
import astropy.io.fits as fits
from desiutil.log import get_logger
from desispec.util import runcmd
from desispec.io import read_raw,read_image,read_fibermap,write_image,write_fiberflat,read_fiberflat
from desispec.io.fluxcalibration import read_average_flux_calibration
from desispec.io.xytraceset import read_xytraceset
from desispec.calibfinder import CalibFinder
from desispec.qproc.io import read_qframe,write_qframe
from desispec.qproc.qextract import qproc_boxcar_extraction
from desispec.qproc.qfiberflat import qproc_apply_fiberflat,qproc_compute_fiberflat
from desispec.qproc.qsky import qproc_sky_subtraction
from desispec.qproc.util import parse_fibers

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description="""Quick extraction and processing script. The input image file is either a raw data fits file from the ICS (in which case one must specify the camera) or a preprocessed image. The output is a frame (a series of spectra) that can be saved to disk and/or displayed. An approximate fiber flat field correction can be computed if the input is a dome screen exposure. For on-sky images, fiber flat field and a rudimentary sky subtraction can be performed. This script relies on the existence of a DESI_SPECTRO_CALIB (older version used DESI_CCD_CALIBRATION_DATA) environment variable pointing to a local copy of the DESI calibration SVN repository.
""",
                                 epilog="""Example: desi_qproc -i desi-00003577.fits --camera r1 --fibermap fibermap-00003577.fits --fibers 12:15 --plot"""
)
parser.add_argument('-i','--image', type = str, default = None, required = True,
                    help = 'path to a preprocessed image fits file')
parser.add_argument('-p','--psf', type = str, default = None, required = True,
                    help = 'PSF file to get the traces')
parser.add_argument('-f','--fibers', type = str, default = None, required = False, help = 'selection of fibers (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be extracted)')
parser.add_argument('--width', type=int, default=3, required=False,
                    help = 'extraction band width (in pixels)')
parser.add_argument('--dist', type=int, default=7, required=False,
                    help = 'distance of center of band from spectral trace (in pixels)')

t0   = time.time()
log  = get_logger()
args = parser.parse_args()

