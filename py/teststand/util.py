import matplotlib.pyplot as plt
import sys
import numpy as np

from desiutil.log import get_logger
from desispec.interpolation import resample_flux
from desispec.qproc.util import parse_fibers

def plot(frame, fibers, opt_err=False, opt_2d=False, label = None, subplot=None,dwave=None) :
    """Plot graph from a given spectra from a fits file and returns figure
    
    ----------
    Parameters
    ----------

    frame : File Directory
    Where the spectra is collected to be plot.

    fibers : fibers to show
    """

    log         = get_logger()
    flux     = frame["FLUX"].data
    ivar     = frame["IVAR"].data

    nfibers     = flux.shape[0]
    if np.max(fibers) >= nfibers :
        log.warning("requested fiber numbers %s exceed number of fibers in file %d"%(str(fibers),nfibers))
        fibers=fibers[fibers<nfibers]
    nfibers=len(fibers)
    flux=flux[fibers]
    ivar=ivar[fibers]
    if ("MASK" in frame) :
        ivar *= (frame["MASK"].data[fibers]==0)
        
    wave        = frame["WAVELENGTH"].data
    
    if dwave is not None :
        minwave=np.min(wave)
        maxwave=np.max(wave)
        nwave = np.linspace(minwave,maxwave,int((maxwave-minwave)/dwave))
        nflux = np.zeros((flux.shape[0],nwave.size))
        nivar = np.zeros((flux.shape[0],nwave.size))
        
        for fiber in range(nfibers) :
            nflux[fiber],nivar[fiber] = resample_flux(nwave,wave,flux[fiber],ivar=ivar[fiber])
        wave = nwave
        flux = nflux
        ivar = nivar

        
    if subplot is None :
        subplot  = plt.subplot(1,1,1)
    
    for fiber in range(nfibers) :
        ok = ivar[fiber]>0
        if label :
            fiber_label = "%s Fiber #%d"%(label,fiber)
        else :
            fiber_label="Fiber #%d"%fiber
            
        log.debug("Plotting fiber %03d" % fiber)
        if opt_err :
            err = np.sqrt(1./ (ivar[fiber] + (ivar[fiber] == 0))) * (ivar[fiber] > 0)
            if len(wave.shape) > 1 :
                subplot.errorbar(wave[fiber][ok], flux[fiber][ok], err[ok], fmt="o-", label=fiber_label)
            else :
                subplot.errorbar(wave[ok], flux[fiber][ok], err[ok], fmt="o-",label=fiber_label)
        else :
            #color="b"
            #if fiber>=10 :
                #color="r"
            if len(wave.shape) > 1 :
                subplot.plot(wave[fiber][ok], flux[fiber][ok], "-",label=fiber_label)
            else :
                subplot.plot(wave[ok], flux[fiber][ok], "-",label=fiber_label)
    
    subplot.set_xlabel("Wavelength [A]")
    
    if opt_2d :
        title="flux"
        if label is not None:
            title = label
        plt.figure(title)
        if len(wave.shape) == 1 :
            plt.imshow(flux[fibers].T,
                       aspect = 'auto',
                       extent = (fibers[0] - 0.5, fibers[-1] + 0.5, wave[0], wave[-1]),
                       origin = 0.,
                       interpolation = "nearest")
            plt.ylabel("Wavelength [A]")
            plt.xlabel("Fiber #")
        else :
            plt.imshow(flux[fibers].T,
                       aspect = 'auto',
                       extent = (fibers[0]-0.5, fibers[-1]+0.5, 0,flux.shape[1]),
                       origin = 0.,
                       interpolation = "nearest")
            plt.ylabel("Y CCD")
            plt.xlabel("Fiber #")
        plt.colorbar()
    

