import matplotlib.pyplot as plt
import sys
import numpy as np

from desiutil.log import get_logger

def parse_fibers(fiber_string) :
    
    if fiber_string is None :
        return None
    
    fibers=[]
    for sub in fiber_string.split(',') :
        if sub.isdigit() :
            fibers.append(int(sub))
            continue
        
        tmp = sub.split(':')
        if ((len(tmp) is 2) and tmp[0].isdigit() == True and tmp[1].isdigit() == True) :
            for f in range(int(tmp[0]),int(tmp[1])) :
                fibers.append(f)
        else :
            log.error("--fibers parsing error.\nCorrect format is either  : --fibers=begin,end (excluded)\nand/or  : --fibers=begin:end (excluded)\nYou can use : --fibers=2,5,6:8,3,10")
            sys.exit(1)
    return np.array(fibers)



def plot_graph(frame, fibers, opt_err=False, opt_2d=False, label = None, subplot=None) :
    """Plot graph from a given spectra from a fits file and returns figure
    
    ----------
    Parameters
    ----------

    frame : File Directory
    Where the spectra is collected to be plot.

    fibers : fibers to show
    """

    log         = get_logger()
    spectra     = frame["FLUX"].data
    ivar        = frame["IVAR"].data
    
    if "MASK" in frame :
        ivar *= (frame["MASK"].data==0)

    wave        = frame["WAVELENGTH"].data
    nfibers     = spectra.shape[0]

    if np.max(fibers) >= nfibers :
        log.warning("requested fiber numbers %s exceed number of fibers in file %d"%(str(fibers),nfibers))
        fibers=fibers[fibers<nfibers]
    
    if subplot is None :
        subplot  = plt.subplot(1,1,1)
    
    for fiber in fibers :
        ok = ivar[fiber]>0
        if label :
            fiber_label = "%s Fiber #%d"%(label,fiber)
        else :
            fiber_label="Fiber #%d"%fiber
            
        log.debug("Plotting fiber %03d" % fiber)
        if opt_err :
            err = np.sqrt(1./ (ivar[fiber] + (ivar[fiber] == 0))) * (ivar[fiber] > 0)
            if len(wave.shape) > 1 :
                subplot.errorbar(wave[fiber][ok], spectra[fiber][ok], err[ok], fmt="o-", label=fiber_label)
            else :
                subplot.errorbar(wave[ok], spectra[fiber][ok], err[ok], fmt="o-",label=fiber_label)
        else :
            #color="b"
            #if fiber>=10 :
                #color="r"
            if len(wave.shape) > 1 :
                subplot.plot(wave[fiber][ok], spectra[fiber][ok], "-",label=fiber_label)
            else :
                subplot.plot(wave[ok], spectra[fiber][ok], "-",label=fiber_label)
    
    subplot.set_xlabel("Wavelength [A]")
    
    if opt_2d :
        title="spectra"
        if label is not None:
            title = label
        plt.figure(title)
        if len(wave.shape) == 1 :
            plt.imshow(spectra[fibers].T,
                       aspect = 'auto',
                       extent = (fibers[0] - 0.5, fibers[-1] + 0.5, wave[0], wave[-1]),
                       origin = 0.,
                       interpolation = "nearest")
            plt.ylabel("Wavelength [A]")
            plt.xlabel("Fiber #")
        else :
            plt.imshow(spectra[fibers].T,
                       aspect = 'auto',
                       extent = (fibers[0]-0.5, fibers[-1]+0.5, 0,spectra.shape[1]),
                       origin = 0.,
                       interpolation = "nearest")
            plt.ylabel("Y CCD")
            plt.xlabel("Fiber #")
        plt.colorbar()
    

