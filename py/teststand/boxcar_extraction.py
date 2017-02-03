import numpy as np
from numpy.polynomial.legendre import legval, legfit

from desispec.log import get_logger

################
#   RETURNS FITS FILE INCLUDING ELECTRONS QUANTITY
################

def boxcar(psf, image_file, fibers=None, width=7, side_bands=False) :
    """Find and returns  wavelength  spectra and inverse variance

        ----------
        Parameters
        ----------

        psf     : File Descriptor
        Where the wavelength is collected.

        pix     : File Descriptor
        Interpreted photons using the wavelength.

        fibers : Optional. If left empty, will extract all fibers.

        -------
        Returns
        -------

        spectra, ivar, wavelength

        """
    log=get_logger()
    log.info("Starting boxcar extraction...")

    # it is a boot or specex psf ?
    psftype=psf[0].header["PSFTYPE"]
    log.info("psf is a '%s'"%psftype)
    if psftype == "bootcalib" :    
        wavemin = psf[0].header["WAVEMIN"]
        wavemax = psf[0].header["WAVEMAX"]
        xcoef   = psf[0].data
        ycoef   = psf[1].data
        xsig    = psf[2].data
    elif psftype == "GAUSS-HERMITE" :
        table=psf[1].data        
        i=np.where(table["PARAM"]=="X")[0][0]
        wavemin=table["WAVEMIN"][i]
        wavemax=table["WAVEMAX"][i]
        xcoef=table["COEFF"][i]
        i=np.where(table["PARAM"]=="Y")[0][0]
        ycoef=table["COEFF"][i]
        i=np.where(table["PARAM"]=="GHSIGX")[0][0]
        xsig=table["COEFF"][i]

    if fibers is None :
        fibers = np.arange(xcoef.shape[0])
    
    log.info("wavelength range : [%f,%f]"%(wavemin,wavemax))
    
    
    log.debug("xcoef.shape",xcoef.shape)
    log.debug("ycoef.shape",ycoef.shape)
    

    flux        = image_file[0].data
    #   Inverse variance of the image's value
    flux_ivar   = image_file["IVAR"].data
    #   Use masked pixels 
    flux_mask   = image_file["MASK"].data
    flux_ivar   *= (flux_mask==0)
    
    #   Variance based on inverse variance's size
    flux_var    = np.zeros(flux_ivar.shape)

    #   Applying a mask that keeps positive value to get the Variance by inversing the inverse variance.
    mask        = (flux_ivar > 0)
    flux_var[mask] = 1./flux_ivar[mask]

    #   Number of pixels in an image 
    #   We are going to extract one flux per fiber per Y pixel (total = nfibers x npix_y)
    npix_y  = flux.shape[0]
    npix_x  = flux.shape[1]
    
    nfibers = xcoef.shape[0]
    if np.max(fibers) >= nfibers :
        log.warning("requested fiber numbers %s exceed number of fibers in file %d"%(str(fibers),nfibers))
        ii=np.where(fibers<nfibers)
        fibers=fibers[ii]
    
    #   Flux as a function of wavelength
    spectra             = np.zeros((fibers.size,npix_y))
    #   Inverse-variance of spectrum
    spectra_ivar        = np.zeros((fibers.size,npix_y))
    #   Wavelength
    wave_of_y           = np.zeros((fibers.size, npix_y))

###
# Using legendre's polynomial to get a spectrum per fiber
###

    for f,fiber in enumerate(fibers) :
        log.info("extracting fiber #%03d"%fiber)
        x1_of_y, x2_of_y, tmp_wave = invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, fiber, npix_y, width)
        wave_of_y[f] = tmp_wave
        
        hw=width//2
        if not side_bands :
            for y in range(npix_y) :
                #   Checking if there's a dead pixel
                nb_invalidPix   = np.sum(flux_ivar[y, x1_of_y[y]:x2_of_y[y]] <= 0)
                if nb_invalidPix == 0 :
                    #   Sum of flux
                    spectra[f, y]   = np.sum(flux[y, x1_of_y[y]:x2_of_y[y]])
                    #   Sum of variance
                    var                 = np.sum(flux_var[y, x1_of_y[y]:x2_of_y[y]])
                    #   Spectrum of inverse variance
                    spectra_ivar[f, y] = 1./var
        else :
           for y in range(npix_y) :
               #   Checking if there's a dead pixel
               
               nb_invalidPix   = np.sum(flux_ivar[y, x1_of_y[y]-hw:x2_of_y[y]+hw+1] <= 0)
               if nb_invalidPix == 0 :
                   #   Sum of flux
                   spectra[f, y]   = np.sum(flux[y, x1_of_y[y]:x2_of_y[y]]) -  np.sum(flux[y, x1_of_y[y]-hw:x1_of_y[y]]) - np.sum(flux[y, x2_of_y[y]:x2_of_y[y]+hw+1])
                   #   Sum of variance
                   var                 = np.sum(flux_var[y, x1_of_y[y]-hw:x2_of_y[y]+hw+1])
                   #   Spectrum of inverse variance
                   spectra_ivar[f, y] = 1./var 


    log.info("Boxcar extraction complete")
    return spectra, spectra_ivar, wave_of_y

def u(wave, wavemin, wavemax) :
    return 2. * (wave - wavemin)/(wavemax - wavemin) - 1.

def invert_legendre_polynomial(wavemin, wavemax, ycoef, xcoef, fiber, npix_y, width=7) :
 
    #   Wavelength array used in 'invert_legendre_polynomial'
    wave                = np.linspace(wavemin, wavemax, 100)
    #   Determines value of Y, so we can know its coeficient and then its position
    y_of_wave           = legval(u(wave, wavemin, wavemax), ycoef[fiber])
    coef                = legfit(u(y_of_wave, 0, npix_y), wave, deg=ycoef[fiber].size)
    wave_of_y           = legval(u(np.arange(npix_y).astype(float), 0, npix_y), coef)
    #   Determines wavelength intensity (x) based on Y
    x_of_y              = legval(u(wave_of_y, wavemin, wavemax), xcoef[fiber])
    #   Ascertain X by using low and high uncertainty
    x1_of_y             = (np.floor(x_of_y).astype(int) - width//2).astype(int)
    x2_of_y             = (np.floor(x_of_y).astype(int) + width//2 + 1).astype(int)
    return (x1_of_y, x2_of_y, wave_of_y)
