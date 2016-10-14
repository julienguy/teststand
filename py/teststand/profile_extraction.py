from desispec.log import get_logger

spectrum_ivar = np.sum(ivar[y,x1:x2]*prof**2)
spectrum=np.sum(ivar[y,x1:x2]*image[y,x1:x2]*prof)/spectrum_ivar
x=np.arange(x1,x2)
prof = exp(-(x-xc)**2/2/sigmax**2)
prof /= np.sum(prof)