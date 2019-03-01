import fitsio
import scipy as sp
import argparse
import matplotlib.pyplot as plt

from desispec.io.xytraceset import read_xytraceset
from desispec.calibfinder import CalibFinder

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Plot the PSF image and slice along a constant wavelength')

    parser.add_argument('--b-path', type=str, default=None, required=True,
        help='Path to b-camera preproc file')

    parser.add_argument('--r-path', type=str, default=None, required=True,
        help='Path to r-camera preproc file')

    parser.add_argument('--z-path', type=str, default=None, required=True,
        help='Path to z-camera preproc file')

    parser.add_argument('--b-line-value', type=float, default=4359.56, required=False,
        help='Value of the emission line for the b-camera')

    parser.add_argument('--r-line-value', type=float, default=7247.1631, required=False,
        help='Value of the emission line for the r-camera')

    parser.add_argument('--z-line-value', type=float, default=8379.9093, required=False,
        help='Value of the emission line for the z-camera')

    parser.add_argument('--b-line-name', type=str, default='HgI', required=False,
        help='Name of the emission line for the b-camera')

    parser.add_argument('--r-line-name', type=str, default='NeI', required=False,
        help='Name of the emission line for the r-camera')

    parser.add_argument('--z-line-name', type=str, default='NeI', required=False,
        help='Name of the emission line for the z-camera')

    parser.add_argument('--b-psf-path', type=str, default=None, required=False,
        help='Path to b-camera psf file (go to $DESI_SPECTRO_CALIB if not given)')

    parser.add_argument('--r-psf-path', type=str, default=None, required=False,
        help='Path to r-camera psf file (go to $DESI_SPECTRO_CALIB if not given)')

    parser.add_argument('--z-psf-path', type=str, default=None, required=False,
        help='Path to z-camera psf file (go to $DESI_SPECTRO_CALIB if not given)')

    args = parser.parse_args()

    offset = 5.
    fmin = 4
    fmax = 27

    dic = {'b':{}, 'r':{}, 'z':{}}
    for c in dic.keys():
        dic[c]['PATH'] = getattr(args,'{}_path'.format(c))
        dic[c]['LINE'] = {'ELE':getattr(args,'{}_line_name'.format(c)), 'LINE':getattr(args,'{}_line_value'.format(c))}

    f, ax = plt.subplots(nrows=6, ncols=1, figsize=(10,10))
    plt.subplots_adjust(top=0.97)

    for i,cam in enumerate(dic.keys()):

        h = fitsio.FITS(dic[cam]['PATH'])
        head = h['IMAGE'].read_header()
        camName = head['CAMERA'].strip()
        d = h['IMAGE'].read()
        w = h['MASK'].read()==0.
        td = d.copy()
        td[~w] = sp.nan
        h.close()

        if getattr(args,'{}_psf_path'.format(cam)) is None:
            cfinder = CalibFinder([head])
            p = cfinder.findfile('PSF')
        else:
            p = getattr(args,'{}_psf_path'.format(cam))
        dic[cam]['PSF'] = read_xytraceset(p)

        ### image of the PSF
        xmin = min( dic[cam]['PSF'].x_vs_wave(fmin,dic[cam]['LINE']['LINE']), dic[cam]['PSF'].x_vs_wave(fmax,dic[cam]['LINE']['LINE']) )
        xmax = max( dic[cam]['PSF'].x_vs_wave(fmin,dic[cam]['LINE']['LINE']), dic[cam]['PSF'].x_vs_wave(fmax,dic[cam]['LINE']['LINE']) )
        ymin = min( dic[cam]['PSF'].y_vs_wave(fmin,dic[cam]['LINE']['LINE']), dic[cam]['PSF'].y_vs_wave(fmax,dic[cam]['LINE']['LINE']) )
        ymax = max( dic[cam]['PSF'].y_vs_wave(fmin,dic[cam]['LINE']['LINE']), dic[cam]['PSF'].y_vs_wave(fmax,dic[cam]['LINE']['LINE']) )
        ax[2*i].imshow(td,interpolation='nearest',origin='lower',cmap='hot')
        ax[2*i].set_xlim(sp.floor(xmin-offset),sp.floor(xmax+offset))
        ax[2*i].set_ylim(sp.floor(ymin-offset),sp.floor(ymax+offset))
        #ax[2*i].set_xticklabels([])
        ax[2*i].set_ylabel(r'$\mathrm{y-axis}$')
        #ax[2*i].grid()

        x = sp.array([ dic[cam]['PSF'].x_vs_wave(f,dic[cam]['LINE']['LINE']) for f in range(fmin-1,fmax+2) ])
        y = sp.array([ dic[cam]['PSF'].y_vs_wave(f,dic[cam]['LINE']['LINE']) for f in range(fmin-1,fmax+2) ])
        fxy = sp.interpolate.interp1d(x,y)
        ax[2*i].plot(x,y,color='white',linestyle='--',linewidth=1)

        ### Slice of the PSF
        label = '\mathrm{'+camName+'}, \lambda_{\mathrm{'+dic[cam]['LINE']['ELE']+'}} = '+str(dic[cam]['LINE']['LINE'])+' \, [\mathrm{\AA{}}]'
        z = sp.array([ d[int(sp.floor(fxy(tx))),int(sp.floor(tx))] for tx in sp.arange(sp.floor(xmin-offset),sp.floor(xmax+offset)) ])
        ax[2*i+1].plot(sp.arange(sp.floor(xmin-offset),sp.floor(xmax+offset)), z, label=r'$'+label+'$')
        ax[2*i+1].plot(sp.arange(sp.floor(xmin-offset),sp.floor(xmax+offset)), sp.zeros(z.size), linewidth=1, linestyle='--', color='black')
        ax[2*i+1].set_xlim(sp.floor(xmin-offset),sp.floor(xmax+offset))
        if 2*i+1==5:
            ax[2*i+1].set_xlabel(r'$x-axis$')
        ax[2*i+1].set_ylabel(r'$\mathrm{IMAGE}$')
        #ax[2*i+1].set_xticklabels([])
        ax[2*i+1].legend(loc=1)
        #ax[2*i+1].grid()

    plt.suptitle(r'$\mathrm{SP'+str(camName[-1])+'}$',fontsize=20)
    plt.show()
