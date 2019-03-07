#!/usr/bin/env python

import os
import scipy as sp
import subprocess
import fitsio
import matplotlib.pyplot as plt
import copy
import argparse
from astropy.time import Time

from desispec.io.xytraceset import read_xytraceset


def preproc(indir,outdir,cam):


    cmd = 'desi_preproc'
    cmd += ' -i {}'.format(indir)
    cmd += ' --cameras '
    for i,tc in enumerate(cam):
        cmd += tc
        if i!=len(cam)-1:
            cmd += ','
    cmd += ' --outdir {}'.format(outdir)
    print(cmd)
    subprocess.call(cmd,shell=True)

    return
def psf_fit(idx,outdir,cam):

    if os.path.isfile('{}/psf-{}-{}.fits'.format(outdir,cam,idx)):
        print('INFO: already done')
        return

    cmd = 'desi_psf_fit'
    cmd += ' -a {}/preproc-{}-{}.fits'.format(outdir,cam,idx)
    cmd += ' --in-psf /software/datasystems/desi_spectro_calib/trunk/spec/sp2/psf-{}-20190301.fits'.format(cam)
    cmd += ' --out-psf {}/psf-{}-{}.fits'.format(outdir,cam,idx)
    if 'b' in cam:
        cmd += ' --trace-deg-wave 4 --trace-deg-x 4 --legendre-deg-x 4 --legendre-deg-wave 2 --single-bundle'
    elif ('r' in cam) or ('z' in cam):
        cmd += ' --legendre-deg-x 4 --legendre-deg-wave 3 --single-bundle'
    print(cmd)
    subprocess.call(cmd,shell=True)

    return

def plot_qproc_traces(idx,outdir,cam):

    cmd = 'desi_qproc'
    cmd += ' -i {}/preproc-{}-{}.fits'.format(outdir,cam,idx)
    cmd += ' -p {}/psf-{}-{}.fits --plot'.format(outdir,cam,idx)
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'plot_traces.py'.format(outdir,cam,idx)
    cmd += ' -p {}/psf-{}-{}.fits'.format(outdir,cam,idx)
    print(cmd)
    subprocess.call(cmd,shell=True)

    return
def plot_bundle(outdir,cam,idx):

    cmd = 'plot_psf_bundle.py'
    for c in cam:
        cmd += ' --{}-path {}/preproc-{}-{}.fits'.format(c[0],outdir,c,idx)
        cmd += ' --{}-psf-path {}/psf-{}-{}.fits'.format(c[0],outdir,c,idx)

    print(cmd)
    subprocess.call(cmd,shell=True)

    return
def plot_resolution(outdir,cam,start,end,exptime):

    ###
    for i in exptime.keys():
        idx = str(i).zfill(8)
        cmd = 'desi_compute_resolution'
        cmd += ' -i {}/psf-{}-{}.fits'.format(outdir,cam,idx)
        cmd += ' -o {}/resolution-{}-{}.txt'.format(outdir,cam,idx)
        print(cmd)
        subprocess.call(cmd,shell=True)

    ###
    colormap = plt.cm.Set1
    colors = [colormap(i) for i in sp.linspace(0, 1,len(exptime.keys()))]
    print(colors)
    cmd = 'desi_plot_resolution'
    cmd += ' -i '
    for i in exptime.keys():
        idx = str(i).zfill(8)
        cmd += ' {}/resolution-{}-{}.txt'.format(outdir,cam,idx)
    cmd += ' --mean --only-mean --legend'
    cmd += ' --labels'
    for i,t in exptime.items():
        cmd += ' {}'.format(t)
    subprocess.call(cmd,shell=True)

    return
def plot_psf_vs_time(outdir,cam,exps,fiber,dicLines):

    dic = copy.deepcopy(exps)
    for k,v in dic.items():
        dic[k][cam] = {}
        idx = str(k).zfill(8)
        p = '{}/psf-{}-{}.fits'.format(outdir,cam,idx)
        try:
            psf = read_xytraceset(p)
        except:
            print('\nERROR',p,'\n')
            continue
        dic[k][cam]['SIGMAX'] = psf.xsig_vs_wave(fiber,dicLines[cam[0]]['LINE'])
        dic[k][cam]['SIGMAY'] = psf.ysig_vs_wave(fiber,dicLines[cam[0]]['LINE'])
        dic[k][cam]['X'] = psf.x_vs_wave(fiber,dicLines[cam[0]]['LINE'])
        dic[k][cam]['Y'] = psf.y_vs_wave(fiber,dicLines[cam[0]]['LINE'])

    ###
    for exptime in sorted(set([ dic[tk]['EXPTIME'] for tk in dic.keys() ])):
        f, ax = plt.subplots(nrows=4, ncols=1, figsize=(10,10))
        plt.subplots_adjust(top=0.95,hspace=0.,wspace=0.)
        plt.suptitle(r'$\mathrm{SP'+str(cam[-1])+', cam = '+cam+',\, exptime='+str(exptime)+'}$',fontsize=20)
        for i,k in enumerate(['SIGMAX','SIGMAY','X','Y']):
            x = sp.array([ dic[tk]['DATE_OBS'] for tk in dic.keys() if dic[tk]['EXPTIME']==exptime and k in dic[tk][cam].keys() ])
            y = sp.array([ dic[tk][cam][k] for tk in dic.keys() if dic[tk]['EXPTIME']==exptime and k in dic[tk][cam].keys() ])
            w = sp.argsort(x)
            x = x[w]
            y = y[w]
            x *= 24.*3600.
            x -= x[0]
            y = 100.*(y-y[0])/y[0]
            ax[i].plot(x,y,'o-')
            if i==len(['SIGMAX','SIGMAY','X','Y'])-1:
                ax[i].set_xlabel(r'$\mathrm{\Delta DATE-OBS} \, [\mathrm{s}]$')
            else:
                ax[i].set_xticklabels([])
            ax[i].grid()
            ax[i].set_ylabel(r'$\mathrm{'+k+'} \, [\%]$')

        #plt.show()
        plt.savefig(outdir+'/psf-vs-time-cam-{}-exptime-{}.png'.format(cam,exptime))
        plt.clf()

    return

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Plot the PSF image and slice along a constant wavelength')

    parser.add_argument('--first-exposure', type=int, default=None, required=True,
        help='First exposure')

    parser.add_argument('--last-exposure', type=int, default=None, required=True,
        help='Last exposure')

    parser.add_argument('-i','--in-dir', type=str, default=None, required=True,
        help = 'Input directory like /exposures/desi/<date>/')

    parser.add_argument('-o','--out-dir', type=str, default=None, required=True,
        help = 'Ouput directory')

    parser.add_argument('--cameras', type=str, default=['b0','r0','z0'], required=False, nargs="*",
        help='List of cameras to look at')

    parser.add_argument('--fiber', type=int, default=4, required=False,
        help='Fiber to use')

    parser.add_argument('--run-preproc', action='store_true',
        help="Run pre-proc")

    parser.add_argument('--run-psffit', action='store_true',
        help="Run psf-fit")

    parser.add_argument('--plot', action='store_true',
        help="Plot PSF vs. time")

    parser.add_argument('--more-plot', action='store_true',
        help="A lot of plot of control")
    
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
    
    args = parser.parse_args()
    
    ###
    dicLines = {'b':{}, 'r':{}, 'z':{}}
    for c in dicLines.keys():
        dicLines[c] = {'ELE':getattr(args,'{}_line_name'.format(c)), 'LINE':getattr(args,'{}_line_value'.format(c))}

    ###
    cam = args.cameras
    indir = args.in_dir
    outdir = args.out_dir
    start = args.first_exposure
    end = args.last_exposure

    ### Exposure date
    exps = { i:{} for i in range(start,end+1) }
    for k in exps.keys():
        idx = str(k).zfill(8)
        tindir = indir+'/{0}/desi-{0}.fits.fz'.format(idx)
        h = fitsio.FITS(tindir)
        exps[k]['DATE_OBS'] = Time(h['SPS'].read_header()['DATE-OBS']).mjd
        exps[k]['EXPTIME'] = h['SPS'].read_header()['EXPTIME']
        h.close()

    ###
    if args.run_preproc:
        for i in exps.keys():
            idx = str(i).zfill(8)
            tindir = indir+'/{0}/desi-{0}.fits.fz'.format(idx)
            preproc(tindir,outdir,cam)

    ###
    if args.run_psffit:
        for i in exps.keys():
            for tcam in cam:
                idx = str(i).zfill(8)
                psf_fit(idx,outdir,tcam)

    ###
    if args.plot:
        for tcam in cam:
            plot_psf_vs_time(outdir,tcam,exps,args.fiber,dicLines)

    ###
    if args.more_plot:
        ### Plot qproc and traces
        for tcam in cam:
            for i in exps.keys():
                idx = str(i).zfill(8)
                plot_qproc_traces(idx,outdir,tcam)

        ### Plot bundle
        for i in exps.keys():
            idx = str(i).zfill(8)
            plot_bundle(outdir,cam,idx)

        ### Plot differences
        for tcam in cam:
            plot_resolution(outdir,tcam,start,end,exps)
            print('\n\n\n')

    return

main()


