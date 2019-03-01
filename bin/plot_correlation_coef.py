#!/usr/bin/env python

import fitsio
import scipy as sp
import argparse
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)

from desispec.preproc import _parse_sec_keyword

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i','--infiles', type=str, nargs="*", default=None, required=True,
        help = 'Path to pairs of exposures')

    parser.add_argument('-o','--outfile', type=str, default=None, required=True,
        help = 'Path where to save ascii files with results')

    parser.add_argument('--cameras', type=str, nargs="*", default=['b0','r0','z0'], required=False,
        help = 'Name of the cameras')

    parser.add_argument('--sections', type=str, nargs="*", default=['A','B','C','D'], required=False,
        help = 'Amplifiers sections of the CCD')

    parser.add_argument('--nb-neighbours', type=int, default=7, required=False,
        help = 'Number of Neigbour to look at')

    parser.add_argument('--plot', action='store_true', required=False,
        help='Plot the results, do not compute them')

    args = parser.parse_args()

    if args.plot:
        plot(args)
    else:
        compute(args)

def compute(args):

    nbNeigh = args.nb_neighbours
    cameras = args.cameras
    sections = args.sections

    p1 = args.infiles[0]
    p2 = args.infiles[1]
    h1 = fitsio.FITS(p1)
    h2 = fitsio.FITS(p2)

    for c in cameras:

        ### Read data
        head1 = h1[c].read_header()
        data1 = h1[c].read()
        head2 = h2[c].read_header()
        data2 = h2[c].read()
        dic = {}
        for sec in sections:
            w1 = _parse_sec_keyword(head1['DATASEC'+sec])
            w2 = _parse_sec_keyword(head2['DATASEC'+sec])
            td1 = data1[w1].astype('float64')-data2[w2].astype('float64')
            std1 = sp.std(td1)

            if sec=='B':
                td1 = td1[:,::-1]
            elif sec=='C':
                td1 = td1[::-1,:]
            elif sec=='D':
                td1 = td1[::-1,::-1]

            mask = (td1<-5.*std1) | (td1>5.*std1)
            td1 -= td1[~mask].mean()
            td1[mask] = 0.

            dic[sec] = td1.copy()

        ### Compute correlation coef
        for i1,sec1 in enumerate(dic.keys()):
            for sec2 in list(dic.keys())[:i1+1]:
                print(c, sec1,sec2)
                corr = compute_correlation(dic[sec1],dic[sec2],nbNeigh)
                sp.savetxt(args.outfile+'/correlation_{}_{}_{}.txt'.format(c,sec1,sec2),corr)

    h1.close()
    h2.close()

    return
def compute_correlation(im1,im2,nbNeigh):
    cor = sp.zeros( (2*nbNeigh+1,2*nbNeigh+1) )

    std1 = sp.std(im1[im1!=0.])
    std2 = sp.std(im2[im2!=0.])
    norm = std1*std2

    for i in range(-nbNeigh,nbNeigh+1):
        for j in range(-nbNeigh,nbNeigh+1):
            mult = im1*sp.roll(sp.roll(im2,i,axis=1),j,axis=0)/norm
            cor[nbNeigh+j,nbNeigh+i] = (mult[mult!=0.]).mean()

    return cor
def plot(args):

    nbNeigh = args.nb_neighbours
    cameras = args.cameras
    sections = args.sections

    for c in cameras:
        for i1,sec1 in enumerate(sections):
            for sec2 in sections[:i1+1]:
                print(c, sec1,sec2)
                tc = sp.loadtxt(args.outfile+'/correlation_{}_{}_{}.txt'.format(c,sec1,sec2))
                tc[tc==0.] = sp.nan
                plt.imshow(tc, vmin=max(-0.2,tc.min()), vmax=min(0.2,tc.max()), origin='lower', extent=(-nbNeigh,nbNeigh,-nbNeigh,nbNeigh))
                plt.title(r'$\mathrm{'+c+': \, Amp \, '+sec1+'\, x \, Amp \, \, '+sec2+'}$',fontsize=20)
                plt.xlabel(r'$\Delta x \, [\mathrm{pix}]$',fontsize=20)
                plt.ylabel(r'$\Delta y \, [\mathrm{pix}]$',fontsize=20)
                cbar = plt.colorbar()
                cbar.set_label(r'$Corr(\Delta x, \Delta y)$',size=20)
                cbar.update_ticks()
                plt.grid()
                plt.show()
                plt.savefig(args.outfile+'/correlation_{}_{}_{}.png'.format(c,sec1,sec2))
                plt.clf()
    return

main()
