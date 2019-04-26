#!/usr/bin/env python

import os
import fitsio
import scipy as sp
import argparse
import matplotlib.pyplot as plt

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

    dic = {}
    corr = {}
    nbPairs = int(len(args.infiles)//2)
    for i in range(nbPairs):
        p1 = args.infiles[2*i]
        p2 = args.infiles[2*i+1]
        h1 = fitsio.FITS(p1)
        h2 = fitsio.FITS(p2)

        for c in cameras:

            head1 = h1[c].read_header()
            data1 = h1[c].read()
            head2 = h2[c].read_header()
            data2 = h2[c].read()
            for sec in sections:
                w1 = _parse_sec_keyword(head1['DATASEC'+sec])
                w2 = _parse_sec_keyword(head2['DATASEC'+sec])
                td1 = data1[w1].astype(sp.float64)-data2[w2].astype(sp.float64)
                print(c,sec,td1.shape)
                if sec=='B':
                    td1 = td1[:,::-1]
                elif sec=='C':
                    td1 = td1[::-1,:]
                elif sec=='D':
                    td1 = td1[::-1,::-1]

                std1 = sp.std(td1)
                mask = (td1<-5.*std1) | (td1>5.*std1)
                td1 -= td1[~mask].mean()
                td1[mask] = 0.

                dic[c+'_'+sec] = td1.copy()

        h1.close()
        h2.close()

        ### Compute correlation coef
        for i1,sec1 in enumerate(dic.keys()):
            for sec2 in list(dic.keys())[:i1+1]:
                print(i,sec1,sec2)
                try:
                    tc = compute_correlation(dic[sec1],dic[sec2],nbNeigh)
                except ValueError:
                    print('ERROR: {} {}, shapes {} and {}'.format(sec1,sec2,dic[sec1].shape,dic[sec2].shape))
                    continue
                if i==0:
                    corr['{}_{}'.format(sec1,sec2)] = tc
                else:
                    corr['{}_{}'.format(sec1,sec2)] += tc

    ### Save
    for i1,sec1 in enumerate(dic.keys()):
        for sec2 in list(dic.keys())[:i1+1]:
            print(sec1,sec2)
            try:
                tc = corr['{}_{}'.format(sec1,sec2)]/nbPairs
            except KeyError:
                print('ERROR: {} {}'.format(sec1,sec2))
                continue
            sp.savetxt(args.outfile+'/correlation_{}_{}.txt'.format(sec1,sec2),tc)

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

    lst_sec = [ c+'_'+sec for sec in sections for c in cameras]

    for i1,sec1 in enumerate(lst_sec):
        for sec2 in list(lst_sec)[:i1+1]:
            print(sec1,sec2)
            p = args.outfile+'/correlation_{}_{}.txt'.format(sec1,sec2)
            if not os.path.isfile(p):
                p = args.outfile+'/correlation_{}_{}.txt'.format(sec2,sec1)
            if not os.path.isfile(p):
                print('ERROR: {} {}'.format(sec1,sec2))
                continue
            tc = sp.loadtxt(p)
            tc[tc==0.] = sp.nan
            #vminmax = max( 0.2, max(abs(tc.min()), abs(tc[tc<0.9999].max())) )
            vminmax = max(abs(tc.min()), abs(tc[tc<0.9999].max()))
            plt.imshow(tc, vmin=-vminmax, vmax=vminmax,
                origin='lower', extent=(-nbNeigh,nbNeigh,-nbNeigh,nbNeigh), cmap='seismic')
            plt.title(r'$\mathrm{'+sec1.replace('_','-')+'\, x \, '+sec2.replace('_','-')+'}$',fontsize=20)
            plt.xlabel(r'$\Delta x \, [\mathrm{pix}]$',fontsize=20)
            plt.ylabel(r'$\Delta y \, [\mathrm{pix}]$',fontsize=20)
            cbar = plt.colorbar()
            cbar.set_label(r'$Corr(\Delta x, \Delta y)$',size=20)
            cbar.update_ticks()
            #plt.grid()
            #plt.show()
            plt.savefig(args.outfile+'/correlation_{}_{}.png'.format(sec1,sec2))
            plt.clf()
    return

main()
