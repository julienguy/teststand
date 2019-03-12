#!/usr/bin/env python

import fitsio
import scipy as sp
import argparse
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)

from desispec.preproc import _parse_sec_keyword

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i','--infiles', type=str, nargs="*", default=None, required=False,
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
    assert len(args.infiles)%2. == 0.

    ### Where to save
    dicCor = {}
    for c in cameras:
        dicCor[c] = {}
        for i1,sec1 in enumerate(sections):
            for sec2 in sections[:i1+1]:
                dicCor[c][sec1+'_'+sec2] = sp.zeros( (2*nbNeigh+1,2*nbNeigh+1) )

    ### Compute
    for p in range(int(len(args.infiles)//2.)):
        p1 = args.infiles[2*p]
        p2 = args.infiles[2*p+1]
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
                    print(p, c, sec1,sec2)
                    dicCor[c][sec1+'_'+sec2] += compute_correlation(dic[sec1],dic[sec2],nbNeigh)

        h1.close()
        h2.close()

    ### Write to ascii
    for c in cameras:
        for i1,sec1 in enumerate(sections):
            for sec2 in sections[:i1+1]:
                dicCor[c][sec1+'_'+sec2] /= len(args.infiles)/2.
                sp.savetxt(args.outfile+'/correlation_{}_{}_{}.txt'.format(c,sec1,sec2),dicCor[c][sec1+'_'+sec2])

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

        f, ax = plt.subplots(nrows=len(sections), ncols=len(sections))
        plt.subplots_adjust(top=0.97, bottom=0.07, wspace=0., hspace=0.)
        plt.suptitle(r'$\mathrm{SP2, '+c+'}$',fontsize=20)

        for i1,sec1 in enumerate(sections):
            for i2,sec2 in enumerate(sections):

                if i2>i1:
                    f.delaxes(ax[i1,i2])
                    continue

                try:
                    tc = sp.loadtxt(args.outfile+'/correlation_{}_{}_{}.txt'.format(c,sec1,sec2))
                except:
                    tc = sp.loadtxt(args.outfile+'/correlation_{}_{}_{}.txt'.format(c,sec2,sec1))
                tc[tc==0.] = sp.nan
                ax[i1,i2].imshow(tc, vmin=-0.2, vmax=0.2, origin='lower', extent=(-nbNeigh,nbNeigh,-nbNeigh,nbNeigh))
                if i1==len(sections)-1:
                    ax[i1,i2].set_xlabel(sec2,fontsize=20)
                else:
                    ax[i1,i2].set_xticks([])
                if i2==0:
                    ax[i1,i2].set_ylabel(sec1,fontsize=20)
                else:
                    ax[i1,i2].set_yticks([])

        plt.show()

    return

main()
