#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
from teststand.graph_tools         import plot_graph
from teststand.graph_tools         import show_graph
from desispec.log import get_logger

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--frame', type = str, default = None, required = True,
                    help = 'path of frame fits file')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('--err', action='store_true', default = None, required = False,
                    help = 'Show errorbar in plot.')
parser.add_argument('-d', '--dim', action='store_true', default = None, required = False,
                    help = 'Show a 2d graph via imshow plus the actual plot.')
parser.add_argument('--log', action='store_true', default = False, required = False,
                    help = 'log scale')

log         = get_logger()
args        = parser.parse_args()
frame_file  = pyfits.open(args.frame)


if args.fibers is not None :
    nb = args.fibers.split(',')
    for i in xrange(len(nb)) :
        if nb[i].isdigit() == False :
            tmp = nb[i].split(':')
            if ((len(tmp) is 2) and tmp[0].isdigit() == True and tmp[1].isdigit() == True) :
                f_begin = int(tmp[0])
                f_end   = int(tmp[1])
                plot_graph(frame_file, nfibers=None, start=f_begin, end=f_end, opt_err=args.err, opt_2d=args.dim,logscale=args.log)
            else :
                log.error("--fibers parsing error.\nCorrect format is either  : --fibers=begin,end (excluded)\nand/or  : --fibers=begin:end (excluded)\nYou can use : --fibers=2,5,6:8,3,10")
                sys.exit(1)
        else :
            nb[i] = int(nb[i])
            plot_graph(frame_file, nfibers=None, start=nb[i], end=None, only=True, opt_err=args.err, opt_2d=args.dim,logscale=args.log)
else :
    #   If you did not ask for a specific fiber/group of fiber, they're all plotted.
    plot_graph(frame_file, nfibers=frame_file[0].data.shape[0],logscale=args.log)

show_graph()
log.info("Script's done")
