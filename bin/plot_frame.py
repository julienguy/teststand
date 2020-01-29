#!/usr/bin/env python


import sys
import argparse
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from teststand.graph_tools         import plot_graph,parse_fibers
from desiutil.log                  import get_logger
import os.path
from desispec.io import read_fibermap
               
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*",
                    help = 'path to one or several frame fits files')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('--err', action='store_true', default = None, required = False,
                    help = 'Show errorbar in plot.')
parser.add_argument('--image', action='store_true', default = None, required = False,
                    help = 'Show a 2d graph via imshow plus the actual plot.')
parser.add_argument('--log', action='store_true', default = False, required = False,
                    help = 'log scale')
parser.add_argument('-b', '--batch', action='store_true', default = False, required = False,
                    help = 'batch mode (to save figure and exit')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'save figure in this file')
parser.add_argument('-l','--legend', action='store_true',help="show legend")
parser.add_argument('--xlim',type = str, default=None, help="min,max xlim for plot")
parser.add_argument('--ylim',type = str, default=None, help="min,max ylim for plot")
parser.add_argument('--labels',type = str, default=None, required = False, nargs="*")
parser.add_argument('--objtype',type = str, default=None, required = False, help="display fibers with this OBJTYPE")
parser.add_argument('--focal-plane',action='store_true', help="focal plane image instead of spectra")
parser.add_argument('--vmin',type=float, default=None)
parser.add_argument('--vmax',type=float, default=None)

log         = get_logger()
args        = parser.parse_args()

fig=plt.figure()
subplot=plt.subplot(1,1,1)

if args.fibers is not None :
    fibers = parse_fibers(args.fibers)
else :
    fibers = None



if args.labels == None or len(args.labels)<len(args.infile) :
    args.labels = []
    for filename in args.infile :
        args.labels.append(os.path.basename(filename))

xx=[]
yy=[]
ff=[]

first=True
for filename,label in zip(args.infile,args.labels) :
    frame_file  = pyfits.open(filename)
    
    if args.objtype is not None :
        fmap = read_fibermap(filename)
        fibers = np.where(fmap["OBJTYPE"]==args.objtype)[0]
        print("fibers with OBJTYPE={} : {}".format(args.objtype,fibers))


    if args.focal_plane :
        
        mflux = np.median(frame_file[0].data,axis=1)
        n1=frame_file[0].data.shape[1]
        #mflux = np.mean(frame_file[0].data[:,n1//2-500:n1//2+500],axis=1)
        fmap=frame_file["FIBERMAP"].data
        x=fmap["FIBERASSIGN_X"]
        y=fmap["FIBERASSIGN_Y"]
        xx.append(x)
        yy.append(y)
        ff.append(mflux)
        continue

    if fibers is None :
        fibers = np.arange(frame_file[0].data.shape[0])
    plot_graph(frame=frame_file,fibers=fibers,opt_err=args.err,opt_2d=args.image,label=label)  

if args.focal_plane :
    xx=np.hstack(xx)
    yy=np.hstack(yy)
    ff=np.hstack(ff)
    if args.vmin is None :
        mf=np.median(ff)
        rms=1.4*np.median(np.abs(ff-mf))
        args.vmin = mf-3*rms
        args.vmax = mf+3*rms
        
    plt.scatter(xx,yy,c=ff,vmin=args.vmin,vmax=args.vmax)
    plt.axis('off')
    plt.colorbar()

if args.log :
    subplot.set_yscale("log")
if args.legend :
    subplot.legend(loc="upper left",fontsize="small")
subplot.grid()

if args.xlim is not None :
    try :
        vv=args.xlim.split(",")
        vmin=float(vv[0])
        vmax=float(vv[1])
        plt.xlim([vmin,vmax])
    except :
        print("failed to interpret xlim")
        print(sys.exc_info())
if args.ylim is not None :
    try :
        vv=args.ylim.split(",")
        vmin=float(vv[0])
        vmax=float(vv[1])
        plt.ylim([vmin,vmax])
    except :
        print("failed to interpret ylim")
        print(sys.exc_info())
    

if args.output :
   fig.savefig(args.output)
   log.info("figure saved in %s"%args.output)

if not args.batch :
    plt.show()


