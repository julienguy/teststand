#!/usr/bin/env python


import sys,string
import argparse
import numpy as np
from desispec.io.fibermap import empty_fibermap,write_fibermap

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-o','--out', type = str, default = None, required = True,
                    help = 'path to output fibermap')
parser.add_argument('--nspec', type = int, default = 0, required = True,
                    help = 'number of spectra')
parser.add_argument('--flavor', type = str, default = "bla", required = False,
                    help = 'flavor')
parser.add_argument('--night', type = str, default = "19751016", required = False,
                    help = 'night')
parser.add_argument('--expid', type = str, default = 0, required = False,
                    help = 'exposure id')
parser.add_argument('--spectro', type = int, default = 1, required = False,
                    help = 'spectro')


args        = parser.parse_args()

fmap = empty_fibermap(args.nspec*(args.spectro+1))
fmap.meta["FLAVOR"]=args.flavor
fmap.meta["NIGHT"]=args.night
fmap.meta["EXPID"]=args.expid
write_fibermap(args.out,fmap)
