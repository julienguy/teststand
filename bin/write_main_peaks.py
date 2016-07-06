#!/usr/bin/env python

import numpy as np
import desispec.bootcalib
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-o','--outfile', type = str, default = None, required=True, help="output ASCII file")
parser.add_argument('--air', action = "store_true", help="wavelength in air (default is vacuum)")
parser.add_argument('--minrelint', type = float, default = 10., required=False, help="relative intensity threshold")

args=parser.parse_args()

vacuum=True
if args.air :
    vacuum=False

llist=desispec.bootcalib.load_arcline_list(camera="all", vacuum=vacuum, lamps=None)
ion=np.array(llist["Ion"])
wave=np.array(llist["wave"])
relint=np.array(llist["RelInt"])
maxrelint=np.max(relint)
ok=np.where(relint>args.minrelint)[0]
wave=wave[ok]
ion=ion[ok]
sort=np.argsort(wave)
ofile=open(args.outfile,"w")
if vacuum :
    ofile.write("# WAVE (A, IN VACUUM) RELINT>%f\n"%args.minrelint)
else :
    ofile.write("# WAVE (A, IN AIR) RELINT>%f\n"%args.minrelint)
for line in sort :
    ofile.write("%f %s\n"%(wave[line],ion[line]))
ofile.close()


