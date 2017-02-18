#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--input', type = str, default = None, required = True,
                    help = 'path to psf stability ASCII file')
parser.add_argument('--title', type = str, default = "stability", required = False,
                    help = 'title')
args=parser.parse_args()

x=np.loadtxt(args.input).T
res_fiber=x[0]
res_wave=x[1]
res_emission_line_rms=x[2]
res_continuum_rms=x[3]
fibers=np.unique(res_fiber)

ylim=[0,np.max(res_emission_line_rms)*1.1]
plt.figure()
plt.subplot(2,1,1)
plt.title(args.title)
for fiber in fibers :
    plt.plot(res_wave[res_fiber==fiber],res_emission_line_rms[res_fiber==fiber],"-")
plt.xlabel("Wavelength")
plt.ylabel("flux error for emission lines")
plt.ylim(ylim)
plt.grid()

plt.subplot(2,1,2)
for fiber in fibers :
    plt.plot(res_wave[res_fiber==fiber],res_continuum_rms[res_fiber==fiber],"-")
plt.xlabel("Wavelength")
plt.ylabel("flux error for continuum")
plt.ylim(ylim)
plt.grid()
plt.show()



