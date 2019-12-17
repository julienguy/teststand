#!/usr/bin/env python


import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table


filename=sys.argv[1]
table=Table.read(filename,format="csv")
plt.plot(table["XPIX"],table["YPIX"],".")

if "PIN_ID" in table.dtype.names :
    ii=table["PIN_ID"]>0
    if np.sum(ii)>0 :
        plt.plot(table["XPIX"][ii],table["YPIX"][ii],"o")
    ii=table["PIN_ID"]==1027
    if np.sum(ii)>0 :
        plt.plot(table["XPIX"][ii],table["YPIX"][ii],"o")
   

plt.show()

