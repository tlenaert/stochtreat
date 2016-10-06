#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator 

import plothelpers
import helperfunctions
import sys

plothelpers.latexify()
fig,ax=plt.subplots(1,1)
fig.subplots_adjust(left=0.19,bottom=0.20,right=0.96,top=0.96)
ax.margins(0.05)

if len(sys.argv[:])==2:
    treattimenumber=sys.argv[1]
    filenamer="results/relapsedat"+str(treattimenumber)+".txt"
    filenamen="results/norelapsedat"+str(treattimenumber)+".txt"
else:
    filenamer=sys.argv[1]
    filenamen=sys.argv[2]
relapse=np.loadtxt(filenamer)
norelapse=np.loadtxt(filenamen)

ax.hist(relapse[:,0],bins=50,normed=True,alpha=0.8,label="relapse")
ax.hist(norelapse[:,0],bins=50,normed=True,alpha=0.8,label="no relapse")

ax.xaxis.set_major_locator(MaxNLocator(6))
ax.set_xlabel("initial response slope")
ax.set_ylabel("probability density")

ax.legend()
plt.show()
