#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import plothelpers
import sys

data=np.loadtxt(sys.argv[1])

plothelpers.latexify()
fig,ax=plt.subplots()
fig.subplots_adjust(left=0.2,bottom=0.2,right=0.95,top=0.95)
ax.plot(data[:,0],data[:,1])
ax.set_xlabel("time in years")
ax.set_ylabel("freq. of LSC ext.")

plt.show()
