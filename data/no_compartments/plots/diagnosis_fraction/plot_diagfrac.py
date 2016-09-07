#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import plothelpers
import helperfunctions
import sys

plothelpers.latexify()
fig,ax=plt.subplots(1,1)
fig.subplots_adjust(left=0.2,bottom=0.2,right=0.95,top=0.96)
ax.margins(0.05)

filename=sys.argv[1]
data=np.loadtxt(filename)

ax.plot(data[:,0],1-data[:,1],marker="o",markersize=5.0,linestyle="None")


ax.set_xlabel("\# of stochastic compartments")
ax.set_ylabel("prob. of no diagnosis")
plt.show()
