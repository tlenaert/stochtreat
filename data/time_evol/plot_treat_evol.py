#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import plothelpers
import scipy.stats
import helperfunctions
import sys

plothelpers.latexify()
fig,ax=plt.subplots(1,1)
fig.subplots_adjust(left=0.2,bottom=0.2,right=0.95,top=0.96)
ax.margins(0.05)

filename=sys.argv[1]
data=np.loadtxt(filename,comments="#")

slope, intercept,bla,blub,bla=scipy.stats.linregress(data[:62,0],data[:62,1])
yfit=intercept+slope*data[:62,0]
# ax.plot(data[:,0],data[:,1],marker="o",markeredgewidth=0.1,
        # markeredgecolor="blue",markersize=1.0,linestyle="None")
ax.plot(data[:,0],data[:,1])
ax.plot(data[:62,0],yfit)

# ax.set_yscale('log')
ax.set_xlim(xmax=1800, xmin=1600)


ax.set_xlabel("time in days")
ax.set_ylabel("burden")
plt.show()
