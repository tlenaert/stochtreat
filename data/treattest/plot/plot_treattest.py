#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import plothelpers
import helperfunctions
import itertools 
import sys

plothelpers.latexify()
fig,ax=plt.subplots(1,1)
fig.subplots_adjust(left=0.2,bottom=0.2,right=0.95,top=0.96)
ax.margins(0.05)

filenames=sys.argv[1:]

cmap = plt.get_cmap('viridis')
indices = np.linspace(0, cmap.N, len(filenames))
my_colors = itertools.cycle([cmap(int(i)) for i in indices])
labels=itertools.cycle(['All patients'])

offset=0.
width=1.
lengthening=1. #25.
for filename in filenames:
    data=np.loadtxt(filename)
    color=next(my_colors)
    label=next(labels)
    ax.plot(data[:,0]*lengthening,data[:,1],alpha=0.7,color=color,label=label )#,marker="o",markersize=2.0,linestyle="None")
    # ax.bar((data[:,0]+offset)*lengthening,data[:,1],width*lengthening,alpha=0.7,color=color,edgecolor="none",linewidth=0.,label=label )#,marker="o",markersize=2.0,linestyle="None")
    offset+=width+0.001

# ax.set_xlim(xmin=0.,xmax=10)

ax.set_xlabel("treatment time in years")
ax.set_ylabel("recurrence probability")
ax.legend()
plt.show()
