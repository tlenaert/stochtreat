#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import helperfunctions
import plothelpers
import sys

filenames = helperfunctions.natural_sort(sys.argv[1:])

data=[np.loadtxt(filename)[:,0] for filename in filenames]

means=[np.mean(x[x!=25.]) for x in data]
stddev=[np.std(x[x!=25.]) for x in data]

# print (data)

# plothelpers.latexify()
fig,ax=plt.subplots(2,1,sharex=True)
fig.subplots_adjust(left=0.1,bottom=0.1,right=0.95,top=0.95)
ax[0].margins(0.05)
ax[1].margins(0.05)

ax[0].violinplot(data, range(1,len(data)+1), points=100, widths=0.7,
                      showmeans=True, showextrema=False, showmedians=False)

# ax[0].set_xlabel("# of stochastic compartments")

ax[1].errorbar(range(1,len(means)+1),means,yerr=stddev,marker="o",linestyle="None")
ax[1].set_ylim(ymin=0)
# ax[1].plot(range(1,len(stddev)+1),stddev,marker="o",linestyle="None")


ax[1].set_xlabel("# of stochastic compartments")
ax[0].set_ylabel("diagnosis time dist.")
ax[1].set_ylabel("diag. time mean/stddev")
plt.show()
