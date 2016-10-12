#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats 
import pandas as pd
from matplotlib.ticker import MaxNLocator 

import plothelpers
import helperfunctions
import sys

filename=sys.argv[1]
data = pd.read_csv(filename, delim_whitespace=True, header = 0,index_col=False)

cols=data.columns.tolist()
# print(cols)
# cols= ["$t_{nolsc}$", "$t_{diag}$","$t_{red}$","init resp","lsc at diag", "1y burden"]
# data.columns = cols

cols.insert(0, cols.pop(cols.index('init_resp')))
print(cols)
# cols[3], cols[0] = cols[0], cols[3]
data=data.ix[:,cols]

correl=data.corr().dropna(axis=(0,1),how='all') 
print(correl.columns)

plothelpers.latexify(fig_height=2.5)
fig,ax=plt.subplots(1,1)
fig.subplots_adjust(left=0.20,bottom=0.03,right=0.88,top=0.72)
ax.margins(0.05)

matplot=ax.matshow(correl,cmap="viridis")

ax.set_xticks(np.arange(len(correl.columns)))
ax.set_yticks(np.arange(len(correl.columns)))
ax.set_yticklabels(list(correl.columns))
ax.set_xticklabels(list(correl.columns), rotation=45)
for tick in ax.xaxis.get_majorticklabels():
    tick.set_horizontalalignment("left")

cb=fig.colorbar(matplot)
cb.set_label("Pearson corr. coeff.")


plt.show()
