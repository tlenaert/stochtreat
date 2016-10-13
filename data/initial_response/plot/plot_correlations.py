#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats 
import pandas as pd
from matplotlib.ticker import MaxNLocator 

import plothelpers
import helperfunctions
import sys

# filename=sys.argv[1]
filenames=sys.argv[1:]
data = pd.concat([pd.read_csv(filename, delim_whitespace=True, header = 0,index_col=False) for filename in filenames],axis=0)

cols=data.columns.tolist()
cols=helperfunctions.natural_sort(cols)
cols.insert(0, cols.pop(cols.index('init.response')))
print(cols)

data=data.ix[:,cols]
cols=[col.replace("_","\_") for col in cols]
data.columns=cols

correl=data.corr().dropna(axis=(0,1),how='all')# method="pearson"
# print(correl.columns)

plothelpers.latexify(columns=2,fig_height=4.5)
fig,ax=plt.subplots(1,1)
fig.subplots_adjust(left=0.20,bottom=0.03,right=0.88,top=0.72)
ax.margins(0.05)

matplot=ax.matshow(correl,cmap="inferno")

ax.set_xticks(np.arange(len(correl.columns)))
ax.set_yticks(np.arange(len(correl.columns)))
ax.set_yticklabels(list(correl.columns))
ax.set_xticklabels(list(correl.columns), rotation=45)
for tick in ax.xaxis.get_majorticklabels():
    tick.set_horizontalalignment("left")

cb=fig.colorbar(matplot)
cb.set_label("Pearson corr. coeff.")


plt.show()
