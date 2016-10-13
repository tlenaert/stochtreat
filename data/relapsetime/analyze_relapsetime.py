#!/usr/bin/env python3

import numpy as np
import sys

filenames = sys.argv[1:]
t_max=10.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#")[:,-1] for filename in filenames]

data=np.concatenate(rawdata)
data[data<=-2.]=np.inf
data[data<0.]=np.NaN

# print(data)
# exit(0)
times=np.linspace(0.,t_max,t_steps)

lsc_frac=[]
no_patients=np.nansum((data>=0.))

for t in times:
    no_lscs=np.nansum((data<t))
    lsc_frac.append(1-no_lscs/float(no_patients))

for t,l in zip(times,lsc_frac):
    print(t,l)


