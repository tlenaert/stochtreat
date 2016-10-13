#!/usr/bin/env python3

import numpy as np
import sys

filenames = sys.argv[1:]
t_max=15.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#")[-1] for filename in filenames]

data=np.concatenate(rawdata)
data[data<-1.]=np.inf
data[data<0.]=np.NaN

times=np.linspace(0.,t_max,t_steps)

lsc_frac=[]
no_patients=(data>=0.).nansum()

for t in times:
    no_lscs=(data<t).nansum()
    lsc_frac.append(no_lscs/no_patients)

for t,l in zip(times,lsc_frac):
    print(t,l)


