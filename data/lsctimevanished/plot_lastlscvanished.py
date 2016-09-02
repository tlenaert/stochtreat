#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys

filenames = sys.argv[1:]
t_max=20.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.concatenate(rawdata)
data[data<0.]=np.inf
print(data)

times=np.linspace(0.,t_max,t_steps)

lsc_frac=[]
no_patients=len(data)

for t in times:
    no_lscs=(data<t).sum()
    lsc_frac.append(no_lscs/no_patients)

plt.plot(times,lsc_frac)

plt.show()

