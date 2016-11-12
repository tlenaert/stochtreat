#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys

array_to_find=np.array([0.00266636,0.00119456,0.000851222])
delta=1e-4

epsilon_i_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
epsilon_c_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
for epsilon_c in epsilon_c_range:
    for epsilon_i in epsilon_i_range:
        if epsilon_i<epsilon_c:
            continue

filenames = helperfunctions.natural_sort(sys.argv[1:])
t_max=10.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#")[:,0] for filename in filenames]
data=np.concatenate(rawdata)

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


