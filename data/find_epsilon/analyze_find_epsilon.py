#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys

# array_to_find=np.array([0.00266636,0.00119456,0.000851222])
array_to_find=np.array([0.00261841,0.00115792,0.000802662 ])
delta=5e-6

##get tested epsilon values
epsilon_min=0.5
epsilon_max=1.0
epsilonvalues=20
epsilon_i_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
epsilon_c_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
eps_list=[]
for epsilon_c in epsilon_c_range:
    for epsilon_i in epsilon_i_range:
        if epsilon_i<epsilon_c:
            continue
        eps_list.append([epsilon_c,epsilon_i])
eps_list=np.array(eps_list)

#process file data
filenames = helperfunctions.natural_sort(sys.argv[1:])
rawdata=[[np.loadtxt(filename,comments="#")] for filename in filenames]
data=np.concatenate(rawdata)

a=data <  array_to_find+delta
b=data >  array_to_find-delta
c=np.logical_and(a,b)

resultindeces=np.where(c.all(axis=1))
result=eps_list[resultindeces]
for y in result:
    for x in y:
        print(x,end=" ")
    print()


