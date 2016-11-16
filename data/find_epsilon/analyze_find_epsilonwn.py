#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys

# array_to_find=np.array([0.00266636,0.00119456,0.000851222])
to_find_data=[
[0.119944948121,0.0318476722642,0.0255437284343],#sokal LOW 
[0.0788205700963,0.0373155705618,0.029792798773],#sokal INT 4
[0.248709037978,0.171248812015,0.0696353922543],#sokal HIGH
[0.119117897496,0.0311775891388,0.0391527108148], #euro LOW 0
[0.0933576305015,0.0753268313524,0.0255725420811],#euro INT 0
[0.502049958281,0.204397630799,0.108763321483]]  #euro HIGH 
array_to_find=np.array(to_find_data[int(sys.argv[1])])

delta_start=5e-1
delta_step=1e-5

##get tested epsilon values
epsilon_min=0.5
epsilon_max=1.0
epsilonvalues=30
epsilon_i_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
epsilon_c_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
epsilon_n_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
eps_list=[]
for epsilon_c in epsilon_c_range:
    for epsilon_i in epsilon_i_range:
        for epsilon_n in epsilon_n_range:
            if epsilon_i<epsilon_c or epsilon_n<epsilon_c:
                continue
            eps_list.append([epsilon_c,epsilon_i,epsilon_n])
eps_list=np.array(eps_list)

#process file data
filenames = helperfunctions.natural_sort(sys.argv[2:])
rawdata=[[np.loadtxt(filename,comments="#")] for filename in filenames]
data=np.concatenate(rawdata)
eps_list=eps_list[~np.isnan(data).any(axis=1)]
data= data[~np.isnan(data).any(axis=1)]

def find_closest_values(data,array_to_find,delta):
    a=data <  array_to_find+delta
    b=data >  array_to_find-delta
    c=np.logical_and(a,b)
    
    resultindices=np.where(c.all(axis=1))[0]
    return(resultindices)

number=np.inf
d=delta_start
inx=[]
while (number > 1):
    oldinx=inx
    inx=find_closest_values(data,array_to_find,d)
    number=len(inx)
    d=d-delta_step
    if (number<1):
        inx=oldinx
resultindices=inx

result_eps=eps_list[resultindices]
resultdata=data[resultindices]

print("#finding ",array_to_find,"precision="+str(d))
for epslist,indice in zip(result_eps,resultindices):
    for x in epslist:
        print(x,end=" ")
    for y in data[indice]:
        print(y,end=" ")
    print()


