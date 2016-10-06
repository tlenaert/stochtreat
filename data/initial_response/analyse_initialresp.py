#!/usr/bin/env python3
#concatenate data and split into groups
#Copyright (C) 2016  Marvin A. Böttcher
#This program is free software: you can redistribute it and/or modify    
#it under the terms of the GNU General Public License as published by    
#the Free Software Foundation, either version 3 of the License, or    
#(at your option) any later version.    
#
#This program is distributed in the hope that it will be useful,    
#but WITHOUT ANY WARRANTY; without even the implied warranty of    
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    
#GNU General Public License for more details.    
#
#You should have received a copy of the GNU General Public License    
#along with this program. If not, see <http://www.gnu.org/licenses/>.


#first argument <treattime>, following arguments: filenames
import numpy as np
import sys

treattime=sys.argv[1]
filenames = sys.argv[2:]
t_max=20.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.concatenate(rawdata)

# print(np.where(data[:,1]==1))
relapse=data[np.where(data[:,-1]==1)]
norelapse=data[np.where(data[:,-1]==0)]

np.savetxt("results/relapsedat"+str(treattime)+".txt",relapse)
np.savetxt("results/norelapsedat"+str(treattime)+".txt",norelapse)
# for x in relapse:
#     for y in x:
#         print(y, end=" ")
#     print()
#
# for x in norelapse:
#     for y in x:
#         print(y, end=" ")
#     print()
