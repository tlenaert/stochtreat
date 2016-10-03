#!/usr/bin/env python3
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

"""Scipt takes the timepoints from multiple 
files, and calcultes the corresponding time histogram."""
import numpy as np
import sys

filenames = sys.argv[1:]
t_max=20.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#")[0,:] for filename in filenames]

data=np.concatenate(rawdata)
data[data<0.]=np.inf

times=np.linspace(0.,t_max,t_steps)

time_frac=[]
no_patients=len(data)

for t in times:
    no_timepoints=(data<t).sum()
    time_frac.append(no_timepoints/no_patients)

for t,l in zip(times,time_frac):
    print(t,l)


