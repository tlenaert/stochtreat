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

""" analyses treatment success of several files and writes to stdout."""
import numpy as np
import helperfunctions
import sys

filenames = helperfunctions.natural_sort(sys.argv[1:])
treattime_range=np.linspace(0.,1.,50.)

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.array(rawdata)

probs=data[:,0]
# print(probs)
# exit(0)


for t,p in zip(treattime_range,probs):
    print(t,p)


