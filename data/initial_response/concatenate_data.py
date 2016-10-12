#!/usr/bin/env python3
#concatenate data 
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


#arguments: filenames to concatenate data from. Excluding everything starting with '#'
import numpy as np
import sys
import re

fancy=True
filenames = sys.argv[1:]

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.concatenate(rawdata)

if fancy:
    with open(filenames[0], 'r') as f:
        first_line = f.readline()
        second_line = f.readline()
    
    numbers=re.findall('\d+', filenames[0])
    names=[]
    for name in second_line.split('<')[1:]:
        name=name.split('>')[0]
        if name=="burden" or name=="relapse":
            name=name+str(numbers[0])
        name=name.replace(" ","")
        names.append(name)
        print(name,end=" ")
    print()

for x in data:
    for y in x:
        print(y, end=" ")
    print()
