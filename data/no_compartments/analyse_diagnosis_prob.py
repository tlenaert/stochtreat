#!/usr/bin/env python3

import numpy as np
import sys

filenames = sys.argv[1:]
t_max=20.
t_steps=100

diagfracs=[]
for filename in filenames:
    with open(filename) as f:
        for line in f:
            if line.startswith("#Fraction diagnosed"):
                diagfracs.append(float(line.split()[-1]))
                break

print(np.mean(diagfracs)/100.)

