#!/usr/bin/env python3

import sys
import numpy as np

surface_data = np.fromfile(sys.argv[1], dtype='float64').reshape((-1, 22))

x, _ = np.hsplit(surface_data, [3])
 
tend = x[-1][0]
print(tend)
