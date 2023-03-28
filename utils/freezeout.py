#!/usr/bin/env python3

import numpy as np
import frzout

# data file is binary, not text
surface_data = np.fromfile('surface.dat', dtype='float64').reshape((-1, 16))

# extract usual sub-arrays
x, sigma, v, _ = np.hsplit(surface_data, [3, 6, 8])

# create mapping of pi components
pi = dict(
   xx=surface_data[:, 11],
   xy=surface_data[:, 12],
   yy=surface_data[:, 13]
)

# extract bulk pressure
Pi = surface_data[:, 15]

# create Surface object
surface = frzout.Surface(x, sigma, v, pi=pi, Pi=Pi)

hrg = frzout.HRG(.155)
parts = frzout.sample(surface, hrg)
for ID, x, p in parts:
   print(ID, x, p)