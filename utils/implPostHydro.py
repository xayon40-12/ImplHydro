#!/usr/bin/env python3

import os
import numpy as np
import frzout

# run frsout
# git: Duke-QCD/trento

print("Running frzout:")
surface_data = np.fromfile('surface.dat', dtype='float64').reshape((-1, 16))

# extract usual sub-arrays
x, sigma, v, _ = np.hsplit(surface_data, [3, 6, 8])
tend = x[-1][0]

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

hrg = frzout.HRG(.155, species='urqmd', res_width=True)
parts = frzout.sample(surface, hrg)
nparts = len(parts)
lines = "# {}\n".format(nparts)
for ID, x, p in parts:
   l = "{} {} {} {} {} {} {} {} {}\n".format(ID, x[0], x[1], x[2], x[3], p[0], p[1], p[2], p[3])
   lines += l

with open("particles_in.dat", "w") as f:
   f.write(lines)

print("tend: {}\nnb parts: {}\n".format(tend, nparts))


# run urqmd
# git: jbernhard/urqmd-afterburner
#
# urqmd needs to be installed before. To compile with recent gcc-fortran,
# the following flags might be needed:
#   export FCFLAGS="-w -fallow-argument-mismatch -O2"
#   export FFLAGS="-w -fallow-argument-mismatch -O2"

print("Running UrQMD:")
urqmd = os.popen("afterburner particles_in.dat particles_out.dat")
output = urqmd.read()
with open("urqmd.log", "w") as f:
   f.write(output)
print("done\n")