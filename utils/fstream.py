#!/usr/bin/env python3

import freestream
import sys
import numpy as np

normalization = 16
ic = np.loadtxt("00.dat")*normalization
half_size = 10.0 # fm
tau_fs = 0.5
fs = freestream.FreeStreamer(ic, half_size, tau_fs)

e = fs.energy_density()
ut = fs.flow_velocity(0)
ux = fs.flow_velocity(1)
uy = fs.flow_velocity(2)
def pi(u,v):
    return fs.shear_tensor(u,v)
pi00 = pi(0,0)
bulk = fs.bulk_pressure() # fs.bulk_pressure(eos)

dx = 0.2
etot = e.sum()*197*dx**3/1000 # GeV
em = e.max()*197/1000 # GeV fm^-3
print("e_tot: {} GeV".format(etot))
print("e_max: {} GeV fm^-3".format(em))

