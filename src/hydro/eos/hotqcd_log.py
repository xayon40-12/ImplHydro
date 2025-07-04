#!/usr/bin/env python3

import numpy as np
import scipy

HBARC = 0.1973; # GeV.fm
hotqcd = np.fromfile("hrg_hotqcd_eos_binary.dat", dtype="float64").reshape((-1,4))
e = np.log(hotqcd[:,0] / HBARC)
p = np.log(hotqcd[:,1] / HBARC)
s = np.log(hotqcd[:,2])
t = np.log(hotqcd[:,3] / HBARC)

pe = scipy.interpolate.CubicSpline(e, p)
se = scipy.interpolate.CubicSpline(e, s)
te = scipy.interpolate.CubicSpline(e, t)

n = 1000
le = np.linspace(e[0], e[-1], n)
lp = pe(le)
ls = se(le)
lt = te(le)

log_hotqcd = np.array(range(4*n), dtype="float64").reshape((-1,4))
log_hotqcd[:,0] = le
log_hotqcd[:,1] = lp
log_hotqcd[:,2] = ls
log_hotqcd[:,3] = lt

log_hotqcd.tofile("hrg_hotqcd_eos_binary_log.dat")
