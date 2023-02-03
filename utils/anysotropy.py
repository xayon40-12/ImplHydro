#!/usr/bin/env python3
import sys 
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

nl = 3
lstyles = [(3*i,(3,5,1,5)) for i in range(nl)]

for (n,linestyle) in zip([d for d in os.listdir('.') if os.path.isdir(d)],lstyles):
    anys = []
    anys2 = []
    count = 0
    dir = "{}/results/".format(n)
    for d in os.listdir(dir):
        dird = dir+d
        ts = sorted(os.listdir(dird), key=float)
        any = np.array([np.fromfile(dird+"/"+t+"/momentum_anysotropy.dat", dtype="float64")[0] for t in ts])
        anys += [any]
        anys2 += [any*any]
        count += 1
    
    t0 = float(ts[0])
    ts = [float(t)-t0 for t in ts]
    anys = np.array(anys).sum(axis=0)/count
    anys2 = np.array(anys2).sum(axis=0)/count
    sig = np.sqrt(anys2-anys*anys)
    err = sig/np.sqrt(count)

    plt.errorbar(ts, anys, err, label=n, linestyle=linestyle)
plt.xlabel("t-t0")
plt.ylabel(r"$\frac{<T^{11}-T^{22}>}{<T^{11}+T^{22}>}$")
plt.legend()
plt.title("trento Pb Pb for L=20fm, 100cells, b=7fm, 100 events")
plt.savefig("momentum_anysotropy.pdf")
plt.close()