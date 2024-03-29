#!/usr/bin/env python3
import sys 
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

a = 3
nl = 4
b = a*nl
lstyles = [(a*i,(a,b,1,b)) for i in range(nl)]

for (i, (n,linestyle)) in zip(range(100000),zip([d for d in os.listdir('.') if os.path.isdir(d)],lstyles)):
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

    plt.errorbar(ts, anys, err, label=n, linestyle=linestyle, elinewidth=1, errorevery=(i,nl))
plt.xlabel("t-t0")
plt.ylabel(r"$\frac{<T^{11}-T^{22}>}{<T^{11}+T^{22}>}$")
plt.legend()
plt.title("trento Pb Pb for L=20fm, 100cells, b=7fm, 100 events")
plt.savefig("momentum_anysotropy.pdf")
plt.close()