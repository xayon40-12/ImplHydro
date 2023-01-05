#!/usr/bin/env python3
import sys 
import os
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
sys.path.append(os.path.abspath("%s/.config/" % (os.path.expanduser("~"))))
import plt_setting
import numpy as np
from collections import defaultdict
from math import sqrt

plt.rcParams['axes.grid'] = False

IDx = 0
IDy = 1
IDiter = 2
IDt00 = 3
ID1De = 6
ID2De = 7

def dd(n):
    if n == 1:
        return {}
    else:
        return defaultdict(lambda: dd(n-1))

datas = dd(8)

def convert(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except:
            return v

dir = "results/"
for d in os.listdir(dir):
    p = dir+d+"/"+os.listdir(dir+d)[0]
    data = np.loadtxt(p+"/data.txt")
    info = {k: convert(v.strip()) for [k, v] in np.loadtxt(p+"/info.txt", dtype=object, delimiter=":")}

    t0 = info["t0"]
    tend = info["tend"]
    scheme = info["scheme"]
    maxdt = info["maxdt"]
    dx = info["dx"]
    nx = info["nx"]
    info["l"] = str(dx*nx)
    if info["ny"] == 1:
        dim = "1D"
    else:
        dim = "2D"
    info["dim"] = dim
    name = info["name"]
    
    # print(p, dim, integration, t0, tend, dx, nx, maxdt)
    datas[dim][name][t0][tend][dx][nx][scheme][maxdt] = (info, data)

def compare(i, vss, wss):
    maxerr = 0
    meanerr = 0
    for (vs, ws) in zip(vss,wss):
        a = vs[i]
        b = ws[i]
        err = abs(a-b)/max(abs(a),abs(b))
        maxerr = max(err,maxerr)
        meanerr += err
    meanerr /= len(vss)
    return (maxerr, meanerr)

def convergence(a, ref=None):
    maxdts = sorted([dt for dt in a])
    if ref is None:
        ref = a[maxdts[0]][1]
    rmaxdts = reversed(maxdts)
    all = []
    for i in rmaxdts:
        (info, v) = a[i]
        cost = info["cost"]
        dt = info["maxdt"]
        (maxerr, meanerr) = compare(IDt00, ref, v)
        # all += [(v, info, cost, maxerr, meanerr)]
        all += [(v, info, dt, cost, maxerr, meanerr)]
    
    return np.array(all, dtype=object)

def info2name(info, scheme=True):
    if scheme:
        return "{}_{}_{}_{}_{}_{}_{}".format(info["dim"],info["name"],info["scheme"],info["t0"],info["tend"],info["dx"],info["nx"])
    else:
        return "{}_{}_{}_{}_{}_{}".format(info["dim"],info["name"],info["t0"],info["tend"],info["dx"],info["nx"])

def convall(l, d):
    [dim,name,t0,tend,dx,nx] = l
    scs = sorted(list(d.keys()))
    any = d[scs[0]]
    dtref = sorted(list(any.keys()))[0]
    info = any[dtref][0]
    refs = {s: d[s][dtref][1] for s in d}
    
    for (dtcost, dci) in [("dt", 2), ("cost", 3)]:
        plt.figure()
        for s0 in scs:
            # for s1 in scs:
            for s1 in [scs[0]]:
                c = convergence(d[s0],refs[s1])
                plt.loglog(c[:,dci],c[:,4], label="{} r {}".format(s0, s1))
        plt.xlabel(dtcost)
        plt.ylabel("relative error")
        plt.title("{} {} t0={} tend={} dx={} cells={}".format(dim, name, t0, tend, dx, nx))
        plt.legend()
        plt.savefig("figures/convergence_{}_{}.pdf".format(dtcost, info2name(info, False)))
        plt.close()

def plot1d(datadts):
    maxdts = sorted([dt for dt in datadts])
    (info, data) = datadts[maxdts[0]]
    tend = info["tend"]
    
    x = data[:,IDx]/tend
    iter = data[:,IDiter]
    y = data[:,IDt00]
    _,ax = plt.subplots()
    l1 ,= ax.plot(x,y)
    ax.set_ylabel("e")
    ax.set_xlabel("x/t")
    ax2 = ax.twinx()
    l2 ,= ax2.plot(x,iter, 'o', color="gray")
    plt.legend([l1,l2],["numerics","iterations"])
    plt.savefig("figures/best_{}.pdf".format(info2name(info)))
    plt.close()

def plot2d(datadts):
    maxdts = sorted([dt for dt in datadts])
    (info, data) = datadts[maxdts[0]]

    n = info["nx"]
    x = np.reshape(data[:,IDx], (n,n))
    y = np.reshape(data[:,IDy], (n,n))
    z = np.reshape(data[:,IDt00], (n,n))
    ziter = np.reshape(data[:,IDiter], (n,n))
    l = x[0][0]
    r = x[0][-1]
    d = y[0][0]
    u = y[-1][0]
    for (n, z) in [("t00", z), ("iter", ziter)]:
        fig, ax = plt.subplots()
        pos = ax.imshow(z, extent=[l,r,d,u])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        cbar = plt.colorbar(pos)
        cbar.set_label(n)
        plt.savefig("figures/best_{}_{}.pdf".format(n, info2name(info)))
        plt.close()

def plotall(l, d):
    if l[0] == "1D":
        plot1d(d)
    else:
        plot2d(d)

def alldata(n, data, f):
    def _all(n, l, d):
        if n == 0:
            f(l, d)
        else:
            for p in d:
                _all(n-1, l+[p], d[p])
    _all(n, [], data)
        
try:
    os.mkdir("figures")
except FileExistsError:
    None
    
alldata(6, datas, convall)
alldata(7, datas, plotall)
