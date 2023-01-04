#!/usr/bin/env python3
import sys 
import os
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
sys.path.append(os.path.abspath("%s/.config/" % (os.path.expanduser("~"))))
import plt_setting
import numpy as np
from collections import defaultdict

def dd(n):
    if n == 1:
        return {}
    else:
        return defaultdict(lambda: dd(n-1))

datas = dd(7)

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
    integration = info["integration"]
    maxdt = info["maxdt"]
    dx = info["dx"]
    nx = info["nx"]
    info["l"] = str(dx*nx)
    if info["ny"] == 1:
        dim = "1D"
    else:
        dim = "2D"
    info["dim"] = dim
    
    # print(p, dim, integration, t0, tend, dx, nx, maxdt)
    datas[dim][t0][tend][dx][nx][integration][maxdt] = (info, data)

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
        (maxerr, meanerr) = compare(2, ref, v)
        all += [(v, info, cost, maxerr, meanerr)]
    
    return np.array(all, dtype=object)

def info2name(info, integration=True):
    if integration:
        return "{}_{}_{}_{}_{}_{}".format(info["dim"],info["integration"],info["t0"],info["tend"],info["dx"],info["nx"])
    else:
        return "{}_{}_{}_{}_{}".format(info["dim"],info["t0"],info["tend"],info["dx"],info["nx"])

def convall(l, d):
    [dim,t0,tend,dx,nx] = l
    impl = d["FixPoint"]
    expl = d["Explicit"]
    dt = sorted([dt for dt in impl])[0] # the dt should be the same for expl and impl
    cimpl = convergence(impl)
    cexpl = convergence(expl)
    
    cimplref = cimpl[-1][0]
    cexplref = cexpl[-1][0]
    cimplrexpl = convergence(impl, cexplref)
    cexplrimpl = convergence(expl, cimplref)
    
    plt.figure()
    for (a,lab) in [(cimpl,"impl r impl"),(cimplrexpl,"impl r expl"),(cexpl,"expl r expl"),(cexplrimpl,"expl r impl")]:
        plt.loglog(a[:,2],a[:,3], label=lab)
    plt.xlabel("cost")
    plt.ylabel("relative error")
    plt.title("{} t0={} tend={} dx={} cells={}".format(dim, t0, tend, dx, nx))
    plt.legend()
    plt.savefig("figures/convergence_{}.pdf".format(info2name(cimpl[-1][1], False)))
    plt.close()

def plot1d(datadts):
    maxdts = sorted([dt for dt in datadts])
    (info, data) = datadts[maxdts[0]]
    tend = info["tend"]

    x = data[:,0]/tend
    y = data[:,2]
    plt.figure()
    plt.plot(x,y)
    plt.xlabel("x/t")
    plt.ylabel("t00")
    plt.savefig("figures/best_{}.pdf".format(info2name(info)))
    plt.close()

def plot2d(datadts):
    maxdts = sorted([dt for dt in datadts])
    (info, data) = datadts[maxdts[0]]

    n = info["nx"]
    x = np.reshape(data[:,0], (n,n))
    y = np.reshape(data[:,1], (n,n))
    z = np.reshape(data[:,2], (n,n))
    l = x[0][0]
    r = x[0][-1]
    d = y[0][0]
    u = y[-1][0]
    plt.figure()
    plt.imshow(z, extent=[l,r,d,u])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("figures/best_{}.pdf".format(info2name(info)))
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
    
alldata(5, datas, convall)
alldata(6, datas, plotall)
