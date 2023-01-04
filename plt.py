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

def convert(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except:
            return v
    
datas = dd(7)

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
    datas[dim][integration][t0][tend][dx][nx][maxdt] = (info, data)

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
    print("cost", "maxerr", "meanerr")
    rmaxdts = reversed(maxdts)
    all = []
    for i in rmaxdts:
        (info, v) = a[i]
        cost = info["cost"]
        (maxerr, meanerr) = compare(2, ref, v)
        print("{} {:e} {:e}".format(cost, maxerr, meanerr))
        all += [(v, info, maxerr, meanerr)]
    
    return all

def testconv(dim, t0, tend, dx, n):
    print("convergence "+dim+":")
    impl = datas[dim]["FixPoint"][t0][tend][dx][n]
    expl = datas[dim]["Explicit"][t0][tend][dx][n]
    dt = sorted([dt for dt in impl])[0]
    print("implicit:")
    cimpl = convergence(impl)
    print("explicit:")
    cexpl = convergence(expl)
    (maxerr, meanerr) = compare(2, impl[dt][1],expl[dt][1])
    print("comparison explicit implicit:\n    maxerr {:e}\n    meanerr {:e}".format(maxerr, meanerr))
    
    cibesterr = cimpl[-2][2]
    cebesterr = cexpl[-2][2]
    if cibesterr < cebesterr:
        best = cimpl[-1]
        print("explicit with implicit as ref:")
        cexpl = convergence(expl, best[0])
        print("implicit with explicit as ref:")
        convergence(impl, cexpl[-1][0])
    else:
        best = cexpl[-1]
        print("implicit with explicit as ref:")
        cimpl = convergence(impl, best[0])
    
    print()
    return best, cimpl, cexpl

testconv("1D", 1.0, 4.5, 0.1, 100)
testconv("2D", 1.0, 4.5, 0.1, 100)

def info2name(info):
    return "{}_{}_{}_{}_{}_{}".format(info["dim"],info["integration"],info["t0"],info["tend"],info["dx"],info["nx"])

def plot1d(datadts):
    maxdts = sorted([dt for dt in datadts])
    (info, data) = datadts[maxdts[0]]

    x = data[:,0]
    y = data[:,2]
    plt.figure()
    plt.plot(x,y)
    plt.savefig("figures/{}.pdf".format(info2name(info)))
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
    plt.savefig("figures/{}.pdf".format(info2name(info)))
    plt.close()

def alldata(data, f):
    def _all(l, d):
        if type(d) is dict:
            f(l, d)
        else:
            for p in d:
                _all(l+[p], d[p])
    _all([],data)
        
def plotall(l, d):
    if l[0] == "1D":
        plot1d(d)
    else:
        plot2d(d)
try:
    os.mkdir("figures")
except FileExistsError:
    None
    
alldata(datas, plotall)
