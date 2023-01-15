#!/usr/bin/env python3
import sys 
import os
import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.append(os.path.abspath("%s/.config/" % (os.path.expanduser("~"))))
import plt_setting
import numpy as np
from collections import defaultdict
from math import sqrt
from riemann import riemann
from gubser import gubser

# plt.rcParams['axes.grid'] = False

IDx = 0
IDy = 1
IDiter = 2
IDt00 = 3
ID1De = 5
ID2De = 6

crop = 6

def dd(n):
    if n == 1:
        return {}
    else:
        return defaultdict(lambda: dd(n-1))

datas = dd(10)

def convert(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except:
            return v

def extractCase(n):
    l = n[-1]
    if l in "0123456789":
        (name, v) = extractCase(n[:-1])
        return (name, int(l)+10*v)
    else:
        return (n,0)

dir = "results/"
for d in os.listdir(dir):
    p = dir+d+"/"+os.listdir(dir+d)[0]
    info = {k: convert(v.strip()) for [k, v] in np.loadtxt(p+"/info.txt", dtype=object, delimiter=":")}

    t0 = info["t0"]
    tend = info["tend"]
    scheme = info["scheme"]
    maxdt = info["maxdt"]
    dx = info["dx"]
    nx = info["nx"]
    t = info["t"]
    info["l"] = str(dx*nx)
    if info["ny"] == 1:
        dim = "1D"
        n = nx
    else:
        dim = "2D"
        n = nx*nx
    info["dim"] = dim
    name = info["name"]
    (name, case) = extractCase(name)
    info["name"] = name
    info["case"] = case

    # data = np.loadtxt(p+"/data.txt")
    data = np.fromfile(p+"/data.dat", dtype="float64").reshape((n,-1))
    
    # print(p, dim, integration, t0, tend, dx, nx, maxdt)
    datas[dim][name][t0][tend][dx][nx][t][case][scheme][maxdt] = (info, data)

def compare(i, vss, wss):
    maxerr = 0
    meanerr = 0
    count = 0
    for (vs, ws) in zip(vss,wss):
        if vs[0] >= -crop and vs[0] <= crop and vs[1] >= -crop and vs[1] <= crop:
            count += 1
            a = vs[i]
            b = ws[i]
            err = abs(a-b)/max(abs(a),abs(b))
            maxerr = max(err,maxerr)
            meanerr += err
    meanerr /= count
    return (maxerr, meanerr)

def convergence(a, ref=None):
    maxdts = sorted([dt for dt in a])
    if ref is None:
        ref = a[maxdts[0]][1]
    all = []
    for i in maxdts:
        (info, v) = a[i]
        cost = info["cost"]
        dt = info["maxdt"]
        avdt = (info["tend"]-info["t0"])/info["tsteps"]
        elapsed = info["elapsed"]
        if info["dim"] == "1D":
            id = ID1De
        else:
            id = ID2De
        (maxerr, meanerr) = compare(id, ref, v)
        # all += [(v, info, cost, maxerr, meanerr)]
        all += [(v, info, maxerr, meanerr, dt, cost, avdt, elapsed)]
    
    return np.array(all, dtype=object)

def info2name(info, scheme=True):
    if scheme:
        return "{}_{}{}_{}_{}_{}_{}_{}_{}".format(info["dim"],info["name"],info["case"],info["scheme"],info["t0"],info["tend"],info["dx"],info["nx"],info["t"])
    else:
        return "{}_{}{}_{}_{}_{}_{}_{}".format(info["dim"],info["name"],info["case"],info["t0"],info["tend"],info["dx"],info["nx"],info["t"])

def convall(l, ds):
    [dim,name,t0,tend,dx,nx,t] = l
    for (dtcost, dci) in [("dt", 4), ("cost", 5), ("avdt", 6), ("elapsed", 7)]:
        plt.rcParams["figure.figsize"] = [8, 5]
        plt.figure()
        plt.xlabel(dtcost)
        plt.ylabel("relative error")
        plt.title("{} {} t0={} tend={} dx={} cells={}".format(dim, name, t0, tend, dx, nx))
        d = ds[list(ds)[0]]
        scs = sorted(list(d.keys()))
        scs.reverse()
        if len(scs) <= 1:
            plt.close()
            return 
        for (s0,col) in zip(scs,plt_setting.clist):
            for case in ds: 
                d = ds[case]
                scs = sorted(list(d.keys()))
                scs.reverse()
                if len(scs) <= 1:
                    plt.close()
                    return 
                any = d[scs[0]]
                dtref = sorted(list(any.keys()))[0]
                info = any[dtref][0]
                refs = {s: d[s][sorted(list(d[s].keys()))[0]][1] for s in d}
    
                for s1 in [scs[0]]:
                    c = convergence(d[s0],refs[s1])
                    plt.loglog(c[1:,dci],c[1:,2], 'o', label="{} r {}".format(s0, s1), color=col, linestyle="-.", linewidth=1)
        labels = []
        for p in plt.gca().get_lines():    # this is the loop to change Labels and colors
            label = p.get_label()
            if label in labels:    # check for Name already exists
                p.set_label('_' + label)       # hide label in auto-legend
            else:
                labels += [label]
        plt.legend()
        plt.savefig("figures/convergence_{}_{}.pdf".format(dtcost, info2name(info, False)))
        plt.close()

def plot1d(datadts):
    maxdt = sorted([dt for dt in datadts])[0]
    for dt in datadts:
        timename = "dt{}".format(dt)
        if dt == maxdt:
            timename = "best_"+timename
        (info, data) = datadts[dt]
        t0 = info["t0"]
        tend = info["tend"]
        n = info["nx"]
    
        name = info["name"]
        void = "Void" in name
        e0 = 10
        emin = 1
        cs2 = 1/3
        eps = 1e-10

        e, v = riemann(e0,emin,cs2,eps,void)
        def p(x):
            return cs2*e(x)
        def t00(x):
            vx = v(x)
            if vx == 1:
                return 0
            else:
                ut2 = 1/(1-v(x)*v(x))
                return (e(x)+p(x))*ut2-p(x)

    
        x = data[:,IDx]
        iter = data[:,IDiter]
        y = data[:,ID1De]
        nl = next(i for (i,v) in zip(range(n),x) if v >= -crop)
        nr = n-1-nl
        x = x[nl:nr]/(tend-t0)
        iter = iter[nl:nr]
        y = y[nl:nr]
        ycontinuum = [e(x) for x in x]
        plt.rcParams["figure.figsize"] = [8, 5]
        _,axs = plt.subplots(2, 1, sharex=True)
        iterations ,= axs[0].plot(x,iter, 'o', color="gray", label="iterations")
        axs[0].set_ylabel("iterations")
        axs[0].legend()
        continuum ,= axs[1].plot(x,ycontinuum, color="black", label="continuum")
        numerics ,= axs[1].plot(x,y, label="numerics")
        axs[1].set_ylabel("e")
        axs[1].set_xlabel("x/t")
        axs[1].legend()
        plt.savefig("figures/{}_{}.pdf".format(timename,info2name(info)))
        plt.close()

def plot2d(datadts):
    maxdts = sorted([dt for dt in datadts])
    (info, ref) = datadts[maxdts[0]]
    (_, data) = datadts[maxdts[1]]

    t = info["tend"]
    n = info["nx"]
    x = data[:,IDx]
    y = data[:,IDy]
    z = data[:,ID2De]
    zref = ref[:,ID2De]
    ziter = data[:,IDiter]
    zgubser = [gubser(x,y,t) for (x,y) in zip(x,y)] # this is energy density not t00
    zerr = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zgubser)]
    zerrref = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zref)]
    nl = next(i for (i,v) in zip(range(n*n),x) if v >= -crop)
    nr = n-1-nl
    x = np.reshape(x, (n,n))[nl:nr,nl:nr]
    y = np.reshape(y, (n,n))[nl:nr,nl:nr]
    z = np.reshape(z, (n,n))[nl:nr,nl:nr]
    ziter = np.reshape(ziter, (n,n))[nl:nr,nl:nr]
    zgubser = np.reshape(zgubser, (n,n))[nl:nr,nl:nr]
    zerr = np.reshape(zerr, (n,n))[nl:nr,nl:nr]
    zerrref = np.reshape(zerrref, (n,n))[nl:nr,nl:nr]
    l = x[0][0]
    r = x[0][-1]
    d = y[0][0]
    u = y[-1][0]
    all = [("e", z), ("err ref", zerrref)]
    if "Gubser" in info["name"] and not "Exponential" in info["name"]:
        all += [("err continuum", zerr)]
    if not "Heun" in info["scheme"]:
        all += [("iter", ziter)]
    nb = len(all)
    plt.rcParams["figure.figsize"] = [2+nb*4, 5]
    fig, axs = plt.subplots(1,nb, sharey=True)
    for (i, (n, z)) in zip(range(nb),all):
        im = axs[i].imshow(z, extent=[l,r,d,u], norm=CenteredNorm(0), cmap="terrain")
        axs[i].set_xlabel("x")
        axs[i].xaxis.tick_top()
        axs[i].xaxis.set_label_position('top') 
        if i == 0:
            axs[i].set_ylabel("y")
        divider = make_axes_locatable(axs[i])
        cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
        fig.add_axes(cax)
        cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
        cbar.formatter.set_powerlimits((0, 0))
        cbar.formatter.set_useMathText(True)
        cbar.update_ticks()
        cbar.set_label(n, labelpad=-60)
    plt.savefig("figures/best_e_{}.pdf".format(info2name(info)), dpi=100)
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
    
alldata(7, datas, convall)
alldata(9, datas, plotall)
