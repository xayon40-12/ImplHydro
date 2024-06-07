#!/usr/bin/env python3
import sys 
import os
import warnings
import matplotlib
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plt_setting
import numpy as np
from collections import defaultdict
from math import sqrt
from riemann import riemann
from gubser import gubser
from scipy import signal
from math import ceil,log
from enum import Enum

warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning) # disable matplotlib deprecation warning
np.seterr(invalid='ignore') # disable invalid waring for numpy as NaN are used to discard data in the void
np.seterr(divide='ignore') # disable divide by zero worining


animate = False
def setAnimate():
    global animate
    animate = True
manycases = False
def setManyCases():
    global manycases
    manycases = True
rejectfails = False
def setRejectFails():
    global rejectfails
    rejectfails = True
useschemename = False
def setUseSchemeName():
    global useschemename
    useschemename = True
t_id = -1
def setTimeId(i):
    global t_id
    t_id = int(i)
    

argActions = [(["-a","--animate"], setAnimate),(["-m","--manycases"], setManyCases), (["-r","--rejectfails"], setRejectFails), (["-s","--schemename"], setUseSchemeName), (["-t","--time_id"], setTimeId)]
for arg in sys.argv[1:]:
    found = False
    for (larg, act) in argActions:
        for l in larg:
            if arg.startswith(l):
                found = True
                if len(l) < len(arg):
                    act(arg[len(l)+1:])
                else:
                    act()
    if not found:
        print("invalid arg: ", arg)
        exit(-1)

# crop = 9
crop = 19
crop_eta = crop/4
defaultfromref = 0
maxerr = 1e-6
num2D = 5
compare_type = 0 # 0: absolute, 1: relative
use_square = False
use_average_ref = False
plot_meanmax_dx = False
plot_conv = True
plot_1D = True
plot_2D = False
# plot_2D = True
def refName(scs):
    # return scs[0]
    # return scs[-1]
    # return "RadauIIA2"
    # return "GL2"
    return "GL1"
    # return "Heun"

e0 = 10
emin = 1
cs2 = 1/3
cs = np.sqrt(cs2)
eps = 1e-10

riem_e, riem_v, riem_ep, riem_vs = riemann(e0,emin,cs2,eps,False)
riem_void_e, riem_void_v, riem_void_ep, riem_void_vs = riemann(e0,emin,cs2,eps,True)

def dd(n):
    if n == 1:
        return {}
    else:
        return defaultdict(lambda: dd(n-1))

datas = dd(20)

def convert(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except:
            return v

maxvoidratio = 0
meanvoidratio = 0
countvoidratio = 0
dir = "results/"
for d in filter(lambda d: os.path.isdir(dir+"/"+d), os.listdir(dir)):
    dird = dir+d
    dirs = list(filter(lambda d: os.path.isdir(dird+"/"+d), os.listdir(dird)))
    ts = sorted(dirs, key=float) 
    if len(ts)-1 < t_id:
        continue
    t = ts[t_id]
    # maxt = ts[int(0.3*(len(ts)-1))]

    p = dird+"/"+t
    info = {k: convert(v.strip()) for [k, v] in np.loadtxt(p+"/info.txt", dtype=object, delimiter=":")}

    def round10(x):
        r = 1e10
        return float(round(r*x)/r)
    t0 = info["t0"]
    tend = info["tend"]
    scheme = info["scheme"]
    maxdt = round10(info["maxdt"])
    dx = info["dx"]
    nx = info["nx"]
    ny = info["ny"]
    nz = info["nz"]
    t = info["t"]
    fails = info["fails"]
    t = round10(t)
    info["l"] = str(dx*nx)
    if info["ny"] == 1:
        dim = "1D"
        n = nx
    elif info["nz"] == 1:
        dim = "2D"
        n = nx*ny
    else:
        dim = "3D"
        n = nx*ny*nz
    info["dim"] = dim
    name = info["name"]
    case = info["case"]
    info["variables"] = info["variables"].split(" ")
    vid = {n: i for (i,n) in enumerate(info["variables"])}
    info["ID"] = vid
    visc = info["viscosity"]
    if visc == "Ideal":
        visc = ()
    if visc == "Shear":
        visc = (info["etaovers"])
        cut = info["energycut"]
    if visc == "Bulk":
        visc = (info["zetaovers"])
        cut = info["energycut"]
    if visc == "Both":
        visc = (info["etaovers"],info["zetaovers"])
        cut = info["energycut"]
    info["visc"] = visc

    # data = np.loadtxt(p+"/data.txt")
    data0 = np.fromfile(dird+"/"+ts[0]+"/data.dat", dtype="float64").reshape((n,-1))
    data = np.fromfile(p+"/data.dat", dtype="float64").reshape((n,-1))
    def find_diff(dir):
        try:
            return np.fromfile(dir+"/diff.dat", dtype="float64").reshape((n,-1))
        except FileNotFoundError:
            return None
    diff = find_diff(p)
    
    if case == 0 and plot_2D:
        ld = len(ts)
        # num = num2D
        # if dim == "1D":
        #     nums = np.array(range(ld))
        # else:
        #     nums = np.array([i for i in range(ld) if i%int(ld/(num-1)) == 0]+[ld-1])
        nums = np.array(range(ld))
        datats = [(float(t), np.fromfile(dird+"/"+t+"/data.dat", dtype="float64").reshape((n,-1)),find_diff(dird+"/"+t)) for t in np.array(ts)[nums]]
        info["datats"] = datats

    e0 = data0[:,vid["e"]]
    e = data[:,vid["e"]]
    info["e_max"] = max(e)
    # if the total energy sum(e) is bigger at the end time than the initial time, 
    # it means that explicit failed
    sum_e = sum(e)
    expl_fail = sum_e != sum_e or sum(e) > 1.1*sum(e0)
    # if a fail (decreasing dt) happens in implicit when dt>dx/10, we consider that implicit failed as we want to test large dt and thus do not want dt to be decreased
    impl_fail = False # no dt decrease anymare, so no need for: fails > 0 and maxdt>dx/10
    # impl_fail = fails > 0
    if rejectfails and (expl_fail or impl_fail): 
        continue # skip explicit or implicit that failed

    cut = 1
    err = abs(e.sum()-e[e>cut].sum())/e.sum()
    while err > maxerr: # find the energy cut such that a relative 1e-6 of the total energy is discarded. This is done to discard the erroneous cell due to numerical diffusion in the vacuum
        cut *= 0.5
        err = abs(e.sum()-e[e>cut].sum())/e.sum()
    info["CUT"] = cut
    meanvoidratio += err
    maxvoidratio = max(maxvoidratio, err)
    countvoidratio += 1

    # extract zero rapidity part
    if dim == "3D":
        dim = "2D"
        zs = sorted(set(data[:,vid["z"]]))
        z = zs[int(len(zs)/2)]
        dataz0 = data[:,vid["z"]]==z
        info["data3d"] = data
        if not diff is None:
            diff = diff[dataz0,:]
        data = data[dataz0,:]
        if case == 0 and plot_2D:
            info["datats"] = [(t, d[dataz0,:], diff) for (t,d,diff) in info["datats"]]

    # print(p, dim, t0, tend, dx, nx, maxdt)
    datas[dim][name][visc][t0][tend][t][case][(dx,nx)][scheme][maxdt] = (info, data, diff)

meanvoidratio /= countvoidratio
print("maxvoidratio: ", maxvoidratio, "meanvoidratio: ", meanvoidratio)
with open("voidratio.txt", "w") as fv:
    fv.write("max_void_ratio: {:e}\nmean_void_ratio: {:e}".format(maxvoidratio, meanvoidratio))
print("finished loading")
# sys.exit(0)

def compare(eid, i, cut, vss, wss):
    maxerr = 0
    meanerr = 0
    squareerr = 0
    count = 0
    for (vs, ws) in zip(vss,wss):
        if  vs[eid] > cut and vs[0] >= -crop and vs[0] <= crop and vs[1] >= -crop and vs[1] <= crop and vs[2] >= -crop_eta and vs[2] <= crop_eta:
            count += 1
            a = vs[i]
            b = ws[i]
            if compare_type == 0:
                err = abs(a-b)
                squareerr += (a-b)**2
            if compare_type == 1:
                err = abs(a-b)/(max(abs(a),abs(b))+1e-100)
                squareerr += ((a-b)/((abs(a)+abs(b))/2+1e-100))**2
            maxerr = max(err,maxerr) # or maxerr += err**inf
            meanerr += err
    if count>0:
        meanerr = meanerr/count
        squareerr = sqrt(squareerr/count)
        # as we used max, there is no need for: maxerr = (maxerr/count)**(1/inf) -> where count**(1/inf) = 1
    return (maxerr, meanerr, squareerr)

def convergence(a, ref=None):
    maxdts = sorted(list(a))
    if ref is None:
        ref = a[maxdts[0]]
    all = []
    for i in maxdts:
        (info, v, diff) = a[i]
        if "data3d" in info:
            v = info["data3d"]
            r = ref[0]["data3d"]
        else:
            r = ref[1]
        cost = info["cost"]
        dt = info["maxdt"]
        avdt = (info["tend"]-info["t0"])/info["tsteps"]
        elapsed = info["elapsed"]
        vid = info["ID"]
        eid = vid["e"]
        id = vid["e"]
        cut = info["CUT"]
        (maxerr, meanerr, squareerr) = compare(eid, id, cut, r, v)
        all += [(v, info, maxerr, meanerr, squareerr, dt, cost, avdt, elapsed)]
    
    return np.array(all, dtype=object), ["v", "info", "maxerr", "meanerr", "squareerr", "dt", "cost", "avdt", "elapsed"]

def info2name(info, scheme=True):
    visc = info["viscosity"]
    if visc == "Ideal":
        visc = "Ideal"
    if visc == "Shear":
        visc = "Shear({})".format(info["etaovers"])
    if visc == "Bulk":
        visc = "Bulk({})".format(info["zeta"])
    if visc == "Both":
        visc = "Both({},{})".format(info["etaovers"],info["zetaovers"])
    if scheme:
        return "{}_{}{}_{}_{}_{}_{}_{}_{}_{:.4e}".format(info["dim"],info["name"],info["case"],visc,info["scheme"],info["t0"],info["tend"],info["dx"],info["nx"],info["t"])
    else:
        return "{}_{}{}_{}_{}_{}_{}_{}_{:.4e}".format(info["dim"],info["name"],info["case"],visc,info["t0"],info["tend"],info["dx"],info["nx"],info["t"])

def lighten(c):
    return tuple(min(1,(a+b*0.2)*2) for a, b in zip(c,np.roll(c,1)))

def convall(l, cnds):
    [dim,name,visc,t0,tend,t] = l
    # allx = [("dt", 5), ("cost", 6), ("avdt", 7), ("elapsed", 8)]
    ally = [((1,0), "full", 2)]
    if use_square:
        ally += [((2,2), "none", 4)]
    ally += [((6,6), "none", 3)]
    allx = [("cost", 6)]
    # ally = [("max", 2)]
    fromref = defaultfromref
    e_max = []
    for (dtcost, dci) in allx:
        fig, axs = plt.subplots(2, figsize=(8,9), sharex=True)
        fig.subplots_adjust(wspace=0,hspace=0)
        fig2, ax2 = plt.subplots(figsize=(8,5))
        def setlabel(ax):
            name = dtcost
            if dtcost == "cost":
                name = r"$n_\mathrm{KT}$"
            ax.set_xlabel(name)
            ax.set_ylabel(r"$|\Delta_\mathrm{ref}|$")
        for ax in axs:
            setlabel(ax)
        setlabel(ax2)
        # plt.title("{} {} t0={} tend={} dx={} cells={}".format(dim, name, t0, tend, dx, nx))
        nbcases = len(cnds)
        alpha = sqrt(1/nbcases)
        markers = {}
        for case in cnds:
            nds = cnds[case]
            dxs = sorted(nds,key=lambda x: x[1])
            (maxdx,minn) = dxs[-1]
            for (ax,dxn) in zip(axs,dxs): # make 100 be plotted before 200
                (dx,n) = dxn
                ds = nds[dxn]
                scs = sorted(list(ds.keys()))
                def prio(s):
                    dts = sorted([dt for dt in ds[s]])
                    mindt = dts[0]
                    return integrationPriority(ds[s][mindt][0])
                scs.sort(key=prio)
                if len(scs) <= 1:
                    plt.close()
                    return 
                s1 = refName(scs)
                ax.text(0.03, 0.05, r"$\Delta x = "+str(dx)+"$ fm", color="black", #, bbox={"facecolor": "white", "pad": 10},
                    transform=ax.transAxes, fontsize=22)
                for (linestyle, fillstyle, mmi) in ally:
                    if len(scs) == 3:
                        colors = [plt_setting.clist[i] for i in [0,2,1]]
                    else:
                        colors = plt_setting.clist
                    for (s0,col) in zip(scs,colors): 
                        ds0 = ds[s0]
                        dts = sorted(list(ds0.keys()))
                        dtref = dts[0]
                        info = ds0[dtref][0]
                        ibig = max([i for i, dt in zip(range(1000),dts) if dt < dx/10 + 1e-10])
                        info_big = ds0[dts[ibig]][0]
                        dx = info["dx"]
                        refs = {s: ds[s][sorted(list(ds[s].keys()))[0]] for s in ds}
                        if "FixPoint" in info["integration"]:
                            schemetype = "Implicit"
                        else:
                            schemetype = "Explicit"
                            e_max = e_max + [info_big["e_max"]]
                        if useschemename:
                            schemetype = info["scheme"]
                        def fuse(r1,r2):
                            (a1,b1,c1) = r1
                            (a2,b2,c2) = r2
                            return (a1, (b1+b2)/2, c1)
                        if use_average_ref:
                            c, cname = convergence(ds0, fuse(refs[scs[0]], refs[scs[1]]))
                        else:
                            c, cname = convergence(ds0, refs[s1])
                        ln = np.log
                        if len(c) > 3:
                            ibig = max([i for i, dt in zip(range(1000),c[:,5]) if dt < dx/10 + 1e-10])
                            dt_order = (ln(c[ibig,mmi])-ln(c[1,mmi]))/(ln(c[ibig,5])-ln(c[1,5]))
                            cost_order = (ln(c[ibig,mmi])-ln(c[1,mmi]))/(ln(c[1,6])-ln(c[ibig,6]))
                            print("{: <10}{: <10}\t{:.2}\t{:.2}".format(cname[mmi], s0, dt_order, cost_order))
                        c = np.array(list(filter(lambda v: v[5]+1e-14>=0.1*dx*2**-5, c)), dtype=object) # use dt=0.1dx*2**-5 as reference
                        al = alpha
                        s = 30
                        sizes = [4*s if abs(dt-dx/10)<1e-10 else s for dt in c[fromref:,5]]
                        facecolors = fillstyle
                        def pl(ax,label):
                            nonlocal facecolors
                            if fillstyle == "full":
                                facecolors = col
                            ax.loglog(c[fromref:,dci],c[fromref:,mmi], label=label, color=col, linestyle=(0,linestyle), linewidth=1, alpha=al)
                            markers[label] = pointstyle(info)
                            ax.scatter(c[fromref:,dci],c[fromref:,mmi],sizes, marker=pointstyle(info), facecolors=facecolors, color=col, alpha=al) #, fillstyle=fillstyle, color=col, alpha=al)
                        pl(ax,schemetype)
                        if "FixPoint" in info["integration"]:
                            if n == 100:
                                col = np.roll(col,2) # lighten(col)
                            label = r"$\Delta x = "+str(dx)+"$ fm"
                            pl(ax2,label)

        def line(col,style,marker=None):
            return matplotlib.lines.Line2D([],[], color=col,linestyle=style,marker=marker,markersize=10)
        def clean(ax):
            handles, labels = [], []
            hs, ls = ax.get_legend_handles_labels()
            for (h,l) in zip(hs,ls):    # this is the loop to change Labels and colors
                col = h.get_color()
                marker = markers[l]
                m = max(col)
                if not (l in labels or m == 1):    # check for Name already exists
                    p = line(col,"-",marker)
                    handles += [p]
                    labels += [l]
            leg = ax.legend(handles, labels, loc="upper right")
            for lh in leg.legendHandles:
                lh.set_alpha(1)
        clean(axs[0])
        labels = ["max"]
        if use_square:
            labels += ["rms"]
        labels += ["mean"]
        handles = [line("black",(0,(1,0)))]
        if use_square:
            handles += [line("black",(0,(2,2)))]
        handles += [line("black",(0,(6,6)))]
        axs[1].legend(handles, labels, loc="upper right")
        clean(ax2)
        prefix = ""
        if compare_type == 0:
            prefix += "absolute_"
        if compare_type == 1:
            prefix += "relative_"
        if use_average_ref:
            prefix += "averageref_"
        else:
            prefix += "{}ref_".format(s1)
        fmt_e_max = ""
        last = -1
        for em in e_max:
            if em != last:
                last = em
                fmt_e_max = "{},{:.2f}".format(fmt_e_max, em)
        fmt_e_max = fmt_e_max[1:]
        e_frz = info["freezeoutenergy"]
        fig.suptitle(r"$\tau = {}$   $\epsilon_m = {}$   $\epsilon_f = {:.2f}$".format(t, fmt_e_max, e_frz))
        fig.savefig("figures/{}convergence_{}_meanmax_crop={}_{}.pdf".format(prefix,dtcost, crop, info2name(info, False)))
        if plot_meanmax_dx:
            fig2.savefig("figures/{}convergence_{}_meanmax_dx_crop={}_{}.pdf".format(prefix,dtcost, crop, info2name(info, False)))
        plt.close()

def pointstyle(info):
    markers = ["s","^","o","v","d","*"]
    l = len(markers)
    i = integrationPriority(info)-1
    return markers[i%l]

def integrationPriority(info):
    scheme = info["scheme"]
    integration = info["integration"]
    if "Heun" in scheme:
        return 1
    if "midpoint" in scheme:
        return 2
    elif "GL1" in scheme:
        return 3
    elif "Explicit" in integration:
        return 4
    else:
        return 5

def mask(data,vid,cut):
    m = (data[:,vid["e"]]>cut).astype(int)
    return np.concatenate((data[:,:4], data[:,4:]/m[:,None]), axis=1)

def plot1d(l, nds):
    [dim,name,visc,t0,tend,t,case] = l
    if case > 0 and not manycases:
        return
    nl = 2
    # lstyles = [(nl*i,(nl,(nl-1)*nl)) for i in range(nl)]
    lstyles1 = [(5*i,(3,7)) for i in range(nl)]
    lstyles2 = [(2*i,(1,3)) for i in range(nl)]

    if "Gubser" in name:
        e = lambda x: gubser(x,x,t)
        v = lambda x: gubser_v(x,x,t)
    elif "Void" in name:
        e = riem_void_e
        v = riem_void_v
    else:
        e = riem_e
        v = riem_v

    def diagonal(v):
        v = v.reshape((n,n))
        return np.array([v[i,i] for i in range(n)])

    datas = nds[list(nds)[0]]
    schemes = sorted(list(datas.keys()))
    dtss = []
    for dxn in nds:
        (dx,n) = dxn
        datas = nds[dxn]
        dtss += [sorted([dt for dt in datas[schemes[0]]])]
    def common(a,b):
        if len(a) > len(b):
            tmp = a
            a = b
            b = a
        if len(a) == 0:
            return []
        if a[0] == b[0]:
            return [a[0]]+common(a[1:],b[1:])
        else:
            return common(a,b[1:])
    if len(dtss) > 1:
        dts = common(dtss[0],dtss[1])
    else:
        dts = dtss[0]
    dts = [dts[0],dts[-1]] # only keep best and worst
    mindt = dts[0]
    maxdt = dts[-1]
    schemes.sort(key=lambda s: integrationPriority(datas[s][mindt][0]))
    for dt in dts:
        timename = "dt{:.4e}".format(dt)
        if dt == maxdt:
            timename = "worst_"+timename
        elif dt == mindt:
            timename = "best_"+timename

        fig,axs = plt.subplots(2, 1, figsize=(8,9), sharex=True)
        if "Riemann" in name:
            axin = axs[0].inset_axes([0.06, 0.09, 0.29, 0.47])
            axin.tick_params(labelsize=14)
            axin.locator_params(nbins=2)
            # if not "Void" in name:
            axinv = axs[0].inset_axes([0.55, 0.45, 0.4, 0.5])
            axinv.tick_params(labelsize=14)
            axinv.locator_params(nbins=2)
        l = 0
        for dxn in nds:
            (dx,n) = dxn
            l = max(l,dx*n/2)
        x = np.linspace(-l,l,1000)
        if not "Gubser" in name:
            x = x/(t-t0)
        yexact = [e(x) for x in x]
        exact ,= axs[0].plot(x,yexact, color="black", label="exact", linewidth=1)
        if "Riemann" in name:
            axin.plot(x,yexact, color="black", label="exact", linewidth=1)
            # if not "Void" in name:
            axinv.plot(x,yexact, color="black", label="exact", linewidth=1)
        dxs = []
        for dxn in sorted(nds,key=lambda x: x[1]):
            (dx,n) = dxn
            dxs.append(dx)
            datas = nds[dxn]
            (info,ref,diffref) = datas[schemes[0]][mindt]
            vid = info["ID"]
            cut = info["CUT"]
            ref = mask(ref,vid,cut)

            x = ref[:,vid["x"]]
            nl = next(i for (i,v) in zip(range(n),x) if v >= -crop)
            nr = n-1-nl
            x = x[nl:nr]
            if not "Gubser" in name:
                x = x/(t-t0)
            yexact = [e(x) for x in x]

            if "Gubser" in name:
                yref = diagonal(ref[:,vid["e"]])
            else:
                yref = ref[:,vid["e"]]
            yref = yref[nl:nr]

            # numericsref ,= axs[1].plot(x,yref, color="gray", label="numerics ref", linestyle="-.", linewidth=2 )

            
            for (scheme,col,(ls1,ls2)) in zip(schemes,plt_setting.clist, zip(lstyles1,lstyles2)):
                linestyle = ls1
                if n == 100:
                    linestyle = ls2
                    # col = lighten(col)
                if dt in datas[scheme]:
                    (sinfo, data, diff) = datas[scheme][dt]
                    vid = sinfo["ID"]
                    cut = sinfo["CUT"]
                    if "FixPoint" in sinfo["integration"]:
                        schemetype = "Implicit"
                    else:
                        schemetype = "Explicit"
                    nbStages = 2
                    if scheme == "GL1":
                        nbStages = 1
                    data = mask(data,vid,cut)
                    if "Gubser" in name:
                        iter = diagonal(data[:,vid["iter"]])
                        y = diagonal(data[:,vid["e"]])
                        yut = diagonal(data[:,vid["ut"]])
                        yux = diagonal(data[:,vid["ux"]])
                    else:
                        iter = data[:,vid["iter"]]
                        y = data[:,vid["e"]]
                        yut = data[:,vid["ut"]]
                        yux = data[:,vid["ux"]]
                    iter = iter[nl:nr]
                    cost = iter*nbStages
                    y = y[nl:nr]
                    yut = yut[nl:nr]
                    yux = yux[nl:nr]
                    yvx = yux/yut
                    yerr = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(y,yexact)]
        
                    numerics ,= axs[0].plot(x,y, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )
                    if "Riemann" in name:
                        axin.plot(x,y, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )
                        # if not "Void" in name:
                        axinv.plot(x,y, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )

                    errexact ,= axs[1].plot(x,yerr, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )
                    # pltcost ,= axs[2].plot(x,cost, '.', label=schemetype)
                    # iterations ,= axs[2].plot(x,iter, '.', label=schemetype)
                    # numericsvx ,= axs[3].plot(x,yvx, label=schemetype, linestyle=linestyle, linewidth=3 )

        if "Riemann" in name:
            axin.set_xlim(-0.75, -0.5)
            axin.set_ylim(8.1, 10.1)
            axs[0].indicate_inset_zoom(axin)
            if not "Void" in name:
                axinv.set_xlim(0.5, 1)
                axinv.set_ylim(0.9, 3.5)
            else:
                axinv.set_xlim(0.8, 1.2)
                axinv.set_ylim(-0.05, 0.1)
            axs[0].indicate_inset_zoom(axinv)
        axs[0].set_ylabel(r"$\epsilon$")
        # axs[0].legend()
        handles, labels = [], []
        hs, ls = axs[0].get_legend_handles_labels()
        def line(col,style):
            return matplotlib.lines.Line2D([],[], color=col,linestyle=style)
        for (h,l) in zip(hs,ls):    # this is the loop to change Labels and colors
            col = h.get_color()
            m = max(col)
            if not (l in labels or m == 1):    # check for Name already exists
                # p = Patch(linestyle="-",color=col)
                p = line(col,"-")
                # h.set_linestyle('-')
                handles += [p]
                labels += [l]
        dxs.sort()
        handles = [line("black",lstyles1[0])]+handles
        labels = [r"$\Delta x = {}$ fm".format(dxs[0])]+labels
        if len(dxs) > 1:
            handles = [line("black",lstyles2[0])] + handles
            labels = [r"$\Delta x = {}$ fm".format(dxs[1])] + labels
        axs[1].legend(handles, labels) #, loc="lower left", bbox_to_anchor=(0, 1.02, 1, 0.2), mode="expand",borderaxespad=0,ncol=3)
        axs[1].set_ylabel(r"$\Delta_\mathrm{exact}$")
        # axs[2].set_ylabel("cost")
        axs[len(axs)-1].set_xlabel("$x/t$")
        plt.subplots_adjust(wspace=0,hspace=0)
        plt.savefig("figures/{}_{}.pdf".format(timename,info2name(info, False)))
        plt.close()


    for scheme in ["GL1"]: 
        fig,axs = plt.subplots(2, 2, figsize=(18,7), sharex=True, sharey=True)
        colors = []
        costmax = 1e-100
        costmin = 1e100
        for j, dxn in zip(range(2), sorted(nds,key=lambda x: x[1])):
            (dx,n) = dxn
            datas = nds[dxn]
            dts = list(filter(lambda dt: abs(dt-dx/10) < 1e-5 or abs(dt-dx/80) < 1e-5, datas[schemes[0]]))
            for i, dt in zip(range(2),sorted(dts, reverse=True)):
                (sinfo, data, diff) = datas[scheme][dt]
                nbStages = 2
                if scheme == "GL1":
                    nbStages = 1
                datats = np.array(sinfo["datats"], dtype=object)
                iter = np.array([d[:,vid["iter"]] for d in datats[:,1]])
                cost = (iter*nbStages).astype(int)
                costmin = min(costmin, cost.min())
                costmax = max(costmax, cost.max())
        for j, dxn in zip(range(2), sorted(nds,key=lambda x: x[1])):
            (dx,n) = dxn
            datas = nds[dxn]
            dts = list(filter(lambda dt: abs(dt-dx/10) < 1e-5 or abs(dt-dx/80) < 1e-5, datas[schemes[0]]))
            for i, dt in zip(range(2),sorted(dts, reverse=True)):
                timename = "dt{:.4e}".format(dt)

                (sinfo, data, diff) = datas[scheme][dt]
                vid = sinfo["ID"]
                if "FixPoint" in sinfo["integration"]:
                    schemetype = "Implicit"
                else:
                    continue
                    schemetype = "Explicit"
                nbStages = 2
                if scheme == "GL1":
                    nbStages = 1

            
                datats = np.array(sinfo["datats"], dtype=object)
                tt = datats[:,0]
                xx = data[:,vid["x"]]
                iter = np.array([d[:,vid["iter"]] for d in datats[:,1]])
                cost = (iter*nbStages).astype(int)
                l = xx[0]
                r = xx[-1]
                d = tt[0]
                u = tt[-1]
            
        
                if "Riemann" in name:
                    if "Void" in name:
                        vs = riem_void_vs
                    else:
                        vs = riem_vs
                    axs[j][i].plot(xx, xx/vs, linestyle="--", color="white", label="shock")
                    axs[j][i].plot(xx, -xx/cs, linestyle="-.", color="white", label="rarefaction")
                im = axs[j][i].imshow(cost, extent=[l,r,d,u], origin="lower", label="aoeu", vmin=costmin, vmax=costmax) #, norm=CenteredNorm(0)) # , cmap="terrain"
                colors += list(np.unique(cost))
                if i == 0:
                    axs[j][i].set_ylabel("$t$ (fm)")
                if j == 1:
                    # axs[j][i].xaxis.tick_top()
                    # axs[j][i].xaxis.set_label_position('top') 
                    axs[j][i].set_xlabel("$x$ (fm)")
                # divider = make_axes_locatable(axs[j][i])
                # cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
                # fig.add_axes(cax)
                # cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
                # cbar.formatter.set_powerlimits((0, 0))
                # cbar.formatter.set_useMathText(True)
                # cbar.update_ticks()
                # cbar.set_label("cost", labelpad=-60)

                # if "Riemann" in name:
                    # axs[j][i].legend(labelcolor="white", facecolor=(0.1,0.1,0.3))
                    # axs[j][i].legend([], [], labelcolor="white", facecolor=(0.1,0.1,0.3))
                axs[j][i].text(0.03, 0.25, r"$\Delta x = "+str(dx)+"$", color="white", transform=axs[j][i].transAxes, fontsize=22)
                pow2 = int(np.log2(0.1*dx/dt))
                if pow2 > 0:
                    pow2 = r"\times 2^{-"+str(pow2)+"}"
                else:
                    pow2 = ""
                axs[j][i].text(0.03, 0.1, r"$\Delta t/\Delta x = 0.1"+pow2+"$", color="white", transform=axs[j][i].transAxes, fontsize=22)
        
        patchs = [Patch(color=im.cmap(im.norm(i)), label=str(i)) for i in np.unique(colors)]
        # axs[0][0].legend(handles=patchs, loc="upper left")
        axs[1][1].legend(handles=patchs, loc="lower right", labelcolor="white", facecolor=(0.3,0.3,0.6))
        # axs[0][0].legend(handles=patchs, loc="upper left", labelcolor="white", facecolor=(0.1,0.1,0.3))
        fig.subplots_adjust(wspace=0.1,hspace=0.1)
        plt.savefig("figures/{}_{}_cost-t_{}.pdf".format(timename,scheme,info2name(sinfo)), dpi=100)
        plt.close()

gref = defaultdict(lambda: None)
greft = defaultdict(lambda: defaultdict(lambda: None))
def plot2d(l, datadts):
    # [dim,name,visc,t0,tend,t,case,dxn,scheme] = l
    [dim,name,visc,t0,tend,t,case,dxn] = l
    if case > 0 and not manycases:
        return
    (dx,n) = dxn
    expl = datadts["Heun"]
    datadts = datadts["GL1"]
    maxdts = sorted([dt for dt in datadts])
    fromref = defaultfromref
    if len(maxdts) == 1:
        fromref = 0
    dt = maxdts[fromref]
    # (info, ref, diffref) = datadts[maxdts[0]]
    if len(expl[dt]) == 0:
        return
    (info, ref, diffref) = expl[dt]
    vid = info["ID"]
    cut = info["CUT"]
    ref = mask(ref,vid,cut)
    (einfo, data, diff) = datadts[dt]
    vid = einfo["ID"]
    cut = einfo["CUT"]

    name = einfo["name"]
    case = einfo["case"]
    if (("Trento" in name and case == 0) or "Gubser" in name):
        datats = np.array(einfo["datats"], dtype=object)
        refdatats = np.array(info["datats"], dtype=object)
    else:
        datats = [(t,data,diff)]
        refdatats = [(t,ref,refdiff)]

    ld = len(datats)
    many = ld > 1
    nid = ceil(log(ld)/log(10))

    if many:
        num = num2D
        nums = np.array([i for i in range(ld) if i%int(ld/(num-1)) == 0]+[ld-1])
        # nb = 5
        nb = 2
        if "Gubser" in einfo["name"]:
            nb += 1
        if "FixPoint" in einfo["integration"]:
            nb += 1
        fig, axs = plt.subplots(nb,num, figsize=(2+num*4, nb*5)) #sharey = True, sharex=True)
        if not hasattr(axs[0], "__len__"):
            axs = [axs]
        colors = []
        ma = 1e-100
        mi = 1e100
        zitermax = 1e-100
        zitermin = 1e100
        for (id, (t,data,diff)) in zip(range(num),datats): #[nums]):
            mdata = mask(data,vid,cut)
            n = einfo["nx"]
            x = mdata[:,vid["x"]]
            y = mdata[:,vid["y"]]
            z = mdata[:,vid["e"]]
            ziter = mdata[:,vid["iter"]]
            nl = next(i for (i,v) in zip(range(n*n),x) if v >= -crop)
            nr = n-1-nl
            x = np.reshape(x, (n,n))[nl:nr,nl:nr].reshape((-1))
            y = np.reshape(y, (n,n))[nl:nr,nl:nr].reshape((-1))
            z = np.reshape(z, (n,n))[nl:nr,nl:nr].reshape((-1))
            ziter = np.reshape(ziter, (n,n))[nl:nr,nl:nr]
            zgubser = np.array([gubser(x,y,t) for (x,y) in zip(x,y)])
            zerr = np.array([(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zgubser)])
            mi = min(mi,zerr.min())
            ma = max(ma,zerr.max())
            zitermin = min(zitermin, ziter.min())
            zitermax = max(zitermax, ziter.max())
        for (id, (t,data,diff), (_,ref,refdiff)) in zip(range(num),datats[nums],refdatats[nums]):
            # global greft
            # if greft[case][t] is None:
            #     greft[case][t] = (data,scheme)
            # elif greft[case][t][1] != scheme:
            #     ref = greft[case][t][0]
            mdata = mask(data,vid,cut)
            mref = mask(ref,vid,cut)
            n = einfo["nx"]
            x = mdata[:,vid["x"]]
            y = mdata[:,vid["y"]]
            z = mdata[:,vid["e"]]
            zref = mref[:,vid["e"]]
            ziter = mdata[:,vid["iter"]]
            # print(name, case, scheme, dt, n, ziter.sum())
            # zpi00 = mdata[:,vid["pi00"]]
            # zpi11 = mdata[:,vid["pi11"]]
            # zpi12 = mdata[:,vid["pi12"]]
            # zpi22 = mdata[:,vid["pi22"]]
            zut = mdata[:,vid["ut"]]
            zux = mdata[:,vid["ux"]]
            zvx = zux/zut
            sgn = np.sign(zvx)
            zvx = np.power(zvx*sgn,0.5)*sgn
            zuy = mdata[:,vid["uy"]]
            zvy = zuy/zut
            sgn = np.sign(zvy)
            zvy = np.power(zvy*sgn,0.5)*sgn
            zgubser = [gubser(x,y,t) for (x,y) in zip(x,y)]
            zerr = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zgubser)]
            # zerrref = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zref)]
            zerrref = [(a-b) for (a,b) in zip(z,zref)]
            sgn = np.sign(zerrref)
            # zerrref = np.power(zerrref*sgn,0.5)*sgn
            nl = next(i for (i,v) in zip(range(n*n),x) if v >= -crop)
            nr = n-1-nl
            x = np.reshape(x, (n,n))[nl:nr,nl:nr]
            y = np.reshape(y, (n,n))[nl:nr,nl:nr]
            z = np.reshape(z, (n,n))[nl:nr,nl:nr]
            zref = np.reshape(zref, (n,n))[nl:nr,nl:nr]
            ziter = np.reshape(ziter, (n,n))[nl:nr,nl:nr]
            zgubser = np.reshape(zgubser, (n,n))[nl:nr,nl:nr]
            zerr = np.reshape(zerr, (n,n))[nl:nr,nl:nr]
            zerrref = np.reshape(zerrref, (n,n))[nl:nr,nl:nr]
            zvx = np.reshape(zvx, (n,n))[nl:nr,nl:nr]
            zvy = np.reshape(zvy, (n,n))[nl:nr,nl:nr]
            zux = np.reshape(zux, (n,n))[nl:nr,nl:nr]
            zuy = np.reshape(zuy, (n,n))[nl:nr,nl:nr]
            zut = np.reshape(zut, (n,n))[nl:nr,nl:nr]
            zur = np.sqrt(zux*zux+zuy*zuy)
            # zpi00 = np.reshape(zpi00, (n,n))[nl:nr,nl:nr]
            # zpi11 = np.reshape(zpi11, (n,n))[nl:nr,nl:nr]
            # zpi12 = np.reshape(zpi12, (n,n))[nl:nr,nl:nr]
            # zpi22 = np.reshape(zpi22, (n,n))[nl:nr,nl:nr]
            l = x[0][0]
            r = x[0][-1]
            d = y[0][0]
            u = y[-1][0]
            # all = [(r"$\epsilon$", z),("vy",zvy), ("err ref",zerrref), ("vx", zvx)]
            # all = [(r"$\epsilon$", z), ("ur", zur), ("ut", zut)]
            # all = [(r"$\epsilon$", z),(r"pi00", zpi00),(r"pi11", zpi11),(r"pi12", zpi12),(r"pi22", zpi22)]
            # all = [(r"$\epsilon$", z), ("ref", zref), ("err ref",zerrref)]
            all = [(r"$\epsilon$", z), (r"$\Delta_\mathrm{Impl-Expl}$",zerrref)]
            # all = [(r"$\epsilon$", z)]
            if "Gubser" in einfo["name"]:
                all += [(r"$\Delta_\mathrm{exact}$", zerr)]
            if "FixPoint" in einfo["integration"]:
                all += [(r"$\mathrm{iterations}$", ziter)]
            for (i, (n, z)) in zip(range(nb),all):
                if "iterations" in n:
                    imit = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower", vmin=zitermin, vmax=zitermax)
                    colors += list(np.unique(z))
                elif n == r"$\Delta_\mathrm{exact}$":
                    im = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower", vmin=mi, vmax=ma)
                else:
                    im = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower")

                axs[i][id].yaxis.tick_right()
                if i == 0:
                    axs[i][id].xaxis.set_label_position('top') 
                    axs[i][id].set_xlabel(r"$\tau$ = {} (fm)".format(t))
                elif i == nb-1:
                    axs[i][id].set_xlabel("$x$ (fm)")
                if id == num-1:
                    axs[i][id].yaxis.set_label_position('right')
                    axs[i][id].set_ylabel("$y$ (fm)")
                if id == 0:
                    axs[i][id].set_ylabel(n, fontsize = 36)
                if i<nb-1:
                    axs[i][id].tick_params(axis='x', which='both', labelbottom=False)
                if id<num-1:
                    axs[i][id].tick_params(axis='y', which='both', labelright=False)

                if not "iterations" in n:
                    divider = make_axes_locatable(axs[i][id])
                    p = list(divider.get_position())
                    p[1] -= p[3]*0.09
                    p[3] *= 0.05
                    cax = fig.add_axes(p)
                    # cax = divider.new_vertical(size="5%", pad=0.2, pack_start=True)
                    fig.add_axes(cax)
                    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
                    if "Gubser" in einfo["name"] and "epsilon" in n:
                        # cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation=20)
                        cbar.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=2))
                    # cbar.formatter.set_powerlimits((0, 0))
                    # cbar.formatter.set_useMathText(True)
                    # cbar.update_ticks()
                    # cbar.set_label(n, labelpad=-60)

        ax = axs[nb-1][0]
        patchs = [Patch(color=imit.cmap(imit.norm(i)), label=str(int(i))) for i in np.unique(colors)]
        ax.legend(handles=patchs, loc="upper left", labelcolor="white", facecolor=(0.3,0.3,0.6))
        # ax.text(0.075, 0.1, r"$\Delta x = "+str(dx)+"$ fm", color="white", transform=ax.transAxes, fontsize=22)
        # plt.title(r"$\Delta x = "+str(dx)+"$ fm")
        # plt.subplots_adjust(wspace=0)
        plt.savefig("figures/many_best_e_{}.pdf".format(info2name(einfo)), dpi=100)
        plt.close()
    
    if not animate:
        datats = [datats[-1]]

    for (id, (t,data,diff)) in zip(range(1000000), datats):
        # global greft
        # if greft[case][t] is None:
        #     greft[case][t] = (data,scheme)
        # elif greft[case][t][1] != scheme:
        #     ref = greft[case][t][0]
        mdata = mask(data,vid,cut)
        n = einfo["nx"]
        x = mdata[:,vid["x"]]
        y = mdata[:,vid["y"]]
        z = mdata[:,vid["e"]]
        zref = ref[:,vid["e"]]
        ziter = mdata[:,vid["iter"]]
        zut = mdata[:,vid["ut"]]
        zux = mdata[:,vid["ux"]]
        zvx = zux/zut
        sgn = np.sign(zvx)
        zvx = np.power(zvx*sgn,0.5)*sgn
        zuy = mdata[:,vid["uy"]]
        zvy = zuy/zut
        sgn = np.sign(zvy)
        zvy = np.power(zvy*sgn,0.5)*sgn
        zgubser = [gubser(x,y,t) for (x,y) in zip(x,y)]
        zerr = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zgubser)]
        zerrref = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(z,zref)]
        sgn = np.sign(zerrref)
        zerrref = np.power(zerrref*sgn,0.5)*sgn
        nl = next(i for (i,v) in zip(range(n*n),x) if v >= -crop)
        nr = n-1-nl
        x = np.reshape(x, (n,n))[nl:nr,nl:nr]
        y = np.reshape(y, (n,n))[nl:nr,nl:nr]
        z = np.reshape(z, (n,n))[nl:nr,nl:nr]
        ziter = np.reshape(ziter, (n,n))[nl:nr,nl:nr]
        zgubser = np.reshape(zgubser, (n,n))[nl:nr,nl:nr]
        zerr = np.reshape(zerr, (n,n))[nl:nr,nl:nr]
        zerrref = np.reshape(zerrref, (n,n))[nl:nr,nl:nr]
        zvx = np.reshape(zvx, (n,n))[nl:nr,nl:nr]
        zvy = np.reshape(zvy, (n,n))[nl:nr,nl:nr]
        l = x[0][0]
        r = x[0][-1]
        d = y[0][0]
        u = y[-1][0]
        # all = [("e", z), ("vy", zvy), ("err ref", zerrref),("vx",zvx)]
        all = [("e", z)]
        # all = [("vx", zvx), ("e", z), ("err ref", zerrref)]
        if "Gubser" in einfo["name"] and not "Exponential" in einfo["name"]:
            all += [("err", zerr)]
        if "FixPoint" in einfo["integration"]:
            all += [("iter", ziter)]
        nb = len(all)
        fig, axs = plt.subplots(1,nb, figsize=(2+nb*4, 5), sharey=True)
        if not hasattr(axs, "__len__"):
            axs = [axs]
        for (i, (n, z)) in zip(range(nb),all):
            im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower")
            # if n == "iter":
            #     im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower", vmin=0, vmax=3)
            # else:
            #     im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower") #, norm=CenteredNorm(0)) # , cmap="terrain"
            axs[i].set_xlabel("$x$ (fm)")
            axs[i].xaxis.tick_top()
            axs[i].xaxis.set_label_position('top') 
            if i == 0:
                axs[i].set_ylabel("$y$ (fm)")
            divider = make_axes_locatable(axs[i])
            cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
            fig.add_axes(cax)
            cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
            cbar.formatter.set_powerlimits((0, 0))
            cbar.formatter.set_useMathText(True)
            cbar.update_ticks()
            cbar.set_label(n, labelpad=-60)
            ax = axs[0]
            ax.text(0.075, 0.1, r"$\Delta x = "+str(dx)+"$ fm", color="white", #, bbox={"facecolor": "white", "pad": 10},
                transform=ax.transAxes, fontsize=22)
        if many and animate:
            figname = "figures/best_e_{}".format(info2name(einfo))
            try:
                os.mkdir(figname)
            except FileExistsError:
                None
            plt.savefig(("{}/{:0>"+str(nid)+"}.pdf").format(figname, id), dpi=100)
        else:
            plt.savefig("figures/best_e_{}.pdf".format(info2name(einfo)), dpi=100)
        plt.close()

        if not diff is None:
            n = einfo["nx"]
            dt00 = np.reshape(diff[:,2], (n,n))
            dt01 = np.reshape(diff[:,3], (n,n))
            dt02 = np.reshape(diff[:,4], (n,n))
            totdt00 = dt00.sum()/abs(data[:,2]).sum()
            totdt01 = dt01.sum()/abs(data[:,3]).sum()
            totdt02 = dt02.sum()/abs(data[:,4]).sum()
            all = [("diff T00",dt00,totdt00),("diff T01",dt01,totdt01),("diff T02",dt02,totdt02)]
            nb = len(all)
            fig, axs = plt.subplots(1,nb, figsize=(2+nb*4, 5), sharey=True)
            if not hasattr(axs, "__len__"):
                axs = [axs]
            for (i, (n, z, tot)) in zip(range(nb),all):
                im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower") #, norm=CenteredNorm(0)) # , cmap="terrain"
                axs[i].set_xlabel("$x$ (fm)")
                axs[i].xaxis.tick_top()
                axs[i].xaxis.set_label_position('top') 
                axs[i].title.set_text("relative: {:.2e}".format(tot))
                if i == 0:
                    axs[i].set_ylabel("$y$ (fm)")
                divider = make_axes_locatable(axs[i])
                cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
                fig.add_axes(cax)
                cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
                cbar.formatter.set_powerlimits((0, 0))
                cbar.formatter.set_useMathText(True)
                cbar.update_ticks()
                cbar.set_label(n, labelpad=-60)
            if many and animate:
                figname = "figures/best_diff_{}".format(info2name(info))
                try:
                    os.mkdir(figname)
                except FileExistsError:
                    None
                plt.savefig(("{}/{:0>"+str(nid)+"}.pdf").format(figname, id), dpi=100)
            else:
                plt.savefig("figures/best_diff_{}.pdf".format(info2name(info)), dpi=100)
            plt.close()

def plotall1D(l, d):
    if l[0] == "1D" or "Gubser" in l[1]:
        plot1d(l, d)
def plotall2D(l, d):
    if l[0] == "2D":
        plot2d(l, d)

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
    
if plot_conv:
    alldata(6, datas, convall)
if plot_1D:
    alldata(7, datas, plotall1D)
if plot_2D:
    alldata(8, datas, plotall2D)
    # alldata(9, datas, plotall2D)
