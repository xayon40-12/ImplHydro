#!/usr/bin/env python3
import sys 
import os
import warnings
import matplotlib
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
from scipy import signal
from math import ceil,log

warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning) # disable matplotlib deprecation warning
np.seterr(invalid='ignore') # disable invalid waring for numpy as NaN are used to discard data in the void
np.seterr(divide='ignore') # disable divide by zero worining
# plt.rcParams['axes.grid'] = False


animate = False
animate_arg = "-animate="
viscous = False
viscous_arg = "-viscous="
for arg in sys.argv:
    if animate_arg in arg:
        animate = bool(arg[len(animate_arg):])
    if viscous_arg in arg:
        viscous = bool(arg[len(viscous_arg):])

IDx = 0
IDy = 1
IDiter = 2
IDt00 = 3
ID1De = 5
ID1Dut = 8
ID1Dux = 9
ID2De = 6
ID2Dut = 9
ID2Dux = 10
ID2Duy = 11
if viscous:
    ID2De  += 7
    ID2Dut += 7 
    ID2Dux += 7 
    ID2Duy += 7 

CUT = 1e-9

crop = 9
fromref = 1

e0 = 10
emin = 1
cs2 = 1/3
eps = 1e-10

riem_e, riem_v = riemann(e0,emin,cs2,eps,False)
riem_void_e, riem_void_v = riemann(e0,emin,cs2,eps,True)

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

maxvoidratio = 0
meanvoidratio = 0
countvoidratio = 0
dir = "results/"
for d in os.listdir(dir):
    dird = dir+d
    ts = sorted(os.listdir(dird), key=float) 
    p = dird+"/"+ts[-1]
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
    def find_diff(dir):
        try:
            return np.fromfile(dir+"/diff.dat", dtype="float64").reshape((n,-1))
        except FileNotFoundError:
            return None
    diff = find_diff(p)
    if ("Trento" in name and case == 0) or "Gubser" in name:
        datats = [(float(t), np.fromfile(dird+"/"+t+"/data.dat", dtype="float64").reshape((n,-1)),find_diff(dird+"/"+t)) for t in ts]
        info["datats"] = datats

    t00 = data[:,IDt00]
    m = t00>CUT
    err = abs(t00.sum()-t00[m].sum())/t00.sum()
    meanvoidratio += err
    maxvoidratio = max(maxvoidratio, err)
    countvoidratio += 1

    # print(p, dim, integration, t0, tend, dx, nx, maxdt)
    datas[dim][name][t0][tend][dx][nx][t][case][scheme][maxdt] = (info, data, diff)

meanvoidratio /= countvoidratio
print("maxvoidratio: ", maxvoidratio, "meanvoidratio: ", meanvoidratio)
with open("voidratio.txt", "w") as fv:
    fv.write("max_void_ratio: {:e}\nmean_void_ratio: {:e}".format(maxvoidratio, meanvoidratio))
print("finished loading")
# sys.exit(0)

def compare(i, vss, wss):
    maxerr = 0
    meanerr = 0
    count = 0
    for (vs, ws) in zip(vss,wss):
        if  vs[IDt00] > CUT and vs[0] >= -crop and vs[0] <= crop and vs[1] >= -crop and vs[1] <= crop:
            count += 1
            a = vs[i]
            b = ws[i]
            err = abs(a-b)/max(abs(a),abs(b))
            maxerr = max(err,maxerr)
            meanerr += err
    if count>0:
        meanerr /= count
    return (maxerr, meanerr)

def convergence(a, ref=None):
    maxdts = sorted([dt for dt in a])
    if ref is None:
        ref = a[maxdts[0]][1]
    all = []
    for i in maxdts:
        (info, v, diff) = a[i]
        cost = info["cost"]
        dt = info["maxdt"]
        avdt = (info["tend"]-info["t0"])/info["tsteps"]
        elapsed = info["elapsed"]
        if info["dim"] == "1D":
            id = ID1De
        else:
            id = ID2De
        (maxerr, meanerr) = compare(id, ref, v)
        all += [(v, info, maxerr, meanerr, dt, cost, avdt, elapsed)]
    
    return np.array(all, dtype=object)

def info2name(info, scheme=True):
    if scheme:
        return "{}_{}{}_{}_{}_{}_{}_{}_{}".format(info["dim"],info["name"],info["case"],info["scheme"],info["t0"],info["tend"],info["dx"],info["nx"],info["t"])
    else:
        return "{}_{}{}_{}_{}_{}_{}_{}".format(info["dim"],info["name"],info["case"],info["t0"],info["tend"],info["dx"],info["nx"],info["t"])

def convall(l, ds):
    [dim,name,t0,tend,dx,nx,t] = l
    # allx = [("dt", 4), ("cost", 5), ("avdt", 6), ("elapsed", 7)]
    ally = [("max", 2), ("mean", 3)]
    allx = [("cost", 5)]
    # ally = [("max", 2)]
    for (dtcost, dci) in allx:
        for (meanmax, mmi) in ally:
            plt.rcParams["figure.figsize"] = [8, 5]
            plt.figure()
            plt.xlabel(dtcost)
            plt.ylabel(meanmax+" error")
            # plt.title("{} {} t0={} tend={} dx={} cells={}".format(dim, name, t0, tend, dx, nx))
            d = ds[list(ds)[0]]
            scs = sorted(list(d.keys()))
            dts = sorted([dt for dt in d[scs[0]]])
            mindt = dts[0]
            scs.sort(key=lambda s: integrationPriority(d[s][mindt][0]["integration"]))
            if len(scs) <= 1:
                plt.close()
                return 
            s1 = scs[0]
            for (s0,col) in zip(scs,plt_setting.clist):
                for case in ds: 
                    d = ds[case]
                    ds0 = d[s0]
                    dtref = sorted(list(ds0.keys()))[0]
                    info = ds0[dtref][0]
                    refs = {s: d[s][sorted(list(d[s].keys()))[0]][1] for s in d}
                    if info["integration"] == "FixPoint":
                        schemetype = "Implicit"
                    else:
                        schemetype = "Explicit"
                    c = convergence(d[s0],refs[s1])
                    plt.loglog(c[fromref:,dci],c[fromref:,mmi], 'o', label=schemetype, color=col, linestyle="-.", linewidth=1, alpha=0.5)
            labels = []
            for p in plt.gca().get_lines():    # this is the loop to change Labels and colors
                label = p.get_label()
                if label in labels:    # check for Name already exists
                    p.set_label('_' + label)       # hide label in auto-legend
                else:
                    labels += [label]
            plt.legend()
            plt.savefig("figures/convergence_{}_{}_crop:{}_{}.pdf".format(dtcost, meanmax, crop, info2name(info, False)))
            plt.close()

def integrationPriority(integration):
    if integration == "Explicit":
        return 1
    else:
        return 2

# def mask(m,data):
#     return np.concatenate((data[:,:3], data[:,3:]/m[:,None]), axis=1)
def mask(data):
    m = (data[:,IDt00]>CUT).astype(int)
    return np.concatenate((data[:,:3], data[:,3:]/m[:,None]), axis=1)

def plot1d(l, datas):
    [dim,name,t0,tend,dx,n,t,case] = l
    schemes = sorted(list(datas.keys()))
    dts = sorted([dt for dt in datas[schemes[0]]])
    dts = [dts[0],dts[-1]] # only keep best and worst
    mindt = dts[0]
    schemes.sort(key=lambda s: integrationPriority(datas[s][mindt][0]["integration"]))
    (info,ref,diffref) = datas["Heun"][mindt]
    ref = mask(ref)
    nl = 2
    # lstyles = [(nl*i,(nl,(nl-1)*nl)) for i in range(nl)]
    lstyles = [(3*i,(3,5,1,5)) for i in range(nl)]
    
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

    x = ref[:,IDx]
    nl = next(i for (i,v) in zip(range(n),x) if v >= -crop)
    nr = n-1-nl
    x = x[nl:nr]
    if not "Gubser" in name:
        x = x/(t-t0)
    ycontinuum = [e(x) for x in x]
    if "Gubser" in name:
        yref = diagonal(ref[:,ID2De])
    else:
        yref = ref[:,ID1De]
    yref = yref[nl:nr]

    for dt in dts:
        timename = "dt{}".format(dt)
        if dt == mindt:
            timename = "best_"+timename

    
        # plt.rcParams["figure.figsize"] = [8, 12]
        # _,axs = plt.subplots(4, 1, sharex=True)
        plt.rcParams["figure.figsize"] = [8, 9]
        _,axs = plt.subplots(3, 1, sharex=True)
        continuum ,= axs[1].plot(x,ycontinuum, color="grey", label="continuum", linewidth=2)
        # numericsref ,= axs[1].plot(x,yref, color="gray", label="numerics ref", linestyle="-.", linewidth=2 )

        for (scheme,linestyle) in zip(schemes, lstyles):
            (sinfo, data, diff) = datas[scheme][dt]
            data = mask(data)
            if "Gubser" in name:
                iter = diagonal(data[:,IDiter])
                y = diagonal(data[:,ID2De])
                yut = diagonal(data[:,ID2Dut])
                yux = diagonal(data[:,ID2Dux])
            else:
                iter = data[:,IDiter]
                y = data[:,ID1De]
                yut = data[:,ID1Dut]
                yux = data[:,ID1Dux]
            iter = iter[nl:nr]
            y = y[nl:nr]
            yut = yut[nl:nr]
            yux = yux[nl:nr]
            yvx = yux/yut
            yerr = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(y,ycontinuum)]
        
            if sinfo["integration"] == "FixPoint":
                schemetype = "Implicit"
            else:
                schemetype = "Explicit"
            iterations ,= axs[0].plot(x,iter, '.', label=schemetype)
            numerics ,= axs[1].plot(x,y, label=schemetype, linestyle=linestyle, linewidth=3 )
            errcontinuum ,= axs[2].plot(x,yerr, label=schemetype, linestyle=linestyle, linewidth=3 )
            # numericsvx ,= axs[3].plot(x,yvx, label=schemetype, linestyle=linestyle, linewidth=3 )

        axs[0].set_ylabel("iterations")
        # axs[0].legend()
        axs[1].set_ylabel("e")
        axs[1].legend()
        axs[2].set_ylabel("continuum err")
        axs[2].set_xlabel("x/t")
        # axs[2].legend()
        # axs[3].set_ylabel("vx")
        # axs[3].set_xlabel("x/t")
        # axs[3].legend()
        plt.savefig("figures/{}_{}.pdf".format(timename,info2name(info, False)))
        plt.close()

def plot2d(l, datadts):
    [dim,name,t0,tend,dx,n,t,case,scheme] = l
    maxdts = sorted([dt for dt in datadts])
    (info, ref, diffref) = datadts[maxdts[0]]
    ref = mask(ref)
    if len(maxdts) == 1:
        fromref = 0
    (einfo, data, diff) = datadts[maxdts[fromref]]

    name = info["name"]
    case = info["case"]
    if (("Trento" in name and case == 0) or "Gubser" in name):
        datats = np.array(einfo["datats"], dtype=object)
    else:
        datats = [(t,data,diff)]

    ld = len(datats)
    many = ld > 1
    nid = ceil(log(ld)/log(10))

    if many:
        num = 5
        nums = np.array([i for i in range(ld) if i%(ld/(num-1)) == 0]+[ld-1])
        nb = 1
        if "Gubser" in info["name"]:
            nb += 1
        if info["integration"] == "FixPoint":
            nb += 1
        plt.rcParams["figure.figsize"] = [2+num*4, nb*5]
        fig, axs = plt.subplots(nb,num, sharey=True) #, sharex=True)
        if not hasattr(axs[0], "__len__"):
            axs = [axs]
        for (id, (t,data,diff)) in zip(range(num),datats[nums]):
            mdata = mask(data)
            n = info["nx"]
            x = mdata[:,IDx]
            y = mdata[:,IDy]
            z = mdata[:,ID2De]
            zref = ref[:,ID2De]
            ziter = mdata[:,IDiter]
            zut = mdata[:,ID2Dut]
            zux = mdata[:,ID2Dux]
            zvx = zux/zut
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
            zvx = np.reshape(zvx, (n,n))[nl:nr,nl:nr]
            l = x[0][0]
            r = x[0][-1]
            d = y[0][0]
            u = y[-1][0]
            all = [("e", z)]
            if "Gubser" in info["name"]:
                all += [("err", zerr)]
            if info["integration"] == "FixPoint":
                all += [("iter", ziter)]
            for (i, (n, z)) in zip(range(nb),all):
                if n == "iter":
                    im = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower", vmin=0, vmax=3)
                else:
                    im = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower") #, norm=CenteredNorm(0)) # , cmap="terrain"
                axs[i][id].xaxis.tick_top()
                axs[i][id].xaxis.set_label_position('top') 
                if i == 0:
                    axs[i][id].set_xlabel("x")
                if i>0:
                    axs[i][id].tick_params(axis='x', which='both', labelbottom=False, labeltop=False)
                if id == 0:
                    axs[i][id].set_ylabel("y")
                divider = make_axes_locatable(axs[i][id])
                cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
                fig.add_axes(cax)
                cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
                cbar.formatter.set_powerlimits((0, 0))
                cbar.formatter.set_useMathText(True)
                cbar.update_ticks()
                cbar.set_label("{} (t = {:.2} fm)".format(n, t), labelpad=-60)

        plt.savefig("figures/many_best_e_{}.pdf".format(info2name(info)), dpi=100)
        plt.close()
    
    if not animate:
        datats = [datats[-1]]

    for (id, (t,data,diff)) in zip(range(1000000), datats):
        mdata = mask(data)
        n = info["nx"]
        x = mdata[:,IDx]
        y = mdata[:,IDy]
        z = mdata[:,ID2De]
        zref = ref[:,ID2De]
        ziter = mdata[:,IDiter]
        zut = mdata[:,ID2Dut]
        zux = mdata[:,ID2Dux]
        zvx = zux/zut
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
        zvx = np.reshape(zvx, (n,n))[nl:nr,nl:nr]
        l = x[0][0]
        r = x[0][-1]
        d = y[0][0]
        u = y[-1][0]
        # all = [("vx", zvx), ("e", z), ("err ref", zerrref)]
        all = [("e", z)]
        # all = [("vx", zvx), ("e", z), ("err ref", zerrref)]
        if "Gubser" in info["name"] and not "Exponential" in info["name"]:
            all += [("err", zerr)]
        if info["integration"] == "FixPoint":
            all += [("iter", ziter)]
        nb = len(all)
        plt.rcParams["figure.figsize"] = [2+nb*4, 5]
        fig, axs = plt.subplots(1,nb, sharey=True)
        if not hasattr(axs, "__len__"):
            axs = [axs]
        for (i, (n, z)) in zip(range(nb),all):
            if n == "iter":
                im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower", vmin=0, vmax=3)
            else:
                im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower") #, norm=CenteredNorm(0)) # , cmap="terrain"
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
        if many and animate:
            figname = "figures/best_e_{}".format(info2name(info))
            try:
                os.mkdir(figname)
            except FileExistsError:
                None
            plt.savefig(("{}/{:0>"+str(nid)+"}.pdf").format(figname, id), dpi=100)
        else:
            plt.savefig("figures/best_e_{}.pdf".format(info2name(info)), dpi=100)
        plt.close()

        if not diff is None:
            n = info["nx"]
            dt00 = np.reshape(diff[:,2], (n,n))
            dt01 = np.reshape(diff[:,3], (n,n))
            dt02 = np.reshape(diff[:,4], (n,n))
            totdt00 = dt00.sum()/abs(data[:,IDt00]).sum()
            totdt01 = dt01.sum()/abs(data[:,IDt00+1]).sum()
            totdt02 = dt02.sum()/abs(data[:,IDt00+2]).sum()
            all = [("diff T00",dt00,totdt00),("diff T01",dt01,totdt01),("diff T02",dt02,totdt02)]
            nb = len(all)
            plt.rcParams["figure.figsize"] = [2+nb*4, 5]
            fig, axs = plt.subplots(1,nb, sharey=True)
            if not hasattr(axs, "__len__"):
                axs = [axs]
            for (i, (n, z, tot)) in zip(range(nb),all):
                im = axs[i].imshow(z, extent=[l,r,d,u], origin="lower") #, norm=CenteredNorm(0)) # , cmap="terrain"
                axs[i].set_xlabel("x")
                axs[i].xaxis.tick_top()
                axs[i].xaxis.set_label_position('top') 
                axs[i].title.set_text("relative: {:.2e}".format(tot))
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
    
alldata(7, datas, convall)
alldata(8, datas, plotall1D)
alldata(9, datas, plotall2D)
