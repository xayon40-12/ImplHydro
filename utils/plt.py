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
sys.path.append(os.path.abspath("%s/.config/" % (os.path.expanduser("~"))))
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
# plt.rcParams['axes.grid'] = False


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

argActions = [(["-a","--animate"], setAnimate),(["-m","--manycases"], setManyCases), (["-r","--rejectfails"], setRejectFails)]
for arg in sys.argv:
    for (larg, act) in argActions:
        if arg in larg:
            act()

CUT = 1e-6

crop = 9
defaultfromref = 1

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
    maxt = ts[-1]

    p = dird+"/"+maxt
    info = {k: convert(v.strip()) for [k, v] in np.loadtxt(p+"/info.txt", dtype=object, delimiter=":")}

    t0 = info["t0"]
    tend = info["tend"]
    scheme = info["scheme"]
    maxdt = info["maxdt"]
    dx = info["dx"]
    nx = info["nx"]
    t = info["t"]
    fails = info["fails"]
    r = 1e10
    t = float(round(r*t)/r)
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
    info["variables"] = info["variables"].split(" ")
    vid = {n: i for (i,n) in enumerate(info["variables"])}
    info["ID"] = vid
    visc = info["viscosity"]
    cut = CUT
    match visc:
        case "Ideal":
            visc = ()
        case "Shear":
            visc = (info["etaovers"])
            cut = info["energycut"]
        case "Bulk":
            visc = (info["zeta"])
            cut = info["energycut"]
        case "Both":
            visc = (info["zeta"],info["etaovers"])
            cut = info["energycut"]
    info["visc"] = visc
    cut = CUT # do not consider energycut when viscosity is disabled
    info["CUT"] = cut

    # data = np.loadtxt(p+"/data.txt")
    data0 = np.fromfile(dird+"/"+ts[0]+"/data.dat", dtype="float64").reshape((n,-1))
    data = np.fromfile(p+"/data.dat", dtype="float64").reshape((n,-1))
    def find_diff(dir):
        try:
            return np.fromfile(dir+"/diff.dat", dtype="float64").reshape((n,-1))
        except FileNotFoundError:
            return None
    diff = find_diff(p)
    if case == 0:
        datats = [(float(t), np.fromfile(dird+"/"+t+"/data.dat", dtype="float64").reshape((n,-1)),find_diff(dird+"/"+t)) for t in ts]
        info["datats"] = datats

    e0 = data0[:,vid["e"]]
    e = data[:,vid["e"]]
    # if the total energy sum(e) is bigger at the end time than the initial time, 
    # it means that explicit failed
    expl_fail = sum(e) > sum(e0)
    # if a fail (decreasing dt) happens in implicit when dt>dx/10, we consider that implicit failed as we want to test large dt and thus do not want dt to be decreased
    impl_fail = fails > 0 and maxdt>dx/10
    if rejectfails and (expl_fail or impl_fail): 
        continue # skip explicit or implicit that failed
    m = e>cut
    err = abs(e.sum()-e[m].sum())/e.sum()
    meanvoidratio += err
    maxvoidratio = max(maxvoidratio, err)
    countvoidratio += 1

    # print(p, dim, t0, tend, dx, nx, maxdt)
    datas[dim][name][visc][t0][tend][t][case][(dx,nx)][scheme][maxdt] = (info, data, diff)

meanvoidratio /= countvoidratio
print("maxvoidratio: ", maxvoidratio, "meanvoidratio: ", meanvoidratio)
with open("voidratio.txt", "w") as fv:
    fv.write("max_void_ratio: {:e}\nmean_void_ratio: {:e}".format(maxvoidratio, meanvoidratio))
print("finished loading")
# sys.exit(0)

def compare(i, cut, vss, wss):
    maxerr = 0
    meanerr = 0
    count = 0
    for (vs, ws) in zip(vss,wss):
        if  vs[i] > cut and vs[0] >= -crop and vs[0] <= crop and vs[1] >= -crop and vs[1] <= crop:
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
    maxdts = sorted(list(a))
    if ref is None:
        ref = a[maxdts[0]][1]
    all = []
    for i in maxdts:
        (info, v, diff) = a[i]
        cost = info["cost"]
        dt = info["maxdt"]
        avdt = (info["tend"]-info["t0"])/info["tsteps"]
        elapsed = info["elapsed"]
        vid = info["ID"]
        id = vid["e"]
        cut = info["CUT"]
        (maxerr, meanerr) = compare(id, cut, ref, v)
        all += [(v, info, maxerr, meanerr, dt, cost, avdt, elapsed)]
    
    return np.array(all, dtype=object)

def info2name(info, scheme=True):
    visc = info["viscosity"]
    match visc:
        case "Ideal":
            visc = "Ideal"
        case "Shear":
            visc = "Shear({})".format(info["etaovers"])
        case "Bulk":
            visc = "Bulk({})".format(info["zeta"])
        case "Both":
            visc = "Both({},{})".format(info["zeta"],info["etaovers"])
    if scheme:
        return "{}_{}{}_{}_{}_{}_{}_{}_{}_{:.4e}".format(info["dim"],info["name"],info["case"],visc,info["scheme"],info["t0"],info["tend"],info["dx"],info["nx"],info["t"])
    else:
        return "{}_{}{}_{}_{}_{}_{}_{}_{:.4e}".format(info["dim"],info["name"],info["case"],visc,info["t0"],info["tend"],info["dx"],info["nx"],info["t"])

def lighten(c):
    return tuple(min(1,(a+b*0.2)*2) for a, b in zip(c,np.roll(c,1)))

def convall(l, cnds):
    [dim,name,visc,t0,tend,t] = l
    # allx = [("dt", 4), ("cost", 5), ("avdt", 6), ("elapsed", 7)]
    ally = [((1,0), "full", 2), ((5,5), "none", 3)]
    allx = [("cost", 5)]
    # ally = [("max", 2)]
    fromref = defaultfromref
    for (dtcost, dci) in allx:
        plt.rcParams["figure.figsize"] = [8, 5]
        fig, ax = plt.subplots()
        fig2, ax2 = plt.subplots()
        def setlabel(ax):
            ax.set_xlabel(dtcost)
            ax.set_ylabel(r"$\Delta\epsilon$")
        setlabel(ax)
        setlabel(ax2)
        # plt.title("{} {} t0={} tend={} dx={} cells={}".format(dim, name, t0, tend, dx, nx))
        nbcases = len(cnds)
        alpha = sqrt(1/nbcases)
        for case in cnds:
            nds = cnds[case]
            dxs = sorted(nds,key=lambda x: x[1])
            (maxdx,minn) = dxs[-1]
            for dxn in dxs: # make 100 be plotted before 200
                (dx,n) = dxn
                ds = nds[dxn]
                scs = sorted(list(ds.keys()))
                dts = sorted([dt for dt in ds[scs[0]]])
                mindt = dts[0]
                scs.sort(key=lambda s: integrationPriority(ds[s][mindt][0]["integration"]))
                def pointstyle(s):
                    if "Explicit" in ds[s][mindt][0]["integration"]:
                        return "s"
                    else:
                        return "o"

                if len(scs) <= 1:
                    plt.close()
                    return 
                s1 = scs[0]
                for (linestyle, fillstyle, mmi) in ally:
                    for (s0,col) in zip(scs,plt_setting.clist):
                        ds0 = ds[s0]
                        dtref = sorted(list(ds0.keys()))[0]
                        info = ds0[dtref][0]
                        dx = info["dx"]
                        refs = {s: ds[s][sorted(list(ds[s].keys()))[0]][1] for s in ds}
                        if "FixPoint" in info["integration"]:
                            schemetype = "Implicit"
                        else:
                            schemetype = "Explicit"
                        c = convergence(ds0,refs[s1])
                        al = alpha
                        if n == 100:
                            col = np.roll(col,2) # lighten(col)
                        s = 30
                        sizes = [4*s if abs(dt-dx/10)<1e-10 else s for dt in c[fromref:,4]]
                        facecolors = fillstyle
                        if fillstyle == "full":
                            facecolors = col
                        def pl(ax,label):
                            ax.loglog(c[fromref:,dci],c[fromref:,mmi], label=label, color=col, linestyle=(0,linestyle), linewidth=1, alpha=al)
                            ax.scatter(c[fromref:,dci],c[fromref:,mmi],sizes, marker=pointstyle(s0), facecolors=facecolors, color=col, alpha=al) #, fillstyle=fillstyle, color=col, alpha=al)
                        if dx == maxdx:
                            pl(ax,schemetype)
                        if schemetype == "Implicit":
                            label = r"$\Delta x = "+str(dx)+"$ fm"
                            pl(ax2,label)

        def line(col,style):
            return matplotlib.lines.Line2D([],[], color=col,linestyle=style)
        def clean(ax):
            handles, labels = [], []
            hs, ls = ax.get_legend_handles_labels()
            for (h,l) in zip(hs,ls):    # this is the loop to change Labels and colors
                col = h.get_color()
                m = max(col)
                if not (l in labels or m == 1):    # check for Name already exists
                    p = line(col,"-")
                    handles += [p]
                    labels += [l]
            handles = [line("black",(0,(1,0))),line("black",(0,(5,5)))]+handles
            labels = ["max", "mean"]+labels
            leg = ax.legend(handles, labels)
            for lh in leg.legendHandles:
                lh.set_alpha(1)
        clean(ax)
        clean(ax2)
        fig.savefig("figures/convergence_{}_meanmax_crop:{}_{}.pdf".format(dtcost, crop, info2name(info, False)))
        fig2.savefig("figures/convergence_{}_meanmax_dx_crop:{}_{}.pdf".format(dtcost, crop, info2name(info, False)))
        plt.close()

def integrationPriority(integration):
    if "Explicit" in integration:
        return 1
    else:
        return 2

def mask(data,vid,cut):
    m = (data[:,vid["e"]]>cut).astype(int)
    return np.concatenate((data[:,:3], data[:,3:]/m[:,None]), axis=1)

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
    dts = common(dtss[0],dtss[1])
    dts = [dts[0],dts[-1]] # only keep best and worst
    mindt = dts[0]
    maxdt = dts[-1]
    schemes.sort(key=lambda s: integrationPriority(datas[s][mindt][0]["integration"]))
    for dt in dts:
        timename = "dt{:.4e}".format(dt)
        if dt == maxdt:
            timename = "worst_"+timename
        elif dt == mindt:
            timename = "best_"+timename

        # plt.rcParams["figure.figsize"] = [8, 12]
        # _,axs = plt.subplots(4, 1, sharex=True)
        plt.rcParams["figure.figsize"] = [8, 9]
        fig,axs = plt.subplots(2, 1, sharex=True)
        if "Riemann" in name:
            axin = axs[0].inset_axes([0.03, 0.05, 0.3, 0.5])
            axin.set_xticklabels([])
            axin.set_yticklabels([])
            # if not "Void" in name:
            axinv = axs[0].inset_axes([0.55, 0.45, 0.4, 0.5])
            axinv.set_xticklabels([])
            axinv.set_yticklabels([])
        l = 0
        for dxn in nds:
            (dx,n) = dxn
            l = max(l,dx*n/2)
        x = np.linspace(-l,l,1000)
        if not "Gubser" in name:
            x = x/(t-t0)
        ycontinuum = [e(x) for x in x]
        continuum ,= axs[0].plot(x,ycontinuum, color="black", label="continuum", linewidth=1)
        if "Riemann" in name:
            axin.plot(x,ycontinuum, color="black", label="continuum", linewidth=1)
            # if not "Void" in name:
            axinv.plot(x,ycontinuum, color="black", label="continuum", linewidth=1)
        for dxn in sorted(nds,key=lambda x: x[1]):
            (dx,n) = dxn
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
            ycontinuum = [e(x) for x in x]

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
                yerr = [(a-b)/max(abs(a),abs(b)) for (a,b) in zip(y,ycontinuum)]
        
                numerics ,= axs[0].plot(x,y, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )
                if "Riemann" in name:
                    axin.plot(x,y, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )
                    # if not "Void" in name:
                    axinv.plot(x,y, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )

                errcontinuum ,= axs[1].plot(x,yerr, label=schemetype, color=col, linestyle=linestyle, linewidth=3 )
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
        handles = [line("black",lstyles2[0]),line("black",lstyles1[0])]+handles
        labels = [r"$\Delta x = 0.2$ fm", r"$\Delta x = 0.1$ fm"]+labels
        axs[1].legend(handles, labels) #, loc="lower left", bbox_to_anchor=(0, 1.02, 1, 0.2), mode="expand",borderaxespad=0,ncol=3)
        axs[1].set_ylabel(r"$\Delta\epsilon$")
        # axs[2].set_ylabel("cost")
        axs[len(axs)-1].set_xlabel("$x/t$")
        plt.subplots_adjust(wspace=0,hspace=0)
        plt.savefig("figures/{}_{}.pdf".format(timename,info2name(info, False)))
        plt.close()


    for dxn in nds:
        (dx,n) = dxn
        datas = nds[dxn]
        for dt in dts:
            timename = "dt{:.4e}".format(dt)
            if dt == maxdt:
                timename = "worst_"+timename
            elif dt == mindt:
                timename = "best_"+timename

            for scheme in schemes: 
                fig,ax = plt.subplots(1, 1)
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
                cost = iter*nbStages
                l = xx[0]
                r = xx[-1]
                d = tt[0]
                u = tt[-1]
            
        
                if "Riemann" in name:
                    if "Void" in name:
                        vs = riem_void_vs
                    else:
                        vs = riem_vs
                    ax.plot(xx, xx/vs, linestyle="--", color="white", label="shock")
                    ax.plot(xx, -xx/cs, linestyle="-.", color="white", label="rarefaction")
                im = ax.imshow(cost, extent=[l,r,d,u], origin="lower", label="aoeu") #, norm=CenteredNorm(0)) # , cmap="terrain"
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position('top') 
                ax.set_xlabel("$x$ (fm)")
                ax.set_ylabel("$t$ (fm)")
                divider = make_axes_locatable(ax)
                cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
                fig.add_axes(cax)
                cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
                cbar.formatter.set_powerlimits((0, 0))
                cbar.formatter.set_useMathText(True)
                cbar.update_ticks()
                cbar.set_label("cost", labelpad=-60)

                if "Riemann" in name:
                    ax.legend()
                ax.text(0.7, 0.1, r"$\Delta x = "+str(dx)+"$ fm", color="white", #, bbox={"facecolor": "white", "pad": 10},
                    transform=ax.transAxes, fontsize=22)
                plt.savefig("figures/{}_{}_cost-t_{}.pdf".format(timename,scheme,info2name(sinfo)), dpi=100)
                plt.close()

gref = defaultdict(lambda: None)
greft = defaultdict(lambda: defaultdict(lambda: None))
def plot2d(l, datadts):
    [dim,name,visc,t0,tend,t,case,dxn,scheme] = l
    if case > 0 and not manycases:
        return
    (dx,n) = dxn
    maxdts = sorted([dt for dt in datadts])
    (info, ref, diffref) = datadts[maxdts[0]]
    vid = info["ID"]
    cut = info["CUT"]
    ref = mask(ref,vid,cut)
    fromref = defaultfromref
    if len(maxdts) == 1:
        fromref = 0
    dt = maxdts[fromref]
    (einfo, data, diff) = datadts[dt]
    vid = einfo["ID"]
    cut = einfo["CUT"]

    name = einfo["name"]
    case = einfo["case"]
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
        # nb = 4
        nb = 1
        if "Gubser" in info["name"]:
            nb += 1
        if "FixPoint" in info["integration"]:
            nb += 1
        plt.rcParams["figure.figsize"] = [2+num*4, nb*5]
        fig, axs = plt.subplots(nb,num, sharey=True) #, sharex=True)
        if not hasattr(axs[0], "__len__"):
            axs = [axs]
        for (id, (t,data,diff)) in zip(range(num),datats[nums]):
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
            # print(name, case, scheme, dt, n, ziter.sum())
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
            # all = [("e", z),("vy",zvy), ("err ref",zerrref), ("vx", zvx)]
            all = [("e", z)]
            if "Gubser" in info["name"]:
                all += [(r"$\Delta\epsilon$", zerr)]
            if "FixPoint" in info["integration"]:
                all += [("iter", ziter)]
            for (i, (n, z)) in zip(range(nb),all):
                if n == "iter":
                    im = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower", vmin=0, vmax=3)
                else:
                    im = axs[i][id].imshow(z, extent=[l,r,d,u], origin="lower") #, norm=CenteredNorm(0)) # , cmap="terrain"
                axs[i][id].xaxis.tick_top()
                axs[i][id].xaxis.set_label_position('top') 
                if i == 0:
                    axs[i][id].set_xlabel("$x$ (fm)")
                if i>0:
                    axs[i][id].tick_params(axis='x', which='both', labelbottom=False, labeltop=False)
                if id == 0:
                    axs[i][id].set_ylabel("$y$ (fm)")
                divider = make_axes_locatable(axs[i][id])
                cax = divider.new_vertical(size="5%", pad=0.6, pack_start=True)
                fig.add_axes(cax)
                cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
                cbar.formatter.set_powerlimits((0, 0))
                cbar.formatter.set_useMathText(True)
                cbar.update_ticks()
                cbar.set_label("{} (t = {:.2} fm)".format(n, t), labelpad=-60)

        ax = axs[nb-1][0]
        ax.text(0.075, 0.1, r"$\Delta x = "+str(dx)+"$ fm", color="white", #, bbox={"facecolor": "white", "pad": 10},
            transform=ax.transAxes, fontsize=22)
        plt.subplots_adjust(wspace=0)
        plt.savefig("figures/many_best_e_{}.pdf".format(info2name(info)), dpi=100)
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
        n = info["nx"]
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
        if "Gubser" in info["name"] and not "Exponential" in info["name"]:
            all += [("err", zerr)]
        if "FixPoint" in info["integration"]:
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
            totdt00 = dt00.sum()/abs(data[:,2]).sum()
            totdt01 = dt01.sum()/abs(data[:,3]).sum()
            totdt02 = dt02.sum()/abs(data[:,4]).sum()
            all = [("diff T00",dt00,totdt00),("diff T01",dt01,totdt01),("diff T02",dt02,totdt02)]
            nb = len(all)
            plt.rcParams["figure.figsize"] = [2+nb*4, 5]
            fig, axs = plt.subplots(1,nb, sharey=True)
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
    
alldata(6, datas, convall)
alldata(7, datas, plotall1D)
alldata(9, datas, plotall2D)
