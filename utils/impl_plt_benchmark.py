#!/usr/bin/env python3
import sys 
import os
import warnings
import matplotlib
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.colors import CenteredNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plt_setting
import numpy as np
from math import sqrt
from scipy import signal
from math import ceil,log
from enum import Enum
import yaml

warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning) # disable matplotlib deprecation warning
np.seterr(invalid='ignore') # disable invalid waring for numpy as NaN are used to discard data in the void
np.seterr(divide='ignore') # disable divide by zero worining

def convall(data):
    ally = [((1,0), "full", 0), ((6,6), "none", 1)]
    allx = [("cost", 0), ("elapsed", 1)]
    for (dtcost, xid) in allx:
        fig, axs = plt.subplots(2, figsize=(8,9), sharex=True)
        fig.subplots_adjust(wspace=0,hspace=0)
        def setlabel(ax):
            name = dtcost
            if dtcost == "cost":
                name = r"$n_\mathrm{KT}$"
            ax.set_xlabel(name)
            ax.set_ylabel(r"$|\Delta_\mathrm{ref}|$")
        for ax in axs:
            setlabel(ax)
        nbcases = len(data['dxs'][0]['cases'])
        alpha = sqrt(1/nbcases)
        markers = {}
        for (ax, dxn) in zip(axs, data['dxs']):
            dx = dxn['dx']
            for case in dxn['cases']:
                ax.text(0.03, 0.05, r"$\Delta x = "+str(dx)+"$ fm", color="black", #, bbox={"facecolor": "white", "pad": 10},
                    transform=ax.transAxes, fontsize=22)
                for (linestyle, fillstyle, yid) in ally:
                    if len(case['integrations']) == 3:
                        colors = [plt_setting.clist[i] for i in [1,2,0]]
                    else:
                        colors = [plt_setting.clist[i] for i in [1,0]]
                    for (inte,col) in zip(case['integrations'],colors): 
                        s = 30
                        sizes = [4*s if abs(a['dt']-dx/10)<1e-10 else s for a in inte['analysis']]
                        facecolors = fillstyle
                        relative_tot_energys = [a['relative_tot_energy'] for a in inte['analysis']]
                        dts = [a['dt'] for a in inte['analysis']]
                        costs = [a['cost'] for a in inte['analysis']]
                        elapseds = [a['elapsed'] for a in inte['analysis']]
                        max_errors = [a['max_error'] for a in inte['analysis']]
                        mean_errors = [a['mean_error'] for a in inte['analysis']]
                        xs = [costs, elapseds][xid]
                        ys = [max_errors, mean_errors][yid]
                        
                        enable = [r < 1e-2 for r in relative_tot_energys]
                        def mask(a,b):
                            return [x for (x,b) in zip(a,b) if b]
                        xs = mask(xs, enable)
                        ys = mask(ys, enable)
                        sizes = mask(sizes, enable)
                        def pl(ax,label):
                            nonlocal facecolors
                            if fillstyle == "full":
                                facecolors = col
                            ax.loglog(xs,ys, label=label, color=col, linestyle=(0,linestyle), linewidth=1, alpha=alpha)
                            if inte['integration'] == 'Explicit':
                                markers[label] = 's'
                            else:
                                markers[label] = 'o'
                            ax.scatter(xs,ys,sizes, marker=markers[label], facecolors=facecolors, color=col, alpha=alpha)
                        if 'Explicit' in inte['integration']:
                            schemetype = 'Heun'
                        else:
                            schemetype = 'GL1'
                        pl(ax,schemetype)

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
            for lh in leg.legend_handles:
                lh.set_alpha(1)
        clean(axs[0])
        labels = ["max"]
        labels += ["mean"]
        handles = [line("black",(0,(1,0)))]
        handles += [line("black",(0,(6,6)))]
        axs[1].legend(handles, labels, loc="upper right")
        fig.savefig("figures/{}_{}.pdf".format(dtcost, data['name']))
        plt.close()

try:
    os.mkdir("figures")
except FileExistsError:
    None
    
with open("benchmark.txt") as f:
    data = yaml.load(f, yaml.Loader)
    convall(data)
