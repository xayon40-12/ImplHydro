#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def parse_value(vs):
    return [float(v) for v in vs.split(":")]

def convert(line):
    [name, centrality, values] = line
    centrality = [float(c) for c in centrality.split("-")]
    values = [parse_value(v) for v in values.strip().split(" ")]
    if len(values) == 1:
        [value, err] = values[0]
        return (name, centrality, (value, err))
    else:            
        poss = [v[0] for v in values]
        errs = [v[2] for v in values]
        values = [v[1] for v in values]
        return (name, centrality, (poss, values, errs))
        

def parse(lines):
    tmp = [convert(l.strip().split("|")) for l in lines.strip().split("\n")]
    name = tmp[0][0]
    centralities = [c for (_, c, _) in tmp]
    values = [v for (_, _, v) in tmp]
    return name, centralities, values

def plot_dndeta_eta(dndeta_eta):
    plt.close()
    name, centralities, valuess = dndeta_eta
    for [cl,cr], (pos, value, err) in zip(centralities, valuess):
        plt.plot(pos, value, label="{}%-{}%".format(cl,cr))
    plt.xlabel("$\eta$")
    plt.ylabel("$\mathrm{dN}_{\mathrm{ch}}/\mathrm{d}\eta$")
    plt.legend()
    plt.savefig("dndeta-eta.pdf")

def plot_dndeta_mid(dndeta_mid):
    plt.close()
    alice5TeV_centralities = [[0,2.5],[2.5,5.0],[5.0,7.5],[7.5,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80]]
    alice5TeV_values = [[2035, 52],[1850, 55],[1666, 48],[1505, 44],[1180, 31],[786, 20],[512, 15],[318, 12],[183, 8],[96.3, 5.8],[44.9, 3.4]]

    name, centralities, valuess = dndeta_mid

    for name, shape, centralities, valuess in [("numerics", "-", centralities, valuess), ("Alice 5.02TeV", "o", alice5TeV_centralities, alice5TeV_values)]:
        centralities_c = [(l+r)/2 for [l,r] in centralities]
        values = [v[0] for v in valuess]
        errs = [v[1] for v in valuess]
        plt.plot(centralities_c, values, shape, label=name)
    plt.xlabel("centrality (%)")
    plt.ylabel("$\mathrm{dN}_{\mathrm{ch}}/\mathrm{d}\eta$")
    plt.legend()
    plt.savefig("dndeta-mid.pdf")

def plot_vn(vns):
    alice5TeV_centralities = [[0,5],[5.0,10],[10,20],[20,30],[30,40],[40,50],[50,60],[60,70],[70,80]]
    alice5TeV_values = [[0.029, 0], [0.047, 0], [0.065, 0], [0.088, 0], [0.099, 0], [0.105, 0], [0.103, 0], [0.098, 0], [0.089, 0]]
    alice = [("Alice 5.02TeV", "o", alice5TeV_centralities, alice5TeV_values)]
    vns = [("$v_{{{}}}${{2}}".format(n[1]), "-", c, v) for (n, c, v) in vns]
    plt.close()
    for name, shape, centralities, valuess in vns+alice:
        centralities_c = [(l+r)/2 for [l,r] in centralities]
        values = [v[0] for v in valuess]
        errs = [v[1] for v in valuess]
        plt.plot(centralities_c, values, shape, label=name)
    plt.xlabel("centrality (%)")
    plt.ylabel("$v_n${2}")
    plt.legend()
    plt.savefig("vn.pdf")

with open("observables.txt") as f:
    obs = [parse(lines) for lines in f.read().split("\n\n")]
    vns = [o for o in obs if o[0].startswith("v")]
    plot_vn(vns)
    dndeta_mid = [o for o in obs if o[0] == "dn/deta-mid"][0]
    plot_dndeta_mid(dndeta_mid)
    dndeta_eta = [o for o in obs if o[0] == "dn/deta"][0]
    plot_dndeta_eta(dndeta_eta)
