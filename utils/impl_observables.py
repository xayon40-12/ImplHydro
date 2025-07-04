#!/usr/bin/env python3

import os
import numpy as np
import yaml
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
    max_c = 70
    for [cl,cr], (pos, value, err) in zip(centralities, valuess):
        if cl < max_c:
            plt.errorbar(pos, value, err, label="{}%-{}%".format(cl,cr))
    fst = True
    for [cl,cr], (pos, value, err) in alice5020["dndeta-eta"]:
        if cl < max_c:
            pi = len([p for p in pos if p < 0])
            name = fst and "ALICE 5.02TeV" or ""
            fst = False
            plt.errorbar(pos[pi:], value[pi:], err[pi:], linestyle="", marker="o", color="black", label=name)
    plt.xlabel("$\eta$")
    plt.ylabel("$\mathrm{dN}_{\mathrm{ch}}/\mathrm{d}\eta$")
    plt.legend()
    plt.savefig("dndeta-eta.pdf")

def plot_d2ndetadpt_pt(d2ndetadpt_pt):
    plt.close()
    name, centralities, valuess = d2ndetadpt_pt
    n = len(valuess)
    max_c = 80
    for [cl,cr], (pos, value, err) in zip(centralities, valuess):
        if cl < max_c:
            plt.errorbar(pos, np.array(value)*10**n, np.array(err)*10**n, label="{}%-{}% $\\times 10^{}$".format(cl,cr,n))
        n -= 1
    fst = True
    n = len(valuess)
    for [cl,cr], (pos, value, err) in alice5020["d2ndetadpt-pt"]:
        if cl < max_c:
            pi = len([p for p in pos if p < 0])
            name = fst and "ALICE 5.02TeV" or ""
            fst = False
            plt.errorbar(pos[pi:], np.array(value[pi:])*10**n, np.array(err[pi:])*10**n, linestyle="", marker="o", markersize = 2, color="black", label=name)
        n -= 1
    plt.xlabel("$p_T$")
    plt.ylabel("$1/\mathrm{N_{evt}}\mathrm{d^2N}_{\mathrm{ch}}/\mathrm{d}\eta\mathrm{dp_T}$")
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("d2ndetadpt-pt.pdf")

def plot_dndeta_mid(dndeta_mid):
    plt.close()

    name, centralities, valuess = dndeta_mid
    centralities_c = [(l+r)/2 for [l,r] in centralities]
    values = [v[0] for v in valuess]
    errs = [v[1] for v in valuess]
    # plt.yscale('log')
    plt.errorbar(centralities_c, values, errs, label="numerics")

    centralities, valuess = alice5020["dndeta"]
    centralities_c = [(l+r)/2 for [l,r] in centralities]
    values = [v[0] for v in valuess]
    errs = [v[1] for v in valuess]
    plt.errorbar(centralities_c, values, errs, linestyle="", marker="o", color="black", label="ALICE 5.02TeV")

    plt.xlabel("centrality (%)")
    plt.ylabel("$\mathrm{dN}_{\mathrm{ch}}/\mathrm{d}\eta$")
    plt.legend()
    plt.savefig("dndeta-mid.pdf")

def plot_vn(vns):
    alice = [alice5020[n[0:2]] for (n, _, _) in vns]
    def name(n):
        return "$v_{{{}}}${{2}}".format(n[1])
    vns = [(name(n), c, v) for (n, c, v) in vns]
    plt.close()
    for name, centralities, valuess in vns:
        centralities_c = [(l+r)/2 for [l,r] in centralities]
        values = [v[0] for v in valuess]
        errs = [v[1] for v in valuess]
        plt.errorbar(centralities_c, values, errs, label=name)
    fst = True
    for centralities, valuess in alice:
        centralities_c = [(l+r)/2 for [l,r] in centralities]
        values = [v[0] for v in valuess]
        errs = [v[1] for v in valuess]
        name = fst and "ALICE 5.02TeV" or ""
        fst = False
        plt.errorbar(centralities_c, values, errs, linestyle="", marker="o", color="black", label=name)
    plt.xlabel("centrality (%)")
    plt.ylabel("$v_n${2}")
    plt.legend()
    plt.savefig("vn.pdf")

def tosym(err):
    if "asymerror" in err:
        e = err["asymerror"]
        return max(-e["minus"],e["plus"])
    else:
        return err["symerror"]

alice5020_path = os.environ.get("ALICE5020",".")

alice5020 = {"v2": None, "v3": None, "v4": None, "dndeta": None, "dndeta-eta": None, "d2ndetadpt-pt": None}
with open(alice5020_path+"/ALICE_5.02TeV-v2.yaml", "r") as f:
    y = yaml.safe_load(f)
    ivs = y["independent_variables"][0]["values"]
    cs = [[v[d] for d in ["low","high"]] for v in ivs]
    dvs = y["dependent_variables"][0]["values"]
    vs = [v["value"] for v in dvs]
    errs = [sum(v["errors"][i]["symerror"] for i in [0,1]) for v in dvs]
    alice5020["v2"] = (cs, list(zip(vs, errs)))

with open(alice5020_path+"/ALICE_5.02TeV-v3-v4.yaml", "r") as f:
    y = yaml.safe_load(f)
    ivs = y["independent_variables"][0]["values"]
    cs = [[v[d] for d in ["low","high"]] for v in ivs]
    for v, n in [(0, "v3"), (1, "v4")]:
        dvs = y["dependent_variables"][v]["values"]
        vs = [v["value"] for v in dvs]
        errs = [sum(v["errors"][i]["symerror"] for i in [0,1]) for v in dvs]
        alice5020[n] = (cs, list(zip(vs, errs)))

with open(alice5020_path+"/ALICE_5.02TeV-dndeta.yaml", "r") as f:
    y = yaml.safe_load(f)
    ivs = y["independent_variables"][0]["values"]
    cs = [[v[d] for d in ["low","high"]] for v in ivs]
    dvs = y["dependent_variables"][0]["values"]
    vs = [v["value"] for v in dvs]
    errs = [sum(v["errors"][i]["symerror"] for i in [0,1]) for v in dvs]
    alice5020["dndeta"] = (cs, list(zip(vs, errs)))
    
with open(alice5020_path+"/ALICE_5.02TeV-dndeta-eta.yaml", "r") as f:
    y = yaml.safe_load(f)
    ivs = y["independent_variables"][0]["values"]
    etas = [sum(v[d] for d in ["low","high"])/2 for v in ivs]
    dndeta_eta = []
    for ci in range(len(y["dependent_variables"])):
        dvs = y["dependent_variables"][ci]["values"]
        c = [float(v) for v in y["dependent_variables"][ci]["qualifiers"][0]["value"].split("-")]
        vs = [v["value"] for v in dvs]
        errs = [v["errors"][0]["symerror"]+tosym(v["errors"][1]) for v in dvs]
        dndeta_eta += [(c, (etas, vs, errs))]
    alice5020["dndeta-eta"] = dndeta_eta
    
with open(alice5020_path+"/ALICE_5.02TeV-pT-spectra.yaml", "r") as f:
    y = yaml.safe_load(f)
    ivs = y["independent_variables"][0]["values"]
    pts = [sum(v[d] for d in ["low","high"])/2 for v in ivs]
    d2ndetadpt_pt = []
    for ci in range(len(y["dependent_variables"])):
        dvs = y["dependent_variables"][ci]["values"]
        c = [float(v) for v in y["dependent_variables"][ci]["qualifiers"][0]["value"].split(" ")[0].split("-")]
        vs = [v["value"] for v in dvs]
        errs = [v["errors"][0]["symerror"]+tosym(v["errors"][1]) for v in dvs]
        d2ndetadpt_pt += [(c, (pts, vs, errs))]
    alice5020["d2ndetadpt-pt"] = d2ndetadpt_pt
    
with open("observables.txt") as f:
    obs = [parse(lines) for lines in f.read().split("\n\n")]
    vns = [o for o in obs if o[0].startswith("v")]
    plot_vn(vns)
    dndeta_mid = [o for o in obs if o[0] == "dn/deta-mid"][0]
    plot_dndeta_mid(dndeta_mid)
    dndeta_eta = [o for o in obs if o[0] == "dn/deta"][0]
    plot_dndeta_eta(dndeta_eta)
    d2ndetadpt_pt = [o for o in obs if o[0] == "d2n/detadpt"][0]
    plot_d2ndetadpt_pt(d2ndetadpt_pt)
