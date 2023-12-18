#!/usr/bin/env python3

import freestream
import sys
import os
import shutil
import numpy as np

random_seed = 1
num = 10
pn = int(1+np.log10(num-1))
half_size = 20 # fm
half_size_eta = half_size/4 # fm


# coefficients from "PHYSICAL REVIEW C 101, 024911 (2020)"
p = 0
k = 0.19
rcp = 0.81
nc = 6
wc = 0.43
dmin = 0.81
tau_fs = 0.37 # fm

# handle arguments
cells, args = [], []
for arg in sys.argv[1:]:
    (args if arg.startswith("-") else cells).append(arg)

sigs = [7.0]
energies = [5020] # GeV
norms = [20] # GeV
bs = [10.5]
random_b = False
usefreestream = False
use3d = False
Nscale = 2.6

for arg in args:
    if "-n=" in arg:
        num = int(arg[3:])

    if "-norm=" in arg:
        norms = [float(arg[6:])]

    if "-Nscale=" in arg:
        Nscale = float(arg[8:])

    # if '-mb' is in the argument list, generate 'many' impact parameters
    if "-mb" == arg:
        bs = [3.3,10.5,13.6] # centrality [5%, 45%, 75%]

    if "-rb" == arg:
        random_b = True
        
    if "-f" == arg:
        usefreestream = True

    if "-3d" == arg:
        use3d = True
        usefreestream = False

for c in cells:
    dx = 2*half_size / float(c)
    for sig, energy, normi in zip(sigs,energies,norms):
        for b in bs:
            if random_b:
                dir_b = "_random"
                bminmax = ""
            else:
                dir_b = b
                bminmax = " --b-min {bmin} --b-max {bmax} ".format(bmin=b, bmax=b)
            dir = "TeV{}/b{}/s{}".format(energy, dir_b, c)
            try:
                os.makedirs(dir)
            except FileExistsError:
                shutil.rmtree(dir)
                os.makedirs(dir)

            # run trento
            # git: Duke-QCD/frzout

            eta_c = int(float(c)/2)
            eta_c += 1-eta_c%2;
            norm = normi/0.1973 # fm^-1

            if use3d:
                # Duke-QCD/trento3d-2.0
                # executable: trento-3
                #
                # u [--form-width]: form width
                # w [-w]: nucleon width
                # nc [-m]: constituent number
                # ki: Structure
                # v [-v]: constituent width ki*(w*nc**0.25)
                # kt,min [-t]: transverse mom scale
                # alpha [--shape-alpha]: shape parameter
                # betha [--shape-beta]: shape parameter
                # Nmid [--mid-norm]: mid norm
                # pmid [--mid-power]: mid power
                # Nfb: fireball normalization Nfb = Nmid*mp*(sqrt(sNN)/mp)**pmid
                # sqrt(sNN) [-s]: collision energy GeV
                # mp: nucleon mass
                # k [-k]: fluctuation
                # f [--flatness]: flatness
                # tau0,Pb: hydrodynamisation time
                # Nscale [--overall-norm]: overall scale
                #
                # From arxiv 2306.08665
                # u -> 0.88
                # w -> 1.3
                # nc -> 16.4
                # ki -> 0.50     ---->  v: 1.308 however v<=w so v: 1.3
                # kt,min -> 0.33
                # alpha -> 4.6
                # beta -> 0.19
                # N200 -> 9.6    ----\  Nmid: 1.713
                # N5020 -> 28.1  ----/  pmid: 0.333
                # k -> 0.104
                # f -> 1.0
                # tau0,Pb -> 1.3
                # N -> 2.0       ----> Nscale: 2.0

                # trento_cmd = "trento-3 Pb Pb -s {s} --random-seed {random_seed} --b-min {bmin} --b-max {bmax} --grid-max {grid_max} --grid-step {grid_step} --nsteps-etas {nsteps_etas} --form-width {u} -w {w} -m {nc} -v {v} -t {ktmin} --shape-alpha {alpha} --shape-beta {beta} --mid-norm {Nmid} --mid-power {pmid} -k {k} --flatness {f} --overall-norm {Nscale} {num} -o {dir}"\
                #     .format(u=0.88, w=1.3, nc=16.4, v=1.3, ktmin=0.33, alpha=4.6, beta=0.19, Nmid=1.713, pmid=0.333, k=0.104, f=1.0, Nscale=2.6, nsteps_etas=eta_c, grid_step=dx, grid_max=half_size, s=energy, random_seed=random_seed, bmin=b, bmax=b, num=num, dir=dir)
                trento_cmd = "trento-3 Pb Pb -s {s} --random-seed {random_seed} --grid-max {grid_max} --grid-step {grid_step} --nsteps-etas {nsteps_etas} {bminmax} --form-width {u} -d {d} -w {w} -m {nc} -t {ktmin} --shape-alpha {alpha} --shape-beta {beta} --mid-norm {Nmid} --mid-power {pmid} -k {k} --flatness {f} --overall-norm {Nscale} {num} -o {dir}"\
                    .format(u=rcp, d=dmin, w=wc, nc=nc, ktmin=0.33, alpha=4.6, beta=0.19, Nmid=0.113, pmid=0.615, k=k, f=1.0, Nscale=12, nsteps_etas=eta_c, grid_step=dx, grid_max=half_size, s=energy, random_seed=random_seed, bminmax=bminmax, num=num, dir=dir)
            else:
                trento_cmd = \
                    "trento Pb Pb --random-seed {r} -p {p} -k {k} -w {w} -m {m} -v {v} -d {d} {bminmax} --cross-section {sig} --normalization {n} --grid-max {l} --grid-step {dx} {num} -o {dir}" \
                    .format(r=random_seed, p=p, k=k, w=rcp, m=nc, v=wc, d=dmin, sig=sig, bminmax=bminmax, l=half_size, dx=dx, n=norm, num=num, dir=dir)
            trento = os.popen(trento_cmd)
            output = trento.read()
            print("Trento:\n{}".format(output))
            if use3d:
                etas_len = 2.0 * np.arccosh(0.5 * energy / 0.2); # from github Duke-QCD/trento3d-2.0, commit 503a7d744df29aefbd2e31b7b623aff07abe2bf2, event.cxx line 63
                detas = etas_len/eta_c
                with open("{}/info.txt".format(dir), "w") as fv:
                   fv.write("dx: {}\ndy: {}\ndetas: {}\n".format(dx,dx,detas))
