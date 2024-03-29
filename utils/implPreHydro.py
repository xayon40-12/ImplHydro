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
half_size_eta = half_size/2 # fm


# coefficients from "PHYSICAL REVIEW C 101, 024911 (2020)"
p = 0.002
sigmafluct = 0.90
rcp = 0.88
nc = 6
wc = 0.53
dmin = 1.12
norm = 20/0.1973 # fm^-1
tau_fs = 0.48 # fm

# handle arguments
cells, args = [], []
for arg in sys.argv[1:]:
    (args if arg.startswith("-") else cells).append(arg)

sigs = [6.4]
bs = [3]
# if '-m' is in the argument list, generate 'many' cases
if "-m" in args:
    sigs = [4.23, 6.4, 7.32]
    bs = [0,3,7]

usefreestream = False
if "-f" in args:
    usefreestream = True

use3d = False
if "-3d" in args:
    use3d = True

for c in cells:
    dx = 2*half_size / float(c)
    for sig in sigs:
        for b in bs:
            dir = "sig{}/b{}/s{}".format(sig, b, c)
            try:
                os.makedirs(dir)
            except FileExistsError:
                shutil.rmtree(dir)
                os.makedirs(dir)

            # run trento
            # git: Duke-QCD/frzout

            if use3d:
                trento_cmd = \
                    "trento3d Pb Pb --eta-max {etam} --eta-step {dx} --random-seed {r} -p {p} -k {k} -w {w} -d {d} --b-min {b} --b-max {b} --cross-section {sig} --normalization {n} --xy-max {l} --xy-step {dx} {num} -o {dir}" \
                    .format(etam=half_size_eta, r=random_seed, p=p, k=sigmafluct, w=rcp, d=dmin, sig=sig, b=b, l=half_size, dx=dx, n=norm, num=num, dir=dir)
            else:
                trento_cmd = \
                    "trento Pb Pb --random-seed {r} -p {p} -k {k} -w {w} -m {m} -v {v} -d {d} --b-min {b} --b-max {b} --cross-section {sig} --normalization {n} --grid-max {l} --grid-step {dx} {num} -o {dir}" \
                    .format(r=random_seed, p=p, k=sigmafluct, w=rcp, m=nc, v=wc, d=dmin, sig=sig, b=b, l=half_size, dx=dx, n=norm, num=num, dir=dir)
            trento = os.popen(trento_cmd)
            output = trento.read()
            print("Trento:\n{}".format(output))

            # run free streaming
            # git: Duke-QCD/freestream
            if usefreestream:
                print("freestream")
                for i in range(num):
                    name = ("{}/{:0<"+str(pn)+"}.dat").format(dir,i)
                    print(name)
                    ic = np.loadtxt(name, dtype=np.float64)
                    fs = freestream.FreeStreamer(ic, half_size, tau_fs)

                    e = fs.energy_density()
                    ut = fs.flow_velocity(0)
                    ux = fs.flow_velocity(1)
                    uy = fs.flow_velocity(2)
                    def pi(u,v):
                        return fs.shear_tensor(u,v)
                    bulk = fs.bulk_pressure()  # fs.bulk_pressure(eos)
                    arr = np.transpose(np.array([e, ut, ux, uy, pi(0,0), pi(0,1), pi(0,2), pi(1,1), pi(1,2), pi(2,2), bulk]).reshape((11,-1)))
                    # np.set_printoptions(threshold=sys.maxsize)
                    # print(arr)
                    data = ("{}/data{:0<"+str(pn)+"}.dat").format(dir,i)
                    arr.tofile(data)

                    # etot = e.sum()*197*dx**3/1000 # GeV
                    # em = e.max()*197/1000 # GeV fm^-3
                    # print("e_tot: {} GeV".format(etot))
                    # print("e_max: {} GeV fm^-3".format(em))

