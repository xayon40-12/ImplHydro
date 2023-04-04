#!/usr/bin/env python3

import freestream
import sys
import os
import shutil
import numpy as np

random_seed = 1
num = 10
pn = int(1+np.log10(num-1))
half_size = 10 # fm
b = 3 # fm
norm = 20000/197 # fm^-1
tau_fs = 0.5 # fm
for c in sys.argv[1:]:
    dx = 2*half_size / float(c)
    for sig in [4.23, 6.4, 7.32]:
        dir = "sig{}/b{}/s{}".format(sig, b, c)
        try:
            os.makedirs(dir)
        except FileExistsError:
            shutil.rmtree(dir)
            os.makedirs(dir)
        trento_cmd = \
            "trento Pb Pb --random-seed {r} --b-min {b} --b-max {b} --cross-section {sig} --normalization {n} --grid-max {l} --grid-step {dx} {num} -o {dir}" \
            .format(r=random_seed, sig=sig, b=b, l=half_size, dx=dx, n=norm, num=num, dir=dir)
        trento = os.popen(trento_cmd)
        output = trento.read()
        print("Trento:\n{}".format(output))

        # free streaming
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

