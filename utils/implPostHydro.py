#!/usr/bin/env python3

import os
import sys
import numpy as np
import frzout
from cmath import exp

# run frsout
# git: Duke-QCD/trento

nb_frzout = 10
ymax = 0.5
dim = 2
particlize = False
use_urqmd = False
for arg in sys.argv[1:]:
   if "-ymax=" in arg:
      ymax = arg[6:]
   if arg == "-3d":
      dim = 3
   if arg == "-p":
      particlize = True
   if arg == "-u":
      use_urqmd = True
      

def particlization():
   print("Running frzout:")
   if dim == 2:
      # t x y sig_t sig_x sig_y vx vy pi00 pi01 pi02 pi11 pi12 pi22 pi33 Pi
      # 0 1 2 3     4     5     6  7  8    9    10   11   12   13   14   15
      surface_data = np.fromfile('surface.dat', dtype='float64').reshape((-1, 16))

      # extract usual sub-arrays
      x, sigma, v, _ = np.hsplit(surface_data, [3, 6, 8])
      tend = x[-1][0]

      # create mapping of pi components
      pi = dict(
         xx=surface_data[:, 11],
         xy=surface_data[:, 12],
         yy=surface_data[:, 13]
      )

      # extract bulk pressure
      Pi = surface_data[:, 15]
   elif dim == 3:
      # t x y z sig_t sig_x sig_y zig_z vx vy vz pi00 pi01 pi02 pi03 pi11 pi12 pi13 pi22 pi23 pi33 Pi
      # 0 1 2 3 4     5     6     7     8  9  10 11   12   13   14   15   16   17   18   19   20   21
      surface_data = np.fromfile('surface.dat', dtype='float64').reshape((-1, 22))

      # extract usual sub-arrays
      x, sigma, v, _ = np.hsplit(surface_data, [4, 8, 11])
      tend = x[-1][0]

      # create mapping of pi components
      pi = dict(
         xx=surface_data[:, 15],
         xy=surface_data[:, 16],
         xz=surface_data[:, 17],
         yy=surface_data[:, 18],
         yz=surface_data[:, 19]
      )

      # extract bulk pressure
      Pi = surface_data[:, 21]

   print("tend: {}".format(tend))
   # create Surface object
   surface = frzout.Surface(x, sigma, v, ymax=ymax, pi=pi, Pi=Pi)

   hrg = frzout.HRG(.148, species='urqmd', res_width=True)
   lines = ""
   for i in range(nb_frzout):
      parts = frzout.sample(surface, hrg)
      nparts = len(parts)
      print("nb parts: {}".format(nparts))
      lines += "# {}\n".format(nparts)
      for ID, x, p in parts:
         l = "{} {} {} {} {} {} {} {} {}\n".format(ID, x[0], x[1], x[2], x[3], p[0], p[1], p[2], p[3])
         lines += l
      lines += "\n"

   with open("particles_in.dat", "w") as f:
      f.write(lines)

   if use_urqmd:
      # run urqmd
      # git: jbernhard/urqmd-afterburner
      #
      # urqmd needs to be installed before. To compile with recent gcc-fortran,
      # the following flags might be needed:
      #   export FCFLAGS="-w -fallow-argument-mismatch -O2"
      #   export FFLAGS="-w -fallow-argument-mismatch -O2"

      print("Running UrQMD:")
      urqmd = os.popen("afterburner particles_in.dat particles_out.dat")
      output = urqmd.read()
      with open("urqmd.log", "w") as f:
         f.write(output)

   print("done\n")

def vn(n, events):
   l = len(events)
   ms = [len(e) for e in events]

   qn = [sum(exp(1j*n*v["phi"]) for v in p) for p in events]
   qn2 = [abs(z)**2 for z in qn]
   
   w2 = [m*(m-1) for m in ms]
   
   m2 = [(qn2[i]-ms[i])/w2[i] for i in range(l)]
   
   cn2_n = 0
   cn2sq_n = 0
   cn2_d = 0
   for i in range(l):
      cn2_n += w2[i]*m2[i]
      cn2sq_n += w2[i]*m2[i]**2
      cn2_d += w2[i]
   cn2 = cn2_n/cn2_d
   cn2sq = cn2sq_n/cn2_d
   sig = np.sqrt(cn2sq-cn2**2)
   err = sig/np.sqrt(l)

   vn = np.sqrt(cn2)
   errvn = err/(2*vn)

   return vn, errvn

totiter = 0
cwd = os.getcwd()
dirs = os.listdir(cwd)
if "results" in dirs:
   res = "{}/results".format(cwd)
   dirs = os.listdir(res)
   allevents = {}
   allcounts = {}
   allinfos = {}
   iter = 0
   for d in dirs:
      if iter >= totiter and totiter != 0:
         break
      iter += 1
      d = "{}/{}".format(res,d)
      os.chdir(d)
      if particlize:
         print(d)
         particlization()

      info = {k: v.strip() for [k, v] in np.loadtxt("info.txt", dtype=object, delimiter=":")}
      info.pop("case")
      events = []
      parts = []
      # each line: ID charge pT ET mT phi y eta
      pnames = ["ID", "charge", "pT", "ET", "mT", "phi", "y", "eta"]
      with open("particles_out.dat", "r") as f:
         for p in f.read().split("\n"):
            if "#" in p:
               if len(parts) > 0:
                  events += [np.array(parts)]
               parts = []
            else:
               if len(p) > 0:
                  vals = {name: np.float64(i) for name, i in zip(pnames, p.split())}
                  parts += [vals]
         events += [np.array(parts)]
      count = len(events)
      id = str(info)
      if id in allevents:
         allevents[id] += events
         allcounts[id] += count
      else:
         allevents[id] = events
         allcounts[id] = count
         allinfos[id] = info
   os.chdir(cwd)
   for k in allevents:
      counts = allcounts[k]
      events = allevents[k]
      info = allinfos[k]
      print(counts, info)

      v2, errv2 = vn(2, events)
      v3, errv3 = vn(3, events)
      v4, errv4 = vn(4, events)
      print(v2)
      with open("v2.txt", "w") as fv:
          fv.write("v2: {:e}\nerr: {:e}\n".format(v2, errv2))
          fv.write("v3: {:e}\nerr: {:e}\n".format(v3, errv3))
          fv.write("v4: {:e}\nerr: {:e}\n".format(v4, errv4))
      
else:
   print("No \"results\" folder found.")