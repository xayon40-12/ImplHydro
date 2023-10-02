#!/usr/bin/env python3

import os
import sys
import numpy as np
import frzout

# run frsout
# git: Duke-QCD/trento

nb_frzout = 10

dim = 2
particlize = False
use_urqmd = False
for arg in sys.argv[1:]:
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
   surface = frzout.Surface(x, sigma, v, pi=pi, Pi=Pi)

   hrg = frzout.HRG(.155, species='urqmd', res_width=True)
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

cwd = os.getcwd()
dirs = os.listdir(cwd)
if "results" in dirs:
   res = "{}/results".format(cwd)
   dirs = os.listdir(res)
   allparts = {}
   allcounts = {}
   allinfos = {}
   for d in dirs:
      d = "{}/{}".format(res,d)
      os.chdir(d)
      if particlize:
         print(d)
         particlization()

      info = {k: v.strip() for [k, v] in np.loadtxt("info.txt", dtype=object, delimiter=":")}
      info.pop("case")
      parts = []
      current = []
      # each line: ID charge pT ET mT phi y eta
      with open("particles_out.dat", "r") as f:
         for p in f.read().split("\n"):
            if "#" in p:
               if len(current) > 0:
                  parts += [np.array(current)]
               current = []
            else:
               if len(p) > 0:
                  vals = np.array([np.float64(i) for i in p.split()])
                  current += [vals]
         parts += [np.array(current)]
      count = len(parts)
      id = str(info)
      if id in allparts:
         allparts[id] += parts
         allcounts[id] += count
      else:
         allparts[id] = parts
         allcounts[id] = count
         allinfos[id] = info
   os.chdir(res)
   for k in allparts:
      counts = allcounts[k]
      parts = allparts[k]
      info = allinfos[k]

      # TODO: compute observables from the particles informations
      print(counts, info)
      
else:
   print("No \"results\" folder found.")