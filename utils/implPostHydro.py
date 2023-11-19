#!/usr/bin/env python3

import os
import sys
import numpy as np
import frzout
from cmath import exp
import math
from multiprocessing import Pool
import matplotlib.pyplot as plt

# run frsout
# git: Duke-QCD/trento


HBARC = 0.1973 # GeV.fm

nb_frzout = 10
nb_threads = 15
ymax = 0.5
dim = 3
particlize = False
use_urqmd = False
for arg in sys.argv[1:]:
   if "-nbfreezeout=" in arg:
      nb_frzout = int(arg[13:])
   if "-ymax=" in arg:
      ymax = arg[6:]
   if "-threads=" in arg:
      nb_threads = int(arg[9:])
   if arg == "-3d":
      dim = 3
   if arg == "-p":
      particlize = True
   if arg == "-u":
      use_urqmd = True

def cont_cov(v):
   v[1] *= -1
   v[2] *= -1
   v[3] *= -1

def MilnetoCart_coord(v):
   t = v[0]*math.cosh(v[3])
   z = v[0]*math.sinh(v[3])

   v[0] = t
   v[3] = z
   
def MilnetoCart_vector(eta, v):
   sh = math.sinh(eta)
   ch = math.cosh(eta)
   at = v[0]*ch + v[3]*sh
   az = v[0]*sh + v[3]*ch

   v[0] = at
   v[3] = az
   
def MilnetoCart_tensor(eta, v): # simetric tensor
   sh = math.sinh(eta)
   ch = math.cosh(eta)
   vtt = v[0][0]*ch*ch + 2.0*v[0][3]*sh*ch + v[3][3]*sh*sh
   vtx = v[0][1]*ch + v[1][3]*sh
   vty = v[0][2]*ch + v[2][3]*sh
   vtz = v[0][0]*sh*ch + v[0][3]*(sh*sh + ch*ch) + v[3][3]*sh*ch
   vxz = v[0][1]*sh + v[1][3]*ch
   vyz = v[0][2]*sh + v[2][3]*ch
   vzz = v[0][0]*sh*sh + 2.0*v[0][3]*sh*ch + v[3][3]*ch*ch

   v[0][0] = vtt
   v[0][1] = vtx
   v[1][0] = vtx
   v[0][2] = vty
   v[2][0] = vty
   v[0][3] = vtz
   v[3][0] = vtz
   v[1][3] = vxz
   v[3][1] = vxz
   v[2][3] = vyz
   v[3][2] = vyz
   v[3][3] = vzz

def MilnetoCart_velocity(eta, v): # v: [vx,vy,veta]
   sh = math.sinh(eta)
   ch = math.cosh(eta)
   vz = (sh+v[2]*ch)/(ch+v[2]*sh)
   vx = v[0]*(ch-vz*sh)
   vy = v[1]*(ch-vz*sh)

   v[0] = vx
   v[1] = vy
   v[2] = vz

def particlization(info):
   print("Running frzout:")
   if dim == 2:
      # t x y sig_t sig_x sig_y vx vy pi00 pi01 pi02 pi11 pi12 pi22 pi33 Pi
      # 0 1 2 3     4     5     6  7  8    9    10   11   12   13   14   15
      surface_data = np.fromfile('surface.dat', dtype='float64').reshape((-1, 16))

      # extract usual sub-arrays
      x, sigma, v, pi, _ = np.hsplit(surface_data, [3, 6, 8, 15])
      Pi = surface_data[:, 15]

      # convert viscosity to GeV
      Pi *= HBARC
      for i in range(len(pi[0])):
         pi[:,i] *= HBARC
         
      tend = x[-1][0]
      pi = [[pi[:,0], pi[:,1], pi[:,2]], [pi[:,1], pi[:,3], pi[:,4]], [pi[:,2], pi[:,4], pi[:,5]]]

      # create mapping of pi components
      pi = dict(
         xx=pi[1][1],
         xy=pi[1][2],
         yy=pi[2][2],
      )

   elif dim == 3:
      # t x y z sig_t sig_x sig_y zig_z vx vy vz pi00 pi01 pi02 pi03 pi11 pi12 pi13 pi22 pi23 pi33 Pi
      # 0 1 2 3 4     5     6     7     8  9  10 11   12   13   14   15   16   17   18   19   20   21
      surface_data = np.fromfile('surface.dat', dtype='float64').reshape((-1, 22))

      # extract usual sub-arrays
      x, sigma, v, pi, _ = np.hsplit(surface_data, [4, 8, 11, 21])
      Pi = surface_data[:, 21]

      # convert viscosity to GeV
      Pi *= HBARC
      for i in range(len(pi[0])):
         pi[:,i] *= HBARC
      
      tend = x[-1][0]
      pi = np.array([[pi[:,0], pi[:,1], pi[:,2], pi[:,3]]
                    ,[pi[:,1], pi[:,4], pi[:,5], pi[:,6]]
                    ,[pi[:,2], pi[:,5], pi[:,7], pi[:,8]]
                    ,[pi[:,3], pi[:,6], pi[:,8], pi[:,9]]])

      for i in range(len(x)):
         eta = x[i][3]
         MilnetoCart_coord(x[i])
         MilnetoCart_vector(-eta, sigma[i])
         MilnetoCart_velocity(eta, v[i])
         MilnetoCart_tensor(eta, pi[:,:,i])

      # create mapping of pi components
      pi = dict(
         xx=pi[1][1],
         xy=pi[1][2],
         xz=pi[1][3],
         yy=pi[2][2],
         yz=pi[2][3],
      )

   print("tend: {}".format(tend))
   # create Surface object
   surface = frzout.Surface(x, sigma, v, ymax=ymax, pi=pi, Pi=Pi)

   hrg = frzout.HRG(.151, species='urqmd')
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

def jackknife(fn, ns, ds):
   n = sum(ns)
   d = sum(ds)
   value = fn(n, d)
   l = len(ns)
   if l > 1:
      fs = [fn(n-ni, d-di) for ni, di in zip(ns,ds)]
      fh = sum(fs)/l
      vh = sum((f-fh)**2 for f in fs)/l
      error2 = (l-1)*vh
      return value, np.sqrt(error2)
   else:
      return value, 0

def vn(n, events):
   events = [[[v for v in f if abs(v["eta"])<0.8 and 0.2 < v["pT"] and v["pT"] < 5.0 and v["charge"] != 0] for f in e] for e in events]

   def fun(n, d):
      return np.sqrt(n/d)

   mss = [[len(f) for f in e] for e in events]

   qnss = [[sum(exp(1j*n*v["phi"]) for v in f) for f in e] for e in events]
   qn2ss = [[abs(qn)**2 for qn in qns] for qns in qnss]
   
   w2ss = [[m*(m-1) for m in ms] for ms in mss]
   
   m2ss = [[(qn2-m)/w2 for qn2, m, w2 in zip(qn2s, ms, w2s)] for qn2s, ms, w2s in zip(qn2ss, mss, w2ss)]
   
   n_cn2s = [sum(w2*m2 for w2, m2 in zip(w2s, m2s)) for w2s, m2s in zip(w2ss, m2ss)]
   d_cn2s = [sum(w2s) for w2s in w2ss]

   val, err = jackknife(fun, n_cn2s, d_cn2s)

   with open("v{}.txt".format(n), "w") as fv:
       msg = "v{}: {:e} {:e}\n".format(n, val, err)
       print(msg, end="")
       fv.write(msg)
   
def dch_deta(events):
   deta = 0.5
   dchdetas = [len([v for v in f if abs(v["eta"]) <= deta and v["charge"] != 0])/(2*deta) for e in events for f in e]
   dchdeta = sum(dchdetas)/len(dchdetas)
   dchdeta2 = sum(x**2 for x in dchdetas)/len(dchdetas)
   errdchdeta = np.sqrt((dchdeta2-dchdeta**2)/len(events))
   with open("dchdeta.txt", "w") as fv:
       msg = "dchdeta: {:e} {:e}\n".format(dchdeta, errdchdeta)
       print(msg, end="")
       fv.write(msg)

def dch_deta_eta(events):
   deta = 0.25
   ran = np.arange(0, 5, deta)
   etas = [eta+deta/2 for eta in ran]
   dchdetass = [[len([v for v in f if eta <= abs(v["eta"]) and abs(v["eta"]) <= eta+deta and v["charge"] != 0])/(2*deta) for e in events for f in e] for eta in ran] # abs and /2 because it is symetrized
   dchdetas = [sum(dchdetas)/len(dchdetas) for dchdetas in dchdetass]
   dchdetas2 = [sum(x**2 for x in dchdetas)/len(dchdetas) for dchdetas in dchdetass]
   errdchdetas = [np.sqrt((d2-d**2)/len(events)) for d, d2 in zip(dchdetas,dchdetas2)]
   plt.errorbar(etas, dchdetas, yerr=errdchdetas)
   plt.ylim([0,2500])
   plt.savefig("dch_deta_eta.pdf")

totiter = 0
cwd = os.getcwd()
dirs = os.listdir(cwd)
if "results" in dirs:
   res = "{}/results".format(cwd)
   dirs = os.listdir(res)
   allevents = {}
   allinfos = {}
   def the_thing(d):
      d = "{}/{}".format(res,d)
      os.chdir(d)

      info = {k: v.strip() for [k, v] in np.loadtxt("info.txt", dtype=object, delimiter=":")}
      info.pop("case")

      if particlize:
         print(d)
         particlization(info)

      freezeout_events = []
      parts = []
      # each line: ID charge pT ET mT phi y eta
      pnames = ["ID", "charge", "pT", "ET", "mT", "phi", "y", "eta"]
      with open("particles_out.dat", "r") as f:
         for p in f.read().split("\n"):
            if "#" in p:
               if len(parts) > 0:
                  freezeout_events += [np.array(parts)]
               parts = []
            else:
               if len(p) > 0:
                  vals = {name: np.float64(i) for name, i in zip(pnames, p.split())}
                  parts += [vals]
         freezeout_events += [np.array(parts)]
      id = str(info)
      return freezeout_events, info
   with Pool(nb_threads) as p:
      for freezeout_events, info in p.map(the_thing, dirs):
         count = len(freezeout_events)
         if id in allevents:
            allevents[id] += [freezeout_events]
         else:
            allevents[id] = [freezeout_events]
            allinfos[id] = info
   os.chdir(cwd)
   for k in allevents:
      events = allevents[k]
      info = allinfos[k]
      # print(len(events), info)
      # print()

      dch_deta(events)
      dch_deta_eta(events)
      vn(2, events)
else:
   print("No \"results\" folder found.")