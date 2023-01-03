#!/usr/bin/env python3
import sys 
import os
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
sys.path.append(os.path.abspath("%s/.config/" % (os.path.expanduser("~"))))
import plt_setting
import numpy as np
from collections import defaultdict

def dd(n):
    if n == 1:
        return {}
    else:
        return defaultdict(lambda: dd(n-1))

def convert(v):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except:
            return v
    
datas = dd(7)

dir = "results/"
for d in os.listdir(dir):
    p = dir+d+"/"+os.listdir(dir+d)[0]
    data = np.loadtxt(p+"/data.txt")
    info = {k: convert(v.strip()) for [k, v] in np.loadtxt(p+"/info.txt", dtype=object, delimiter=":")}

    t0 = info["t0"]
    tend = info["tend"]
    integration = info["integration"]
    maxdt = info["maxdt"]
    dx = info["dx"]
    nx = info["nx"]
    info["l"] = str(dx*nx)
    if info["ny"] == 1:
        dim = "1D"
    else:
        dim = "2D"
    info["dim"] = dim
    
    # print(p, dim, integration, t0, tend, dx, nx, maxdt)
    datas[dim][integration][t0][tend][dx][nx][maxdt] = (info, data)

def plot1d():
    datadts = datas["1D"]["FixPoint"][1.0][4.5][0.1][100]
    maxdts = sorted([dt for dt in datadts])
    print(maxdts)
    (info, data) = datadts[maxdts[-1]]

    x = data[:,0]
    y = data[:,2]
    plt.plot(x,y)
    plt.show()

def plot2d():
    datadts = datas["2D"]["FixPoint"][1.0][4.5][0.1][100]
    maxdts = sorted([dt for dt in datadts])
    (info, data) = datadts[maxdts[0]]

    x = np.reshape(data[:,0], (100,100))
    y = np.reshape(data[:,1], (100,100))
    z = np.reshape(data[:,2], (100,100))
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.view_init(90,0)
    # ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    print(x)
    print(y)
    l = x[0][0]
    r = x[0][-1]
    d = y[0][0]
    u = y[-1][0]
    plt.imshow(z, extent=[l,r,d,u])
    plt.show()

plot1d()
plot2d()
