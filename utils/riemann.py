#!/usr/bin/env python3

# ref: arXiv:nucl-th/9809044v1 14 Sep 1998, 
# Dirk H. Rischke, Fluid Dynamics for Relativistic Nuclear Collisions

import sys
from math import sqrt

def riemann(e0,emin,cs2,eps,void):
    sys.setrecursionlimit(100000)
    cs = sqrt(cs2)

    def v(e):
        return (1-(e/e0)**(2*cs/(1+cs2)))/(1+(e/e0)**(2*cs/(1+cs2)))
    def xt(e):
        return (v(e) - cs)/(1 - v(e)*cs)

    def cont(el, vv):
        vel = v(el)
        vr = -vv
        vl = (vel+vr)/(1+vel*vr)
        er = emin
        gl = 1/(1-vl**2)
        gr = 1/(1-vr**2)
        c1 = el*gl*vl-er*vr*gr
        c2 = 4*el*gl*vl*vl+el-(4*er*gr*vr*vr+er)

        return sqrt(c1*c1+c2*c2)
                
    def grad(f, a, b):
        fab = f(a,b)
        da = f(a+eps,b)-fab
        db = f(a,b+eps)-fab
        d = sqrt(da*da+db*db)
        return (fab, da/d, db/d)
        
    def next(f, e, _, a, b):
        (fab, da,db) = grad(f, a, b)
        return (e*0.99, fab, a-e*da, b-e*db)

    # search the energy and speed of the chock front
    def search(err, e, v, i):
        if i<eps:
            return (e, v)
        else:
            (nerr, fab, ne, nv) = next(cont, err, 0, e, v)
            return search(nerr, ne, nv, fab)

    (eq, vc) = search(1e-1, (e0+emin)/3, 0.5, 1)

    def f(x):
        if x < xt(e0):
            return e0
        elif x < xt(eq):
            return e0*((1-cs)/(1+cs)*(1-x)/(1+x))**((1+cs2)/(2*cs))
        elif x < vc:
            return eq
        else:
            return emin

    def f2(x):
        if x < xt(e0):
            return e0
        elif x<1:
            return e0*((1-cs)/(1+cs)*(1-x)/(1+x))**((1+cs2)/(2*cs))
        else:
            return 0

    if void:
        return f2
    else:
        return f

def main() -> int:
    e0 = 10
    emin = 1
    cs2 = 1/3
    eps = 1e-10
    if len(sys.argv) == 2:
        void = (sys.argv[1] == "void")
    else:
        void = False

    fun = riemann(e0,emin,cs2,eps,void)

    for l in sys.stdin:
        print(fun(float(l)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
