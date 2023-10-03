#!/usr/bin/env python3

from sympy import *

# $k^\nu = \partial_t T^{t\nu} = \partial_t ((\epsilon + p(\epsilon))u^t u^\nu - p g^{t\nu})$ 

t, k0, k1, k2, k3 = symbols("t k0 k1 k2 k3", real=True)
[e, u1, u2, u3] = [Function(n)(t) for n in "e u1 u2 u3".split()]
u0 = sqrt(1+sum([u*u for u in [u1,u2,u3]]))
p = Function("p")(e)

u = [u0, u1, u2, u3]
g = [1, 0, 0, 0]

vs = [e, u1, u2, u3]
dtvs = [diff(v,t) for v in vs]
ks = [diff((e+p)*u0*u[i]-p*g[i], t).simplify() for i in range(4)]

m = Matrix([[diff(k,dtv).simplify() for dtv in dtvs] for k in ks])
mdtvs = Matrix([[dtv] for dtv in dtvs])
mks = Matrix([[k] for k in ks])
pprint((m*mdtvs-mks).applyfunc(simplify))
# pprint(m)

dpde = diff(p,e)
et, dpdet, u0t, u1t, u2t, u3t = symbols("e dpde u0 u1 u2 u3", real=True, positive=True)
def removet(a):
    return a.subs(u1*u1+1, u0t*u0t - u2*u2 - u3*u3).simplify().subs(2*u1*u1+1, u0t*u0t + u1*u1 - u2*u2 - u3*u3).simplify().subs(dpde,dpdet).subs(u1, u1t).subs(u2, u2t).subs(u3, u3t).subs(e, et).simplify()
mm = m.applyfunc(removet)
pprint(mm)

D0 = removet(dpde*u0**2-dpde+u0**2)
DP0 = removet(dpde*u0**2+dpde+u0**2)
DU0 = removet(u0**2*D0)
DU1 = removet(-dpde*u0**2*u1**2-dpde*u1**2-u0**2*u1**2)
DU2 = removet(-dpde*u0**2*u2**2-dpde*u2**2-u0**2*u2**2)
DU3 = removet(-dpde*u0**2*u3**2-dpde*u3**2-u0**2*u3**2)
D1 = removet(DU0+DU1)
D2 = removet(DU0+DU2)
D3 = removet(DU0+DU3)
D12 = removet(DU0+DU1+DU2)
D123 = removet(DU0+DU1+DU2+DU3)
Ds = [D0, DP0, D1, D2, D3, D12, D123]
[D0e, DP0e, D1e, D2e, D3e, D12e, D123e] = [d.expand() for d in Ds]
D0s = symbols("D0", real=True)
D0ps = symbols("DP0", real=True)
D1s = symbols("D1", real=True)
D2s = symbols("D2", real=True)
D3s = symbols("D3", real=True)
D12s = symbols("D12", real=True)
D123s = symbols("D123", real=True)
for d in Ds:
    if d == D0e or d == DP0e:
        pprint(d)
    else:
        pprint(d.subs(D0e, D0s).subs(DP0e, D0ps))
m2iv = mm[0:2,0:2].inv().applyfunc(factor).subs(D1e, D1s).subs(D0e, D0s)
pprint(m2iv)
m3iv = mm[0:3,0:3].inv().applyfunc(factor).subs(D12e, D12s).subs(D1e, D1s).subs(D2e, D2s).subs(D0e, D0s).subs(DP0e, D0ps)
pprint(m3iv)
m4iv = mm.inv().applyfunc(factor).subs(D123e, D123s).subs(D1e, D1s).subs(D2e, D2s).subs(D3e, D3s).subs(D0e, D0s).subs(DP0e, D0ps)

# m2 = m.subs(u1*u1+1, u0t*u0t - u2*u2 - u3*u3).applyfunc(simplify).subs(2*u1*u1+1, u0t*u0t + u1*u1 - u2*u2 - u3*u3).applyfunc(simplify)
# pprint(m2)
# miv = m2.inv()
# pprint(miv)
