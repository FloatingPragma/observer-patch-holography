# Independent check of the pixel closure under two readings (coordinator, 2026-11-11).
import numpy as np
from scipy.optimize import brentq
pi = np.pi
def ell2(t, N=120):
    n = np.arange(0, N+1); j = n/2.0; d = 2*j+1; C = j*(j+1)
    w = d*np.exp(-t*C); return float((w*np.log(d)).sum()/w.sum())
def ell3(t, N=60):
    p, q = np.meshgrid(np.arange(N+1), np.arange(N+1))
    d = (p+1)*(q+1)*(p+q+2)/2.0; C = (p*p+q*q+p*q+3*p+3*q)/3.0
    w = d*np.exp(-t*C); return float((w*np.log(d)).sum()/w.sum())
EP = 1.220890e19
def chain(P, aU, b=(33/5, 1.0, -3.0), beta=4.0, read='mZ'):
    MU = EP*np.exp(-2*pi)*P**(1/6); Ecell = EP/np.sqrt(P)
    v = Ecell*np.exp(-2*pi/(beta*aU))
    def al(mu):
        return [1.0/(1.0/aU + bi/(2*pi)*np.log(MU/mu)) for bi in b]
    mu = 91.0
    for _ in range(200):
        a1, a2, a3 = al(mu); aY = 0.6*a1
        mz = v*np.sqrt(pi*(a2+aY)); mw = v*np.sqrt(pi*a2)
        new = mz if read == 'mZ' else mw
        if abs(new-mu) < 1e-12*mu: break
        mu = new
    a1, a2, a3 = al(mu); aY = 0.6*a1
    lhs = ell2(4*pi*pi*a2) + ell3(4*pi*pi*a3)
    return dict(MU=MU, Ecell=Ecell, v=v, mu=mu, mz=v*np.sqrt(pi*(a2+aY)), mw=v*np.sqrt(pi*a2), a1=a1, a2=a2, a3=a3,
                s2=aY/(aY+a2), aem_inv=1.0/(a2*aY/(a2+aY)), lhs=lhs)
def solve(P, **kw):
    f = lambda x: chain(P, 1.0/x, **kw)['lhs'] - P/4
    return brentq(f, 18.0, 40.0)
P = 1.6309721
# (a) cut reading: unified coupling at the cell, both factors at the same t
g = lambda x: ell2(4*pi*pi/x) + ell3(4*pi*pi/x) - P/4
xc = brentq(g, 5.0, 60.0)
print('cut reading: 1/alpha =', round(xc, 4), ' t =', round(4*pi*pi/xc, 4), ' v/E* =', np.exp(-2*pi*xc/4)/np.sqrt(P))
# (b) paper reading: couplings at the electroweak anchor, solve alpha_U at M_U
for read in ('mZ', 'mW'):
    x = solve(P, read=read); c = chain(P, 1.0/x, read=read)
    print(f"paper reading at {read}: 1/alpha_U = {x:.4f}  M_U = {c['MU']:.4e}  v = {c['v']:.3f} GeV  v/E* = {c['v']/EP:.4e}  mZ = {c['mz']:.3f}  mW = {c['mw']:.3f}")
    print(f"    1/a1 = {1/c['a1']:.3f} 1/a2 = {1/c['a2']:.3f} 1/a3 = {1/c['a3']:.3f}  sin2 = {c['s2']:.5f}  1/a_em(mu) = {c['aem_inv']:.3f}  alpha_s = {c['a3']:.5f}")
    print(f"    ell2 = {ell2(4*pi*pi*c['a2']):.5f} ell3 = {ell3(4*pi*pi*c['a3']):.5f} P/4 = {P/4:.5f}")
# closure with measured couplings at mZ (NOT BLIND)
a2m, a3m = 1/29.58, 0.1180
print('measured couplings at mZ (NOT BLIND): lhs =', ell2(4*pi*pi*a2m)+ell3(4*pi*pi*a3m), ' P/4 =', P/4)
# left side along the chain's running, as a function of scale
x = solve(P); c = chain(P, 1.0/x)
for mu in (c['Ecell'], c['MU'], 1e12, 1e6, 1e3, c['mu'], 30.0):
    L = np.log(c['MU']/mu); a2 = 1/(x + 1.0/(2*pi)*L); a3 = 1/(x - 3.0/(2*pi)*L)
    print(f"  mu = {mu:.3e}  1/a2 = {1/a2:.3f} 1/a3 = {1/a3:.3f}  lhs = {ell2(4*pi*pi*a2)+ell3(4*pi*pi*a3):.5f}")
