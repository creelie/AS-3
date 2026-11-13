"""Level-one Rademacher expansion of the finite parts psi_m^F of 1/Phi_10.

(-1)^{l+1} c^F_m(n,l) = (-1)^l e^{i pi/4} ( M + K2 + K3 ),  Delta = 4mn - l^2 > 0, with
 M  = sum_{g>=1} (2pi/g) sum_{polar (nt,lt)} c^F_m(nt,lt) Kl_g(Dt, D; lt, l) (|Dt|/D)^{23/4} I_{23/2}(pi sqrt(|Dt| D)/(m g))
 K2 = -sqrt(2m) d(m) sum_g g^{-1/2} Kl_g(-4m, D; 0, l) (4m/D)^6 I_12(2 pi sqrt(D/m)/g)
 K3 = d(m)/(2pi) sum_g g^{-2} sum_{h mod 2mg} Kl_g(-4m-h^2, D; h, l) (4m/D)^{25/4}
        int_{-1/sqrt m}^{1/sqrt m} f_{g,h}(u) (1-mu^2)^{25/4} I_{25/2}(2pi sqrt(D(1-mu^2)/m)/g) du
 Kl_g(A, D; lt, l) = sum_{0<=-delta<g, (delta,g)=1} e((alpha A + delta D)/(4 m g)) W_{g,delta}(lt,l),
 W_{g,delta}(lt,l) = (2 m g i)^{-1/2} sum_{T=0}^{g-1} e((alpha x^2 - 2 l x + delta l^2)/(4 m g)), x = lt - 2mT,
 alpha = delta^{-1} mod g in [0,g) (alpha=0 if g=1).
"""
import mpmath as mp
from math import gcd
from mpmath.calculus.quadrature import TanhSinh

def ee(x):  # e(x) for rational/mpf x
    return mp.expjpi(2*x)

def pairs(g):
    if g == 1:
        return [(0, 0)]
    out = []
    for nd in range(0, g):
        d = -nd
        if gcd(g, d) == 1:
            a = pow(d % g, -1, g)
            out.append((d, a))
    return out

class Rad:
    def __init__(self, m, polar, dm, GM=40, GK=14, dps=40, ts_degree=6):
        mp.mp.dps = dps
        self.m, self.polar, self.dm = m, polar, dm  # polar: dict (nt,lt)->value (nonzero ones)
        self.GM, self.GK = GM, GK
        self.pairs = {g: pairs(g) for g in range(1, max(GM, GK)+1)}
        self.Wcache = {}
        # tanh-sinh nodes on [0, 1/sqrt m] (integrand even in u: combine f(u)+f(-u)=2Re f(u))
        a = 1/mp.sqrt(m)
        h = mp.mpf(2)**(-ts_degree)
        self.qnodes = []
        k = 0
        while True:
            t = k*h
            for tt in ((t,) if k == 0 else (t, -t)):
                sh = mp.pi/2*mp.sinh(tt)
                x = mp.tanh(sh)
                w = mp.pi/2*mp.cosh(tt)/mp.cosh(sh)**2
                u = a*(1+x)/2
                if u > 0 and u < a:
                    self.qnodes.append((u, w*a/2))
            if mp.pi/2*mp.cosh(t)/mp.cosh(mp.pi/2*mp.sinh(t))**2 < mp.mpf(10)**(-mp.mp.dps-5):
                break
            k += 1
        self.h = h
        # precompute f values: F[g][h] = list over nodes of Re f(u)+f(-u) = 2 Re f(u)
        self.fv = {}
        for g in range(1, GK+1):
            for hh in range(0, 2*m*g):
                vals = []
                for (u, w) in self.qnodes:
                    vals.append(2*mp.re(self.f(g, hh, u)))
                self.fv[(g, hh)] = vals

    def f(self, g, hh, u):
        m = self.m
        if hh % (2*m*g) != 0:
            return mp.pi**2/mp.sinh(mp.pi*u/g - 1j*mp.pi*hh/(2*m*g))**2
        y = mp.pi*u/g
        if abs(y) < mp.mpf('1e-6'):
            return mp.pi**2*(-mp.mpf(1)/3 + y**2/15 - 2*y**4/189)
        return mp.pi**2/mp.sinh(y)**2 - mp.mpf(g)**2/u**2

    def W(self, g, d, a, lt, l):
        key = (g, d, lt % (2*self.m), l)
        if key in self.Wcache: return self.Wcache[key]
        m = self.m
        s = mp.mpc(0)
        for T in range(g):
            x = lt - 2*m*T
            s += ee(mp.mpf(a*x*x - 2*l*x + d*l*l)/(4*m*g))
        val = s*mp.expjpi(mp.mpf(-1)/4)/mp.sqrt(2*m*g)
        self.Wcache[key] = val
        return val

    def Kl(self, g, A, D, lt, l):
        m = self.m
        s = mp.mpc(0)
        for (d, a) in self.pairs[g]:
            s += ee(mp.mpf(a*A + d*D)/(4*m*g))*self.W(g, d, a, lt, l)
        return s

    def terms(self, n, l):
        m = self.m; D = 4*m*n - l*l; assert D > 0
        Mg = []; K2g = []; K3g = []
        for g in range(1, self.GM+1):
            s = mp.mpc(0)
            for (nt, lt), c in self.polar.items():
                Dt = 4*m*nt - lt*lt
                s += c*self.Kl(g, Dt, D, lt, l)*(mp.mpf(-Dt)/D)**(mp.mpf(23)/4)*mp.besseli(mp.mpf(23)/2, mp.pi*mp.sqrt(-Dt*D)/(m*g))
            Mg.append(2*mp.pi/g*s)
        for g in range(1, self.GM+1):
            K2g.append(-mp.sqrt(2*m)*self.dm*mp.mpf(g)**(-0.5)*self.Kl(g, -4*m, D, 0, l)*(mp.mpf(4*m)/D)**6*mp.besseli(12, 2*mp.pi*mp.sqrt(mp.mpf(D)/m)/g))
        for g in range(1, self.GK+1):
            kl = [self.Kl(g, -4*m - hh*hh, D, hh, l) for hh in range(2*m*g)]
            tot = mp.mpc(0)
            for i, (u, w) in enumerate(self.qnodes):
                bes = (1 - m*u*u)**(mp.mpf(25)/4)*mp.besseli(mp.mpf(25)/2, 2*mp.pi*mp.sqrt(D*(1-m*u*u)/m)/g)
                fs = mp.mpc(0)
                for hh in range(2*m*g):
                    fs += kl[hh]*self.fv[(g, hh)][i]
                tot += w*fs*bes
            tot *= self.h
            K3g.append(self.dm/(2*mp.pi)*mp.mpf(g)**(-2)*(mp.mpf(4*m)/D)**(mp.mpf(25)/4)*tot)
        return Mg, K2g, K3g

    def value(self, n, l):
        Mg, K2g, K3g = self.terms(n, l)
        tot = sum(Mg) + sum(K2g) + sum(K3g)
        return (-1)**l*mp.expjpi(mp.mpf(1)/4)*tot, (Mg, K2g, K3g)

def tail_bound(m, D, polar, dm, NM, NK):
    """Explicit bound for the tails gamma>NM (M, K2) and gamma>NK (K3) of the three series,
    multiplied by |e^{i pi/4}|=1 (see the tail proposition in expansion.tex)."""
    mp.mp.dps = 30
    TM = mp.mpf(0)
    for (nt, lt), c in polar.items():
        Dt = -(4*m*nt - lt*lt)
        TM += abs(c)*(mp.pi*Dt/(2*m))**(mp.mpf(23)/2)*mp.exp(mp.pi**2*Dt*D/(4*m*m*NM*NM))
    TM *= 2*mp.pi/mp.gamma(mp.mpf(25)/2)*mp.mpf(NM)**(-mp.mpf(21)/2)/(mp.mpf(21)/2)
    TK2 = mp.sqrt(2*m)*dm*(2*mp.pi)**12/mp.gamma(13)*mp.exp(mp.pi**2*D/(m*NM*NM))*mp.mpf(NM)**(-mp.mpf(21)/2)/(mp.mpf(21)/2)
    TK3 = dm/(2*mp.pi)*(2/mp.sqrt(m))*(mp.pi**4*m*m/2)*(2*mp.pi)**(mp.mpf(25)/2)/mp.gamma(mp.mpf(27)/2)*mp.exp(mp.pi**2*D/(m*NK*NK))*mp.mpf(NK)**(-mp.mpf(21)/2)/(mp.mpf(21)/2)
    return TM + TK2 + TK3
