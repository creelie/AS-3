"""Numerical confirmation of the exact expansion (Theorem R) against exact coefficients.
Truncation gamma<=GM in M and K2, gamma<=GK in K3; mpmath working precision."""
import sys, math, time
import mpmath as mp
from exact_psi import dcoeffs, compute_psiF
from polar_sb import polar_list, polar_sb

def e(x): return mp.expjpi(2*x)

def cosets(g):
    if g == 1: return [(0, 0)]
    return [(-nd, pow((-nd) % g, -1, g)) for nd in range(g) if math.gcd(g, nd) == 1]

def Wsum(m, g, dl, al, lt, l):
    s = mp.mpc(0)
    for T in range(g):
        x = lt - 2*m*T
        s += e(mp.mpf(al*x*x - 2*l*x + dl*l*l)/(4*m*g))
    return s/mp.sqrt(2*m*g)          # = sqrt(i) * W_{g,delta}

def Kl(m, g, A, D, lt, l):
    return sum(e(mp.mpf(al*A + dl*D)/(4*m*g))*Wsum(m, g, dl, al, lt, l) for (dl, al) in cosets(g))

def fgh(m, g, h, u):
    if h % (2*m*g):
        return mp.pi**2/mp.sinh(mp.pi*u/g - 1j*mp.pi*h/(2*m*g))**2
    y = mp.pi*u/g
    if abs(y) < mp.mpf('1e-8'): return -mp.pi**2/3
    return mp.pi**2/mp.sinh(y)**2 - mp.mpf(g)**2/u**2

def rad(m, n, l, polar, dm, GM=10, GK=4):
    D = 4*m*n - l*l
    M = mp.mpc(0); K2 = mp.mpc(0); K3 = mp.mpc(0)
    for g in range(1, GM+1):
        s = 0
        for (nt, lt), c in polar.items():
            Dt = 4*m*nt - lt*lt
            s += c*Kl(m, g, Dt, D, lt, l)*(mp.mpf(-Dt)/D)**(mp.mpf(23)/4)*mp.besseli(mp.mpf(23)/2, mp.pi*mp.sqrt(-Dt*D)/(m*g))
        M += 2*mp.pi/g*s
        K2 += -mp.sqrt(2*m)*dm*mp.mpf(g)**-0.5*Kl(m, g, -4*m, D, 0, l)*(mp.mpf(4*m)/D)**6*mp.besseli(12, 2*mp.pi*mp.sqrt(mp.mpf(D)/m)/g)
    for g in range(1, GK+1):
        for h in range(2*m*g):
            kl = Kl(m, g, -4*m - h*h, D, h, l)
            I = mp.quad(lambda u: fgh(m, g, h, u)*(1-m*u*u)**(mp.mpf(25)/4)*mp.besseli(mp.mpf(25)/2, 2*mp.pi*mp.sqrt(D*(1-m*u*u)/m)/g), [-1/mp.sqrt(m), 0, 1/mp.sqrt(m)])
            K3 += dm/(2*mp.pi)*mp.mpf(g)**-2*kl*(mp.mpf(4*m)/D)**(mp.mpf(25)/4)*I
    return (-1)**l*(M + K2 + K3), M, K2, K3

if __name__ == '__main__':
    mp.mp.dps = 30
    d = dcoeffs(400)
    cases = [(1, 1, 0), (1, 1, 1), (1, 3, 1), (1, 10, 1), (2, 2, 0), (2, 2, 1), (2, 2, 2), (2, 5, 1),
             (3, 3, 0), (3, 3, 2), (3, 3, 3), (3, 6, 1), (4, 4, 1), (4, 4, 4), (5, 5, 3), (5, 7, 0),
             (6, 6, 5), (6, 6, 6), (7, 7, 2), (8, 8, 8), (9, 9, 4), (10, 10, 7), (12, 12, 12)]
    print('m n l Delta exact (-1)^(l+1)c^F | rel.err of truncated series | |Im|/exact | M_1(leading)/exact')
    for (m, n, l) in cases:
        t0 = time.time()
        ex = compute_psiF(m, n)[(n, l)]
        val = (-1)**(l+1)*ex
        polar = {(nt, lt): polar_sb(m, nt, lt, d) for (nt, lt) in polar_list(m)}
        # extend lt to the full residue system -m<lt<=m
        full = {}
        for (nt, lt), c in polar.items():
            full[(nt, lt)] = c
            if 0 < lt < m: full[(nt, -lt)] = c
        r, M, K2, K3 = rad(m, n, l, full, d[m])
        D = 4*m*n - l*l
        lead = 2*mp.pi/mp.sqrt(2*m)*(m+2)*(mp.mpf(m*m+4*m)/D)**(mp.mpf(23)/4)*mp.besseli(mp.mpf(23)/2, mp.pi*mp.sqrt(D*(1+mp.mpf(4)/m)))
        print(m, n, l, D, val, mp.nstr(abs(r.real-val)/val, 3), mp.nstr(abs(r.imag)/val, 3), mp.nstr(lead/val, 6), f'{time.time()-t0:.1f}s', flush=True)
