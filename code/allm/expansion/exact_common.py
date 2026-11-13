"""Exact Fourier coefficients of the Fourier--Jacobi coefficients psi_m of 1/Phi_10,
computed from the Borcherds/Gritsenko--Nikulin product formula
   Phi_10 = q y p prod_{(n,l,m)>0} (1 - q^n y^l p^m)^{c(4nm-l^2)},  c = coeffs of 2 phi_{0,1},
i.e. 1/Phi_10 = p^{-1} (Delta(tau) phi_{-2,1}(tau,z))^{-1} exp(sum_{l>=1} p^l (2phi_{0,1})|V_l).
Bivariate truncated series are dicts  {(n, r): int}  (power of q, power of y)."""
from fractions import Fraction
from functools import lru_cache
import flint

def qseries_eta_power(power, N):
    """coefficients a_0..a_N of prod_{n>=1}(1-q^n)^power, via a(n) = (-power/n) sum sigma(j) a(n-j)."""
    sig = [0]*(N+1)
    for dd in range(1, N+1):
        for k in range(dd, N+1, dd):
            sig[k] += dd
    a = [1] + [0]*N
    for n in range(1, N+1):
        s = sum(sig[j]*a[n-j] for j in range(1, n+1))
        v = -power*s
        assert v % n == 0
        a[n] = v//n
    return a

@lru_cache(None)
def dcoef(N):
    """d(-1..N): 1/Delta = sum_{n>=-1} d(n) q^n ; returns dict n->d(n)."""
    a = qseries_eta_power(-24, N+1)
    return {n: a[n+1] for n in range(-1, N+1)}

def bmul(A, B, qmax):
    C = {}
    for (n1, r1), a in A.items():
        if a == 0: continue
        for (n2, r2), b in B.items():
            n = n1+n2
            if n > qmax: continue
            k = (n, r1+r2)
            C[k] = C.get(k, 0) + a*b
    return {k: v for k, v in C.items() if v != 0}

def phi_m21_parts(qmax):
    """F(q,y) = prod (1-q^n y)^2 (1-q^n/y)^2/(1-q^n)^4, so phi_{-2,1} = (y-2+1/y) F."""
    F = {(0, 0): 1}
    for n in range(1, qmax+1):
        fac = {(0, 0): 1, (n, 1): -2, (2*n, 2): 1}
        F = bmul(F, fac, qmax)
        fac = {(0, 0): 1, (n, -1): -2, (2*n, -2): 1}
        F = bmul(F, fac, qmax)
    e = qseries_eta_power(-4, qmax)
    E = {(n, 0): e[n] for n in range(qmax+1) if e[n]}
    return bmul(F, E, qmax)

@lru_cache(None)
def phi01_c(Dmax):
    """c0(D) = coefficient of phi_{0,1} at discriminant D=4n-r^2, for -1<=D<=Dmax.
    Uses phi_{0,1} = 12 phi_{-2,1} wp/(2 pi i)^2 with
    wp/(2pi i)^2 = 1/12 + y/(1-y)^2 + sum_{n>=1} sum_{d|n} d (y^d-2+y^{-d}) q^n,
    so phi_{0,1} = F*[ (y-2+1/y) + 12 + 12 (y-2+1/y) S ]."""
    qmax = Dmax//4 + 2
    F = phi_m21_parts(qmax)
    S = {}
    for n in range(1, qmax+1):
        for d in range(1, n+1):
            if n % d == 0:
                for (r, c) in ((d, d), (0, -2*d), (-d, d)):
                    S[(n, r)] = S.get((n, r), 0) + c
    P = {(0, 1): 1, (0, 0): -2, (0, -1): 1}
    inner = dict(P)
    inner[(0, 0)] = inner.get((0, 0), 0) + 12
    for k, v in bmul(P, S, qmax).items():
        inner[k] = inner.get(k, 0) + 12*v
    phi = bmul(F, inner, qmax)
    c = {}
    for (n, r), v in phi.items():
        D = 4*n - r*r
        if D <= Dmax:
            if D in c and c[D] != v:
                raise ValueError("phi01 not a function of D")
            c[D] = v
    for D in range(-1, Dmax+1):
        if D % 4 in (0, 3):
            c.setdefault(D, 0)
    return c

def cc(D, table):
    """coefficient c(D) of 2 phi_{0,1}"""
    if D < -1: return 0
    if D % 4 in (1, 2): return 0
    return 2*table[D]

def jL(j, qmax, table):
    """the integral series j*L_j, L_j = (2phi_{0,1})|V_j  (weight 0):
    coefficient at q^n y^r = sum_{a | gcd(n,r,j)} (j/a) c((n j/a^2)*4 ... ) """
    out = {}
    for n in range(0, qmax+1):
        rmax = int((4*n*j + j*j)**0.5) + 1
        for r in range(-rmax, rmax+1):
            s = 0
            for a in range(1, j+1):
                if j % a == 0 and n % a == 0 and r % a == 0:
                    s += (j//a) * cc(4*(n//a)*(j//a) - (r//a)**2, table)
            if s:
                out[(n, r)] = s
    return out

def ek_list(K, qmax):
    """e_0..e_K with exp(sum L_l X^l) = sum e_k X^k (integer coefficients)."""
    table = phi01_c(4*qmax*K + 4)
    JL = [None] + [jL(j, qmax, table) for j in range(1, K+1)]
    e = [{(0, 0): 1}]
    for k in range(1, K+1):
        acc = {}
        for j in range(1, k+1):
            for kk, v in bmul(JL[j], e[k-j], qmax).items():
                acc[kk] = acc.get(kk, 0) + v
        ek = {}
        for kk, v in acc.items():
            if v % k: raise ValueError("non-integral e_k")
            if v: ek[kk] = v//k
        e.append(ek)
    return e

def finite_part_coeffs(m, nmax):
    """c^F_m(n,l) for -1<=n<=nmax, 0<=l<=2m, read off from the expansion of psi_m in the
    region 1<|y|<|q|^{-1} (where the Appell--Lerch part has no y^l, 0<=l<=2m)."""
    qmax = nmax + 2
    e = ek_list(m+1, qmax)[m+1]
    F = phi_m21_parts(qmax)
    # 1/phi_{-2,1} = [sum_{k>=1} k y^{-k}] * 1/F ; 1/F = prod (1-q^n)^4/((1-q^ny)^2(1-q^n/y)^2)
    invF = {(0, 0): 1}
    for n in range(1, qmax+1):
        g1 = {}
        for k in range(0, qmax//n + 1):
            g1[(n*k, k)] = k+1
        g2 = {}
        for k in range(0, qmax//n + 1):
            g2[(n*k, -k)] = k+1
        invF = bmul(invF, g1, qmax)
        invF = bmul(invF, g2, qmax)
    e4 = qseries_eta_power(4, qmax)
    invF = bmul(invF, {(n, 0): e4[n] for n in range(qmax+1) if e4[n]}, qmax)
    H = bmul(e, invF, qmax)
    # c'(n,l) = sum_{k>=1} k [q^n y^{l+k}] H
    cp = {}
    for (n, r), v in H.items():
        for l in range(0, 2*m+1):
            k = r - l
            if k >= 1:
                cp[(n, l)] = cp.get((n, l), 0) + k*v
    d = dcoef(qmax+1)
    out = {}
    for n in range(-1, nmax+1):
        for l in range(0, 2*m+1):
            s = 0
            for j in range(-1, n+1):
                s += d[j]*cp.get((n-j, l), 0)
            out[(n, l)] = s
    return out
