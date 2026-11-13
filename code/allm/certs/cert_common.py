"""Common routines for the interval certificates (python-flint arb; every real is a ball)."""
import math
from flint import arb, ctx
ctx.prec = 160
from polar_sb import polar_list, polar_sb

def dtable(N):
    """exact d(n), n=-1..N:  p_24(n+1), via  k a_k = 24 sum_{j=1}^k sigma(j) a_{k-j}."""
    sig = [0]*(N+3)
    for i in range(1, N+3):
        for j in range(i, N+3, i): sig[j] += i
    a = [1] + [0]*(N+1)
    for k in range(1, N+2):
        s = 0
        for j in range(1, k+1): s += sig[j]*a[k-j]
        assert (24*s) % k == 0
        a[k] = 24*s//k
    return {n: a[n+1] for n in range(-1, N+1)}

PI = arb.pi()
SQ2 = arb(2).sqrt()
NU = arb(23)/2
def zeta11m1():
    # zeta(11)-1 = sum_{g>=2} g^-11 <= sum_{g=2}^{50} g^-11 + int_{50}^inf x^-11 dx
    s = sum(arb(1)/arb(g)**11 for g in range(2, 51)) + arb(1)/(10*arb(50)**10)
    return s
Z = zeta11m1()
CG = 1 + SQ2*Z           # factor for gamma>=2 in the mock terms

def eps_bar(m, D, d, polar=None):
    """upper bound (arb) for |eps| at magnetic charge m, discriminant >= D (decreasing in D).
    returns (total, eps_P, eps_gamma, eps_K2, eps_K3)."""
    if polar is None:
        polar = {k: polar_sb(m, *k, d) for k in polar_list(m)}
    D = arb(D)
    t0 = m*m + 4*m
    X = PI*(D*t0).sqrt()/m
    IX = X.bessel_i(NU)
    eP = arb(0)
    for (nt, lt), c in polar.items():
        if (nt, lt) == (-1, m): continue
        t = lt*lt - 4*m*nt
        w = 2 if 0 < lt < m else 1
        r = arb(t)/t0
        eP += w*arb(c)/(m+2)*r**(arb(23)/4)*(X*r.sqrt()).bessel_i(NU)/IX
    eg = SQ2*Z*(1+eP)
    E = 2*PI/(arb(2*m).sqrt())*(m+2)*(arb(t0)/D)**(arb(23)/4)*IX
    y = 2*PI*(D/m).sqrt()
    K2 = CG*d[m]*(arb(4*m)/D)**6*y.bessel_i(arb(12))
    mn = 2*PI/y
    mn = mn.sqrt()
    fac = arb(2) if (mn > 2) else (mn if mn < 2 else arb(2))   # min(2, sqrt(2 pi/y)); ball compare
    if not (mn < 2 or mn > 2): fac = arb(2)
    K3 = CG*d[m]*SQ2*PI*m/3*(arb(4*m)/D)**(arb(25)/4)*y.bessel_i(arb(25)/2)*fac
    tot = eP + eg + K2/E + K3/E
    return tot, eP, eg, K2/E, K3/E

def up(x):
    """rigorous upper bound of a ball as a float (rounded up)."""
    return float(x.mid() + x.rad()) * (1 + 1e-15) + 1e-300

# ---- outward rounding of printed bounds ----------------------------------------------
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
getcontext().prec = 80
def ub(x):
    """float upper bound of a ball (or of a float that is already an upper bound)."""
    v = float(x.mid() + x.rad()) if isinstance(x, arb) else float(x)
    return v + abs(v)*1e-15 + 1e-300
def lb(x):
    """float lower bound of a ball (or of a float that is already a lower bound)."""
    v = float(x.mid() - x.rad()) if isinstance(x, arb) else float(x)
    return v - abs(v)*1e-15 - 1e-300
def fup(x, nd=6):
    """decimal string >= x with nd decimals (rounded upward)."""
    return str(Decimal(ub(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_CEILING))
def fdown(x, nd=6):
    """decimal string <= x with nd decimals (rounded downward)."""
    return str(Decimal(lb(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_FLOOR))
def fupe(x, sig=3):
    """scientific string >= x > 0 with sig significant digits (rounded upward)."""
    v = Decimal(ub(x)); e = v.adjusted()
    q = v.scaleb(-e).quantize(Decimal(1).scaleb(-(sig-1)), rounding=ROUND_CEILING)
    if q >= 10: q = (q/10).quantize(Decimal(1).scaleb(-(sig-1)), rounding=ROUND_CEILING); e += 1
    return f'{q}e{e:+03d}'
