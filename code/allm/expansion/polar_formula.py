"""Wall (continued-fraction) formula for the polar coefficients of psi_m^F:
 c^F_m(nt, lt) = sum_{G in W(m,nt,lt)} l_G d(m_G) d(n_G),
 W = {G=(p q; r s) in SL2(Z): r,s>=1, q/s <= lt/2m < p/r},
 m_G = r^2 nt + p^2 m - p r lt, n_G = s^2 nt + q^2 m - q s lt, l_G = -2 r s nt - 2 p q m + (ps+qr) lt."""
from fractions import Fraction
from math import floor
from exact_common import dcoef

def walls(m, nt, lt, RS):
    w0 = Fraction(lt, 2*m)
    out = []
    for r in range(1, RS+1):
        for s in range(1, RS+1):
            q = floor(w0*s)
            if not (Fraction(1+q*r, r*s) > w0):
                continue
            if (1+q*r) % s: continue
            p = (1+q*r)//s
            mG = r*r*nt + p*p*m - p*r*lt
            nG = s*s*nt + q*q*m - q*s*lt
            lG = -2*r*s*nt - 2*p*q*m + (p*s+q*r)*lt
            out.append((p, q, r, s, mG, nG, lG))
    return out

def polar_value(m, nt, lt, RS=None, with_walls=False):
    RS = RS or (2*m + 2)
    D = dcoef(4*m*m + 10*m + 10)
    tot = 0; used = []
    for (p, q, r, s, mG, nG, lG) in walls(m, nt, lt, RS):
        if mG >= -1 and nG >= -1:
            tot += lG*D[mG]*D[nG]; used.append((p, q, r, s, mG, nG, lG))
    return (tot, used) if with_walls else tot
