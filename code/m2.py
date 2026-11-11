import sys; sys.path.insert(0,'.')
from series import *
from fractions import Fraction
def m2_sequences(N):
    """returns dict L[j][k] = (2/3)*(-1)^{j+1}*[ (xi_j - 4800 F_j)/Delta ] at u-exponent k (u=q^{1/8}),
    as Fractions; d*(n,l,2) = L[l mod 4][8n - l^2]."""
    U = 8*N + 64
    def up(p, s=8):  # series in q -> series in u
        c = p.coeffs(); out = [0]*(len(c)*s)
        for i, v in enumerate(c): out[i*s] = v
        return trunc(P(out), U)
    def ser(d):  # dict exp->coef in u
        c = [0]*(U+1)
        for e, v in d.items():
            if 0 <= e <= U: c[e] += v
        return P(c)
    Nq = U//8 + 2
    E2, E4, E6 = up(eis(2, Nq)), up(eis(4, Nq)), up(eis(6, Nq))
    th0 = ser({8*n*n: (1 if n == 0 else 2) for n in range(0, isqrt(U//8)+2)})
    th1 = ser({2*r*r: 2 for r in range(1, isqrt(U//2)+2, 2)})
    th0p = ser({8*n*n: 2*n*n for n in range(1, isqrt(U//8)+2)})       # q d/dq
    # th1' = sum (r^2/4) q^{r^2/4}: 264*th1' = sum 2*66 r^2 u^{2r^2}
    th1p264 = ser({2*r*r: 2*66*r*r for r in range(1, isqrt(U//2)+2, 2)})
    th0_2 = ser({16*n*n: (1 if n == 0 else 2) for n in range(0, isqrt(U//16)+2)})
    th1_2 = ser({4*r*r: 2 for r in range(1, isqrt(U//4)+2, 2)})
    th1_h = ser({r*r: 2 for r in range(1, isqrt(U)+2, 2)})
    M = lambda a, b: mul(a, b, U)
    g = P([5])*E6 - P([11])*M(E2, E4)
    f0 = M(th1p264, E4) + M(g, th1)
    f1 = P([264])*M(th0p, E4) + M(g, th0)
    # xi * eta^12 (times 2 for xi1 to avoid 1/2)
    X0 = -(M(M(th0_2, th1), f0) + M(M(th1_2, th0), f1))
    X1x2 = M(th1_h, M(f0, th0) + M(f1, th1))
    X2 = -(M(M(th1_2, th1), f0) + M(M(th0_2, th0), f1))
    # 1/(eta^12 Delta) = u^{-12} prod (1-u^{8n})^{-36}
    i36 = up(qpoch_pow(-36, Nq))
    i24 = up(qpoch_pow(-24, Nq))
    H12 = hurwitz(U + 8)
    def Hh(n): return Fraction(H12[n], 12) if n >= 0 else Fraction(0)
    def cnr(n, r):
        if n < 0 or 8*n - r*r < 0: return Fraction(0)
        s = Hh(8*n - r*r)
        if n % 2 == 0 and r % 2 == 0 and (8*n - r*r) % 4 == 0: s += 2*Hh((8*n - r*r)//4)
        return s
    L = {}
    for j, (X, half) in ((0, (X0, 1)), (1, (X1x2, 2)), (2, (X2, 1))):
        a = M(X, i36).coeffs() + [0]*(U+1)     # (xi_j eta^12)/(eta^12 Delta) at u^{k+12}  (times half)
        # F_j = sum_{D>=0} c((D+j^2)/8, j) q^{D/8} = u^D ; F_j/Delta = u^{-8} F_j * prod^{-24}
        Fd = {}
        for D in range(0, U+1):
            if (D + j*j) % 8 == 0:
                v = cnr((D + j*j)//8, j)
                if v: Fd[D] = v
        i24c = i24.coeffs() + [0]*(U+1)
        out = {}
        sgn = (-1)**(j+1)
        for k in range(-12, U - 64):
            xi = Fraction(int(a[k+12]), half) if k + 12 >= 0 else 0
            fF = sum(v*int(i24c[k+8-D]) for D, v in Fd.items() if 0 <= k+8-D) if True else 0
            out[k] = Fraction(2, 3)*sgn*(xi - 4800*fF)
        L[j] = out
    return L
if __name__ == '__main__':
    L = m2_sequences(6)
    for l in (0, 1, 2):
        print(l, [L[l % 4][8*n - l*l] for n in range(-1, 6) if 8*n - l*l in L[l%4]])
