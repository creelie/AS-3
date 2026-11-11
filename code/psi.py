# Fourier coefficients c^F_m(n,l) of the finite parts psi^F_m of 1/Phi10 (DMZ)
import sys, pickle
from math import isqrt
import flint

def c0_table(Dmax):
    """coefficients of 2*phi_{0,1}, as function of D=4n-l^2 (D>=-1)."""
    # phi_{0,1} = 4 sum_{i=2,3,4} theta_i(z)^2/theta_i(0)^2 ; compute via phi_{0,1}= (theta1^2/eta^6)*... use
    # phi_{0,1} = 12 phi_{-2,1} * wp-type? simpler: phi_{0,1} = 4 sum theta_i^2(z)/theta_i^2(0)
    X = Dmax//4 + 3          # q-degree needed
    W = 4*X + 4              # in w = q^{1/4}
    R = isqrt(W) + 3
    def th(kind):            # dict (w-exp, y-exp) -> coeff  of theta_kind(tau,z), q^{1/8} units handled by w=q^(1/8)? use q^(1/8)
        pass
    # Use exact formula via Jacobi theta in variable q^(1/8): theta2 = sum q^{(r+1/2)^2/2} y^{r+1/2}, theta3 = sum q^{r^2/2} y^r, theta4 = sum (-1)^r q^{r^2/2} y^r
    # Work with u = q^(1/8): theta2: exponent 4(r+1/2)^2 = (2r+1)^2 ; y exponent (2r+1)/2 -> store 2*yexp
    from collections import defaultdict
    U = 8*X + 8
    def theta(kind, withz=True):
        d = defaultdict(int)
        for r in range(-R*2, R*2+1):
            if kind == 2:
                e = (2*r+1)**2; ye = (2*r+1) if withz else 0; s = 1
            else:
                e = 4*r*r; ye = 2*r if withz else 0; s = (-1 if r % 2 else 1) if kind == 4 else 1
            if e <= U: d[(e, ye)] += s
        return d
    def mul(A, B):
        C = defaultdict(int)
        for (e1, y1), v1 in A.items():
            for (e2, y2), v2 in B.items():
                if e1+e2 <= U: C[(e1+e2, y1+y2)] += v1*v2
        return C
    from fractions import Fraction as Fr
    total = defaultdict(Fr)
    for kind in (2, 3, 4):
        num = mul(theta(kind), theta(kind))
        den = mul(theta(kind, False), theta(kind, False))
        # den is series in u with y^0; invert
        dl = [0]*(U+1)
        for (e, _), v in den.items(): dl[e] += v
        e0 = min(i for i in range(U+1) if dl[i]); lead = dl[e0]
        inv = [Fr(0)]*(U+1); inv[0] = Fr(1, lead)
        for n in range(1, U+1):
            sacc = Fr(0)
            for k in range(1, n+1):
                if e0+k <= U and dl[e0+k]: sacc += dl[e0+k]*inv[n-k]
            inv[n] = -sacc/lead
        for (e, y), v in num.items():
            for k in range(0, U+1-e+e0):
                if inv[k]: total[(e-e0+k, y)] += 4*v*inv[k]
    c0 = {}
    for (e, y), v in total.items():
        if v == 0: continue
        assert e % 8 == 0 and y % 2 == 0, (e, y)
        n = e//8; l = y//2; D = 4*n - l*l
        if n > X - 2: continue
        assert v.denominator == 1
        if D in c0: assert c0[D] == 2*v, (D, c0[D], v)
        else: c0[D] = 2*int(v)
    return c0

if __name__ == '__main__':
    c0 = c0_table(40)
    print(sorted(c0.items())[:12])
