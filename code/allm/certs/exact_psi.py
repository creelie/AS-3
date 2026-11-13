"""Exact Fourier coefficients of the finite parts psi_m^F of 1/Phi_10.

Method (independent of any Rademacher input):
  1/Phi_10 = (p*phi_{10,1})^{-1} * exp( sum_{N>=1} p^N Phi_N ),  Phi_N = 2 phi_{0,1}|V_N,
  phi_{10,1} = q y (1-y^{-1})^2 prod (1-q^b)^20 (1-q^b y)^2 (1-q^b/y)^2,
which is the product formula for Phi_10 grouped by the power of p.  Hence
  psi_m = F_{m+1} / phi_{10,1},   F = exp(sum p^N Phi_N) = sum F_N p^N,
with N F_N = sum_{j=1}^N (j Phi_j) F_{N-j} (all coefficients integers).
In the region |q|<|y|<1 the polar part d(m) eta^{-24} A_{2,m} has no monomial y^{-j}, 0<=j<=2m,
so the coefficient of q^n y^{-j} (0<=j<=2m) of psi_m equals c^F_m(n,j)=c^F_m(n,-j).
"""
import sys, json, math, time
from flint import fmpz_poly

class BS:
    """truncated series sum_{n=0}^{NQ} q^n y^{off_n} P_n(y)."""
    def __init__(self, NQ, rows=None):
        self.NQ = NQ
        self.rows = rows if rows is not None else [(0, fmpz_poly([]))]*(NQ+1)

def shift(p, k):
    if k == 0 or p == 0: return p
    return fmpz_poly([0]*k + [int(c) for c in p.coeffs()])

def radd(a, b):
    (oa, pa), (ob, pb) = a, b
    if pa == 0: return b
    if pb == 0: return a
    if oa <= ob: return (oa, pa + shift(pb, ob-oa))
    return (ob, pb + shift(pa, oa-ob))

def rtrunc(r, yhi):
    o, p = r
    if p == 0 or yhi is None: return r
    L = yhi - o + 1
    if L <= 0: return (0, fmpz_poly([]))
    c = p.coeffs()
    if len(c) <= L: return r
    return (o, fmpz_poly([int(x) for x in c[:L]]))

def bmul(A, B, NQ, yhi=None):
    out = [(0, fmpz_poly([]))]*(NQ+1)
    for i, (oa, pa) in enumerate(A.rows):
        if pa == 0: continue
        for j, (ob, pb) in enumerate(B.rows):
            if i+j > NQ: break
            if pb == 0: continue
            out[i+j] = radd(out[i+j], (oa+ob, pa*pb))
    if yhi is not None:
        out = [rtrunc(r, yhi) for r in out]
    return BS(NQ, out)

def phi01_coeffs(Dmax):
    """C0(D) for D=4n-l^2 <= Dmax, via phi_{0,1} = 12*phi_{-2,1}*wp/(2 pi i)^2."""
    NQ = Dmax//4 + 2
    # Pi = prod (1-q^n y)^2 (1-q^n/y)^2 (1-q^n)^-4 ; rows as dict n -> {e: c}
    import collections
    def mul(a, b):
        c = collections.defaultdict(lambda: collections.defaultdict(int))
        for i, ra in a.items():
            for j, rb in b.items():
                if i+j > NQ: continue
                for e1, c1 in ra.items():
                    for e2, c2 in rb.items():
                        c[i+j][e1+e2] += c1*c2
        return c
    P = {0: {0: 1}}
    for n in range(1, NQ+1):
        for e, pw in ((1, 2), (-1, 2)):
            f = {}
            for k in range(0, pw+1):
                if n*k > NQ: break
                f[n*k] = {e*k: (-1)**k*math.comb(pw, k)}
            P = mul(P, f)
        f = {}
        k = 0
        while n*k <= NQ:
            f[n*k] = {0: math.comb(k+3, 3)}
            k += 1
        P = mul(P, f)
    # B = (y-2+1/y)/12 + 1 + (y-2+1/y) sum_{n>=1} sum_{d|n} d (y^d+y^-d-2) q^n ; phi01 = 12*B*Pi
    Bs = {0: {1: 1, 0: 22, -1: 1}}   # 12*((y-2+1/y)/12+1) = y+10+1/y
    for n in range(1, NQ+1):
        s = collections.defaultdict(int)
        for d in range(1, n+1):
            if n % d == 0:
                s[d] += d; s[-d] += d; s[0] -= 2*d
        row = collections.defaultdict(int)
        for e, c in s.items():
            for e2, c2 in ((1, 1), (0, -2), (-1, 1)):
                row[e+e2] += 12*c*c2
        Bs[n] = dict(row)
    Bs[0] = {1: 1, 0: 10, -1: 1}
    R = mul(Bs, P)
    C0 = {}
    for n, row in R.items():
        for e, c in row.items():
            if c == 0: continue
            D = 4*n - e*e
            if D <= Dmax and n <= Dmax//4:
                if D in C0: assert C0[D] == c, (D, C0[D], c)
                C0[D] = c
    return C0

def dcoeffs(N):
    """d(n), n=-1..N, coefficients of 1/eta^24 = sum d(n) q^n."""
    a = [0]*(N+2); a[0] = 1   # index n+1
    # prod (1-q^k)^{-24}
    for k in range(1, N+2):
        for _ in range(24):
            for i in range(k, N+2):
                a[i] += a[i-k]
    return {n: a[n+1] for n in range(-1, N+1)}

def build_G(NQ, yhi):
    # E_n with offset -n
    jL = {}
    for j in range(1, NQ+1):
        c = [0]*(2*j+1)
        for k in range(1, j+1):
            if j % k == 0:
                c[j+k] += 2*(j//k); c[j-k] += 2*(j//k)
        jL[j] = fmpz_poly(c)
    E = [fmpz_poly([1])]
    for n in range(1, NQ+1):
        acc = fmpz_poly([])
        for j in range(1, n+1):
            acc += jL[j]*E[n-j]          # offsets -j + -(n-j) = -n
        c = [int(x) for x in acc.coeffs()]
        assert all(x % n == 0 for x in c)
        E.append(fmpz_poly([x//n for x in c]))
    h = [0]*(NQ+1); h[0] = 1
    for b in range(1, NQ+1):
        for _ in range(20):
            for i in range(b, NQ+1):
                h[i] += h[i-b]
    rows = []
    for n in range(0, NQ+1):
        acc = fmpz_poly([])
        for i in range(0, n+1):
            if h[i]: acc += E[n-i]*fmpz_poly([0]*(i)+[1]) * h[i] if False else shift(E[n-i], i)*h[i]
        # acc has offset -n (E_{n-i} offset -(n-i), shifted by i)
        c = [int(x) for x in acc.coeffs()]
        # multiply by y (1-y)^{-2}: offset -n+1, then two prefix sums, truncated at exponent yhi
        off = -n + 1
        L = yhi - off + 1
        c = c + [0]*(max(0, L-len(c)))
        c = c[:L]
        for _ in range(2):
            t = 0
            for i in range(L):
                t += c[i]; c[i] = t
        rows.append((off, fmpz_poly(c)))
    return BS(NQ, rows)

def compute_psiF(m, nmax, C0=None):
    """return dict (n,j) -> c^F_m(n,j) for -1<=n<=nmax, 0<=j<=2m."""
    NQ = nmax + 2
    M1 = m + 1
    emaxF = lambda n: int(math.isqrt(4*M1*n + M1*M1))
    Dneed = 4*M1*NQ
    if C0 is None: C0 = phi01_coeffs(Dneed)
    def c0(D):
        if D < -1: return 0
        return C0.get(D, 0)
    # j*Phi_j
    jPhi = {}
    for j in range(1, M1+1):
        rows = []
        for n in range(0, NQ+1):
            E = emaxF(n) if j == M1 else int(math.isqrt(4*j*n + j*j))
            coeffs = []
            for e in range(-E, E+1):
                g = math.gcd(math.gcd(j, n), abs(e))
                s = 0
                for k in range(1, g+1):
                    if g % k == 0:
                        s += (j//k)*2*c0((4*j*n - e*e)//(k*k)) if (4*j*n-e*e) % (k*k) == 0 else 0
                coeffs.append(s)
            rows.append((-E, fmpz_poly(coeffs)))
        jPhi[j] = BS(NQ, rows)
    F = {0: BS(NQ, [(0, fmpz_poly([1]))] + [(0, fmpz_poly([]))]*NQ)}
    for N in range(1, M1+1):
        acc = BS(NQ)
        for j in range(1, N+1):
            pr = bmul(jPhi[j], F[N-j], NQ)
            acc = BS(NQ, [radd(a, b) for a, b in zip(acc.rows, pr.rows)])
        rows = []
        for (o, p) in acc.rows:
            c = [int(x) for x in p.coeffs()]
            assert all(x % N == 0 for x in c)
            rows.append((o, fmpz_poly([x//N for x in c])))
        F[N] = BS(NQ, rows)
    FM = F[M1]
    # G = q*(1/phi_{10,1}) = y (1-y)^{-2} H(q) E(q,y),  H = prod (1-q^b)^{-20},
    # E = prod (1-q^b y)^{-2}(1-q^b/y)^{-2} = exp(2 sum_n q^n sum_{k|n} (y^k+y^-k)/k)
    yhi = emaxF(NQ) + 1
    G = build_G(NQ, yhi)
    P = bmul(FM, G, NQ)   # = q * psi_m
    out = {}
    for n in range(-1, nmax+1):
        o, p = P.rows[n+1]
        c = [int(x) for x in p.coeffs()]
        for j in range(0, 2*m+1):
            idx = -j - o
            out[(n, j)] = c[idx] if 0 <= idx < len(c) else 0
    return out

if __name__ == '__main__':
    t0 = time.time()
    m = int(sys.argv[1]); nmax = int(sys.argv[2])
    r = compute_psiF(m, nmax)
    print(m, nmax, time.time()-t0)
    for n in range(-1, min(nmax, 4)+1):
        print(n, [r[(n, j)] for j in range(0, m+1)])
