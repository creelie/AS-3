"""Exact certificate for small discriminants.
(a) Exact psi_m^F (exact_psi.compute_psiF) for 1<=m<=8; d*(n,l,m)=(-1)^(l+1) c^F_m(n,l) at every reduced
    triple 0<=l<=m<=n with 0<Delta<=240 (every reduced triple with Delta<=240 has m<=8 since Delta>=3m^2).
    Checks d*>0, and computes eps = d*/E - 1 with E in ball arithmetic; compares with eps_bar of cert_common.
(b) Doubling for every reduced f0=(n0,l0,a), l0 odd, Delta'<=60, and the three index-2 sublattice forms
    f0(2x,y), f0(x,2y), f0(x+y,2y), reduced.
(c) Cross-checks: doubling over all (n',l',m), l' odd, |l'|<=m<=199, 4Delta'<=240; torsion sum d*_r>0 for
    2<=r<=15, n in r^2 Z, l in r Z, |l|<=m<=199, 0<Delta<=240."""
import time, math, sys
from flint import arb
from exact_psi import compute_psiF
from cert_common import PI, NU, eps_bar, dtable, up, fup, fdown
from polar_sb import polar_list, polar_sb

t00 = time.time()
DMAX = 240
MMAX = 8
assert 3*(MMAX+1)**2 > DMAX

def reduce_form(A, B, C):
    """SL2(Z)-reduce A x^2 + B x y + C y^2 (positive definite) to |B|<=A<=C; return (A,B,C)."""
    assert A > 0 and C > 0 and B*B - 4*A*C < 0
    while True:
        if C < A:
            A, B, C = C, -B, A
        k = (A - B) // (2*A)          # brings B into (-A, A]
        if k:
            B, C = B + 2*k*A, A*k*k + B*k + C
        if C < A: continue
        return A, B, C

def reduced_triple(n, l, m):
    """reduced triple (N, |L|, M) with |L|<=M<=N equivalent to n x^2 + l x y + m y^2."""
    A, B, C = reduce_form(m, l, n)     # m x^2 + l x y + n y^2 is SL2-equivalent to n x^2 - l x y + m y^2
    return (C, abs(B), A)

# ---------------- (a) exact values
DS = {}
psi = {}
dd = dtable(400)
for m in range(1, MMAX+1):
    t0 = time.time()
    nmax = (DMAX + m*m)//(4*m)
    c = compute_psiF(m, nmax)
    psi[m] = c
    cnt = 0
    for n in range(m, nmax+1):
        for l in range(0, m+1):
            D = 4*m*n - l*l
            if 0 < D <= DMAX:
                v = (-1)**(l+1)*c[(n, l)]
                DS[(n, l, m)] = v; cnt += 1
    print(f'm={m}: nmax={nmax}, {cnt} reduced triples (0<=l<=m<=n, 0<Delta<=240), psi^F in {time.time()-t0:.1f}s', flush=True)
# anchors
assert DS[(1, 1, 1)] == 25353 and DS[(1, 0, 1)] == 50064 and DS[(2, 0, 1)] == 1127472
assert DS[(10, 7, 10)] == 2544996616618246559487086544 if (10, 7, 10) in DS else True
neg = [k for k, v in DS.items() if v <= 0]
print(f'(a) {len(DS)} reduced triples with Delta<=240; nonpositive values: {len(neg)}')
assert not neg

def E_ball(D, m):
    t0 = m*m + 4*m
    X = PI*(arb(D)*t0).sqrt()/m
    return 2*PI/arb(2*m).sqrt()*(m+2)*(arb(t0)/D)**(arb(23)/4)*X.bessel_i(NU)

polar = {m: {k: polar_sb(m, *k, dd) for k in polar_list(m)} for m in range(1, MMAX+1)}
epsb_cache = {}
print('(a) eps = d*/E - 1 over the reduced triples with Delta<=240:')
print(' m | #triples | min eps | max eps | max|eps| (Delta<=60) | max|eps| (Delta>=61) | max |eps|/eps_bar')
worst = {}
for m in range(1, MMAX+1):
    es = []
    for (n, l, mm), v in DS.items():
        if mm != m: continue
        D = 4*m*n - l*l
        e = arb(v)/E_ball(D, m) - 1
        if (m, D) not in epsb_cache:
            epsb_cache[(m, D)] = eps_bar(m, D, dd, polar[m])[0]
        eb = epsb_cache[(m, D)]
        assert abs(e) < eb, (m, n, l, e, eb)   # sanity: the analytic bound dominates
        es.append((D, e, eb, (n, l)))
    lo = min(es, key=lambda x: x[1].mid())
    hi = max(es, key=lambda x: x[1].mid())
    a60 = max([up(abs(x[1])) for x in es if x[0] <= 60] or [0])
    a61 = max([up(abs(x[1])) for x in es if x[0] >= 61] or [0])
    rr = max(float((abs(x[1])/x[2]).mid()) for x in es)
    worst[m] = (lo, hi)
    print(f'{m:2d} | {len(es):4d} | {fdown(lo[1],6)} at (n,l)={lo[3]} D={lo[0]} | {fup(hi[1],6)} at (n,l)={hi[3]} D={hi[0]} | {fup(a60,6) if a60 else '-'} | {fup(a61,6) if a61 else '-'} | {fup(rr,4)}', flush=True)

def dstar(n, l, m):
    k = reduced_triple(n, l, m)
    return DS[k], k

# ---------------- (b) doubling, reduced f0 with Delta'<=60
print('(b) doubling d*(g) > 2 d*(f0), g in {f0(2x,y), f0(x,2y), f0(x+y,2y)}, f0 reduced, l0 odd, Delta\'<=60')
cnt = 0; minr = None; byD = {}
for a in range(1, MMAX+1):
    for n0 in range(a, 61):
        for l0 in range(-a, a+1):
            if l0 % 2 == 0: continue
            Dp = 4*a*n0 - l0*l0
            if not (0 < Dp <= 60): continue
            assert Dp % 4 == 3
            v0, k0 = dstar(n0, l0, a)
            for (gn, gl, gm) in ((4*n0, 2*l0, a), (n0, 2*l0, 4*a), (n0, 2*n0 + 2*l0, n0 + 2*l0 + 4*a)):
                assert 4*gm*gn - gl*gl == 4*Dp
                vg, kg = dstar(gn, gl, gm)
                assert a <= kg[2] <= 4*a
                r = arb(vg)/(2*v0)
                assert r > 1
                cnt += 1
                if minr is None or r < minr[0]: minr = (r, (n0, l0, a), kg)
                if Dp not in byD or r < byD[Dp]: byD[Dp] = r
print(f'    {cnt} pairs (f0,g) checked; minimal ratio d*(g)/(2d*(f0)) = {fdown(minr[0],4)} at f0={minr[1]}, g~{minr[2]}')
print('    minimal ratio by Delta\': ' + ', '.join(f'{D}:{fdown(r,1)}' for D, r in sorted(byD.items())))

# ---------------- (c) cross-checks with the published counts
cnt = 0; minr = None
for m in range(1, 200):
    for l in range(-m, m+1):
        if l % 2 == 0: continue
        n = (l*l)//(4*m) + 1
        while 4*(4*m*n - l*l) <= DMAX:
            if 4*m*n - l*l > 0:
                v1, _ = dstar(n, l, m); v2, _ = dstar(4*n, 2*l, m)
                r = arb(v2)/(2*v1); assert r > 1; cnt += 1
                if minr is None or r < minr[0]: minr = (r, (n, l, m))
            n += 1
print(f'(c) doubling over (n\',l\',m), l\' odd, |l\'|<=m<=199, 4Delta\'<=240: {cnt} triples, min ratio {fdown(minr[0],2)} at {minr[1]}')
cnt = 0; minv = None
for r in range(2, 16):
    divs = [s for s in range(1, r+1) if r % s == 0]
    for m in range(1, 200):
        for l in range(-m, m+1):
            if l % r: continue
            n = ((l*l)//(4*m) // (r*r)) * r*r
            while 4*m*n - l*l <= DMAX:
                D = 4*m*n - l*l
                if D > 0:
                    tot = 0
                    for s in divs:
                        v, _ = dstar(n//(s*s), l//s, m)
                        tot += s*(-1)**((l - l//s) % 2)*v
                    assert tot > 0, (r, n, l, m, tot)
                    cnt += 1
                    q = arb(tot)
                    if minv is None or tot < minv[0]: minv = (tot, (r, n, l, m))
                n += r*r
print(f'(c) torsion sums d*_r>0 for 2<=r<=15, n in r^2Z, l in rZ, |l|<=m<=199, 0<Delta<=240: {cnt} triples, all positive; smallest {minv[0]} at (r,n,l,m)={minv[1]}')
print(f'total time {time.time()-t00:.1f}s')
