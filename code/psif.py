import sys, pickle, time
from math import isqrt
import flint
from psi import c0_table

LRANGE = lambda m, lmax: range(-lmax, lmax+1)
def compute(M, NQ):
    """return dict cF[m][(n,l)] for 1<=m<=M, -1<=n<=NQ, all l, of psi^F_m."""
    c0 = c0_table(4*(M+1)*NQ + 8)
    Y = 2*(M+1) + NQ + 2          # bound on |y-exponent| of F_n
    W = 4*Y + 3
    def enc(d):                    # d: dict (qdeg,c)->int, offset Y
        L = [0]*((NQ+1)*W)
        for (b, c), v in d.items():
            if b <= NQ and v:
                assert abs(c) <= Y, (b, c)
                L[b*W + c + Y] = v
        return flint.fmpz_poly(L)
    def dec(P, off):               # off = total offset in product
        co = P.coeffs(); d = {}
        for i, v in enumerate(co):
            if v:
                b, r = divmod(i, W); 
                if b <= NQ: d[(b, r - off)] = int(v)
        return d
    # G_j = j * L_j,  L_a = sum_{k|a} (1/k) sum_{b,c} c0(4(a/k) b - c^2) q^{bk} y^{ck}
    G = [None]
    for j in range(1, M+2):
        d = {}
        for k in range(1, j+1):
            if j % k: continue
            a = j//k
            for b in range(0, NQ//k + 1):
                cm = isqrt(4*a*b + 1)
                for c in range(-cm, cm+1):
                    D = 4*a*b - c*c
                    if D < -1: continue
                    e = c0.get(D, 0)
                    if e: d[(b*k, c*k)] = d.get((b*k, c*k), 0) + (j//k)*e
        G.append(d)
    Genc = [None] + [enc(g) for g in G[1:]]
    F = [{(0, 0): 1}]
    for n in range(1, M+2):
        acc = {}
        for j in range(1, n+1):
            prod = dec(Genc[j]*enc(F[n-j]), 2*Y)
            for k, v in prod.items(): acc[k] = acc.get(k, 0) + v
        Fn = {}
        for k, v in acc.items():
            if v:
                assert v % n == 0
                Fn[k] = v//n
        F.append(Fn)
    # P2 = q^{-1} prod (1-q^n)^{-24} * prod 1/((1-q^n y)^2 (1-q^n/y)^2); store qdeg shifted by +1
    P2 = {(0, 0): 1}
    def mulser(A, B):
        C = {}
        for (b1, c1), v1 in A.items():
            for (b2, c2), v2 in B.items():
                if b1+b2 <= NQ+1:
                    C[(b1+b2, c1+c2)] = C.get((b1+b2, c1+c2), 0) + v1*v2
        return C
    for n in range(1, NQ+2):
        for (sgn, cnt) in ((0, 20), (1, 2), (-1, 2)):
            for _ in range(cnt):
                geo = {(n*t, sgn*t): 1 for t in range(0, (NQ+1)//n + 1)}
                P2 = mulser(P2, geo)
    # p24
    p24 = {}
    for (b, c), v in P2.items(): pass
    eta = {0: 1}
    EL = max(M, NQ) + 5; E = [0]*EL; E[0] = 1
    for n in range(1, EL):
        for _ in range(24):
            for k in range(n, EL): E[k] += E[k-n]
    # E[k] = coefficient of q^{k-1} in 1/eta^24 i.e. p24(k)
    out = {}
    for m in range(1, M+1):
        # psi_m = P2 * (y/(1-y)^2) * F_{m+1}  (q-degree of P2 shifted by 1)
        H = mulser(P2, F[m+1])        # qdeg shifted +1
        # apply y/(1-y)^2: coefficient at l: sum_{c<l} (l-c) H_c
        lmax = m + 2*NQ + 10
        psi = {}
        byq = {}
        for (b, c), v in H.items():
            byq.setdefault(b, {})[c] = v
        for b, row in byq.items():
            n = b - 1
            if n > NQ: continue
            for l in LRANGE(m, lmax):
                s = 0
                for c, v in row.items():
                    if c < l: s += (l - c)*v
                if s: psi[(n, l)] = s
        # polar part p24(m+1) * (1/eta^24) * A_{2,m}
        A = {}
        def add(b, c, v):
            if b <= NQ+1: A[(b, c)] = A.get((b, c), 0) + v
        # s = 0: sum_{k>=1} k y^k   (we keep l within lmax)
        for k in range(1, lmax+1): add(0, k, k)
        s = 1
        while m*s*s + s <= NQ+1:
            base = m*s*s + s
            for k in range(0, (NQ+1-base)//s + 1):
                add(base + s*k, 2*m*s + 1 + k, k+1)
            s += 1
        s = -1
        while m*s*s + s <= NQ+1:
            base = m*s*s + s  # q^{base} y^{2ms+1} * sum_j (j+1) X^{-j-2}, X^{-1} = q^{|s|} y^{-1}
            for j in range(0, NQ+2):
                b = base + (-s)*(j+2); c = 2*m*s + 1 - (j+2)
                if b > NQ+1: break
                add(b, c, j+1)
            s -= 1
        # multiply by 1/eta^24 (with q^{-1}): coefficients E[k] at q^{k-1}
        Pm = {}
        for (b, c), v in A.items():
            for k in range(0, NQ+2-b+1):
                if b + k - 1 <= NQ:
                    Pm[(b+k-1, c)] = Pm.get((b+k-1, c), 0) + v*E[k]*E[m+1]
        cF = {}
        for key in set(psi) | set(Pm):
            n, l = key
            if abs(l) > lmax - 5: continue
            v = psi.get(key, 0) - Pm.get(key, 0)
            if v: cF[key] = v
        out[m] = cF
    return out

if __name__ == '__main__':
    M, NQ = int(sys.argv[1]), int(sys.argv[2])
    t = time.time(); out = compute(M, NQ); print('time', time.time()-t)
    pickle.dump(out, open(f'cF_{M}_{NQ}.pkl', 'wb'))
    for m in range(1, min(M, 4)+1):
        print(m, sorted((k, v) for k, v in out[m].items() if k[0] <= 2 and abs(k[1]) <= m+1))
