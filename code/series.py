# one-variable q-series for the m=1 and m=2 single-centred indices (Bringmann-Murthy forms)
import flint
from fractions import Fraction
from math import isqrt
P = flint.fmpz_poly
def trunc(p, N):
    c = p.coeffs()[:N+1]; return P(c)
def mul(a, b, N): return trunc(a*b, N)
def inv(a, N):
    # power series inverse of a with a(0)=+-1
    c = a.coeffs(); assert abs(c[0]) == 1
    res = [0]*(N+1); res[0] = c[0]
    # Newton would be faster; use flint's series via fmpq? simple: use doubling
    r = P([c[0]]); k = 1
    while k <= N:
        k2 = min(2*k, N+1)
        e = trunc(a, k2-1)
        r = trunc(r*(P([2]) - trunc(e*r, k2-1)), k2-1)
        k = k2
    return trunc(r, N)
def qpoch_pow(k, N, step=1):
    # prod_{n>=1} (1-q^{n*step})^k truncated: via euler pentagonal for k=1 then power
    c = [0]*(N+1)
    j = 0
    while True:
        done = True
        for jj in (j, -j) if j else (0,):
            e = step*jj*(3*jj-1)//2
            if e <= N: c[e] += (-1)**(jj % 2); done = False
        if done and j > 0: break
        j += 1
    e1 = P(c)
    if k >= 0:
        r = P([1]); b = e1; kk = k
        while kk:
            if kk & 1: r = mul(r, b, N)
            b = mul(b, b, N); kk >>= 1
        return r
    else:
        return inv(qpoch_pow(-k, N, step), N)
def sigma(n, k):
    s = 0; d = 1
    while d*d <= n:
        if n % d == 0:
            s += d**k
            if d*d != n: s += (n//d)**k
        d += 1
    return s
def eis(k, N, scale=1):
    # E_k normalized: 1 + c*sum sigma_{k-1} q^n
    cst = {2: -24, 4: 240, 6: -504}[k]
    return P([1] + [cst*sigma(n, k-1) for n in range(1, N+1)])
def hurwitz(Nmax):
    # 12*H(N) as integers for 0..Nmax (H(0)=-1/12)
    H12 = [0]*(Nmax+1); H12[0] = -1
    # count forms ax^2+bxy+cy^2, disc b^2-4ac=-N, reduced: |b|<=a<=c, b>=0 if |b|=a or a=c; weights 1/2 for a=b=c, 1/3 for a=c,b=0
    for N in range(1, Nmax+1):
        if N % 4 not in (0, 3): continue
        tot = 0  # in units of 1/6
        b = N % 2
        while b*b <= N//3:
            m = (b*b + N)//4
            a = max(b, 1)
            while a*a <= m:
                if m % a == 0:
                    c = m//a
                    if a == b or a == c or b == 0:
                        w = 6
                        if b == 0 and a == c: w = 3
                        if a == b == c: w = 2
                        tot += w
                    else:
                        tot += 12
                a += 1
            b += 2
        H12[N] = 2*tot
    return H12
