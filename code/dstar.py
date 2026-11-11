import sys, pickle, time
import psif
from math import isqrt, gcd
psif.LRANGE = lambda m, lmax: range(0, m+1)
D = int(sys.argv[1])
M = isqrt(D//3) + 1
NQ = D//3 + 2
t=time.time()
cF = psif.compute(M, NQ)
print('computed', time.time()-t, M, NQ)
def reduce(n, l, m):
    # form n x^2 + l x y + m y^2 (n=Q^2/2, l=Q.P, m=P^2/2) -> reduced (a,b,c) |b|<=a<=c
    a, b, c = m, l, n
    while True:
        if c < a: a, c = c, a; b = -b
        if abs(b) > a:
            k = (b + a) // (2*a) if b > 0 else -((-b + a) // (2*a))
            # b -> b - 2 a k
            c = c - b*k + a*k*k; b = b - 2*a*k
            continue
        if c < a: continue
        return a, b, c
def dstar(n, l, m):
    a, b, c = reduce(n, l, m)
    assert 4*a*c - b*b == 4*m*n - l*l
    v = cF[a].get((c, abs(b)))
    if v is None: raise KeyError((n,l,m,a,b,c))
    return (-1)**(abs(b)+1)*v if False else (1 if b % 2 else -1)*v
pickle.dump(cF, open('cF_%d.pkl'%D,'wb'))
