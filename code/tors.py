import sys, pickle
sys.argv=['x','0']
exec(open('dstar.py').read().split("D = int")[0])
cF = pickle.load(open('cF.pkl','rb'))
from math import isqrt
def reduce(n, l, m):
    a, b, c = m, l, n
    while True:
        if c < a: a, c = c, a; b = -b
        if abs(b) > a:
            k = round(b/(2*a))
            c = c - b*k + a*k*k; b = b - 2*a*k
            continue
        if c < a: continue
        return a, b, c
def dstar(n, l, m):
    a, b, c = reduce(n, l, m)
    v = cF[a][(c, abs(b))]
    return (1 if b % 2 else -1)*v
Dmax = int(sys.argv_[0]) if False else 120
bad = []; cnt = 0; worst = None
for r in (2,3,4,5,6):
    for m in range(1, 60):
        for l0 in range(-60, 61):
            for n0 in range(1, 60):
                n, l = r*r*n0, r*l0
                Dl = 4*m*n - l*l
                if Dl <= 0 or Dl > Dmax: continue
                try:
                    tot = 0; terms=[]
                    for s in range(1, r+1):
                        if r % s: continue
                        t = s*(-1)**((l - l//s) % 2)*dstar(n//(s*s), l//s, m)
                        tot += t; terms.append(t)
                except KeyError:
                    continue
                cnt += 1
                ratio = tot/terms[0]
                if worst is None or ratio < worst[0]: worst = (ratio, r, n, l, m, Dl, terms)
                if tot <= 0: bad.append((r, n, l, m, Dl, terms))
print(cnt, 'checked; bad', len(bad)); print(bad[:10]); print('worst ratio', worst)
