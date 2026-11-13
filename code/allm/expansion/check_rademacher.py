"""Check the level-one Rademacher expansion against exact coefficients from the product formula,
for 1<=m<=6, 0<=l<=m, 0<Delta=4mn-l^2<=200. Polar data from the wall formula (check_polar.py)."""
import json, time, sys
import mpmath as mp
from rademacher import Rad, tail_bound
from polar_formula import polar_value
from exact_common import dcoef
ex = json.load(open("exact_cF.json"))
GM, GK, DPS = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
DMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 200
out = open(f"checks_rademacher_GM{GM}_GK{GK}.log", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.write(s+"\n"); out.flush()
P(f"== Rademacher series (truncated: gamma<={GM} in M,K2; gamma<={GK} in K3; mp.dps={DPS}) vs exact ==")
P("columns: m n l Delta exact  rel.err  |Im|/|exact|  explicit tail bound/|exact|")
worst = {}
t0 = time.time()
for m in range(1, 7):
    polar = {}
    for lt in range(2*m):
        for nt in range(-1, m+1):
            if 4*m*nt - lt*lt < 0:
                v = polar_value(m, nt, lt)
                if v: polar[(nt, lt)] = v
    dm = dcoef(m)[m]
    R = Rad(m, polar, dm, GM=GM, GK=GK, dps=DPS, ts_degree=6)
    P(f"m={m}: d(m)={dm}, nonzero polar coefficients {sorted(polar.items())}")
    cnt = 0; wr = mp.mpf(0); wi = mp.mpf(0)
    for l in range(0, m+1):
        for n in range(0, 60):
            D = 4*m*n - l*l
            if D <= 0 or D > DMAX: continue
            key = f"{n},{l}"
            if key not in ex[str(m)]: continue
            exact = (-1)**(l+1)*ex[str(m)][key]
            v, _ = R.value(n, l)
            rel = abs(mp.re(v) - exact)/abs(exact)
            im = abs(mp.im(v))/abs(exact)
            tb = tail_bound(m, D, polar, dm, GM, GK)/abs(exact)
            mp.mp.dps = DPS
            P(f"{m} {n} {l} {D} {exact} {mp.nstr(rel,3)} {mp.nstr(im,3)} {mp.nstr(tb,3)}")
            cnt += 1; wr = max(wr, rel); wi = max(wi, im)
    P(f"m={m}: {cnt} coefficients, max rel.err {mp.nstr(wr,3)}, max |Im|/|exact| {mp.nstr(wi,3)}, t={time.time()-t0:.0f}s")
P("done")
