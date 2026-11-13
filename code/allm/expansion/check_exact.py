import json, sys, time
from exact_common import *
t0 = time.time()
log = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); log.append(s)
P("== exact coefficients of psi_m^F from the product formula ==")
tab = phi01_c(40)
P("phi_{0,1} c(D), D=-1..12:", [tab[D] for D in range(-1, 13) if D % 4 in (0, 3)])
d = dcoef(12)
P("d(n), n=-1..6:", [d[n] for n in range(-1, 7)])
res = {}
targets = {1: 52, 2: 28, 3: 19, 4: 15, 5: 13, 6: 11, 7: 11, 8: 11, 9: 11, 10: 11, 11: 12, 12: 13}
for m, nmax in targets.items():
    cf = finite_part_coeffs(m, nmax)
    res[m] = {f"{n},{l}": v for (n, l), v in cf.items()}
    # consistency: c(n,l) = c(n+m-l, 2m-l)  (evenness + elliptic invariance)
    bad = 0; tot = 0
    for (n, l), v in cf.items():
        n2 = n + m - l
        if (n2, 2*m-l) in cf:
            tot += 1
            if cf[(n2, 2*m-l)] != v: bad += 1
    # depends only on Delta and l mod 2m (l=0 vs l=2m)
    P(f"m={m}: computed n<=" + str(nmax), f"symmetry checks {tot}, failures {bad}",
      f"[c^F(n,l)=0 for n<=-2 by construction]  t={time.time()-t0:.1f}s")
    sys.stdout.flush()
P("anchors:")
for l in range(3):
    P(f"  (-1)^(l+1) c^F_1(2,{l}) =", (-1)**(l+1)*res[1][f"2,{l}"])
for l in range(3):
    P(f"  (-1)^(l+1) c^F_2(1,{l}) =", (-1)**(l+1)*res[2][f"1,{l}"])
P("  c^F_10(1,6) =", res[10]["1,6"], " (Rossello: 1393564656)")
P("  c^F_m(-1,m), m=1..12:", [res[m][f"-1,{m}"] for m in range(1, 13)], "(expected m+2)")
json.dump({str(m): v for m, v in res.items()}, open("exact_cF.json", "w"))
open("checks_exact.log", "w").write("\n".join(log) + "\n")
