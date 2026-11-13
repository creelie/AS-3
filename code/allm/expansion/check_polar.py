import json
from polar_formula import polar_value, walls
ex = json.load(open("exact_cF.json"))
log = []
def P(*a):
    s=" ".join(str(x) for x in a); print(s); log.append(s)
P("== polar coefficients: wall formula (R2) vs exact coefficients ==")
total = 0; fails = 0; nonzero = 0; maxrs = 0; viol_t0 = 0
for m in range(1, 13):
    cnt = 0
    for lt in range(0, 2*m):
        for nt in range(-4, m+1):
            if 4*m*nt - lt*lt >= 0: continue
            val, used = polar_value(m, nt, lt, RS=4*m+6, with_walls=True)
            for w in used:
                maxrs = max(maxrs, max(w[2], w[3]) / (2*m))
            exact = ex[str(m)].get(f"{nt},{lt}", 0 if nt <= -2 else None)
            assert exact is not None
            total += 1; cnt += 1
            if val != exact:
                fails += 1; P("  MISMATCH", m, nt, lt, val, exact)
            if exact != 0:
                nonzero += 1
                if lt*lt - 4*m*nt > m*m + 4*m: viol_t0 += 1
    P(f"m={m}: {cnt} polar pairs (nt>=-4) checked")
P(f"total {total} polar pairs, {nonzero} nonzero, mismatches {fails}")
P(f"max over contributing walls of max(r,s)/(2m) = {maxrs:.3f} (Lemma: <=1)")
P(f"nonzero polar coefficients with |Delta~| > m^2+4m: {viol_t0} (Lemma: 0)")
P("examples m=10:", [(nt, lt, polar_value(10, nt, lt)) for (nt, lt) in [(-1,10),(0,1),(0,7),(1,9),(2,19)]])
open("checks_polar.log","w").write("\n".join(log)+"\n")
