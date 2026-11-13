"""Assemble checks.log from the three check logs and summarise the Rademacher run per m."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
rad = os.path.join(here, "checks_rademacher_GM40_GK16.log")
out = []
for f in ("checks_exact.log", "checks_polar.log"):
    out.append(open(os.path.join(here, f)).read().rstrip())
lines = open(rad).read().rstrip().splitlines()
out.append("\n".join(lines))
summ = {}
for L in lines:
    p = L.split()
    if len(p) == 8 and p[0].isdigit():
        m = int(p[0]); rel, im, tb = map(float, p[5:8])
        s = summ.setdefault(m, [0, 0.0, 0.0, 0.0, 0])
        s[0] += 1; s[1] = max(s[1], rel); s[2] = max(s[2], im); s[3] = max(s[3], tb)
        s[4] += rel <= tb
out.append("== summary of the Rademacher check (relative errors; tail = explicit bound of Prop. 4.8 / |exact|) ==")
for m in sorted(summ):
    c, r, i, t, ok = summ[m]
    out.append(f"m={m}: {c} coefficients, max rel.err {r:.3g}, max |Im|/|exact| {i:.3g}, max tail/|exact| {t:.3g}, rel.err <= tail bound in {ok}/{c}")
if not lines[-1].startswith("done"):
    sys.exit("Rademacher log incomplete")
open(os.path.join(here, "checks.log"), "w").write("\n\n".join(out) + "\n")
print("\n".join(out[-1 - len(summ):]))
