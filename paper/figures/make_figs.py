import sys, os, pickle, math
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..', '..', 'code')
sys.path.insert(0, CODE)
os.chdir(HERE)
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from series import *
from m1 import m1_sequences
from m2q import m2q
from fractions import Fraction as F
from math import comb, log10
plt.rcParams.update({'font.size': 10, 'font.family': 'serif'})
cF = pickle.load(open(os.path.join(CODE, 'data', 'cF_240.pkl'), 'rb'))
def reduce(n, l, m):
    a, b, c = m, l, n
    while True:
        if c < a: a, c = c, a; b = -b
        if abs(b) > a:
            k = round(b/(2*a)); c = c-b*k+a*k*k; b = b-2*a*k; continue
        if c < a: continue
        return a, b, c
def ds(n, l, m):
    a, b, c = reduce(n, l, m); return (1 if b % 2 else -1)*cF[a][(c, abs(b))]

# Figure 1: divisor lattice of r = 60 with the pairing s -> s/2
fig = plt.figure(figsize=(6.2, 5.0)); ax = fig.add_subplot(111, projection='3d')
for a in range(3):
    for b in range(2):
        for c in range(2):
            neg = (a == 2)
            col = '#b2182b' if neg else '#2166ac'
            ax.scatter(a, b, c, s=90, color=col, depthshade=True, edgecolor='k')
            ax.text(a+0.06, b+0.06, c+0.06, str(2**a*3**b*5**c), fontsize=8)
            for (da, db, dc) in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                if a+da <= 2 and b+db <= 1 and c+dc <= 1:
                    ax.plot([a, a+da], [b, b+db], [c, c+dc], color='0.75', lw=0.8)
            if neg:
                ax.quiver(a, b, c, -0.85, 0, 0, color='#1b7837', lw=1.6, arrow_length_ratio=0.25)
ax.set_xlabel('exponent of 2'); ax.set_ylabel('exponent of 3'); ax.set_zlabel('exponent of 5')
ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1]); ax.set_zticks([0, 1])
ax.view_init(elev=22, azim=-58)
ax.set_title('Divisors $s$ of $r=60$ for $\\ell/r$ odd: negative terms (red) and $s\\mapsto s/2$', fontsize=9)
plt.tight_layout(); plt.savefig('fig_pairing.png', dpi=200); plt.close()

# Figure 2: doubling margin, 3D scatter (m, D', log10 ratio) from exact data 4D' <= 240
xs, ys, zs = [], [], []
for m in range(1, 200):
    for l in range(-m, m+1):
        if l % 2 == 0: continue
        for n in range(1, 400):
            D = 4*m*n-l*l
            if D <= 0: continue
            if 4*D > 240: break
            f = ds(4*n, 2*l, m); g = ds(n, l, m)
            xs.append(m); ys.append(D); zs.append(log10(f/(2*g)))
fig = plt.figure(figsize=(6.2, 5.0)); ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(xs, ys, zs, c=zs, cmap='viridis', s=8)
ax.set_xlabel("$m=P^2/2$"); ax.set_ylabel("$D'$"); ax.set_zlabel('$\\log_{10}$ of ratio')
ax.set_title("$d^*(4n',2\\ell',m)\\,/\\,2d^*(n',\\ell',m)$ for $\\ell'$ odd, $4D'\\leq 240$", fontsize=9)
ax.view_init(elev=20, azim=-130)
fig.colorbar(sc, shrink=0.6, pad=0.1)
plt.tight_layout(); plt.savefig('fig_doubling.png', dpi=200); plt.close()

# Figure 3: m=1 and m=2 doubling margins to large D'
K0, K1 = m1_sequences(4*120+8)
x1 = [4*n-1 for n in range(1, 121)]; y1 = [log10(int(K0[4*n])/(2*int(K1[n+1]))) for n in range(1, 121)]
RP, VSP, C2P, H1P, R, S = m2q(8*60+20)
x2, y2 = [], []
for j in range(0, 58):
    Dp = 8*j+7; K = (Dp+3)//2; M = (Dp+9)//8
    L2 = F(2, 3)*(RP[K]+F(4800, 12)*C2P[K-1]); L1 = F(2, 3)*(VSP[M]-F(4800, 12)*H1P[M])
    x2.append(Dp); y2.append(log10(L2/(2*L1)))
fig, ax = plt.subplots(figsize=(6.0, 3.6))
ax.plot(x1, y1, 'o-', ms=2.5, lw=1, color='#2166ac', label="$m=1$: $\\kappa_0(4D')/2\\kappa_1(D')$")
ax.plot(x2, y2, 's-', ms=2.5, lw=1, color='#b2182b', label="$m=2$: $L_2(4D')/2L_1(D')$")
ax.set_xlabel("$D'$"); ax.set_ylabel('$\\log_{10}$ of the ratio'); ax.legend(frameon=False)
ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig('fig_margins.png', dpi=200); plt.close()

# Figure 4: the signed terms of the torsion sum, 3D bars, for r=12, m=1,2
fig = plt.figure(figsize=(6.4, 4.8)); ax = fig.add_subplot(111, projection='3d')
divs = [1, 2, 3, 4, 6, 12]
from m2 import m2_sequences
L2s = m2_sequences(140)
def ds12(n, l, m):
    D = 4*m*n-l*l
    if m == 1:
        return int(K1[(D+1)//4+1]) if l % 2 else int(K0[D//4+1])
    j = abs(l) % 4; j = 1 if j == 3 else j
    return L2s[j][D]
rows = []
for k, (m, l0, n0) in enumerate(((1, 1, 1), (2, 1, 1), (1, 1, 2))):
    r = 12; n, l = r*r*n0, r*l0
    if 4*m*n-l*l <= 0: continue
    terms = [s*(-1)**((l-l//s) % 2)*ds12(n//(s*s), l//s, m) for s in divs]
    rows.append(((m, l0, n0), terms))
for k, (lab, terms) in enumerate(rows):
    for i, t in enumerate(terms):
        h = math.copysign(log10(abs(t)), t)
        ax.bar3d(i, k, 0, 0.6, 0.6, h, color='#b2182b' if t < 0 else '#2166ac', alpha=0.85, shade=True)
ax.set_xticks([i+0.3 for i in range(len(divs))]); ax.set_xticklabels([str(s) for s in divs])
ax.set_yticks([k+0.3 for k in range(len(rows))])
ax.set_yticklabels(["$(%d,%d,%d)$" % (144*a[2], 12*a[1], a[0]) for a, _ in rows], fontsize=7)
ax.set_xlabel('$s$'); ax.set_zlabel('signed $\\log_{10}|$term$|$')
ax.set_title('Terms $s(-1)^{\\ell-\\ell/s}d^*(n/s^2,\\ell/s,m)$ for $r=12$', fontsize=9)
ax.view_init(elev=24, azim=-60)
plt.tight_layout(); plt.savefig('fig_terms.png', dpi=200); plt.close()
for lab, t in rows: print(lab, t, sum(t))
