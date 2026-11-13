"""Polar coefficients of psi_m^F from the wall formula (Theorem R, part (ii)),
evaluated along the Stern-Brocot chain of x0 = lt/(2m)."""
import sys, json, math, time
from exact_psi import dcoeffs, compute_psiF

def polar_list(m):
    """all (nt, lt) with 0<=lt<=m, nt>=-1, 4 m nt - lt^2 < 0"""
    out = []
    for lt in range(0, m+1):
        nt = -1
        while 4*m*nt - lt*lt < 0:
            out.append((nt, lt)); nt += 1
    return out

def polar_sb(m, nt, lt, d, chain=False):
    Q = lambda x, y: m*x*x - lt*x*y + nt*y*y
    B = lambda u, v: 2*m*u[0]*v[0] - lt*(u[0]*v[1] + v[0]*u[1]) + 2*nt*u[1]*v[1]
    v, u = (0, 1), (1, 1)        # interval [q/s, p/r) = [0,1)
    tot = 0; steps = []
    while True:
        A, Bv = Q(*u), Q(*v)
        if A < -1 or Bv < -1: break
        lg = -B(u, v)
        assert lg > 0 and lg*lg == 4*A*Bv + (lt*lt - 4*m*nt)
        tot += lg*d[A]*d[Bv]; steps.append((u, v, A, Bv, lg))
        w = (u[0]+v[0], u[1]+v[1])
        if lt*w[1] < 2*m*w[0]: u = w
        else: v = w
    return (tot, steps) if chain else tot

if __name__ == '__main__':
    M = int(sys.argv[1])
    d = dcoeffs(4*M*M)
    bad = 0; cnt = 0; nz = 0; maxsteps = 0
    for m in range(1, M+1):
        ex = compute_psiF(m, max(0, (m*m)//(4*m) + 1))
        for (nt, lt) in polar_list(m):
            v, st = polar_sb(m, nt, lt, d, chain=True)
            cnt += 1; nz += (v != 0); maxsteps = max(maxsteps, len(st))
            if v != ex[(nt, lt)]:
                bad += 1; print('MISMATCH', m, nt, lt, v, ex[(nt, lt)])
            # check: c^F vanishes for nt=-2 (outside list) and positivity
            assert v > 0
        print(f'm={m}: {len(polar_list(m))} polar pairs, all match exact' if bad == 0 else f'm={m} bad={bad}', flush=True)
    print(f'TOTAL polar pairs {cnt}, nonzero {nz}, mismatches {bad}, max chain length {maxsteps}')
