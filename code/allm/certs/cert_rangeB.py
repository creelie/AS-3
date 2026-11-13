"""Certificate, range B: magnetic charges 1 <= m <= 99.
For each m it bounds |eps(f)| for all reduced forms f of magnetic charge m and discriminant >= D,
at the values of D needed in the proofs (the bound is decreasing in D)."""
import time, sys
from cert_common import *
t00 = time.time()
M1 = 100
d = dtable(4*M1)
rows = []
maxP = (0, None); maxF = (0, None); maxG = (0, None); maxAll3 = (0, None)
print('m | Dmin=3m^2: eps_bar (eps_P, eps_gamma, eps_K2, eps_K3) | D=max(61,3m^2) | D=max(244,3m^2) | D=241 (m<=2)')
for m in range(1, M1):
    t0 = time.time()
    polar = {k: polar_sb(m, *k, d) for k in polar_list(m)}
    assert all(v > 0 for v in polar.values())
    e3 = eps_bar(m, 3*m*m, d, polar)
    eF = eps_bar(m, max(61, 3*m*m), d, polar)
    eG = eps_bar(m, max(244, 3*m*m), d, polar)
    s = f'{m:3d} | {fup(e3[0],6)} ({fup(e3[1],6)}, {fupe(e3[2],3)}, {fupe(e3[3],3)}, {fupe(e3[4],3)}) | {fup(eF[0],6)} | {fup(eG[0],6)}'
    if m <= 2:
        e241 = eps_bar(m, 241, d, polar)
        s += f' | {fupe(e241[0],4)}'
    else:
        if up(e3[0]) > maxAll3[0]: maxAll3 = (up(e3[0]), m)
        if m >= 4 and up(e3[0]) > maxP[0]: maxP = (up(e3[0]), m)
    if up(eF[0]) > maxF[0]: maxF = (up(eF[0]), m)
    if up(eG[0]) > maxG[0]: maxG = (up(eG[0]), m)
    print(s + f'   [{len(polar)} polar terms, {time.time()-t0:.1f}s]', flush=True)
print(f'SUMMARY range B (1<=m<=99):')
print(f'  max_{{3<=m<=99}} eps_bar(m,3m^2) = {fup(maxAll3[0],6)} at m={maxAll3[1]}')
print(f'  max_{{4<=m<=99}} eps_bar(m,3m^2) = {fup(maxP[0],6)} at m={maxP[1]}')
print(f'  max_{{1<=m<=99}} eps_bar(m,max(61,3m^2))  = {fup(maxF[0],6)} at m={maxF[1]}')
print(f'  max_{{1<=m<=99}} eps_bar(m,max(244,3m^2)) = {fup(maxG[0],6)} at m={maxG[1]}')
print(f'  total time {time.time()-t00:.1f}s')
