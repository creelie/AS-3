import sys; sys.path.insert(0,'.')
from series import *
from m2q import m2q
from math import comb, sqrt
N = 400
RP, VSP, C2P, H1P, R, S = m2q(N)
Rc = [int(x) for x in R.coeffs()]; Sc = [int(x) for x in S.coeffs()]
print('R[0..6]', Rc[:7]); print('S[0..6]', Sc[:7])
assert all(x > 0 for x in Rc[2:N]), 'R>=0 from 2'
assert all(x > 0 for x in Sc[:N]), 'S>0'
p36 = [int(x) for x in qpoch_pow(-36, N).coeffs()]
# step: bracket
for K in range(3, N):
    br = Rc[2]*p36[K-2] + Rc[3]*p36[K-3] - 6*p36[K] - 3912*p36[K-1]
    assert br >= 0
Psi = trunc(qpoch_pow(1, N, 4)*qpoch_pow(-1, N), N)
psi36 = P([1])
for _ in range(36): psi36 = mul(psi36, Psi, N)
R4 = P([0,0,0,0] + Rc[4:N+1])
rs = [int(x) for x in mul(R4, psi36, N).coeffs()]
rho = [rs[4*i+5] for i in range((N-5)//4)]
V2 = P([1 if isqrt(8*n+1)**2 == 8*n+1 else 0 for n in range(N+1)])
sig = [int(x) for x in mul(V2, S, N).coeffs()]
C1 = [rho[i] >= 2738*sig[i] for i in range(len(rho))]
print('C1p all', all(C1), 'min ratio', min(rho[i]/(2738*sig[i]) for i in range(len(rho))))
fails = [j for j in range(0, len(rho)) if not (rho[0]-2738*sig[0])*p36[j] > 2*(36*sig[j+1] + sig[j+2])]
print('C2p fails', fails)
# S bound check
assert all(abs(Sc[n]) <= 545694*(1+2*sqrt(n))**2*n**5 for n in range(1, N))
print('max S(n)/((1+2sqrt n)^2 n^5)', max(Sc[n]/((1+2*sqrt(n))**2*n**5) for n in range(1,N)))
# exact I2 small
ok = True
for j in range(0, 45):
    Dp = 8*j+7; K=(Dp+3)//2; M=(Dp+9)//8
    L2 = 12*RP[K] + 4800*C2P[K-1]; L1 = 12*VSP[M] - 4800*H1P[M]   # times 18
    if not (L2 > 2*L1): ok=False
    if j < 6: print(j, Dp, L2//18 if L2%18==0 else L2/18, L1/18)
print('I2 exact j<45', ok)
