import sys; sys.path.insert(0,'.')
from series import *
from m1 import m1_sequences
from m2q import m2q
from math import comb, sqrt
N = 400
# ratio lemma checks
for k in (24, 30, 36):
    p = qpoch_pow(-k, N).coeffs()
    assert all(int(p[i+1]) <= (k+1)*int(p[i]) for i in range(N-1)), k
Psi = trunc(qpoch_pow(1, N, 4)*qpoch_pow(-1, N), N)   # (q^4;q^4)/(q;q)
ps = [int(c) for c in Psi.coeffs()]
assert all(c >= 1 for c in ps), 'Psi>=1/(1-q)'
for k in (30, 36):
    pk = qpoch_pow(-k, N)
    lhs = trunc(mul(P(ps), P([1]), N), N)
    psik = P([1])
    for _ in range(k): psik = mul(psik, Psi, N)
    # dilation identity check
    pk4 = qpoch_pow(-k, N, 4)
    assert mul(psik, pk4, N) == pk, 'dilation'
    psk = [int(c) for c in psik.coeffs()]
    assert all(psk[i+1] <= (k+1)*psk[i] for i in range(len(psk)-1)), 'ratio Psi^k'
    assert all(psk[i] >= comb(i+k-1, k-1) for i in range(len(psk)))
print('ratio lemma, dilation, Psi bounds ok up to', N)
# (I1)
E4 = eis(4, N)
T = P([1 if isqrt(4*n+1)**2 == 4*n+1 else 0 for n in range(N+1)])
th0 = P([(1 if n == 0 else 2) if isqrt(n)**2 == n else 0 for n in range(N+1)])
E4T = mul(E4, T, N)
E4c = E4T - P([0, 225])
psi30 = P([1])
for _ in range(30): psi30 = mul(psi30, Psi, N)
bser = [int(c) for c in mul(E4c, psi30, N).coeffs()]
b = [bser[4*i] for i in range(N//4)]
c = [int(x) for x in mul(E4, th0, N).coeffs()]
C1 = [b[i+1] >= 961*c[i] for i in range(N//4 - 1)]
print('C1 holds for i <', N//4-1, all(C1), 'min ratio', min(b[i+1]/(961*c[i]) for i in range(N//4-1)))
p30 = [int(x) for x in qpoch_pow(-30, N).coeffs()]
fails = [n for n in range(1, N-1) if not p30[n] > 30*c[n] + c[n+1]]
print('C2 fails at n =', fails)
# check analytic bound for c(i)
assert all(c[i] <= (1+2*sqrt(i))*288.5*i**3 for i in range(1, N))
K0, K1 = m1_sequences(4*40+8)
print('I1 exact n=1..40:', all(int(K0[4*n]) > 2*int(K1[n+1]) for n in range(1, 41)))
for n in range(1, 8): print(n, int(K0[4*n]), int(K1[n+1]))
