"""Certificate, range A: magnetic charges m >= 100 (all reduced forms, Delta >= 3m^2).
  eps_P <= T(M1) = sum_{(nt,k) != (-1,0)} 2 S(k+nt) S(nt) exp(-pi sqrt3 h(nt,k,max(M1,k))),
  h(nt,k,m) = (k(2m-k) + 4m(nt+1))/(2m+4),   S(j) = sum_{i=-1}^j d(i).
  Tail (k+nt > L): S(j) <= 1.05 j^-6 exp(4 pi sqrt j) (j>=1) and h >= theta (k+4nt+4), theta=M1/(2M1+4).
  Mock terms: explicit bound decreasing in m, evaluated at m=M1."""
from cert_common import *
from decimal import Decimal
from flint import arb
import time
t00 = time.time()
M1 = 100
L = 400            # exact S up to L
d = dtable(L+5)
S = {}; acc = 0
for j in range(-1, L+1): acc += d[j]; S[j] = acc
c = PI*arb(3).sqrt()
theta = arb(M1)/(2*M1+4)
# --- check of the analytic bound for S(j) on the exact range 1<=j<=L
P0 = (arb(24)*(-2*PI).exp()/(1-(-2*PI).exp())**2).exp()   # prod (1-z^n)^-24 <= exp(24 z/(1-z)^2), z=e^{-2pi}
assert P0 < arb('1.05')
Sbound = lambda j: arb('1.05')*arb(j)**(-6)*(4*PI*arb(j).sqrt()).exp()
for j in range(1, L+1):
    assert arb(S[j]) < Sbound(j), j
# --- main finite part
T = arb(0)
for nt in range(-1, L+1):
    for k in range(0, L+1-nt):
        if (nt, k) == (-1, 0): continue
        mm = max(M1, k)
        h = arb(k*(2*mm-k) + 4*mm*(nt+1))/(2*mm+4)
        T += 2*arb(S[k+nt])*arb(S[nt])*(-c*h).exp()
Tfinite = T
# --- tail j = k+nt > L:  <= 2 e^{-4 c theta} C1 sum_{j>L} S(j) e^{-c theta j}
C1 = arb(0)
for nt in range(-1, L+1):
    C1 += arb(S[nt])*(-3*c*theta*nt).exp()
# tail of C1 beyond L: terms <= 1.05 exp(4pi sqrt(nt) - 3 c theta nt), ratio <= exp(2pi/sqrt(L) - 3 c theta)
r1 = (2*PI/arb(L).sqrt() - 3*c*theta).exp(); assert r1 < 1
C1 += Sbound(L+1)*(-3*c*theta*(L+1)).exp()/(1-r1)
r2 = (2*PI/arb(L).sqrt() - c*theta).exp(); assert r2 < 1
tailsum = Sbound(L+1)*(-c*theta*(L+1)).exp()/(1-r2)
Ttail = 2*(-4*c*theta).exp()*C1*tailsum
T = Tfinite + Ttail
print(f'M1={M1}, exact S(j) for j<={L}; analytic bound S(j) <= 1.05 j^-6 e^(4 pi sqrt j) verified for 1<=j<={L}')
print(f'T finite part  = {Tfinite.str(10)}')
print(f'T tail bound   = {Ttail.str(5)}')
print(f'eps_P bound T  <= {fup(T,6)}')
eg = SQ2*Z*(1+T)
print(f'eps_gamma <= sqrt2 (zeta(11)-1)(1+T) <= {fupe(eg,4)}   [zeta(11)-1 <= {fupe(Z,7)}]')
# --- mock terms for m >= M1 at Delta=3m^2 (eps_K decreasing in Delta)
X0 = PI*arb(3).sqrt()*M1
kappa = (2*PI*X0).sqrt()*(-X0).exp()*X0.bessel_i(NU)
print(f'kappa = sqrt(2 pi X0) e^-X0 I_23/2(X0) at X0 = pi sqrt3 * {M1}: {kappa.str(8)}')
def logB2(m):
    m = arb(m)
    A = CG*arb('1.05')*(arb(4)/3)**6*arb(3)**(arb(23)/4)*arb(2)**(arb(1)/4)/(2*PI*kappa)
    return A, A*m**(-arb(49)/4)*((4*PI+2*PI*arb(3).sqrt())*m.sqrt() - c*(m*m+4*m).sqrt()).exp()
A2, B2 = logB2(M1)
B3 = B2*SQ2*PI*M1/3*(arb(4)/(3*M1))**(arb(1)/4)*(2*PI/(2*PI*(arb(3*M1)).sqrt())).sqrt()
print(f'eps_K2 <= {A2.str(6)} m^(-49/4) exp((4pi+2pi sqrt3) sqrt m - pi sqrt3 sqrt(m^2+4m)) <= {fupe(B2,4)} for m>={M1}')
print(f'eps_K3 <= eps_K2 bound * (sqrt2 pi m/3)(4/(3m))^(1/4) sqrt(2pi/y) <= {fupe(B3,4)} for m>={M1}')
tot = T + eg + B2 + B3
print(f'RANGE A: for every m >= {M1} and every reduced form, |eps| <= {fup(tot,6)}')
print(f'sum for the paper: eps_bar <= {fup(T,7)} + {fup(eg,7)} + {fupe(B2+B3,2)} <= {fup(T,7)}+{fup(eg,7)}+{fupe(B2+B3,2)}')
print(f'  exact sum of these printed upper bounds <= {Decimal(fup(T,7))+Decimal(fup(eg,7))+Decimal(fupe(B2+B3,2))}')
assert tot < arb('0.298921')
print('assert: |eps| < 0.298921 for m >= 100  (OK)')
print(f'time {time.time()-t00:.1f}s')
