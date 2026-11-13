"""Certificate for the analytic part of the doubling inequality (Delta' >= 61).
rho = E(4D,b)/E(D,a) >= (2 sqrt2/3) 2^(-23/2) 2^(-1/2) exp(h0 pi sqrt D),  h0 = 2 sqrt2 - sqrt5,
for 1<=a<=b<=4a; the right side increases with D.  Needed: rho > 2(1+e0)/(1-e0) with
e0 = uniform bound for |eps| (Theorem E).  Also a direct ball check of the four factor bounds on a grid."""
from flint import arb, ctx
from cert_common import fup, fdown
ctx.prec = 160
PI = arb.pi(); NU = arb(23)/2
e0 = arb('0.298921')            # Theorem E: cert_rangeA (m>=100); range B (4<=m<=99); exact values (m<=3)
need = 2*(1+e0)/(1-e0)
h0 = 2*arb(2).sqrt() - arb(5).sqrt()
def rho_low(D):
    return 2*arb(2).sqrt()/3*arb(2)**(-12)*(h0*PI*arb(D).sqrt()).exp()
print(f'h0 = 2 sqrt2 - sqrt5 = {h0.str(10)}')
print(f'required ratio 2(1+e0)/(1-e0) with e0=0.298921: {need.str(8)}')
for D in (61, 63, 100, 1000):
    print(f'rho_low({D}) = {rho_low(D).str(8)}')
assert rho_low(61) > need
print(f'rho(61) >= {fdown(rho_low(61),2)};  2(1+e0)/(1-e0) <= {fup(need,4)};  rho(61)(1-e0)/(2(1+e0)) >= {fdown(rho_low(61)/need,3)}')
# direct check of rho >= rho_low on a grid (sanity of the factor bounds)
def E(D, m):
    t0 = m*m + 4*m
    return 2*PI/arb(2*m).sqrt()*(m+2)*(arb(t0)/D)**(arb(23)/4)*(PI*(arb(D)*t0).sqrt()/m).bessel_i(NU)
worst = None
for a in range(1, 41):
    for b in range(a, 4*a+1):
        for D in (3*a*a, 3*a*a+50, 3*a*a+500):
            r = E(4*D, b)/E(D, a)/rho_low(D)
            assert r > 1, (a, b, D)
            if worst is None or r < worst[0]: worst = (r, a, b, D)
print(f'grid 1<=a<=40, a<=b<=4a, D in {{3a^2, 3a^2+50, 3a^2+500}}: min rho/rho_low = {worst[0].str(6)} at (a,b,D)={worst[1:]}')
print('DOUBLING (analytic part): for Delta\'>=61, d*(g) >= (1-e0)E(g) > 2(1+e0)E(f) >= 2 d*(f).')
