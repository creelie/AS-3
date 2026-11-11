import sys; sys.path.insert(0,'.')
from series import *
from fractions import Fraction as F
from math import comb, sqrt, prod, e
from m1 import m1_sequences
from m2q import m2q
Pi1=F(prod(range(30,38)),prod(range(5,13))); print('Pi1',Pi1,float(Pi1))
Pi2=F(prod(range(13,41)),prod(range(8,36))); print('Pi2',Pi2,float(Pi2))
z3=1.2020569031595942; z5=1.0369277551433699
print('240z3',240*z3,'3024z5',3024*z5,'7920z3',7920*z3)
N=300
E2,E4,E6=eis(2,N),eis(4,N),eis(6,N)
g=[int(x) for x in (P([5])*E6-P([11])*mul(E2,E4,N)).coeffs()]
print('g first',g[:5]); assert all(abs(g[a])<=9520.3*a**5 for a in range(1,N))
print('max |g|/a^5',max(abs(g[a])/a**5 for a in range(1,N)))
e4=[int(x) for x in E4.coeffs()]; assert all(e4[a]<=288.5*a**3 for a in range(1,N))
RP,VSP,C2P,H1P,R,S=m2q(N)
T = P([1 if isqrt(4*n+1)**2 == 4*n+1 else 0 for n in range(N+1)])
th0 = P([(1 if n == 0 else 2) if isqrt(n)**2 == n else 0 for n in range(N+1)])
DT = P([n*c for n, c in enumerate(T.coeffs())]); Dth0 = P([n*c for n, c in enumerate(th0.coeffs())])
gP=P(g)
phi0=[int(x) for x in (mul(P([66])*T+P([264])*DT,E4,N)+mul(gP,T,N)).coeffs()]
phi1=[int(x) for x in (P([264])*mul(Dth0,E4,N)+mul(gP,th0,N)).coeffs()]
print('phi0',phi0[:4],'phi1',phi1[:4])
assert all(abs(phi0[n])<=104726*(1+sqrt(n))*n**5 for n in range(1,N))
assert all(abs(phi1[n])<=85685*(1+2*sqrt(n))*n**5 for n in range(1,N))
Sc=[int(x) for x in S.coeffs()]
assert all(abs(Sc[n])<=380822*(1+sqrt(n))*(1+2*sqrt(n))*n**5 for n in range(1,N))
V2 = P([1 if isqrt(8*n+1)**2 == 8*n+1 else 0 for n in range(N+1)])
sig=[int(x) for x in mul(V2,S,N).coeffs()]
assert all(abs(sig[i])<=5.5164e6*i**6.5 for i in range(1,N))
print('sigma0..3',sig[:4])
Rc=[int(x) for x in R.coeffs()]; print('R0..5',Rc[:6])
print('rho0',36*Rc[4]+Rc[5], 'rho0-2738*108', 36*Rc[4]+Rc[5]-2738*108)
# C1 analytic constants
print('C1: need 256/24*Pi1 >= 961*288.5*3:', 256/24*float(Pi1), 961*288.5*3)
print('C1p: R4*16384/5040*Pi2 =',Rc[4]*16384/5040*float(Pi2),' vs 2738*5.5164e6=',2738*5.5164e6)
# C2 m=1: f(n)=C(n+29,29)/(n+1)^3.5 > 26830.5 from n=9
for n in range(7,11): print('m1 f',n,comb(n+29,29)/(n+1)**3.5)
for j in range(0,7): print('m2 f',j,comb(j+35,35),0.822*(j+2)**6.5, (36*Rc[4]+Rc[5]-295704)*comb(j+35,35) > 2*37*5.5164e6*(j+2)**6.5)
print('2*37*5.5164e6/ (rho0-295704)=',2*37*5.5164e6/(36*Rc[4]+Rc[5]-295704))
print('monotone consts', 3.5*e**1.75, 6.5*e**(6.5/7))
K0,K1=m1_sequences(4*12+8)
print('I1 table'); 
for n in range(1,10): print(n,4*n-1,int(K0[4*n]),int(K1[n+1]), int(K0[4*n])>2*int(K1[n+1]))
RP,VSP,C2P,H1P,R,S=m2q(120)
print('I2 table')
for j in range(0,6):
    Dp=8*j+7;K=(Dp+3)//2;M=(Dp+9)//8
    L2=F(2,3)*(RP[K]+F(4800,12)*C2P[K-1]); L1=F(2,3)*(VSP[M]-F(4800,12)*H1P[M])
    print(j,Dp,L2,L1,L2>2*L1)
