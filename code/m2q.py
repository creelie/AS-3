import sys; sys.path.insert(0,'.')
from series import *
def m2q(N):
    """q-series for m=2. Returns (L2num, L1num) where
    3/2*L2(4D') = RP[K] + 4800*C2P[K-1],  K=(D'+3)/2
    3/2*L1(D')  = VSP[M] - 4800*H1P[M],   M=(D'+9)/8  ; class numbers scaled: we return 12x integer parts"""
    H12 = hurwitz(8*N + 16)
    E2, E4, E6 = eis(2, N), eis(4, N), eis(6, N)
    T = P([1 if isqrt(4*n+1)**2 == 4*n+1 else 0 for n in range(N+1)])
    T2 = P([1 if (n % 2 == 0 and isqrt(2*n+1)**2 == 2*n+1) else 0 for n in range(N+1)])  # T(q^2): q^{2k(k+1)}
    V2 = P([1 if isqrt(8*n+1)**2 == 8*n+1 else 0 for n in range(N+1)])   # q^{k(k+1)/2}
    th0 = P([(1 if n == 0 else 2) if isqrt(n)**2 == n else 0 for n in range(N+1)])
    th0q2 = P([(1 if n == 0 else 2) if (n % 2 == 0 and isqrt(n//2)**2 == n//2) else 0 for n in range(N+1)])
    DT = P([n*c for n, c in enumerate(T.coeffs())])
    Dth0 = P([n*c for n, c in enumerate(th0.coeffs())])
    M = lambda a, b: mul(a, b, N)
    g = P([5])*E6 - P([11])*M(E2, E4)
    phi0 = M(P([66])*T + P([264])*DT, E4) + M(g, T)
    phi1 = P([264])*M(Dth0, E4) + M(g, th0)
    R = P([0, 8])*M(M(T2, T), phi0) + M(M(th0q2, th0), phi1)
    R = trunc(R, N)
    S = P([2])*M(phi0, th0) + P([2])*M(T, phi1)
    i36 = qpoch_pow(-36, N); i24 = qpoch_pow(-24, N)
    RP = M(R, i36).coeffs() + [0]*(N+2)
    VSP = M(M(V2, S), i36).coeffs() + [0]*(N+2)
    def c2_12(n):  # 12*c(n,2)
        v = H12[8*n - 4] if 8*n-4 >= 0 else 0
        if n % 2 == 0: v += 2*H12[2*n - 1] if 2*n-1 >= 0 else 0
        return v
    C2 = P([c2_12(k+1) for k in range(N+1)])
    H1 = P([0] + [H12[8*l - 1] for l in range(1, N+1)])
    C2P = M(C2, i24).coeffs() + [0]*(N+2)
    H1P = M(H1, i24).coeffs() + [0]*(N+2)
    return [int(x) for x in RP], [int(x) for x in VSP], [int(x) for x in C2P], [int(x) for x in H1P], R, S
if __name__ == '__main__':
    from fractions import Fraction as Fr
    RP, VSP, C2P, H1P, R, S = m2q(40)
    print('R', R.coeffs()[:5], 'S', S.coeffs()[:5])
    from m2 import m2_sequences
    L = m2_sequences(30)
    for Dp in range(7, 120, 8):
        K = (Dp+3)//2; Mm = (Dp+9)//8
        L2 = Fr(2,3)*(RP[K] + Fr(4800,12)*C2P[K-1])
        L1 = Fr(2,3)*(VSP[Mm] - Fr(4800,12)*H1P[Mm])
        print(Dp, L2 == L[2].get(4*Dp), L1 == L[1].get(Dp))
