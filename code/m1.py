import sys; sys.path.insert(0,'.')
from series import *
def m1_sequences(N):
    """kappa0(4D) for D+1<=N via index, kappa1 at index n+1<=N.
    returns K0[k] = [6E4T/(q)^30 + 648 H0/(q)^24]_k, K1[k] = 3[E4 th0/(q)^30]_{k} - 648[H1*/(q)^24]_{k-1}"""
    H12 = hurwitz(4*N+4)
    E4 = eis(4, N); i30 = qpoch_pow(-30, N); i24 = qpoch_pow(-24, N)
    T = P([1 if isqrt(4*n+1)**2 == 4*n+1 else 0 for n in range(N+1)])   # q^{k(k+1)}
    th0 = P([ (1 if n == 0 else 2) if isqrt(n)**2 == n else 0 for n in range(N+1)])
    A = mul(mul(E4, T, N), i30, N)
    # 648 H0 = 54 * (12 H(4n))
    H0 = P([H12[4*n] for n in range(N+1)])
    B = mul(H0, i24, N)
    K0 = [6*a + 54*b for a, b in zip(A.coeffs()+[0]*(N+1), B.coeffs()+[0]*(N+1))][:N+1]
    C = mul(mul(E4, th0, N), i30, N).coeffs() + [0]*(N+1)
    H1 = P([H12[4*n+3] for n in range(N+1)])
    Dd = mul(H1, i24, N).coeffs() + [0]*(N+1)
    K1 = [0] + [3*C[k] - 54*Dd[k-1] for k in range(1, N+1)]
    return K0, K1
if __name__ == '__main__':
    K0, K1 = m1_sequences(12)
    # d*(n, l=1, m=1): D=4n-1 -> K1[n+1]; d*(n,l=0,1): D=4n -> kappa0(4D') index D'+1 with 4D'=4n: K0[n+1]
    print('l=1 (D=4n-1), n=0..:', K1[1:10])
    print('l=0 (D=4n),   n=0..:', K0[1:10])
