import sys, time; sys.path.insert(0,'.')
from m1 import m1_sequences
from m2 import m2_sequences
from math import log
NM1 = int(sys.argv[1]); NM2 = int(sys.argv[2])
t=time.time()
K0, K1 = m1_sequences(NM1)
worst=None; bad=0
for n in range(1, NM1//4+1):
    big, small = int(K0[4*n]), int(K1[n+1])
    assert big > 0 and small > 0
    if big <= 2*small: bad+=1
    r = log(big)-log(2*small)
    if worst is None or r < worst[0]: worst=(r, n)
print('m=1: D\'=4n-1 up to', 4*(NM1//4)-1, 'bad', bad, 'min log-ratio', worst, time.time()-t)
t=time.time()
L = m2_sequences(NM2)
mx = max(L[2].keys())
worst=None; bad=0; cnt=0
D=7
while 4*D <= mx:
    big, small = L[2][4*D], L[1][D]
    assert big > 0 and small > 0
    cnt+=1
    if big <= 2*small: bad+=1
    r = log(big)-log(2*small)
    if worst is None or r < worst[0]: worst=(r, D)
    D += 8
print('m=2: D\' up to', D-8, 'count', cnt, 'bad', bad, 'min log-ratio', worst, time.time()-t)
