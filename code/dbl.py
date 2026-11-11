import pickle, sys
cF = pickle.load(open('cF.pkl','rb'))
def reduce(n, l, m):
    a, b, c = m, l, n
    while True:
        if c < a: a, c = c, a; b = -b
        if abs(b) > a:
            k = round(b/(2*a)); c = c - b*k + a*k*k; b = b - 2*a*k; continue
        if c < a: continue
        return a, b, c
def dstar(n,l,m):
    a,b,c = reduce(n,l,m); return (1 if b%2 else -1)*cF[a][(c,abs(b))]
Dmax=120
worst=None; cnt=0
for a in range(1,8):
  for b in range(-a,a+1):
    for c in range(a, 200):
      D=4*a*c-b*b
      if D<=0 or 4*D>Dmax: continue
      g=dstar(c,b,a)
      for (n2,l2,m2) in ((4*c,2*b,a),(c,2*b,4*a),(c+b+a,2*(b+2*a),4*a)):
        try: f=dstar(n2,l2,m2)
        except KeyError: continue
        cnt+=1
        r=f/(2*abs(g))
        if worst is None or r<worst[0]: worst=(r,(a,b,c),(n2,l2,m2),f,g)
print(cnt,worst)
