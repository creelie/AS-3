import pickle,sys
cF=pickle.load(open(sys.argv[1],'rb')); Dmax=int(sys.argv[2])
def reduce(n,l,m):
    a,b,c=m,l,n
    while True:
        if c<a: a,c=c,a; b=-b
        if abs(b)>a:
            k=round(b/(2*a)); c=c-b*k+a*k*k; b=b-2*a*k; continue
        if c<a: continue
        return a,b,c
def ds(n,l,m):
    a,b,c=reduce(n,l,m); return (1 if b%2 else -1)*cF[a][(c,abs(b))]
# doubling: all classes of forms (n',l',m) with l' odd and 4D'<=Dmax, reduced reps
cnt=0;worst=None;miss=0;mn=set()
for a in range(1,80):
  for b in range(-a,a+1):
    if b%2==0: continue
    for c in range(a,400):
      D=4*a*c-b*b
      if D<=0 or 4*D>Dmax: continue
      # m can be any of the represented values; the pair (n',l',m) and its double (4n',2l',m)
      # depends on which coefficient is m; use all three SL2 reps with the given m: (c,b,a),(a,-b,c),(a+b+c,b+2a?)
      for (n1,l1,m1) in ((c,b,a),(a,-b,c),(a+b+c,b+2*c,c)):
        if l1%2==0: continue
        try: f=ds(4*n1,2*l1,m1); g=ds(n1,l1,m1)
        except KeyError: miss+=1; continue
        cnt+=1; mn.add(m1); r=f/(2*g)
        if worst is None or r<worst[0]: worst=(r,(n1,l1,m1),f,g)
print('doubling checks',cnt,'missing',miss,'worst f/(2g)',worst, 'max m',max(mn))
cnt=0;worst=None;miss=0;ms=set()
for m in range(1,200):
  for l in range(-m,m+1):
    if l%2==0: continue
    for n in range(1,400):
      D=4*m*n-l*l
      if D<=0: continue
      if 4*D>Dmax: break
      try: f=ds(4*n,2*l,m); g=ds(n,l,m)
      except KeyError: miss+=1; continue
      cnt+=1; ms.add(m); r=f/(2*g)
      if worst is None or r<worst[0]: worst=(r,(n,l,m),f,g)
print('full doubling checks',cnt,'missing',miss,'worst',worst,'m range',min(ms),max(ms))
