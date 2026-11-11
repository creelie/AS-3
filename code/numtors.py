import pickle,sys,os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'numgen.py')).read().split('# doubling')[0])
cnt=0;bad=0;worst=None
for r in range(2,16):
  for m in range(1,400):
    for l0 in range(-400,401):
      l=r*l0
      if abs(l)>m: continue   # T-reduced representative
      for n0 in range(1,400):
        n=r*r*n0; D=4*m*n-l*l
        if D<=0: continue
        if D>Dmax: break
        tot=0; first=None
        for s in range(1,r+1):
          if r%s: continue
          t=s*(-1)**((l-l//s)%2)*ds(n//(s*s),l//s,m); tot+=t
        cnt+=1
        if tot<=0: bad+=1
        print if False else None
print('torsion charges',cnt,'nonpositive',bad)
