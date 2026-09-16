import itertools, random, time, sys, numpy as np
M=4
verts=[v for v in itertools.product([-1,0,1],repeat=M) if any(v)]
def neg(v): return tuple(-x for x in v)
def leq(x,y): return all(xi==0 or xi==yi for xi,yi in zip(x,y))
edges=[(x,y) for x in verts for y in verts if x!=y and leq(x,y)]
free=[v for v in verts if next(x for x in v if x)==1]; fidx={v:i for i,v in enumerate(free)}; n=len(free)
def rep(v): return (fidx[v],1) if v in fidx else (fidx[neg(v)],-1)
print(f"S^3 complex: |V|={len(verts)}, free={n}, |E|={len(edges)}")
labels=[1,-1,2,-2,3,-3]; mags=[1,2,3]
V=set()
for x,y in edges:
    (i,si),(j,sj)=rep(x),rep(y)
    if i==j: continue
    for a in labels:
        for b in labels:
            if si*a==-sj*b: V.add(frozenset([(i,a),(j,b)]))
def is_viol(al): return any(frozenset(p) in V for p in itertools.combinations(al,2))
def act(g,al):
    pi,eps,lpi,leps=g; out=[]
    for i,l in al:
        v=free[i]; w=tuple(eps[k]*v[pi[k]] for k in range(M)); j,s=rep(w)
        out.append((j, s*leps[abs(l)]*(1 if l>0 else -1)*lpi[abs(l)]))
    return tuple(sorted(out))
ID=tuple(range(M)); I1={1:1,2:2,3:3}; E1={1:1,2:1,3:1}
gens=[((1,0,2,3),(1,)*M,I1,E1),((0,2,1,3),(1,)*M,I1,E1),((0,1,3,2),(1,)*M,I1,E1),(ID,(-1,1,1,1),I1,E1),
      (ID,(1,)*M,{1:2,2:1,3:3},E1),(ID,(1,)*M,{1:1,2:3,3:2},E1),(ID,(1,)*M,I1,{1:-1,2:1,3:1})]
d=3; t=time.time(); orb={}; norb=0
for S in itertools.combinations(range(n),d):
    for ls in itertools.product(labels,repeat=d):
        a=tuple(zip(S,ls))
        if a in orb or not is_viol(a): continue
        stack=[a]; orb[a]=norb
        while stack:
            b=stack.pop()
            for g in gens:
                c=act(g,b)
                if c not in orb: orb[c]=norb; stack.append(c)
        norb+=1
print(f"degree 3: {len(orb)} violating partial labelings in {norb} orbits  [{time.time()-t:.0f}s]"); sys.stdout.flush()
rng=random.Random(5)
def row(L):
    r=np.zeros(norb,dtype=np.int64)
    for S in itertools.combinations(range(n),d):
        a=tuple((i,L[i]) for i in S)
        if a in orb: r[orb[a]]+=1
    return r
Ls=[[rng.choice(labels) for _ in range(n)] for _ in range(4*norb+50)]
A=np.array([row(L) for L in Ls]); b=np.ones(len(Ls),dtype=np.int64)
# F_2 solve
def gf2(A,b):
    rows=[int(''.join('1' if x&1 else '0' for x in r[::-1]),2) for r in A]; piv={}
    for r,bb in zip(rows,b):
        bb&=1
        while r:
            c=r.bit_length()-1
            if c in piv: r^=piv[c][0]; bb^=piv[c][1]
            else: piv[c]=(r,bb); break
        else:
            if bb: return None
    x=0
    for c in sorted(piv):
        r,bb=piv[c]; val=bb^(bin(r&x&~(1<<c)).count('1')&1)
        if val: x|=1<<c
    return x
x=gf2(A,b)
if x is None: print("F_2 degree 3 (symmetric): NO symmetric certificate")
else:
    xs=np.array([(x>>k)&1 for k in range(norb)])
    fresh=[[rng.choice(labels) for _ in range(n)] for _ in range(3000)]
    bad=sum(1 for L in fresh if (row(L)@xs)%2!=1)
    print(f"F_2 degree 3 (symmetric): certificate FOUND, {xs.sum()} of {norb} orbits used, failures on 3000 fresh = {bad}")
# Q solve
xq,res,rk,_=np.linalg.lstsq(A.astype(float),b.astype(float),rcond=None)
print(f"Q   degree 3 (symmetric): relative residual {np.linalg.norm(A@xq-b)/np.sqrt(len(b)):.3e}  ({'certificate' if np.linalg.norm(A@xq-b)/np.sqrt(len(b))<1e-9 else 'NO certificate'})")
