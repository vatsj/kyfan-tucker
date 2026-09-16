import itertools, sys
import numpy as np
from collections import defaultdict
exec(open('build.py').read().split('# ---------- octahedron')[0])
exec(open('primal.py').read().split('# ======================= TUCKER')[0].split("exec(open('build.py')")[0])
# (reuse violating_pairs, is_viol, group, act from primal.py by re-exec of its defs)
src=open('primal.py').read(); defs=src.split("def violating_pairs")[1].split("# ======================= TUCKER")[0]
exec("def violating_pairs"+defs)
n=len(free); labels=[1,-1,2,-2]; V=violating_pairs(labels); G=group(labels); d=2

# --- describe a size-2 partial labeling by an invariant "type" ---
def vtype(i,j):
    x,y=free[i],free[j]; nx=neg(x); 
    if leq(x,y) or leq(y,x): rel='comparable'
    elif leq(nx,y) or leq(y,nx): rel='anti-comparable'
    else: rel='incomparable'
    return rel,(size(x),size(y))
def ptype(a):
    (i,li),(j,lj)=a
    rel,sz=vtype(i,j)
    # sign relative to comparability direction is messy; record magnitudes and whether signs agree
    return (rel, tuple(sorted(sz)), 'same-mag' if abs(li)==abs(lj) else 'diff-mag', 'same-sign' if (li>0)==(lj>0) else 'opp-sign')

# --- invariant unknowns: orbits of NON-violating size-2 partial labelings under full group ---
orb={}; orbs=[]
for S in itertools.combinations(range(n),d):
    for ls in itertools.product(labels,repeat=d):
        a=tuple(zip(S,ls))
        if a in orb or is_viol(a,V): continue
        imgs={act(g,a) for g in G}
        for b in imgs: orb[b]=len(orbs)
        orbs.append(sorted(imgs))
print(f"non-violating size-2 partial labelings: {len(orb)} in {len(orbs)} full-group orbits")

# --- SA consistency system on invariant E ---
def form(beta):   # returns dict orbit->coeff (mod 2) for E[beta], beta a partial labeling of size<=2
    beta=tuple(sorted(beta))
    if len(beta)==d: return ({orb[beta]:1} if beta in orb else {})
    used={v for v,_ in beta}; vmin=min(v for v in range(n) if v not in used); f=defaultdict(int)
    for l in labels:
        for c,x in form(beta+((vmin,l),)).items(): f[c]^=x
    return {c:x for c,x in f.items() if x}
rows=[]; rhs=[]
for k in range(d):
    for S in itertools.combinations(range(n),k):
        for ls in itertools.product(labels,repeat=k):
            beta=tuple(zip(S,ls)); used=set(S); vmin=min(v for v in range(n) if v not in used)
            for w in range(n):
                if w in used or w==vmin: continue
                f=defaultdict(int)
                for c,x in form(beta).items(): f[c]^=x
                for l in labels:
                    for c,x in form(beta+((w,l),)).items(): f[c]^=x
                f={c:x for c,x in f.items() if x}
                if f: rows.append(f); rhs.append(0)
rows.append(form(())); rhs.append(1)
# solve over F2 (small)
m=len(orbs); A=np.zeros((len(rows),m),dtype=np.uint8)
for r,f in enumerate(rows):
    for c in f: A[r,c]=1
b=np.array(rhs,dtype=np.uint8)
def gf2_solve_all(A,b):
    A=A.copy(); b=b.copy(); R,C=A.shape; r=0; piv=[]
    for c in range(C):
        nz=np.nonzero(A[r:,c])[0]
        if len(nz)==0: continue
        p=r+nz[0]; A[[r,p]]=A[[p,r]]; b[[r,p]]=b[[p,r]]
        m_=A[:,c].astype(bool); m_[r]=False
        A[m_]^=A[r]; b[m_]^=b[r]; piv.append(c); r+=1
        if r==R: break
    if b[r:].any(): return None,None
    x=np.zeros(C,dtype=np.uint8)
    for k_,c in enumerate(piv): x[c]=b[k_]
    freecols=[c for c in range(C) if c not in piv]
    return x,freecols
x,freecols=gf2_solve_all(A,b)
if x is None: print("no FULLY-invariant degree-2 pseudo-solution over F_2"); sys.exit()
print(f"fully-invariant degree-2 pseudo-solution EXISTS; free parameters: {len(freecols)}")
print("\nE on non-violating size-2 orbits (particular solution):")
for k_,o in enumerate(orbs):
    a=o[0]; print(f"  E={int(x[k_])}  size={len(o):3d}  type={ptype(a)}  rep={a}")
