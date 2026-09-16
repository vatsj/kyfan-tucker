import itertools, sys
import numpy as np
from collections import defaultdict
exec(open('pseudo2.py').read().split("# --- invariant unknowns")[0])


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
    return x,[c for c in range(C) if c not in piv]

def solve_restricted(gens, support, label):
    # orbits under subgroup <gens> of non-violating size-2 partial labelings satisfying support()
    orb={}; orbs=[]
    for S in itertools.combinations(range(n),d):
        for ls in itertools.product(labels,repeat=d):
            a=tuple(zip(S,ls))
            if a in orb or is_viol(a,V) or not support(a): continue
            stack=[a]; orb[a]=len(orbs); members=[a]
            while stack:
                b_=stack.pop()
                for g in gens:
                    c=act(g,b_)
                    if c not in orb: orb[c]=len(orbs); members.append(c); stack.append(c)
            orbs.append(members)
    def form(beta):
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
    m=len(orbs); A=np.zeros((len(rows),m),dtype=np.uint8)
    for r,f in enumerate(rows):
        for c in f: A[r,c]=1
    x,freecols=gf2_solve_all(A,np.array(rhs,dtype=np.uint8))
    print(f"{label:55s} unknowns={m:5d}  ->  {'NO solution' if x is None else f'solution, {len(freecols)} free params'}")
    return x,orbs

ID=(0,1,2); I1={1:1,2:2}; E1={1:1,2:1}
full=[((1,0,2),(1,1,1),I1,E1),((0,2,1),(1,1,1),I1,E1),(ID,(-1,1,1),I1,E1),(ID,(1,1,1),{1:2,2:1},E1),(ID,(1,1,1),I1,{1:-1,2:1})]
z3z3=[((1,2,0),(1,1,1),I1,E1)]  # only one Z3 available on labels with 2 magnitudes; coordinate Z3
hemi=[((1,0,2),(1,1,1),I1,E1),(ID,(-1,1,1),I1,E1),(ID,(1,1,1),{1:2,2:1},E1),(ID,(1,1,1),I1,{1:-1,2:1})]  # stabilizer of coordinate 3
labelsym=[(ID,(1,1,1),{1:2,2:1},E1),(ID,(1,1,1),I1,{1:-1,2:1})]
none=[]
distinct=lambda a: abs(a[0][1])!=abs(a[1][1])
anyS=lambda a: True
solve_restricted(none,anyS,"unrestricted")
solve_restricted(none,distinct,"support: distinct magnitudes only")
solve_restricted(z3z3,anyS,"Z3 (coordinate cycle) invariant")
solve_restricted(hemi,anyS,"hemisphere-stabilizer invariant")
solve_restricted(labelsym,anyS,"label-symmetry invariant (signed perms of labels)")
solve_restricted(labelsym,distinct,"label-sym invariant + distinct magnitudes")
x,orbs=solve_restricted(hemi,distinct,"hemisphere-stabilizer invariant + distinct magnitudes")
if x is not None:
    print("\nsolution by orbit (hemisphere-stabilizer, distinct magnitudes):")
    for k_,o in enumerate(orbs):
        if x[k_]:
            a=o[0]; (i,li),(j,lj)=a
            print(f"  E=1 size={len(o):3d} rep: {free[i]}->{li:+d}, {free[j]}->{lj:+d}   type={ptype(a)}")
