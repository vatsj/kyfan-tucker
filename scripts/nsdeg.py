import itertools, sys, time
import numpy as np
from collections import defaultdict
exec(open('build.py').read().split('# ---------- octahedron')[0])   # reuse complex + rep()

def build_violating(labels):
    q=len(labels); li={a:k for k,a in enumerate(labels)}
    V=set()
    for x,y in edges:
        (i,si),(j,sj)=rep(x),rep(y)
        if i==j: continue
        for a in labels:
            for b in labels:
                if si*a==-sj*b: V.add(frozenset([(i,li[a]),(j,li[b])]))
    return V

def ns_dual_system(labels, d):
    """Degree-d Nullstellensatz dual (Sherali-Adams-style consistency) in the one-hot quotient.
       Returns (rows, ncols) with rows as list of (dict col->coeff, rhs)."""
    n=len(free); q=len(labels); V=build_violating(labels)
    def violating(beta):
        return any(frozenset(p) in V for p in itertools.combinations(beta,2))
    col={}
    for S in itertools.combinations(range(n),d):
        for ls in itertools.product(range(q),repeat=d):
            beta=tuple(zip(S,ls))
            if not violating(beta): col[beta]=len(col)
    memo={}
    def form(beta):
        beta=tuple(sorted(beta))
        if beta in memo: return memo[beta]
        if len(beta)==d:
            f={col[beta]:1} if beta in col else {}
        else:
            used={v for v,_ in beta}; vmin=min(v for v in range(n) if v not in used)
            f=defaultdict(int)
            for l in range(q):
                for c,x in form(beta+((vmin,l),)).items(): f[c]=(f[c]+x)
            f=dict(f)
        memo[beta]=f; return f
    rows=[]
    for k in range(d):
        for S in itertools.combinations(range(n),k):
            for ls in itertools.product(range(q),repeat=k):
                beta=tuple(zip(S,ls)); used=set(S)
                vmin=min(v for v in range(n) if v not in used)
                for w in range(n):
                    if w in used or w==vmin: continue
                    f=defaultdict(int)
                    for c,x in form(beta).items(): f[c]+=x
                    for l in range(q):
                        for c,x in form(beta+((w,l),)).items(): f[c]-=x
                    f={c:x for c,x in f.items() if x!=0}
                    if f: rows.append((f,0))
    rows.append((dict(form(())),1))
    return rows,len(col)

def refutable_F2(rows,ncols):
    R=len(rows); W=(ncols+1+63)//64
    M=np.zeros((R,W),dtype=np.uint64)
    for r,(f,rhs) in enumerate(rows):
        for c,x in f.items():
            if x%2: M[r,c//64]^=np.uint64(1)<<np.uint64(c%64)
        if rhs%2: M[r,ncols//64]^=np.uint64(1)<<np.uint64(ncols%64)
    r=0
    for c in range(ncols):
        w,b=c//64,c%64
        colbits=(M[r:,w]>>np.uint64(b))&np.uint64(1)
        nz=np.nonzero(colbits)[0]
        if len(nz)==0: continue
        p=r+nz[0]
        if p!=r: M[[r,p]]=M[[p,r]]
        mask=((M[:,w]>>np.uint64(b))&np.uint64(1)).astype(bool); mask[r]=False
        M[mask]^=M[r]
        r+=1
        if r==R: break
    rhsbit=(M[r:,ncols//64]>>np.uint64(ncols%64))&np.uint64(1)
    return bool(rhsbit.any()), r

def refutable_Fp(rows,ncols,p=1000003):
    # sparse Gaussian elimination mod p (dict rows)
    rows=[({c:x%p for c,x in f.items() if x%p},rhs%p) for f,rhs in rows]
    rows=[(f,rhs) for f,rhs in rows if f or rhs]
    pivots={}   # col -> (row dict, rhs) with leading col
    for f,rhs in rows:
        f=dict(f)
        while f:
            c=min(f)
            if c in pivots:
                pf,prhs=pivots[c]; x=f[c]
                for cc,y in pf.items():
                    v=(f.get(cc,0)-x*y)%p
                    if v: f[cc]=v
                    else: f.pop(cc,None)
                rhs=(rhs-x*prhs)%p
            else:
                inv=pow(f[c],p-2,p)
                f={cc:(y*inv)%p for cc,y in f.items()}; rhs=(rhs*inv)%p
                pivots[c]=(f,rhs); break
        else:
            if rhs: return True, len(pivots)
    return False, len(pivots)

labels=[1,-1,2,-2]
for d in [2,3]:
    t=time.time(); rows,ncols=ns_dual_system(labels,d)
    print(f"d={d}: unknowns (non-violating size-{d} partial labelings)={ncols}, equations={len(rows)}  [build {time.time()-t:.1f}s]")
    t=time.time(); ref,rk=refutable_F2(rows,ncols)
    print(f"   F_2 : degree-{d} refutation exists = {ref}   (rank {rk})  [{time.time()-t:.1f}s]")
    t=time.time(); ref,rk=refutable_Fp(rows,ncols)
    print(f"   F_p : degree-{d} refutation exists = {ref}   (rank {rk})  [{time.time()-t:.1f}s]")
    sys.stdout.flush()
