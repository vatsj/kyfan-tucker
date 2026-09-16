import itertools, sys, time, random
import numpy as np
from collections import defaultdict
exec(open('build.py').read().split('# ---------- octahedron')[0])
n=len(free)

def violating_pairs(labels):
    V=set()
    for x,y in edges:
        (i,si),(j,sj)=rep(x),rep(y)
        if i==j: continue
        for a in labels:
            for b in labels:
                if si*a==-sj*b: V.add(frozenset([(i,a),(j,b)]))
    return V

def is_viol(alpha,V):
    return any(frozenset(p) in V for p in itertools.combinations(alpha,2))

# ---- symmetry group: B_3 on the complex x signed perms of magnitudes ----
def group(labels):
    mags=sorted({abs(l) for l in labels})
    G=[]
    for pi in itertools.permutations(range(3)):
        for eps in itertools.product([1,-1],repeat=3):
            for lpi in itertools.permutations(mags):
                for leps in itertools.product([1,-1],repeat=len(mags)):
                    G.append((pi,eps,dict(zip(mags,lpi)),dict(zip(mags,leps))))
    return G
def act(g,alpha):
    pi,eps,lpi,leps=g; out=[]
    for i,l in alpha:
        v=free[i]; w=tuple(eps[k]*v[pi[k]] for k in range(3))
        j,s=rep(w); tl=leps[abs(l)]*(1 if l>0 else -1)*lpi[abs(l)]
        out.append((j,s*tl))
    return tuple(sorted(out))

def random_labeling(labels,rng): return [rng.choice(labels) for _ in range(n)]
def restrict(L,S): return tuple((i,L[i]) for i in S)

def gf2_solve(rows_bits, ncols, rhs):
    """rows_bits: list of python ints (bitsets over ncols). Solve x with A x = rhs over F2. Returns x (int bitset) or None."""
    piv={}  # leading col -> (row, rhs)
    for r,b in zip(rows_bits,rhs):
        while r:
            c=r.bit_length()-1
            if c in piv: r^=piv[c][0]; b^=piv[c][1]
            else: piv[c]=(r,b); break
        else:
            if b: return None
    # back-substitute
    x=0
    for c in sorted(piv):
        r,b=piv[c]
        # r has leading bit c; other bits are lower and (after sorting ascending) already solved
        val=b ^ (bin(r & x & ~(1<<c)).count('1')&1)
        if val: x|=1<<c
    return x

def fp_solve(A,b,p=1000003):
    """Dense GE mod p. A: (R,C) int64. Returns solution vector or None."""
    A=A.copy()%p; b=b.copy()%p; R,C=A.shape; r=0; pivcols=[]
    for c in range(C):
        nz=np.nonzero(A[r:,c])[0]
        if len(nz)==0: continue
        pi=r+nz[0]
        if pi!=r: A[[r,pi]]=A[[pi,r]]; b[[r,pi]]=b[[pi,r]]
        inv=pow(int(A[r,c]),p-2,p); A[r]=(A[r]*inv)%p; b[r]=(b[r]*inv)%p
        f=A[:,c].copy(); f[r]=0
        A=(A-np.outer(f,A[r]))%p; b=(b-f*b[r])%p
        pivcols.append(c); r+=1
        if r==R: break
    if np.any(b[r:]!=0): return None
    x=np.zeros(C,dtype=np.int64)
    for k,c in enumerate(pivcols): x[c]=b[k]
    return x

# ======================= TUCKER, labels ±1,±2 =======================
labels=[1,-1,2,-2]; V=violating_pairs(labels); rng=random.Random(1)
G=group(labels)
# sanity: group preserves violating set
Vimg={frozenset(act(G[rng.randrange(len(G))],tuple(sorted(p)))) for p in V}
print("group preserves violating pairs:", Vimg==V, " |G| =",len(G))

def primal_unrestricted_F2(labels,V,d,target,nsamp,nverify,rng):
    cols={}
    for S in itertools.combinations(range(n),d):
        for ls in itertools.product(labels,repeat=d):
            a=tuple(zip(S,ls))
            if is_viol(a,V): cols[a]=len(cols)
    def row(L):
        r=0
        for S in itertools.combinations(range(n),d):
            a=restrict(L,S)
            if a in cols: r|=1<<cols[a]
        return r
    Ls=[random_labeling(labels,rng) for _ in range(nsamp)]
    x=gf2_solve([row(L) for L in Ls],len(cols),[target(L) for L in Ls])
    if x is None: return len(cols),None,None
    bad=0
    for _ in range(nverify):
        L=random_labeling(labels,rng)
        if (bin(row(L)&x).count('1')&1)!=target(L): bad+=1
    return len(cols),bin(x).count('1'),bad

t=time.time()
nc,wt,bad=primal_unrestricted_F2(labels,V,3,lambda L:1,9000,20000,rng)
print(f"Tucker F_2 d=3 primal: {nc} violating size-3 unknowns; certificate found = {wt is not None}"
      + (f", support {wt}, failures on 20000 fresh labelings = {bad}" if wt is not None else "") + f"  [{time.time()-t:.0f}s]")

def primal_symmetric_Fp(labels,V,d,nsamp,nverify,rng,p=1000003):
    # orbits of violating size-d partial labelings
    orb={}; norb=0
    for S in itertools.combinations(range(n),d):
        for ls in itertools.product(labels,repeat=d):
            a=tuple(zip(S,ls))
            if a in orb or not is_viol(a,V): continue
            imgs={act(g,a) for g in G}
            for b in imgs: orb[b]=norb
            norb+=1
    def row(L):
        r=np.zeros(norb,dtype=np.int64)
        for S in itertools.combinations(range(n),d):
            a=restrict(L,S)
            if a in orb: r[orb[a]]+=1
        return r
    A=np.array([row(random_labeling(labels,rng)) for _ in range(nsamp)]); b=np.ones(nsamp,dtype=np.int64)
    x=fp_solve(A,b,p)
    if x is None: return len(orb),norb,None,None
    bad=sum(1 for _ in range(nverify) if int(row(random_labeling(labels,rng))@x)%p!=1)
    return len(orb),norb,int(np.count_nonzero(x)),bad

for d in [3,4]:
    t=time.time()
    nv,norb,supp,bad=primal_symmetric_Fp(labels,V,d,4000,3000,rng)
    print(f"Tucker F_p d={d} symmetric primal: {nv} violating unknowns in {norb} orbits; certificate found = {supp is not None}"
          + (f", support {supp} orbits, failures on 3000 fresh = {bad}" if supp is not None else "") + f"  [{time.time()-t:.0f}s]")
    sys.stdout.flush()

# ======================= KY FAN, labels ±1..±3, target A+ + 1 mod 2 =======================
labels3=[1,-1,2,-2,3,-3]; V3=violating_pairs(labels3)
t=time.time()
nc,wt,bad=primal_unrestricted_F2(labels3,V3,3,lambda L:(pos_alt(L)+1)&1,22000,20000,rng)
print(f"Ky Fan F_2 d=3 primal: {nc} violating size-3 unknowns; certificate found = {wt is not None}"
      + (f", support {wt}, failures on 20000 fresh labelings = {bad}" if wt is not None else "") + f"  [{time.time()-t:.0f}s]")
