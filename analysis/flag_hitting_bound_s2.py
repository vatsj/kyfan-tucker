import sys, itertools, time; sys.path.insert(0,'.'); import numpy as np
from kyfan.complex import SignedComplex, leq
from kyfan.abstract_gadget import build_dual
from kyfan.labels import label_set
from kyfan import linalg
cx=SignedComplex(3); free=cx.free; N=len(free); labels=[1,-1,2,-2]
# pick the flag class: free reps of the pole chain U0 = e3 < e3-e1 < e3-e1-e2  (ranks 1,2,3)
U=[(0,0,1),(-1,0,1),(-1,-1,1)]
Ur=[cx.rep(v) for v in U]; Uidx=[i for i,_ in Ur]; others=[i for i in range(N) if i not in Uidx]
# edges among free reps with signs
E=[]
for x,y in cx.edges:
    (i,si),(j,sj)=cx.rep(x),cx.rep(y)
    if i!=j: E.append((i,si,j,sj))
Eoff=[(i,si,j,sj) for i,si,j,sj in E if i in others and j in others]
Emix=[(i,si,j,sj) for i,si,j,sj in E if (i in others)!=(j in others)]
Ein=[(i,si,j,sj) for i,si,j,sj in E if i in Uidx and j in Uidx]
print(f"flag free reps {Uidx}; other free vertices {len(others)}; in-U constraints {len(Ein)} (all pairs)")
t=time.time(); domain_types={}
for ls in itertools.product(labels,repeat=len(others)):
    L=dict(zip(others,ls))
    if any(si*L[i]==-sj*L[j] for i,si,j,sj in Eoff): continue
    dom={u:set(labels) for u in Uidx}
    for i,si,j,sj in Emix:
        if i in L: dom[j].discard(-si*L[i]*sj)    # sj*lab(j) != -si*L[i]  => lab(j) != -si*L[i]*sj... careful: sj*y != -si*x => y != -si*x*sj
        else:      dom[i].discard(-sj*L[j]*si)
    key=tuple(tuple(sorted(dom[u])) for u in Uidx)
    domain_types[key]=domain_types.get(key,0)+1
print(f"non-violating restrictions: {sum(domain_types.values())}, distinct domain triples: {len(domain_types)}  [{time.time()-t:.0f}s]")
# residual CSP on the 3 free reps: constraints from Ein with signs; build the degree-2 dual per domain type, keep the consistent ones (= gadgets)
full=list(itertools.product(labels,repeat=3)); fidx={f:k for k,f in enumerate(full)}
def viol_full(f):
    L=dict(zip(Uidx,f)); return any(si*L[i]==-sj*L[j] for i,si,j,sj in Ein)
rows_all=[]; ngad=0; nsat=0
for key,cnt in domain_types.items():
    doms=[list(k) for k in key]
    # SA-dual at degree 2 on 3 variables with these domains; violating = domain or in-U complementary (with signs!)
    def viol(alpha):
        L=dict(alpha)
        return any(l not in doms[Uidx.index(u)] for u,l in alpha) or any(u in L and w in L and su*L[u]==-sw*L[w] for u,su,w,sw in Ein)
    col={}
    for size in range(3):
        for S in itertools.combinations(Uidx,size):
            for ls in itertools.product(labels,repeat=size):
                a=tuple(zip(S,ls))
                if not viol(a): col[a]=len(col)
    rows=[]
    for a,c in col.items():
        if len(a)==2: continue
        used={u for u,_ in a}
        for v in Uidx:
            if v in used: continue
            f={c:1}
            for l in labels:
                b=tuple(sorted(a+((v,l),)))
                if b in col: f[col[b]]=f.get(col[b],0)^1
            rows.append(({k:x for k,x in f.items() if x},0))
    rows.append(({col[()]:1},1))
    res=linalg.gf2_dense(rows,len(col),want_solution=True)
    sat=any(not viol(tuple(zip(Uidx,f))) for f in itertools.product(*doms))
    if sat: nsat+=cnt; continue
    if not res.consistent: continue
    ngad+=cnt
    # pseudo-solution space: E on full assignments = derived from size-2 values by extension (any vertex). Get affine space: particular + nullspace of the homogeneous system
    # extract full-assignment values E(f) = sum over extensions... simpler: E(1_f) for full f = E[f restricted]: we need E on size-3 (full) = not in col (degree 2 only). Use consistency: E[f] is not defined at degree 2!
    # Instead: the constraint sum_{alpha in T} E(1_alpha) uses full alphas -> need E at degree 3 extension? No: restricted certificate has size-3 terms 1_alpha; E is a degree-2 functional... 
    rows_all.append((key,cnt))
print("gadget domain types:",len(rows_all),"| restrictions that are gadgets:",ngad,"| satisfiable residuals:",nsat)

# ---- constraints on T (subset of the 64 full labelings of the 3 free reps) from every gadget type ----
def Zspace(key):
    doms=[list(k) for k in key]; pts=[f for f in itertools.product(*doms)]; pidx={f:k for k,f in enumerate(pts)}
    # violating beta of size <=2 inside the domains: complementary pairs (with signs) among U's reps
    cons=[]
    for (u,su,w,sw) in Ein:
        iu,iw=Uidx.index(u),Uidx.index(w)
        for a in doms[iu]:
            for b in doms[iw]:
                if su*a==-sw*b:
                    row=[k for k,f in enumerate(pts) if f[iu]==a and f[iw]==b]
                    cons.append(row)
    # Z = {e on pts : sum over each constraint row = 0}; nullspace over F2
    A=np.zeros((len(cons),len(pts)),dtype=np.uint8)
    for r,row in enumerate(cons): A[r,row]=1
    # nullspace via RREF
    M=A.copy(); R,C=M.shape; piv=[]; r=0
    for c in range(C):
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        p=r+nz[0]; M[[r,p]]=M[[p,r]]
        m_=M[:,c].astype(bool); m_[r]=False; M[m_]^=M[r]; piv.append(c); r+=1
        if r==R: break
    freec=[c for c in range(C) if c not in piv]; basis=[]
    for fc in freec:
        v=np.zeros(C,dtype=np.uint8); v[fc]=1
        for k,pc in enumerate(piv):
            if M[k,fc]: v[pc]=1
        basis.append(v)
    return pts,basis
constraints=[]   # (vector over 64 full labelings, rhs)
for key,cnt in rows_all:
    pts,basis=Zspace(key)
    for e in basis:
        vec=np.zeros(64,dtype=np.uint8)
        for k,f in enumerate(pts):
            if e[k]: vec[fidx[f]]=1
        constraints.append((vec,int(e.sum()%2)))
A=np.array([v for v,_ in constraints]); b=np.array([r for _,r in constraints],dtype=np.uint8)
print(f"constraints on T: {len(constraints)} from {len(rows_all)} gadget types; Z-dimensions: {[len(Zspace(k)[1]) for k,_ in rows_all]}")
# min-weight solution over 64 variables: solve, then enumerate the solution space if small
M=np.concatenate([A,b[:,None]],axis=1).astype(np.uint8); R=M.shape[0]; piv=[]; r=0
for c in range(64):
    nz=np.nonzero(M[r:,c])[0]
    if len(nz)==0: continue
    p=r+nz[0]; M[[r,p]]=M[[p,r]]; m_=M[:,c].astype(bool); m_[r]=False; M[m_]^=M[r]; piv.append(c); r+=1
    if r==R: break
assert not M[r:,64].any(), "inconsistent"
freec=[c for c in range(64) if c not in piv]; print(f"rank {r}, free variables {len(freec)}")
x0=np.zeros(64,dtype=np.uint8)
for k,c in enumerate(piv): x0[c]=M[k,64]
null=[]
for fc in freec:
    v=np.zeros(64,dtype=np.uint8); v[fc]=1
    for k,pc in enumerate(piv):
        if M[k,fc]: v[pc]=1
    null.append(v)
best=None
if len(freec)<=22:
    for mask in range(1<<len(freec)):
        x=x0.copy()
        for k in range(len(freec)):
            if (mask>>k)&1: x^=null[k]
        w=int(x.sum())
        if best is None or w<best: best=w
    print(f"minimum number of full-flag terms per antipodal flag class on S^2: {best}   (minimum certificate uses 12)")
else: print("solution space too large to enumerate here")

# exact minimum odd-hitting set for the 16 supports
supports=[]
for key,cnt in rows_all:
    pts,basis=Zspace(key); e=basis[0]
    supports.append(frozenset(fidx[f] for k,f in enumerate(pts) if e[k]))
verts=sorted(set().union(*supports)); print(f"supports: {len(supports)} of sizes {sorted(len(s) for s in supports)}, touching {len(verts)} full labelings; max degree {max(sum(1 for s in supports if v in s) for v in verts)}")
best=[len(verts)+1]; bestT=[None]
def bb(chosen, parity):
    if len(chosen)>=best[0]: return
    bad=[k for k,s in enumerate(supports) if parity[k]%2==0]
    if not bad: best[0]=len(chosen); bestT[0]=set(chosen); return
    s=supports[bad[0]]
    for v in sorted(s):
        if v in chosen: continue
        chosen.add(v)
        for k,t in enumerate(supports):
            if v in t: parity[k]+=1
        bb(chosen,parity)
        for k,t in enumerate(supports):
            if v in t: parity[k]-=1
        chosen.remove(v)
bb(set(),[0]*len(supports))
print(f"EXACT minimum full-flag terms per antipodal flag class: {best[0]}  ->  degree-3 certificates need >= 24 x {best[0]} = {24*best[0]} size-3 terms (minimum certificate has 288)")
# check the actual minimum certificate's T for this flag satisfies all parities
terms=[tuple(sorted((int(a),int(b)) for a,b in (p.split(':') for p in line.split()))) for line in open('analysis/isd_s2_d3_cert.txt') if not line.startswith('#')]
Tcert={fidx[tuple(dict(t)[u] for u in Uidx)] for t in terms if len(t)==3 and set(u for u,_ in t)==set(Uidx)}
print(f"certificate's terms on this flag: {len(Tcert)}; parities ok: {all(len(Tcert&s)%2==1 for s in supports)}")
