import itertools, random, time
from collections import Counter

def posalt(ls):   # labels positively alternating: distinct magnitudes, signs +,-,+,... in magnitude order
    s=sorted(ls,key=abs)
    return len({abs(l) for l in s})==len(s) and all((l>0)==(i%2==0) for i,l in enumerate(s))
def negalt(ls): return posalt([-l for l in ls])
def g(ls):        # local lemma value on a label tuple of an n-simplex
    return (sum(posalt(ls[:i]+ls[i+1:]) for i in range(len(ls))) + posalt(ls) + negalt(ls)) % 2
def complementary(ls): return any(a==-b for a,b in itertools.combinations(ls,2))

# ---- (1) local lemma, exhaustive: g == 0 on every non-complementary tuple ----
for n in range(0,4):
    for k in range(1,6):
        labs=[s*i for i in range(1,k+1) for s in (1,-1)]
        bad=sum(1 for ls in itertools.product(labs,repeat=n+1) if not complementary(ls) and g(ls))
        nz =sum(1 for ls in itertools.product(labs,repeat=n+1) if complementary(ls) and g(ls))
        print(f"n={n} k={k}: violations of local lemma = {bad:>2}   (g nonzero on {nz} complementary tuples)")
print()

# ---- (2) the telescoped identity on random labelings, m = 3,4,5 ----
def complex_(m):
    verts=[v for v in itertools.product([-1,0,1],repeat=m) if any(v)]
    def leq(x,y): return all(xi==0 or xi==yi for xi,yi in zip(x,y))
    def top(vs):  # maximal chains within vertex set vs, by rank
        byrank={}; 
        for v in vs: byrank.setdefault(sum(1 for t in v if t),[]).append(v)
        chains=[[v] for v in byrank.get(1,[])]
        for r in range(2,max(byrank)+1):
            chains=[c+[w] for c in chains for w in byrank[r] if leq(c[-1],w)]
        return chains
    return verts,top
for m in [3,4,5]:
    k=m+1; labs=[s*i for i in range(1,k+1) for s in (1,-1)]
    verts,top=complex_(m); rng=random.Random(m)
    free=[v for v in verts if next(t for t in v if t)==1]; fidx={v:i for i,v in enumerate(free)}
    def lab(v,L): return L[fidx[v]] if v in fidx else -L[fidx[tuple(-t for t in v)]]
    # levels j=1..m: hemisphere H^{(j)} = top simplices of the [j]-complex (padded with zeros) whose top element has x_j=+1
    levels=[]
    for j in range(1,m+1):
        vj=[v+(0,)*(m-j) for v in complex_(j)[0]]
        Hj=[c for c in complex_(j)[1](complex_(j)[0]) if c[-1][j-1]==1]
        levels.append([[v+(0,)*(m-j) for v in c] for c in Hj])
    Tm=top(verts)
    stats=Counter(); t=time.time()
    for trial in range(400):
        L=[rng.choice(labs) for _ in free]
        lhs=(sum(posalt([lab(v,L) for v in c]) for c in Tm)+1)%2
        rhs=sum(g(tuple(lab(v,L) for v in c)) for lev in levels for c in lev)%2
        stats['identity holds']+=(lhs==rhs)
        if all(lab(x,L)!=-lab(y,L) for c in Tm for x,y in itertools.combinations(c,2)): stats['no-complementary labelings seen']+=1
    ncert=sum(len(lev) for lev in levels)
    print(f"m={m} (S^{m-1}, {len(verts)} vertices, {len(Tm)} top simplices, labels ±1..±{k}): identity holds on {stats['identity holds']}/400 random labelings; certificate = {ncert} local-lemma instances, degree {m}  [{time.time()-t:.0f}s]")
