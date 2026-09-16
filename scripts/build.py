import itertools, random, sys, time
import numpy as np
from collections import defaultdict, Counter

# ---------- signed-subset complex (barycentric subdivision of cross-polytope boundary), m=3 ----------
m = 3
verts = [v for v in itertools.product([-1,0,1], repeat=m) if any(v)]
def neg(v): return tuple(-x for x in v)
def leq(x,y): return all(xi==0 or xi==yi for xi,yi in zip(x,y))
def size(v): return sum(1 for x in v if x)
edges = [(x,y) for x in verts for y in verts if x!=y and leq(x,y)]
tris  = [(x,y,z) for x in verts for y in verts for z in verts
         if size(x)==1 and size(y)==2 and size(z)==3 and leq(x,y) and leq(y,z)]
print(f"complex: |V|={len(verts)} |E|={len(edges)} |T|={len(tris)}  Euler={len(verts)-len(edges)+len(tris)}")

# free vertices = one rep per antipodal pair
free = [v for v in verts if next(x for x in v if x)==1]
fidx = {v:i for i,v in enumerate(free)}
def rep(v):           # (free index, sign)
    return (fidx[v],+1) if v in fidx else (fidx[neg(v)],-1)
print(f"free vertices: {len(free)}")

def lab(v, L):        # L = list of labels on free vertices
    i,s = rep(v); return s*L[i]

def complementary_edges(L):
    return sum(1 for x,y in edges if lab(x,L) == -lab(y,L))

def pos_alt(L):       # number of positively alternating triangles
    c = 0
    for t in tris:
        ls = sorted((lab(v,L) for v in t), key=abs)
        if len({abs(l) for l in ls})==3 and ls[0]>0 and ls[1]<0 and ls[2]>0: c += 1
    return c

# ---------- octahedron (unsubdivided) degeneracy demo ----------
oct_tris = [(s1,s2,s3) for s1 in (1,-1) for s2 in (1,-1) for s3 in (1,-1)]
def oct_posalt(l1,l2,l3):
    c=0
    for s in oct_tris:
        ls = sorted((s[0]*l1,s[1]*l2,s[2]*l3), key=abs)
        if len({abs(l) for l in ls})==3 and ls[0]>0 and ls[1]<0 and ls[2]>0: c+=1
    return c
labs4 = [l for k in range(1,5) for l in (k,-k)]
vals = Counter()
for l1,l2,l3 in itertools.product(labs4, repeat=3):
    if len({abs(l1),abs(l2),abs(l3)})==3:          # = no complementary edge on the octahedron
        vals[oct_posalt(l1,l2,l3)] += 1
print(f"octahedron, labels ±1..±4, satisfying labelings: {sum(vals.values())}, A+ histogram: {dict(vals)}")

# ---------- DFS over labelings on the subdivided complex ----------
# precompute pairwise forbidden combos: for free i<j, set of (a,b) with a on i, b on j forbidden
def forbidden(labels):
    F = defaultdict(set)
    for x,y in edges:
        (i,si),(j,sj) = rep(x),rep(y)
        if i==j: continue
        for a in labels:
            for b in labels:
                if si*a == -sj*b:
                    F[(min(i,j),max(i,j))].add((a,b) if i<j else (b,a))
    return F

def dfs_solutions(labels, cap=200000, seed=0):
    F = forbidden(labels); n=len(free); L=[None]*n; out=[]
    rng = random.Random(seed)
    def ok(i,a):
        for j in range(i):
            if (L[j],a) in F.get((j,i),()): return False
        return True
    def rec(i):
        if len(out)>=cap: return
        if i==n: out.append(tuple(L)); return
        ls = labels[:]; rng.shuffle(ls)
        for a in ls:
            if ok(i,a):
                L[i]=a; rec(i+1); L[i]=None
    rec(0); return out

t=time.time()
sols2 = dfs_solutions([1,-1,2,-2])
print(f"Tucker (labels ±1,±2): solutions found = {len(sols2)}   [{time.time()-t:.1f}s]")

t=time.time()
sols3 = dfs_solutions([1,-1,2,-2,3,-3], cap=200000)
hist = Counter(pos_alt(L) for L in sols3)
assert all(complementary_edges(L)==0 for L in sols3[:2000])
print(f"Ky Fan (labels ±1..±3): sampled {len(sols3)} solutions, A+ histogram: {dict(sorted(hist.items()))}  [{time.time()-t:.1f}s]")
print("all odd:", all(k%2==1 for k in hist))
