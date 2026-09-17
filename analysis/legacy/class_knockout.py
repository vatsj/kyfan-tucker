"""Class-level knockout / greedy on S^{m-1} at degree m-1: classes = (simplex?, #edge pairs, magnitude multiset)."""
import itertools, time, sys
from collections import defaultdict
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.group import z3z3_generators
from kyfan.complex import leq, neg, rank
m=int(sys.argv[1]); k=m-1
cx=SignedComplex(m); labels=label_set(k); V=violating_pairs(cx,labels); d=m-1; gens=z3z3_generators(m,list(range(1,k+1)))
def is_edge(x,y): return leq(x,y) or leq(y,x) or leq(neg(x),y) or leq(y,neg(x))
def is_simplex(a):
    vs=[cx.free[i] for i,_ in a]
    for signs in itertools.product((1,-1),repeat=len(vs)):
        ws=sorted([tuple(s*t for t in v) for s,v in zip(signs,vs)],key=rank)
        if all(leq(ws[i],ws[i+1]) for i in range(len(ws)-1)): return True
    return False
def cls(a):
    ne=sum(1 for p,q in itertools.combinations(range(len(a)),2) if is_edge(cx.free[a[p][0]],cx.free[a[q][0]]))
    mags=tuple(sorted(abs(l) for _,l in a)); rel={}
    for x in mags:
        rel.setdefault(x,len(rel)+1)
    return ('S' if is_simplex(a) else 'N', ne, tuple(rel[x] for x in mags))
col,members=unknowns(cx,labels,d,V,gens=gens); rows,ncols,_=sa_dual_system(cx,labels,d,V=V,gens=gens,col=col)
bycls=defaultdict(set)
for a,c in col.items(): bycls[cls(a)].add(c)
classes=sorted(bycls); print(f"S^{m-1} degree {d}: {len(classes)} classes")
def solve_without(removed):
    rs=[({c:x for c,x in f.items() if c not in removed},b) for f,b in rows]; rs=[(f,b) for f,b in rs if f or b]
    return solve(rs,ncols,"F2")
print("single-class knockouts:")
necessary=[]
for C in classes:
    res=solve_without(bycls[C]); 
    if not res.consistent: necessary.append(C)
    print(f"  remove {str(C):24s} ({len(bycls[C]):6d} orbit-cols): {'solution, %d free' % res.n_free if res.consistent else 'NO solution'}", flush=True)
print("necessary classes:", necessary)
for order_name,order in (("largest first",sorted(classes,key=lambda C:-len(bycls[C]))),("smallest first",sorted(classes,key=lambda C:len(bycls[C])))):
    removed=set()
    for C in order:
        cand=removed|bycls[C]
        if solve_without(cand).consistent: removed=cand
    kept=[C for C in classes if not bycls[C]<=removed]
    print(f"greedy ({order_name}) minimal sufficient class family: {kept}", flush=True)
