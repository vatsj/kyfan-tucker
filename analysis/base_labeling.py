"""Pseudo-solutions near a base labeling L_0 (one antipodal pair of complementary edges), on S^{m-1} at degree d=m-1,
for several L_0 (different valid equatorial labelings with a single +-1 pair). usage: base_labeling.py m [ntries]"""
import sys, time, random
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.labels import complementary_edges, forbidden_pairs
from kyfan.dual import unknowns, sa_dual_system, solve

m = int(sys.argv[1]); ntries = int(sys.argv[2]) if len(sys.argv) > 2 else 3
cxe, cx = SignedComplex(m - 1), SignedComplex(m)
labels = label_set(m - 1); V = violating_pairs(cx, labels); d = m - 1
F = forbidden_pairs(cxe, labels); ne = cxe.n_free
others = [l for l in labels if abs(l) != 1]
def equatorial(seed):
    rng = random.Random(seed); L = [None] * ne
    def ok(i, a): return all((L[j], a) not in F.get((j, i), ()) for j in range(i))
    def rec(i):
        if i == ne: return True
        ls = others[:]; rng.shuffle(ls)
        for a in ls:
            if ok(i, a):
                L[i] = a
                if rec(i + 1): return True
                L[i] = None
        return False
    for v0 in rng.sample(range(ne), ne):     # which free vertex carries +1
        L = [None] * ne; L[v0] = 1
        order = [i for i in range(ne) if i != v0]
        # simple DFS over the remaining vertices in index order
        def rec2(k):
            if k == len(order): return True
            i = order[k]; ls = others[:]; rng.shuffle(ls)
            for a in ls:
                if all((L[j], a) not in F.get((min(i, j), max(i, j)), ()) if j < i else (a, L[j]) not in F.get((i, j), ()) for j in range(ne) if L[j] is not None and j != i):
                    L[i] = a
                    if rec2(k + 1): return True
                    L[i] = None
            return False
        if rec2(0): return list(L)
    return None
seen = set()
for seed in range(ntries):
    Leq = equatorial(seed)
    if Leq is None or tuple(Leq) in seen: continue
    seen.add(tuple(Leq)); assert complementary_edges(cxe, Leq) == 0
    L0 = [1 if not any(v[:m - 1]) else Leq[cxe.fidx[v[:m - 1]]] for v in cx.free]
    print(f"m={m} seed {seed}: L_eq={Leq}  L_0 complementary edges={complementary_edges(cx, L0)}", flush=True)
    for h in range(d):
        t = time.time()
        support = lambda a, h=h: sum(1 for i, l in a if L0[i] != l) <= h
        col, _ = unknowns(cx, labels, d, V, support=support)
        rows, ncols, _ = sa_dual_system(cx, labels, d, V=V, support=support, col=col)
        res = solve(rows, ncols, "F2")
        print(f"   h={h}: unknowns {ncols:8d} -> {'solution, %d free' % res.n_free if res.consistent else 'NO solution'}  [{time.time()-t:.0f}s]", flush=True)
