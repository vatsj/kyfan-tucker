"""Task 2 (a'): marginal knockout on S^3, degree 3. For each full-group type of size-k partial labelings (k=1,2),
impose E[beta] = 0 for all beta of that type (extra rows sum_l E[beta+(v0,l)] = 0) on the Z3xZ3-invariant system."""
import time, pickle
from collections import defaultdict
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.group import z3z3_generators, full_group_generators, orbits, act
from kyfan.violating import partial_labelings, is_violating
from kyfan.analysis import describe

cx = SignedComplex(4); labels = label_set(3); mags = [1, 2, 3]; V = violating_pairs(cx, labels); d = 3; n = cx.n_free
z3 = z3z3_generators(4, mags); full = full_group_generators(4, mags)
col, members = unknowns(cx, labels, d, V, gens=z3)
rows, ncols, _ = sa_dual_system(cx, labels, d, V=V, gens=z3, col=col)
base = solve(rows, ncols, "F2"); print("base: consistent", base.consistent, "rank", base.rank, flush=True)
for k in (1, 2):
    items = [b for b in partial_labelings(n, labels, k) if not is_violating(b, V)]
    _, ktypes = orbits(items, full, lambda g, a: act(cx, g, a))
    print(f"\nsize-{k} marginals: {len(items)} partial labelings in {len(ktypes)} full-group types", flush=True)
    necessary = []
    for t, mem in enumerate(ktypes):
        _, zorbs = orbits(mem, z3, lambda g, a: act(cx, g, a))   # one extra row per Z3xZ3 orbit of beta
        extra = []
        for zo in zorbs:
            f = defaultdict(int)
            def expand(beta):   # E[beta] as a sum of size-d unknowns, extending along the smallest unused vertex
                if len(beta) == d:
                    c = col.get(beta)
                    if c is not None: f[c] += 1
                    return
                used = {v for v, _ in beta}; v0 = min(v for v in range(n) if v not in used)
                for l in labels: expand(tuple(sorted(beta + ((v0, l),))))
            expand(zo[0])
            f = {c: x for c, x in f.items() if x % 2}
            if f: extra.append((f, 0))
        res = solve(rows + extra, ncols, "F2")
        tag = 'solution' if res.consistent else 'NO solution'
        if not res.consistent: necessary.append(t)
        print(f"  E=0 on type {t:3d} (size {len(mem):6d}, {len(extra)} rows): {tag:12s} rank {res.rank}  [{describe(cx, mem[0])}]", flush=True)
    print(f"NECESSARY size-{k} marginal types ({len(necessary)} of {len(ktypes)}):")
    for t in necessary: print(f"   {describe(cx, ktypes[t][0])}   rep={ktypes[t][0]}")
