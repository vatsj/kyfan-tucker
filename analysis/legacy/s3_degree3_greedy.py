"""Task 2 (c): greedy minimal *type* support. Starting from all 808 full-group types, try removing types one at a time
(largest first) and keep the removal if a Z3xZ3-invariant pseudo-solution still exists (exact, odd group).
The result is an inclusion-minimal family of types that suffices — every remaining type is necessary relative to the family."""
import time, pickle
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.group import z3z3_generators
from kyfan.analysis import describe

cx = SignedComplex(4); labels = label_set(3); mags = [1, 2, 3]; V = violating_pairs(cx, labels); d = 3
types, _ = pickle.load(open('analysis/legacy/s3_degree3_types.pkl', 'rb'))
z3 = z3z3_generators(4, mags)
col, members = unknowns(cx, labels, d, V, gens=z3)
rows, ncols, _ = sa_dual_system(cx, labels, d, V=V, gens=z3, col=col)
removed = set()
order = sorted(range(len(types)), key=lambda t: -len(types[t]))
t0 = time.time()
for i, t in enumerate(order):
    cand = removed | {col[a] for a in types[t]}
    rs = [({c: x for c, x in f.items() if c not in cand}, b) for f, b in rows]
    rs = [(f, b) for f, b in rs if f or b]
    res = solve(rs, ncols, "F2")
    if res.consistent:
        removed = cand
    keep = 'REMOVED' if res.consistent else 'kept   '
    print(f"[{i:3d}/{len(types)}] type {t:3d} size {len(types[t]):6d}: {keep}  (rank {res.rank})  [{time.time()-t0:.0f}s]", flush=True)
kept = [t for t in range(len(types)) if not ({col[a] for a in types[t]} <= removed)]
print(f"\nminimal sufficient family: {len(kept)} types, {sum(len(types[t]) for t in kept)} partial labelings of {sum(map(len, types))}")
for t in kept:
    print(f"  type {t:3d} size {len(types[t]):6d}  {describe(cx, types[t][0])}   rep={types[t][0]}")
pickle.dump(kept, open('analysis/legacy/s3_degree3_greedy.pkl', 'wb'))
