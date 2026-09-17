"""Task 2 (a): full-group orbit types of non-violating size-3 partial labelings on S^3 and the knockout analysis.
Run: .venv/bin/python analysis/legacy/s3_degree3_types.py  (~15 min). Output also saved to analysis/legacy/s3_degree3_types.out."""
import time, pickle, sys
from collections import Counter
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.group import z3z3_generators
from kyfan.analysis import full_types, knockout, describe

cx = SignedComplex(4); labels = label_set(3); mags = [1, 2, 3]; V = violating_pairs(cx, labels); d = 3
t = time.time(); type_of, types = full_types(cx, labels, d, V, mags)
print(f"{len(type_of)} non-violating size-3 partial labelings in {len(types)} full-group orbit types [{time.time()-t:.0f}s]", flush=True)
print("type sizes:", Counter(len(m) for m in types).most_common())
gens = z3z3_generators(4, mags)
col, members = unknowns(cx, labels, d, V, gens=gens)
rows, ncols, _ = sa_dual_system(cx, labels, d, V=V, gens=gens, col=col)
t = time.time()
res = knockout(cx, labels, d, V, gens, types, base_rows=rows, col=col, verbose=True)
necessary = [t_ for t_, ok, _ in res if not ok]
print(f"\nknockout done [{time.time()-t:.0f}s]. NECESSARY types ({len(necessary)} of {len(types)}):")
for t_ in necessary:
    print(f"  type {t_:4d} size {len(types[t_]):7d}  {describe(cx, types[t_][0])}   rep={types[t_][0]}")
pickle.dump((types, res), open('analysis/legacy/s3_degree3_types.pkl', 'wb'))
