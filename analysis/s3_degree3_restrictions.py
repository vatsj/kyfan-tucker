"""Task 2 (b): restricted degree-3 pseudo-solutions on S^3. Odd-group (Z3xZ3) averaging makes the support tests exact;
the label-group / hemisphere-stabilizer tests are about *invariant* solutions only (even order)."""
import time
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.group import z3z3_generators, label_group_generators, hemisphere_stabilizer_generators, generated_group

cx = SignedComplex(4); labels = label_set(3); mags = [1, 2, 3]; V = violating_pairs(cx, labels); d = 3
z3 = z3z3_generators(4, mags)
distinct = lambda a: len({abs(l) for _, l in a}) == 3
no_equal_pair = distinct
def run(label, gens, support=None):
    t = time.time()
    col, members = unknowns(cx, labels, d, V, gens=gens, support=support)
    rows, ncols, _ = sa_dual_system(cx, labels, d, V=V, gens=gens, support=support, col=col)
    res = solve(rows, ncols, "F2")
    print(f"{label:60s} unknowns={ncols:7d} rows={len(rows):7d} -> {'solution, rank %d, %d free' % (res.rank, res.n_free) if res.consistent else 'NO solution'}  [{time.time()-t:.0f}s]", flush=True)
run("unrestricted (Z3xZ3-averaged, exact)", z3)
run("support: all three magnitudes distinct (Z3xZ3, exact)", z3, distinct)
lg = label_group_generators(4, mags); print("label group order", len(generated_group(lg, 4, mags)))
run("label-group invariant (signed perms of magnitudes)", lg)
hs = hemisphere_stabilizer_generators(4, mags); print("hemisphere stabilizer order", len(generated_group(hs, 4, mags)))
run("hemisphere-stabilizer invariant", hs)
run("label-group invariant + distinct magnitudes", lg, distinct)
run("hemisphere-stabilizer invariant + distinct magnitudes", hs, distinct)
