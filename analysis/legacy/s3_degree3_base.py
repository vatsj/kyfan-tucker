"""Task 2 (c): pseudo-solutions near a base labeling L_0 with exactly one antipodal pair of complementary edges.
L_0 = pullback of a valid equatorial labeling (S^2, labels +-1..+-3) using -1 on exactly one antipodal pair, with +1 on e_4.
Support restriction: E[alpha] = 0 unless alpha disagrees with L_0 in <= h positions (h = 0,1,2; h=3 is unrestricted).
This is implied by (weaker than) mu being supported within Hamming distance h of L_0."""
import time, sys
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.labels import dfs_valid_labelings, complementary_edges, forbidden_pairs
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.complex import neg

cx3, cx4 = SignedComplex(3), SignedComplex(4)
labels = label_set(3); V4 = violating_pairs(cx4, labels)
# equatorial labeling: free vertex 0 gets +1, all others in {+-2,+-3}
F = forbidden_pairs(cx3, labels); n3 = cx3.n_free
L = [None] * n3; L[0] = 1
def ok(i, a): return all((L[j], a) not in F.get((j, i), ()) for j in range(i))
def rec(i):
    if i == n3: return True
    for a in (2, -2, 3, -3):
        if ok(i, a):
            L[i] = a
            if rec(i + 1): return True
            L[i] = None
    return False
assert rec(1), "no equatorial labeling with a single +-1 pair"
Leq = list(L); print("equatorial labeling:", Leq, " complementary edges:", complementary_edges(cx3, Leq))
def pi(v): return v[:3]
L0 = []
for v in cx4.free:
    p = pi(v)
    L0.append(1 if not any(p) else Leq[cx3.fidx[p]])
print("L_0 complementary edges on S^3:", complementary_edges(cx4, L0), "(expect 2 = one antipodal pair)")
for h in (0, 1, 2):
    t = time.time()
    support = lambda a, h=h: sum(1 for i, l in a if L0[i] != l) <= h
    col, members = unknowns(cx4, labels, 3, V4, support=support)
    rows, ncols, _ = sa_dual_system(cx4, labels, 3, V=V4, support=support, col=col)
    res = solve(rows, ncols, "F2", want_solution=True)
    print(f"h={h}: unknowns {ncols:8d} rows {len(rows):8d} -> {'solution, rank %d, %d free, support %d' % (res.rank, res.n_free, int(res.solution.sum())) if res.consistent else 'NO solution'}  [{time.time()-t:.0f}s]", flush=True)
