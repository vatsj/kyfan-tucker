"""Chain gadget for the S^{m-1} lower bound. Equator S^{m-2} labeled with magnitudes 1..m-2 except one antipodal pair
of top simplices +-sigma carrying +-(m-1); U = e_m + lifts of sigma (a top simplex of S^{m-1} through the pole).
If the residual's degree-(m-1) dual is consistent, the sphere's F_2 degree is >= m (restriction argument, ball.Residual).
usage: gadget.py m [ntries]"""
import sys, time, random
from kyfan import SignedComplex, label_set
from kyfan.labels import forbidden_pairs, complementary_edges
from kyfan.ball import Ball, Residual
from kyfan.complex import leq, rank

m = int(sys.argv[1]); ntries = int(sys.argv[2]) if len(sys.argv) > 2 else 20
cxe = SignedComplex(m - 1); small = label_set(m - 2); big = m - 1
F = forbidden_pairs(cxe, label_set(m - 1))
sigma = cxe.top[0]                       # a top simplex of the equator; all are equivalent under the group
print(f"m={m}: equator S^{m-2}, sigma = {sigma}")
sig_idx = {cxe.rep(v)[0] for v in sigma}
from kyfan.labels import random_valid_labeling
def labeling(seed):
    fixed = {}
    for v in sigma:
        i, s = cxe.rep(v); fixed[i] = -s * big          # lambda(v) = -big on sigma, +big on -sigma
    # the rest uses magnitudes 1..m-2 only: pass the small label set but keep the fixed big labels
    return random_valid_labeling(cxe, label_set(m - 1), fixed=fixed, seed=seed, allowed=small)
results = []
for seed in range(ntries):
    t = time.time(); L = labeling(seed)
    if L is None: print(f"seed {seed}: no labeling found"); continue
    assert complementary_edges(cxe, L) == 0
    B = Ball(m, L); em = B.e_m()
    U = [em] + [i for i, y in enumerate(B.cap) if any(y[:-1]) and cxe.label(y[:-1], L) == -big]
    assert len(U) == m
    R = Residual(m, L, U)
    d, log = R.degree(m - 1, verbose=False)
    print(f"seed {seed}: L_eq found [{time.time()-t:.1f}s]; U={[B.cap[i] for i in U]} domains={[R.domains[k] for k in range(R.n)]} -> residual degree {'%d' % d if d else '> %d  ==> sphere degree >= %d' % (m-1, m)}   dual log {[(dd,nc,rk,c) for dd,nc,_,rk,c in log]}", flush=True)
    results.append((L, U, d))
print(f"\nsummary: {sum(1 for *_, d in results if d is None)} of {len(results)} gadgets have residual degree {m} (= conjectured sphere degree)")
