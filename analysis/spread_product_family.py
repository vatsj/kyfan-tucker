"""A rigorous product sub-family of gadget labelings.

Given a base labeling L in G_m, find a set I of vertices and an alternative label alt[i] != L[i] for each i in I such
that EVERY one of the 2^|I| labelings obtained by switching any subset of I to its alternative lies in G_m. (Checked
pairwise: a labeling is in G_m iff each vertex satisfies its static domain condition and no edge is complementary, and
both are conditions on single vertices and pairs.)  Then the uniform distribution on these 2^|I| gadgets satisfies
P[match a on A] <= 2^(-|A ∩ I|) exactly. Greedy, so |I| is a lower bound on the best such set.

Usage:  PYTHONPATH=. python analysis/spread_product_family.py --m 5 6 7 8
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from analysis.spread_mcmc import GadgetSampler


def product_family(g, rng, tries=20):
    best = {}
    for _ in range(tries):
        alt = {}
        for i in rng.permutation(g.movable):
            cand = [int(a) for a in g.allowed(i) if a != g.L[i]]
            rng.shuffle(cand)
            for a in cand:
                # a must not be complementary to the alternative of any chosen neighbour
                if all(int(j) not in alt or a != s * alt[int(j)] for j, s in zip(g.nbr_idx[i], g.nbr_sgn[i])):
                    alt[int(i)] = a
                    break
        if len(alt) > len(best):
            best = alt
    return best


def verify(g, alt, rng, trials=300):
    """Independent check: random subsets of I switched, full validity + static conditions re-tested."""
    L0 = g.L.copy()
    lab_pos = {int(l): k for k, l in enumerate(g.labels)}
    I = list(alt)
    for _ in range(trials):
        L = L0.copy()
        for i in I:
            if rng.random() < 0.5:
                L[i] = alt[i]
        for i in g.movable:
            assert g.static_ok[i, lab_pos[int(L[i])]]
            assert not np.any(L[i] == g.nbr_sgn[i] * L[g.nbr_idx[i]])
    return True


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, nargs="+", default=[5, 6, 7])
    ap.add_argument("--sweeps", type=int, default=300)
    a = ap.parse_args()
    for m in a.m:
        t0 = time.time()
        g = GadgetSampler(m, seed=1)
        rng = np.random.default_rng(7)
        alt0 = product_family(g, rng); verify(g, alt0, rng)
        g.sweep(a.sweeps)
        alt1 = product_family(g, rng); verify(g, alt1, rng)
        M = len(g.movable)
        print(f"m={m}: movable {M};  product family around the explicit gadget: |I| = {len(alt0)} ({len(alt0)/M:.0%});  "
              f"around a sampled gadget: |I| = {len(alt1)} ({len(alt1)/M:.0%})   [verified on 300 random switches each; {time.time()-t0:.0f}s]")
        sys.stdout.flush()
