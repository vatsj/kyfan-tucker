"""Per-vertex marginals of gadget-compatible labelings: is a constant fraction of the sphere forced?

For the Spread Lemma (docs/size_lower_bound.md, section 4) to hold for a distribution on G_m, at most O(m) vertices may
have a (near-)deterministic label. This script measures, on the family G_m of analysis/spread_mcmc.py:

  * RIGOROUS: the static domain size of every movable vertex (labels not excluded by the gadget domain conditions).
    A vertex with static domain size 1 has the same label in every member of G_m.
  * SAMPLED: the per-vertex marginal of the Glauber chain, p_max(v) = max_l P[L(v) = l], pooled over chains, and the
    largest disagreement between chains (a mixing sanity check; ergodicity is not proved, so "forced in the samples"
    is evidence, while "free in the samples" is a fact about G_m).

Usage:  PYTHONPATH=. python analysis/spread_marginals.py --m 4 5 6 7 8 --chains 4 --sweeps 1500
"""
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np

from kyfan.complex import leq, neg, rank
from analysis.spread_mcmc import GadgetSampler


def run(m, chains, sweeps, burn, thin, seed):
    t0 = time.time()
    samplers = [GadgetSampler(m, seed=seed + c) for c in range(chains)]
    g = samplers[0]
    n, N, mv = g.n, g.N, g.movable
    nlab = len(g.labels)
    lab_pos = {int(l): k for k, l in enumerate(g.labels)}
    static_size = g.static_ok[mv].sum(axis=1)
    counts = np.zeros((chains, N, nlab), dtype=np.int64)
    nsamp = 0
    for c, s in enumerate(samplers):
        s.sweep(burn)
    for it in range(sweeps // thin):
        for c, s in enumerate(samplers):
            s.sweep(thin)
            idx = np.fromiter((lab_pos[int(l)] for l in s.L), dtype=np.int64, count=N)
            counts[c, np.arange(N), idx] += 1
        nsamp += 1
    P = counts / nsamp                                   # per chain
    pooled = P.mean(axis=0)
    pmax = pooled[mv].max(axis=1)
    between = np.abs(P[:, mv, :] - pooled[mv][None]).max()
    # where are the near-forced vertices? comparable to the chain sigma (the "rigid zone") or not
    w = [tuple(-1 if k < r else 0 for k in range(n)) for r in range(1, n + 1)]
    def near_sigma(v):
        return any(leq(vv, x) or leq(x, vv) for vv in (v, neg(v)) for x in w)
    near = np.array([near_sigma(g.cxe.free[i]) for i in mv])
    forced95 = pmax >= 0.95
    minent = float(-np.log2(pmax).sum())
    print(f"m={m} (S^{m-1}; equator has {N} free vertices, {len(g.fixed)} fixed on sigma, {len(mv)} movable)  "
          f"[{time.time()-t0:.0f}s, {chains} chains x {nsamp} samples]")
    print(f"   static domain size (rigorous): =1: {(static_size==1).sum()}, =2: {(static_size==2).sum()}, "
          f">=3: {(static_size>=3).sum()};  median {int(np.median(static_size))} of {nlab} labels")
    print(f"   sampled p_max: >=0.95: {forced95.sum()} ({forced95.mean():.1%} of movable), >=0.8: {(pmax>=0.8).sum()}, "
          f">=0.6: {(pmax>=0.6).sum()};  median {np.median(pmax):.2f}, mean {pmax.mean():.2f}")
    print(f"   of the {forced95.sum()} near-forced: {int((forced95 & near).sum())} comparable to sigma, "
          f"{int((forced95 & ~near).sum())} not;  movable vertices comparable to sigma: {int(near.sum())}")
    print(f"   sum of -log2 p_max over movable vertices: {minent:.0f} bits ({minent/len(mv):.2f} per vertex);  "
          f"max |chain marginal - pooled|: {between:.2f}")
    return dict(m=m, movable=len(mv), forced95=int(forced95.sum()), minent=minent)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, nargs="+", default=[4, 5, 6, 7])
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--sweeps", type=int, default=1500)
    ap.add_argument("--burn", type=int, default=300)
    ap.add_argument("--thin", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    for m in a.m:
        run(m, a.chains, a.sweeps, a.burn, a.thin, a.seed)
        sys.stdout.flush()
