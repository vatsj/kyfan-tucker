"""Fast qualitative checks of the Spread-Lemma experiment (analysis/spread_mcmc.py; docs/spread_evidence.md).

Short chains at m = 4, 5 (the stored run uses 4 x 3000 sweeps for m = 4..8).  Asserted with generous tolerances:
  * every recorded state is gadget-compatible (valid, agrees with the explicit gadget on sigma, static domain conditions);
  * the mean flag density of the largest magnitude class is below the explicit gadget's 1/(2n) by a factor >= 1.5;
  * the two-chain match probability on random sets decays geometrically over |A| = 1, 2, 3, 4 (fitted beta < 0.9, and
    strictly decreasing).
No digit of analysis/spread_mcmc.out is asserted here; that file is reproduced by rerunning the command in its header.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
import spread_mcmc as SM
from kyfan.realize_explicit import explicit_Leq


def _assert_compatible(S):
    """The current state of a GadgetSampler is in G_m."""
    L = S.L
    for i in S.fixed:
        assert L[i] == S.L0[i]
    lab_pos = {int(l): k for k, l in enumerate(S.labels)}
    for i in range(S.N):
        assert S.static_ok[i, lab_pos[int(L[i])]]
        assert not np.any(S.nbr_sgn[i] * L[S.nbr_idx[i]] == L[i])       # no complementary edge


@pytest.mark.parametrize("m", [4, 5])
def test_chains_stay_gadget_compatible_and_flag_density_drops(m):
    n = m - 1
    S0 = SM.GadgetSampler(m)
    assert list(S0.L0) == list(explicit_Leq(m))
    e = S0.stats()
    assert abs(e["flag_density"] - 1 / (2 * n)) < 1e-12
    traces, samplers = SM.run_chains(m, chains=2, sweeps=300, record_every=2, log=False)
    for S in samplers:
        _assert_compatible(S)
    fd = np.concatenate(traces["flag_density"])
    assert len(fd) == 2 * 113                      # 2 chains x records at sweeps 75, 77, ..., 299
    assert fd.mean() * 1.5 < 1 / (2 * n)
    summ = SM.summarize(traces, 2)
    assert summ["n_mags"]["mean"] == n
    assert summ["hamming"]["mean"] > 0.05


@pytest.mark.parametrize("m", [4, 5])
def test_match_probability_decays_geometrically(m):
    out = SM.match_probabilities(m, nsamp=400, burn=100, thin=3, sizes=(1, 2, 3, 4), trials=20, ntargets=100, log=False)
    sizes, means, errs = out["random"]
    assert all(p == p and 0 < p <= 1 for p in means)
    assert all(a > b for a, b in zip(means, means[1:]))
    beta = float(np.exp(np.polyfit(sizes, np.log(means), 1)[0]))
    assert beta < 0.9
