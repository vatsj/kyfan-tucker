"""Ky Fan lower bound (docs/kyfan_lower_bound.md): the separating functional
    E(f) = sum_{lambda(x_i) in S_i} f(lambda)   (mod 2)
over a box on one top simplex sigma_0 = (x_1 < ... < x_{n+1}), S_1 = {+1,+2}, S_i = {+i,-i}, other vertices fixed.
Even |S_i| kills every function of <= n vertices; only +-sigma_0 contribute to E(A_+), and the box has exactly one
alternating tuple, so E(A_+ + 1) = 1: no degree-n certificate for k >= n+1.
"""
import itertools
import random

from .labels import label_set, pos_alt_count
from .violating import partial_labelings


def box(cx, S1=(1, 2), sigma0=None, base_label=1):
    """The labelings in the box: sigma_0's vertices range over S_1 x S_2 x ... , all other free vertices = base_label."""
    n = cx.n
    sigma0 = cx.top[0] if sigma0 is None else sigma0
    reps = [cx.rep(v) for v in sigma0]
    S = [list(S1)] + [[i, -i] for i in range(2, n + 2)]
    out = []
    for choice in itertools.product(*S):
        L = [base_label] * cx.n_free
        for (i, s), l in zip(reps, choice):
            L[i] = s * l          # lambda(x_j) = l  =>  lambda(free rep) = s*l
        out.append(L)
    return out


def E(box_labelings, f):
    return sum(f(L) for L in box_labelings) % 2


def check(cx, S1=(1, 2), samples=3000, seed=0):
    """(E(A_+), E(1), #size-n indicators with E != 0 among those tested, #tested).
    Indicator test is exhaustive when the number of size-n partial labelings is small, sampled otherwise."""
    n = cx.n
    labels = label_set(n + 1)
    B = box(cx, S1)
    eA = E(B, lambda L: pos_alt_count(cx, L))
    e1 = E(B, lambda L: 1)
    rng = random.Random(seed)
    total = 0
    for _ in partial_labelings(cx.n_free, labels, n):
        total += 1
        if total > samples:
            break
    if total <= samples:
        alphas = partial_labelings(cx.n_free, labels, n)
    else:
        alphas = (tuple(zip(sorted(rng.sample(range(cx.n_free), n)), [rng.choice(labels) for _ in range(n)]))
                  for _ in range(samples))
    bad = tested = 0
    for a in alphas:
        bad += E(B, lambda L, a=a: int(all(L[i] == l for i, l in a)))
        tested += 1
    return eA, e1, bad, tested
