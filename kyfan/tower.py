"""The tower theorem (CLAUDE.md §2): local lemma g_n and the telescoped identity
    A_+^{(m)} + 1 = sum_{j=1..m} sum_{sigma in H^{(j)}} g_{j-1}(lambda|_sigma)   (mod 2),
an identity on all labelings whose right side has degree m = n+1 and is supported on complementary tuples.
"""
import itertools
import random

from .complex import SignedComplex
from .labels import label_set, posalt, negalt, complementary, pos_alt_count


def g(ls):
    """Local lemma on an (n+1)-tuple of labels: sum_i [ls minus ls_i pos-alt] + [ls pos-alt] + [ls neg-alt] mod 2."""
    return (sum(posalt(ls[:i] + ls[i + 1:]) for i in range(len(ls))) + posalt(ls) + negalt(ls)) % 2


def local_lemma_violations(n, k):
    """(#non-complementary tuples with g != 0, #complementary tuples with g != 0) over labels +-1..+-k, tuples of size n+1."""
    labs = label_set(k)
    bad = nz = 0
    for ls in itertools.product(labs, repeat=n + 1):
        if g(ls):
            if complementary(ls):
                nz += 1
            else:
                bad += 1
    return bad, nz


def tower_levels(cx):
    """H^{(j)} for j = 1..m: hemisphere top simplices of the [j]-complex, padded with zeros to length m."""
    levels = []
    for j in range(1, cx.m + 1):
        cj = SignedComplex(j)
        levels.append([tuple(v + (0,) * (cx.m - j) for v in c) for c in cj.hemisphere_top()])
    return levels


def tower_rhs(cx, levels, L):
    return sum(g(tuple(cx.label(v, L) for v in c)) for lev in levels for c in lev) % 2


def check_identity(cx, k=None, trials=400, seed=None):
    """Number of random labelings (labels +-1..+-k, default k = m+1) on which the telescoped identity holds."""
    k = cx.m + 1 if k is None else k
    labs = label_set(k)
    rng = random.Random(cx.m if seed is None else seed)
    levels = tower_levels(cx)
    holds = 0
    for _ in range(trials):
        L = [rng.choice(labs) for _ in range(cx.n_free)]
        holds += ((pos_alt_count(cx, L) + 1) % 2 == tower_rhs(cx, levels, L))
    return holds, sum(len(lev) for lev in levels)
