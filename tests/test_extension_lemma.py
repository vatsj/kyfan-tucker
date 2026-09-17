"""docs/tree_gadget.md, Lemma (star): for S a proper subset of the leaves with |S| <= n-2, a leaf l not in S, and
any assignment beta on S, the number of picks x making beta u {l -> x} an (S u {l})-pattern is 1 if beta is an S-pattern
and 0 or 2 otherwise. This is exactly what makes E consistent. Verified exhaustively for every unordered binary tree
shape with n <= 6 leaves (1 + 2 + 3 + 6 shapes; fast tier) and all 11 shapes with n = 7 (slow tier, ~40 s).
"""
from collections import Counter

import pytest

from kyfan.tree_rule import leaves_and_paths, rule_patterns, all_tree_shapes

TREES = [(t, n) for n in (3, 4, 5, 6) for t in all_tree_shapes(n)]
TREES_N7 = [(t, 7) for t in all_tree_shapes(7)]


def _check_star(tree, n):
    import itertools
    leaves, paths = leaves_and_paths(tree)
    doms = {j: [s * m for m, s in paths[j - 1]] for j in leaves}     # tree-pick labels only (no P here)
    cnt = Counter()
    violation = None
    for k in range(0, n - 1):                                        # |S| <= n-2
        for S in itertools.combinations(leaves, k):
            S = frozenset(S)
            pats = rule_patterns(tree, S)
            for beta_labels in itertools.product(*[doms[j] for j in sorted(S)]):
                beta = frozenset(zip(sorted(S), beta_labels))
                is_pat = beta in pats
                for l in leaves:
                    if l in S:
                        continue
                    P2 = rule_patterns(tree, S | {l})
                    ext = sum(1 for x in doms[l] if beta | {(l, x)} in P2)
                    cnt[(is_pat, ext)] += 1
                    if (is_pat and ext != 1) or (not is_pat and ext % 2):
                        violation = (S, beta, l, ext)
    return cnt, violation


@pytest.mark.parametrize("tree,n", TREES)
def test_extension_lemma(tree, n):
    cnt, violation = _check_star(tree, n)
    assert violation is None
    assert set(cnt) <= {(True, 1), (False, 0), (False, 2)}
    assert cnt[(True, 1)] > 0


@pytest.mark.slow
@pytest.mark.parametrize("tree,n", TREES_N7)
def test_extension_lemma_n7(tree, n):
    cnt, violation = _check_star(tree, n)
    assert violation is None
    assert set(cnt) <= {(True, 1), (False, 0), (False, 2)}
    assert cnt[(True, 1)] > 0
