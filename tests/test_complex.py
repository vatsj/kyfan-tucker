"""build.py: complex counts, conventions."""
import itertools
import math
from collections import Counter

from kyfan import SignedComplex, label_set, posalt
from kyfan.labels import octahedron_pos_alt, dfs_valid_labelings, pos_alt_count, complementary_edges
from kyfan.group import full_group, act, generated_group, full_group_generators, z3z3_generators, assert_odd_order
from kyfan.violating import violating_pairs


def test_counts_m3(cx3):
    assert (len(cx3.verts), len(cx3.edges), len(cx3.top)) == (26, 72, 48)
    assert cx3.euler() == 2
    assert cx3.n_free == 13


def test_counts_general():
    for m in (2, 3, 4, 5):
        cx = SignedComplex(m)
        assert len(cx.verts) == 3 ** m - 1
        assert cx.n_free == (3 ** m - 1) // 2
        assert len(cx.top) == math.factorial(m) * 2 ** m
        assert len(cx.hemisphere_top()) == len(cx.top) // 2


def test_free_reps_first_nonzero_plus(cx4):
    for v in cx4.free:
        assert next(x for x in v if x) == 1
    for v in cx4.verts:
        i, s = cx4.rep(v)
        assert tuple(s * x for x in cx4.free[i]) == v


def test_pos_alt_convention():
    assert posalt((1, -2, 3)) and posalt((3, -2, 1)) and not posalt((-1, 2, -3))
    assert not posalt((1, -1, 2)) and not posalt((1, 2, 3))
    assert posalt(()) is True


def test_octahedron_degenerate():
    """Unsubdivided octahedron: A_+ == 1 on all 192 valid labelings with labels +-1..+-4."""
    labs = label_set(4)
    vals = Counter(octahedron_pos_alt(*ls) for ls in itertools.product(labs, repeat=3)
                   if len({abs(l) for l in ls}) == 3)
    assert dict(vals) == {1: 192}


def test_tucker_unsat_and_kyfan_parity_m3(cx3):
    assert dfs_valid_labelings(cx3, label_set(2)) == []
    sols = dfs_valid_labelings(cx3, label_set(3), cap=3000)
    assert len(sols) == 3000
    assert all(complementary_edges(cx3, L) == 0 for L in sols[:200])
    assert all(pos_alt_count(cx3, L) % 2 == 1 for L in sols)


def test_group_preserves_violating_pairs_and_edges(cx3, L2, V32):
    """The check that primal.py got wrong: every element of the full group must preserve V."""
    G = full_group(3, [1, 2])
    assert len(G) == 384
    for g in G:
        assert {frozenset(act(cx3, g, tuple(sorted(p)))) for p in V32} == V32


def test_generated_group_orders():
    assert len(generated_group(full_group_generators(3, [1, 2]), 3, [1, 2])) == 384
    assert assert_odd_order(z3z3_generators(4, [1, 2, 3]), 4, [1, 2, 3]) == 9
    assert assert_odd_order(z3z3_generators(3, [1, 2]), 3, [1, 2]) == 3
