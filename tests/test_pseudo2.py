"""pseudo2.py / pseudo2b.py / pseudo2c.py: structure of the degree-2 pseudo-solution on S^2 (results.md, structural facts)."""
from kyfan.dual import unknowns, sa_dual_system, solve
from kyfan.group import (full_group_generators, hemisphere_stabilizer_generators, label_group_generators,
                         z3z3_generators)
from kyfan.complex import leq, neg, rank


def pair_type(cx, a):
    """Coarse full-group invariant of a size-2 partial labeling: ('edge'|'unrelated', sorted ranks, 'same'|'diff' magnitude)."""
    (i, li), (j, lj) = a
    x, y = cx.free[i], cx.free[j]
    edge = any(leq(u, w) or leq(w, u) for u in (x, neg(x)) for w in (y,))
    return ('edge' if edge else 'unrelated', tuple(sorted((rank(x), rank(y)))), 'same' if abs(li) == abs(lj) else 'diff')


def _solve(cx, L2, gens=None, support=None, want=False):
    rows, ncols, col = sa_dual_system(cx, L2, 2, gens=gens, support=support)
    return solve(rows, ncols, "F2", want_solution=want), ncols, col


def test_unrestricted(cx3, L2):
    res, ncols, _ = _solve(cx3, L2)
    assert ncols == 1104 and res.consistent and res.rank == 507 and res.n_free == 597


def test_no_full_group_invariant_solution(cx3, L2):
    res, ncols, _ = _solve(cx3, L2, gens=full_group_generators(3, [1, 2]))
    assert ncols == 20 and not res.consistent


def test_restrictions(cx3, L2):
    distinct = lambda a: abs(a[0][1]) != abs(a[1][1])
    res, ncols, _ = _solve(cx3, L2, support=distinct)
    assert ncols == 624 and not res.consistent
    res, ncols, _ = _solve(cx3, L2, gens=z3z3_generators(3, [1, 2]))          # coordinate Z_3 only (k=2)
    assert ncols == 368 and res.consistent and res.n_free == 201
    res, ncols, _ = _solve(cx3, L2, gens=hemisphere_stabilizer_generators(3, [1, 2]))
    assert ncols == 44 and not res.consistent
    res, ncols, _ = _solve(cx3, L2, gens=label_group_generators(3, [1, 2]))
    assert ncols == 198 and not res.consistent
    res, ncols, _ = _solve(cx3, L2, gens=label_group_generators(3, [1, 2]), support=distinct)
    assert ncols == 78 and not res.consistent
    res, ncols, _ = _solve(cx3, L2, gens=hemisphere_stabilizer_generators(3, [1, 2]), support=distinct)
    assert ncols == 19 and not res.consistent


def test_knockout_necessary_types(cx3, L2, V32):
    """Removing any one full-group orbit type: solution survives except for exactly four types."""
    _, members = unknowns(cx3, L2, 2, V32, gens=full_group_generators(3, [1, 2]))
    assert len(members) == 20
    necessary = []
    for mem in members:
        ms = set(mem)
        res, _, _ = _solve(cx3, L2, support=lambda a, ms=ms: a not in ms)
        if not res.consistent:
            necessary.append(pair_type(cx3, mem[0]))
    assert sorted(necessary) == sorted([('edge', (1, 2), 'diff'), ('edge', (1, 3), 'diff'), ('edge', (2, 3), 'diff'),
                                        ('unrelated', (1, 1), 'same')])
    # each of those coarse types is a single orbit, so the coarse description identifies the orbit
    types = [pair_type(cx3, mem[0]) for mem in members]
    for t in necessary:
        assert types.count(t) == 1
