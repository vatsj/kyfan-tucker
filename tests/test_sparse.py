"""Task 1: the sparse GF(2) solver (gf2solve/) reproduces the dense verdicts and settles S^3 d=3 via the dual."""
import time
import numpy as np
import pytest

from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.dual import sa_dual_system, unknowns, solve
from kyfan.group import z3z3_generators, assert_odd_order
from kyfan import linalg


def _check_solution(rows, sol):
    return all(sum(int(sol[c]) * x for c, x in f.items()) % 2 == b % 2 for f, b in rows)


def test_sparse_matches_dense_m3(cx3, L2):
    for d, rank, consistent in ((2, 507, True), (3, 8400, False)):
        rows, ncols, _ = sa_dual_system(cx3, L2, d)
        rs = linalg.gf2_sparse(rows, ncols, want_solution=True)
        rd = linalg.gf2_dense(rows, ncols, want_solution=True)
        assert (rs.consistent, rs.rank) == (rd.consistent, rd.rank) == (consistent, rank)
        if consistent:
            assert _check_solution(rows, rs.solution) and _check_solution(rows, rd.solution)


def test_sparse_nullspace_m3_d2(cx3, L2):
    rows, ncols, _ = sa_dual_system(cx3, L2, 2)
    rs = linalg.gf2_sparse(rows, ncols, want_solution=True, want_null=True)
    assert len(rs.null_basis) == rs.n_free == 597
    hom = [(f, 0) for f, _ in rows]
    for fc, supp in rs.null_basis[:20]:
        v = np.zeros(ncols, dtype=np.uint8)
        v[supp] = 1
        assert v[fc] == 1 and _check_solution(hom, v)


def test_sparse_random_systems_vs_dense():
    rng = np.random.default_rng(0)
    for trial in range(30):
        R, C, k = rng.integers(5, 60), rng.integers(5, 60), rng.integers(1, 6)
        rows = [({int(c): 1 for c in rng.choice(C, size=min(k, C), replace=False)}, int(rng.integers(0, 2))) for _ in range(R)]
        rs = linalg.gf2_sparse(rows, C, want_solution=True)
        rd = linalg.gf2_dense(rows, C, want_solution=True)
        assert (rs.consistent, rs.rank) == (rd.consistent, rd.rank), trial
        if rs.consistent:
            assert _check_solution(rows, rs.solution)


def test_s3_degree3_pseudo_solution_via_dual(cx4):
    """Tucker S^3 lower bound (results.md #3) through the dual: the Z_3 x Z_3-invariant degree-3 SA-dual is consistent
    (valid over F_2 since |G| = 9 is odd: averaging a pseudo-solution gives an invariant one). The solution is then
    checked against every unreduced consistency equation. ~30 s, replaces the 7-min sampled primal (z9)."""
    labels = label_set(3)
    V = violating_pairs(cx4, labels)
    gens = z3z3_generators(4, [1, 2, 3])
    assert assert_odd_order(gens, 4, [1, 2, 3]) == 9
    col, members = unknowns(cx4, labels, 3, V, gens=gens)
    assert (len(col), len(members)) == (1834800, 204048)
    rows, ncols, _ = sa_dual_system(cx4, labels, 3, V=V, gens=gens, col=col)
    assert len(rows) == 111233
    res = solve(rows, ncols, "F2", want_solution=True)
    assert res.consistent and res.rank == 94064
    full_rows, _, _ = sa_dual_system(cx4, labels, 3, V=V, gens=gens, col=col, reduce_rows=False)
    assert len(full_rows) == 987265
    assert _check_solution(full_rows, res.solution)
