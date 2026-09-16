"""TREE_GADGET_THEOREM.md, statement (i): the explicit rule reproduces the solver's unique pseudo-solution.

For every binary conflict tree, `rule_E` (the closed-form E) equals `unique_E` (build the degree-n SA dual and solve
over F_2), and both have support exactly 3^n. Cases up to n=6 are fast; n=7 (36,409 unknowns) uses the sparse solver.
"""
import pytest
from kyfan.tree_rule import unique_E, rule_E

TREES = [
    ((1, (2, 'x', 'x'), 'x'), 3),
    ((1, (2, (3, 'x', 'x'), 'x'), 'x'), 4),
    ((1, (2, 'x', 'x'), (3, 'x', 'x')), 4),
    ((1, (2, (3, (4, 'x', 'x'), 'x'), 'x'), 'x'), 5),
    ((1, (2, (3, 'x', 'x'), 'x'), (4, 'x', 'x')), 5),
    ((1, (2, (3, 'x', 'x'), (4, 'x', 'x')), 'x'), 5),
    ((1, (2, (3, (4, (5, 'x', 'x'), 'x'), 'x'), 'x'), 'x'), 6),   # caterpillar
    ((1, (2, (3, 'x', 'x'), (4, 'x', 'x')), (5, 'x', 'x')), 6),
]
TREE_N7 = ((1, (2, (3, (4, 'x', 'x'), 'x'), (5, 'x', 'x')), (6, 'x', 'x')), 7)


@pytest.mark.parametrize("tree,n", TREES)
def test_rule_matches_solver(tree, n):
    doms, E, res = unique_E(tree, n, method="sparse")
    assert res.consistent and res.rank == res.ncols          # unique pseudo-solution
    assert E == rule_E(tree, n)
    assert len(E) == 3 ** n


@pytest.mark.slow
def test_rule_matches_solver_n7():
    tree, n = TREE_N7
    doms, E, res = unique_E(tree, n, method="sparse")         # 36,409 unknowns
    assert res.consistent and res.rank == res.ncols
    assert E == rule_E(tree, n)
    assert len(E) == 3 ** n
