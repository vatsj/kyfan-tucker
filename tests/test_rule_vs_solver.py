"""docs/tree_gadget.md, statement (i): the explicit rule reproduces the solver's unique pseudo-solution.

For every unordered binary conflict tree shape, `rule_E` (the closed-form E) equals `unique_E` (build the degree-n SA
dual and solve over F_2), and both have support exactly 3^n. All shapes with n <= 6 (fast tier); all 11 shapes with
n = 7 (36,409 unknowns each, sparse solver; slow tier, ~16 s).
"""
import pytest
from kyfan.tree_rule import unique_E, rule_E, all_tree_shapes

TREES = [(t, n) for n in (3, 4, 5, 6) for t in all_tree_shapes(n)]
TREES_N7 = [(t, 7) for t in all_tree_shapes(7)]


@pytest.mark.parametrize("tree,n", TREES)
def test_rule_matches_solver(tree, n):
    doms, E, res = unique_E(tree, n, method="sparse")
    assert res.consistent and res.rank == res.ncols          # unique pseudo-solution
    assert E == rule_E(tree, n)
    assert len(E) == 3 ** n


@pytest.mark.slow
@pytest.mark.parametrize("tree,n", TREES_N7)
def test_rule_matches_solver_n7(tree, n):
    doms, E, res = unique_E(tree, n, method="sparse")         # 36,409 unknowns
    assert res.consistent and res.rank == res.ncols
    assert E == rule_E(tree, n)
    assert len(E) == 3 ** n
