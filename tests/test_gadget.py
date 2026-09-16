"""Restriction lower bounds (kyfan.gadget): S^3 >= 4 (known) and S^4 >= 5 (NEW, settles the conjecture at n=4)."""
import pytest
from kyfan.gadget import verify, GADGETS


@pytest.mark.parametrize("m,support", [(4, 27), (5, 81)])
def test_chain_gadget_lower_bound(m, support):
    r = verify(m)
    assert len(r["U"]) == m
    assert r["unsat"] and r["proper_subsets_satisfiable"]
    assert r["consistent"] and r["solution_verified"]
    assert r["rank"] == r["unknowns"]            # the degree-(m-1) pseudo-solution is unique
    assert r["support"] == support               # 3^(m-1)


def test_gadget_m5_domains():
    r = verify(5)
    assert r["domains"] == {(-1, -1, -1, -1, 1): [1, -2, -4], (-1, -1, -1, 0, 1): [1, 2, -4],
                            (-1, -1, 0, 0, 1): [-1, 3, -4], (-1, 0, 0, 0, 1): [-1, -3, -4], (0, 0, 0, 0, 1): [4]}


def test_gadget_degree_is_sharp_m4():
    """At degree m the residual dual is inconsistent (the trivial degree-|U| certificate), as it must be."""
    r = verify(4, d=4)
    assert not r["consistent"]


@pytest.mark.slow
@pytest.mark.parametrize("m,support", [(6, 243), (7, 729)])
def test_chain_gadget_lower_bound_high(m, support):
    """S^5 >= 6 (5 s) and S^6 >= 7 (~20 s): tree gadgets realized by analysis/gadget3.py."""
    r = verify(m)
    assert len(r["U"]) == m and r["unsat"] and r["proper_subsets_satisfiable"]
    assert r["consistent"] and r["solution_verified"] and r["rank"] == r["unknowns"] and r["support"] == support


def test_abstract_tree_gadgets():
    """Binary conflict trees + pole: degree = number of variables, unique pseudo-solution of support 3^n."""
    from kyfan.abstract_gadget import degree, unsat
    def tree_domains(tree, path=()):
        if tree == 'x':
            return [set(path)]
        mag, left, right = tree
        return tree_domains(left, path + (mag,)) + tree_domains(right, path + (-mag,))
    for n, tree in ((3, (1, (2, 'x', 'x'), 'x')), (4, (1, (2, 'x', 'x'), (3, 'x', 'x'))),
                    (5, (1, 'x', (2, 'x', (3, 'x', (4, 'x', 'x')))))):
        doms = [{n}] + [D | {-n} for D in tree_domains(tree)]
        assert unsat(doms)
        d, info = degree(doms)
        assert d == n + 1
        assert info[-2][2] == info[-2][1] and info[-2][4] == 3 ** n     # unique, support 3^n
