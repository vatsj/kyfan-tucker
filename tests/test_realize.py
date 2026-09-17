"""Statement (ii), realizability: a deterministic, backtrack-free construction of an equatorial labeling of S^{m-1}
whose pole-chain residual is the caterpillar binary conflict tree — so `kyfan.tree_rule` (statement (i)) applies and
Tucker's F_2 degree on S^{m-1} is exactly m. This gives a second path to the S^4..S^7 lower bounds, independent of the
searched `GADGETS`, and one that is uniform in m. `chain_domains` (no full sphere) is validated against `residual`.
"""
import pytest

from kyfan.gadget import (realize, check_realization, chain_domains, residual, chain_U, verify, GADGETS)
from kyfan.complex import rank


@pytest.mark.parametrize("m", [4, 5, 6, 7, pytest.param(8, marks=pytest.mark.slow)])
def test_chain_domains_matches_residual(m):
    """The fast up/down-set domain computation reproduces the full-sphere residual on the stored gadgets."""
    n = m - 1
    Leq = GADGETS[m]
    fast = chain_domains(m, Leq)
    U = chain_U(m, Leq)
    _, dom, _ = residual(m, Leq, U)
    truth = {rank(y) - 1: set(dom[y]) for y in U}
    assert all(fast[r] == truth[r] for r in range(n + 1))


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_caterpillar_realize(m):
    """Deterministic caterpillar realizer: no backtracking, and the residual is a degree-(m-1) gadget."""
    Leq, backtracks = realize(m)
    assert backtracks == 0                       # the greedy construction is forced
    r = check_realization(m, Leq)
    assert r["unsat"] and r["unique"]            # UNSAT, unique degree-(m-1) pseudo-solution


@pytest.mark.parametrize("m", [4, 5, 6])
def test_caterpillar_full_sphere(m):
    """The realized labeling passes the full-sphere verify (support exactly 3^(m-1))."""
    Leq, _ = realize(m)
    r = verify(m, Leq)
    assert r["unsat"] and r["consistent"] and r["solution_verified"]
    assert r["rank"] == r["unknowns"] and r["support"] == 3 ** (m - 1)


@pytest.mark.slow
def test_caterpillar_realize_m7():
    """S^6: caterpillar realized deterministically (11,743-unknown degree-6 dual, sparse solver)."""
    Leq, backtracks = realize(7)
    assert backtracks == 0
    r = check_realization(7, Leq)
    assert r["unsat"] and r["unique"] and r["unknowns"] == 11743


@pytest.mark.veryslow
def test_caterpillar_realize_m8():
    """S^7: caterpillar realized deterministically (~45 s for the equatorial CSP, then a 95,901-unknown dual)."""
    Leq, backtracks = realize(8)
    assert backtracks == 0
    r = check_realization(8, Leq)
    assert r["unsat"] and r["unique"]
