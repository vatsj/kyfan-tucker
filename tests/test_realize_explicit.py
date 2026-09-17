"""Statement (ii), explicit (REALIZATION_THEOREM.md): the closed-form labeling lambda(u) = s*(n-q+1) is a valid
equatorial labeling whose pole-chain residual is the reverse-caterpillar binary conflict tree. With statement (i)
(test_rule_vs_solver, test_extension_lemma) and the restriction lemma, Tucker's F_2 degree on S^n is exactly n+1 for
every n. No search: the labeling is a formula, so this is a theorem, checked here for m = 3..8.
"""
import pytest

from kyfan.realize_explicit import check, explicit_Leq, reverse_caterpillar_domains
from kyfan.gadget import verify


@pytest.mark.parametrize("m", [3, 4, 5, 6, 7, 8])
def test_explicit_residual_dual(m):
    """Equatorial validity, exact reverse-caterpillar residual domains, and a consistent unique degree-(m-1) dual
    (=> Tucker F_2 degree on S^{m-1} = m). m=8 is a 95,901-unknown sparse solve (~3 s)."""
    r = check(m, solve=True, method="sparse")
    assert r["equatorial_violations"] == 0
    assert r["domains_match"]
    assert r["unsat"] and r["consistent"] and r["unique"]


@pytest.mark.parametrize("m", [4, 5, 6])
def test_explicit_full_sphere(m):
    """The explicit labeling passes the full-sphere verify (support exactly 3^(m-1))."""
    r = verify(m, explicit_Leq(m))
    assert r["unsat"] and r["consistent"] and r["solution_verified"]
    assert r["rank"] == r["unknowns"] and r["support"] == 3 ** (m - 1)


def test_reverse_caterpillar_shape():
    """The target domains are exactly the reverse-caterpillar leaf paths (pole = {n})."""
    n = 5
    t = reverse_caterpillar_domains(n)
    assert t[0] == {5}
    assert t[1] == {4, -5}
    assert t[2] == {3, -4, -5}
    assert t[3] == {2, -3, -4, -5}
    assert t[4] == {1, -2, -3, -4, -5}
    assert t[5] == {-1, -2, -3, -4, -5}
