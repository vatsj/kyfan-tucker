"""tower.py + kyfan_lower.py: Ky Fan degree exactly n+1 (results.md #5)."""
import pytest
from kyfan import SignedComplex
from kyfan.tower import local_lemma_violations, check_identity
from kyfan.lower import check, box, E
from kyfan.labels import pos_alt_count


def test_local_lemma_exhaustive():
    for n in range(4):
        for k in range(1, 6):
            bad, _ = local_lemma_violations(n, k)
            assert bad == 0, (n, k)


@pytest.mark.parametrize("m,ncert", [(3, 29), (4, 221), (5, 2141)])
def test_telescoped_identity(m, ncert):
    holds, count = check_identity(SignedComplex(m), trials=200)
    assert holds == 200 and count == ncert


@pytest.mark.parametrize("m", [3, 4, 5])
def test_kyfan_lower_functional(m):
    eA, e1, bad, tested = check(SignedComplex(m), samples=1000)
    assert (eA, e1, bad) == (1, 0, 0) and tested >= 1000


def test_kyfan_lower_needs_collision():
    """The negation-closed box S_1 = {+1,-1} gives E(A_+) = 0 — the collision at magnitude 2 is essential."""
    cx = SignedComplex(4)
    assert E(box(cx, S1=(1, -1)), lambda L: pos_alt_count(cx, L)) == 0
    assert E(box(cx, S1=(1, 2)), lambda L: pos_alt_count(cx, L)) == 1
