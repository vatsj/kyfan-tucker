"""Explicit certificates from the random-point primal (HISTORY.md table rows 1, 2 (d = 3, 4), 4)."""
import random
import pytest

from kyfan.primal import (sampled_primal, tucker_target, kyfan_target, verify_random, verify_exhaustive_F2,
                          kyfan_target_monomials, verify_exact_F2)
from kyfan.group import full_group_generators


def test_tucker_s2_F2_degree3_certificate_exact(cx3, L2, V32):
    rng = random.Random(1)
    cert, ncols = sampled_primal(cx3, L2, 3, tucker_target, 9000, rng, V=V32)
    assert ncols == 5616
    assert cert is not None
    assert verify_random(cx3, cert, tucker_target, L2, 2000, rng) == 0
    # exact: canonical one-hot basis expansion, target = the constant 1
    assert verify_exact_F2(cert, {(): 1}, L2)


@pytest.mark.slow
def test_tucker_s2_F2_degree3_certificate_exhaustive(cx3, L2, V32):
    rng = random.Random(1)
    cert, _ = sampled_primal(cx3, L2, 3, tucker_target, 9000, rng, V=V32)
    assert verify_exhaustive_F2(cx3, cert, L2, 1)


def test_tucker_s2_Fp_symmetric_degree3_4_none(cx3, L2, V32):
    """Full-group orbit primal over F_p (valid: |G| = 384 is invertible mod p). Inconsistent sampled system is sound."""
    gens = full_group_generators(3, [1, 2])
    rng = random.Random(1)
    cert, norb = sampled_primal(cx3, L2, 3, tucker_target, 4000, rng, field=1000003, gens=gens, V=V32)
    assert norb == 55 and cert is None
    cert, norb = sampled_primal(cx3, L2, 4, tucker_target, 4000, rng, field=1000003, gens=gens, V=V32)
    assert norb == 657 and cert is None


def test_kyfan_s2_F2_degree3_certificate_exact(cx3, L3, V33):
    rng = random.Random(1)
    target = kyfan_target(cx3)
    cert, ncols = sampled_primal(cx3, L3, 3, target, 22000, rng, V=V33)
    assert ncols == 13176
    assert cert is not None
    assert verify_random(cx3, cert, target, L3, 1000, rng) == 0
    assert verify_exact_F2(cert, kyfan_target_monomials(cx3, L3), L3)


@pytest.mark.slow
def test_tucker_s2_F2_full_group_symmetric_none_d3_d4(cx3, L2, V32):
    """Negative/methodological: fully-symmetric F_2 certificates do not exist at d=3,4 although an asymmetric d=3 one does.
    (|G| = 384 is even, so this says nothing about unrestricted certificates.)"""
    import numpy as np
    from kyfan.primal import violating_columns, FastRows
    from kyfan import linalg
    gens = full_group_generators(3, [1, 2])
    for d, norb in ((3, 55), (4, 657)):
        col, members = violating_columns(cx3, L2, d, V32, gens)
        assert len(members) == norb
        M = FastRows(cx3, L2, d, col).packed_F2_system(tucker_target, 4 * norb + 50, np.random.default_rng(5))
        assert not linalg.gf2_dense_packed(M, norb).consistent


@pytest.mark.veryslow
def test_tucker_s2_Q_degree5_symmetric_certificate(cx3, L2, V32):
    """HISTORY.md table row 2, d=5: 924,480 violating unknowns in 5,482 full-group orbits; float rank 1,976; residual ~1e-14;
    fresh-labeling check; and an exact certificate mod p = 1000003. ~5 min (the dense mod-p solve dominates)."""
    import numpy as np
    from kyfan.primal import violating_columns, FastRows
    from kyfan import linalg
    gens = full_group_generators(3, [1, 2])
    col, members = violating_columns(cx3, L2, 5, V32, gens)
    assert (len(col), len(members)) == (924480, 5482)
    fr = FastRows(cx3, L2, 5, col)
    rng = np.random.default_rng(7)
    nrows = int(1.6 * len(members))
    A = np.array([fr.counts(rng.integers(0, 4, cx3.n_free)) for _ in range(nrows)], dtype=np.int64)
    b = np.ones(nrows)
    x, _, rk, _ = np.linalg.lstsq(A.astype(float), b, rcond=None)
    assert rk == 1976
    assert np.linalg.norm(A @ x - b) / np.sqrt(nrows) < 1e-12
    fresh = np.array([fr.counts(rng.integers(0, 4, cx3.n_free)) for _ in range(1000)], dtype=np.int64)
    assert np.abs(fresh @ x - 1).max() < 1e-9
    assert linalg.fp_dense_solve(A, np.ones(nrows, dtype=np.int64)) is not None
