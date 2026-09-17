"""Tucker on S^3 (labels +-1..+-3) at degree 3 (HISTORY.md table row 3 lower bound, and the negative facts)."""
import random
import time
import numpy as np
import pytest

from kyfan import label_set, violating_pairs
from kyfan.primal import violating_columns, sampled_system, FastRows, tucker_target
from kyfan.group import full_group_generators, z3z3_generators, assert_odd_order
from kyfan import linalg


@pytest.mark.slow
def test_s3_full_group_symmetric_degree3_none(cx4):
    """Full-group orbit search at d=3: none over F_2 (uninformative, |G| even) and none over R (residual ~7e-2)."""
    labels = label_set(3)
    V = violating_pairs(cx4, labels)
    gens = full_group_generators(4, [1, 2, 3])
    col, members = violating_columns(cx4, labels, 3, V, gens)
    assert len(col) == 299280 and len(members) == 232
    rng = random.Random(5)
    A, b, _, _ = sampled_system(cx4, labels, 3, tucker_target, 4 * 232 + 50, rng, gens, V)
    rows = [int(''.join('1' if v & 1 else '0' for v in r[::-1]), 2) for r in A]
    assert linalg.gf2_bitset_solve(rows, list(b)) is None
    assert linalg.real_lstsq_residual(A, b) > 1e-3


@pytest.mark.veryslow
def test_s3_degree3_F2_no_certificate_z9(cx4):
    """Z_3 x Z_3-averaged (odd order: valid over F_2) unrestricted sampled primal at d=3 is inconsistent:
    33,312 orbits, rank 27,584. ~7 min with the dense solver."""
    labels = label_set(3)
    V = violating_pairs(cx4, labels)
    gens = z3z3_generators(4, [1, 2, 3])
    assert assert_odd_order(gens, 4, [1, 2, 3]) == 9
    col, members = violating_columns(cx4, labels, 3, V, gens)
    assert len(members) == 33312
    fr = FastRows(cx4, labels, 3, col)
    t = time.time()
    M = fr.packed_F2_system(tucker_target, int(1.15 * len(members)), np.random.default_rng(11))
    res = linalg.gf2_dense_packed(M, len(members))
    print(f"z9: rank {res.rank}, consistent={res.consistent}  [{time.time() - t:.0f}s]")
    assert not res.consistent and res.rank == 27584
