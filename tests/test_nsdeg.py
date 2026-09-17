"""Tucker on S^2 (labels +-1,+-2) — F_2 degree exactly 3, F_p not refutable at 3 (HISTORY.md table rows 1, 2)."""
from kyfan.dual import sa_dual_system, solve


def test_tucker_s2_degree2_consistent(cx3, L2):
    rows, ncols, _ = sa_dual_system(cx3, L2, 2)
    assert ncols == 1104
    r2 = solve(rows, ncols, "F2")
    assert r2.consistent and r2.rank == 507 and r2.n_free == 597
    rp = solve(rows, ncols, 1000003)
    assert rp.consistent and rp.rank == 507


def test_tucker_s2_degree3_F2_refutable_Fp_not(cx3, L2):
    rows, ncols, _ = sa_dual_system(cx3, L2, 3)
    assert ncols == 12688
    r2 = solve(rows, ncols, "F2")
    assert not r2.consistent and r2.rank == 8400
    rp = solve(rows, ncols, 1000003)
    assert rp.consistent and rp.rank == 8401
