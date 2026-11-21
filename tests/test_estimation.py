"""The estimation experiment of analysis/estimation/README.md: exact enumeration on S^2 and the low-degree
deduction-projection estimator. Every number quoted there is reproduced here (d <= 3 fast; d = 4 slow)."""
import itertools
import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "analysis", "estimation")
sys.path.insert(0, HERE)
from lowdeg import estimate            # noqa: E402
from run_s2 import build               # noqa: E402

Z_S2 = 22193664                        # valid labelings of S^2 with labels +-1..+-3, of 6^13


def test_exact_at_full_degree_on_s1():
    N, q, ef, sf = build(2, 3)
    z = n = 0.0
    for idx in itertools.product(range(q), repeat=N):
        v = 1.0
        for ax, t in ef:
            v *= t[tuple(idx[a] for a in ax)]
        s = 1.0
        for ax, t in sf:
            s *= t[tuple(idx[a] for a in ax)]
        z += v
        n += v * s
    z /= q ** N
    n /= q ** N
    assert n == pytest.approx(-z)                               # Ky Fan on S^1
    assert estimate(ef, N) == pytest.approx(z)
    assert estimate(ef + sf, N) == pytest.approx(n)


@pytest.fixture(scope="module")
def s2():
    return build(3, 3)


@pytest.mark.parametrize("d, z_ratio, parity", [(0, 0.83, 0.059), (1, 0.83, 0.059), (2, 1.08, 0.116), (3, 1.03, 0.035)])
def test_low_degree_estimates_on_s2(s2, d, z_ratio, parity):
    N, q, ef, sf = s2
    z = estimate(ef, d)
    n = estimate(ef + sf, d)
    assert z / (Z_S2 / q ** N) == pytest.approx(z_ratio, abs=0.006)
    assert n / z == pytest.approx(parity, abs=0.0006)           # truth is -1


@pytest.mark.slow
def test_low_degree_estimate_on_s2_degree4(s2):
    N, q, ef, sf = s2
    z = estimate(ef, 4)
    assert z / (Z_S2 / q ** N) == pytest.approx(1.00, abs=0.006)
    assert estimate(ef + sf, 4) / z == pytest.approx(-0.072, abs=0.0006)


@pytest.mark.slow
def test_exact_enumeration_on_s2(tmp_path):
    if shutil.which("gcc") is None:
        pytest.skip("needs gcc")
    subprocess.run([sys.executable, os.path.join(HERE, "enum_gen.py"), "3", "3"], check=True, cwd=tmp_path,
                   capture_output=True)
    shutil.copy(os.path.join(HERE, "enum.c"), tmp_path)
    subprocess.run(["gcc", "-O2", "-o", "enum", "enum.c"], check=True, cwd=tmp_path)
    out = subprocess.run(["./enum"], check=True, cwd=tmp_path, capture_output=True, text=True).stdout.split("\n")
    assert out[0] == f"Z {Z_S2}"
    hist = {int(l.split()[1]): int(l.split()[2]) for l in out[1:] if l.startswith("A ")}
    assert all(a % 2 == 1 for a in hist) and sum(hist.values()) == Z_S2          # Ky Fan: A is odd on every valid labeling
    assert hist[1] == 11686752
