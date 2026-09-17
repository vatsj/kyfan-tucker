"""Minimum certificate size on S^2 (SIZE_LOWER_BOUND.md; `analysis/isd.py`).

Degree <= 3: the minimum number of monomials in an F_2 certificate of Tucker on S^2 is exactly 304.
  <= : the certificate stored in analysis/isd_s2_d3_cert.txt (16 size-2 + 288 size-3 monomials) is verified exhaustively.
  >= : CP-SAT proves optimality on the sampled system (true certificates are a subset of its solutions, so the
       sampled minimum is a lower bound; the verified certificate attains it).
"""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
import isd
from kyfan.primal import verify_exhaustive_F2

CERT = os.path.join(os.path.dirname(__file__), '..', 'analysis', 'isd_s2_d3_cert.txt')


def load_cert(path):
    cert = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        cert[tuple((int(i), int(l)) for i, l in (tok.split(':') for tok in line.split()))] = 1
    return cert


@pytest.mark.slow
def test_s2_degree3_min_size_304():
    cert = load_cert(CERT)
    assert len(cert) == 304 and all(len(a) <= 3 for a in cert)
    M, members, labels, cx = isd.build(3, 2, 3, 9500, seed=0)
    assert verify_exhaustive_F2(cx, cert, labels, 1)
    C = len(members)
    rank, piv = isd.reduce(M, C)
    assert rank == 4176
    col = {a: i for i, mem in enumerate(members) for a in mem}
    hint = np.zeros(C, dtype=np.uint8)
    hint[[col[a] for a in cert]] = 1
    assert isd.check_solution(M, hint, C)      # row operations preserve the solution set
    status, w, x, bound = isd.cpsat_min_weight(M, C, rank, piv, seconds=600, workers=8, hint=hint, log=False)
    assert status == "OPTIMAL" and w == 304 and round(bound) == 304
