"""The exact small numbers cited by the paper (docs/size_lower_bound.md; scripts in analysis/).

  * 304  = minimum number of monomials in an F_2 certificate of Tucker on S^2 at degree <= 3.
           Here: the stored certificate (analysis/isd_s2_d3_cert.txt) has 304 violating monomials of size <= 3 and is an
           identity on all 4^13 labelings (slow tier; a 2,000-labeling random check runs in the fast tier).  The matching
           lower bound (CP-SAT optimality) is tests/test_min_size.py.
  * t = 8 on S^2 (Theorem 1'): all 16 gadget domain types at the pole flag, exact odd-hitting number by branch and bound.
  * t >= 12 on S^3 (Theorem 1'): the 32 stored realizing labelings (analysis/flag_hitting_bound_m4_realizations.txt) are
           re-verified without CP-SAT, expanded to 192 gadget types by label orbits, and the parity system is solved
           exactly (CP-SAT, OPTIMAL with proven bound 12).
  * the best degree-<= 4 Z_3-invariant certificate found on S^2 (analysis/isd_s2_d4_z3_cert.txt) has size profile
           {2: 24, 3: 288} -- no size-4 monomial -- and is an identity on all 4^13 labelings (slow; random check fast).
"""
import os
import sys
from collections import Counter

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
import flag_hitting_bound as F
from kyfan import SignedComplex, label_set, violating_pairs
from kyfan.violating import is_violating
from kyfan.primal import verify_exhaustive_F2, verify_random, tucker_target
from test_min_size import load_cert, CERT

HERE = os.path.dirname(__file__)
CERT_D4 = os.path.join(HERE, '..', 'analysis', 'isd_s2_d4_z3_cert.txt')
REALIZATIONS_M4 = os.path.join(HERE, '..', 'analysis', 'flag_hitting_bound_m4_realizations.txt')


def _load_s2_cert(path, size, profile):
    cx = SignedComplex(3)
    labels = label_set(2)
    V = violating_pairs(cx, labels)
    cert = load_cert(path)
    assert len(cert) == size
    assert dict(Counter(len(a) for a in cert)) == profile
    assert all(is_violating(a, V) for a in cert)
    assert all(0 <= i < cx.n_free and l in labels for a in cert for i, l in a)
    return cx, labels, cert


# ---------------------------------------------------------------------------------------------------- 304
def test_min_certificate_304_random():
    cx, labels, cert = _load_s2_cert(CERT, 304, {2: 16, 3: 288})
    rng = np.random.default_rng(0)
    assert verify_random(cx, cert, tucker_target, labels, 2000, rng) == 0


@pytest.mark.slow
def test_min_certificate_304_exhaustive():
    cx, labels, cert = _load_s2_cert(CERT, 304, {2: 16, 3: 288})
    assert verify_exhaustive_F2(cx, cert, labels, 1)


# ---------------------------------------------------------------------------------------------------- t = 8 on S^2
def test_t8_on_s2_branch_and_bound():
    gadgets = F.family_all_s2()
    assert len(gadgets) == 16
    supports = []
    for D in gadgets:
        cons, dim = F.constraints_of(D)
        assert dim == 1
        (S, rhs), = cons
        assert rhs == 1 and len(S) == 3 ** (3 - 2)
        supports.append(S)
    assert F.min_odd_hitting_bb(supports) == 8


# ---------------------------------------------------------------------------------------------------- t >= 12 on S^3
def test_t12_on_s3_from_stored_realizations():
    m, n = 4, 3
    reals = F.load_realizations(REALIZATIONS_M4)
    assert len(reals) == 32
    keys = set()
    for name, perm, target, L in reals:
        assert F.check_realization(m, target, L)          # non-violating off the pole chain, exact residual domains
        assert F.is_gadget(target, n)                     # the residual has a degree-n pseudo-solution
        keys.add(tuple(frozenset(t) for t in target))
    assert len(keys) == 32
    fam = {}
    for key in keys:
        for D2 in F.label_orbit([set(k) for k in key], n):
            fam.setdefault(tuple(frozenset(x) for x in D2), D2)
    gadgets = list(fam.values())
    assert len(gadgets) == 192
    cons = []
    for D in gadgets:
        c, dim = F.constraints_of(D)
        assert dim == 1 and len(c[0][0]) == 3 ** (n - 1) and c[0][1] == 1
        cons += c
    st, t, bound, T = F.min_parity_hitting(cons, seconds=120)
    assert st == "OPTIMAL" and t == 12 and round(bound) == 12
    for S, rhs in cons:                                   # the minimizer really satisfies every parity
        assert len(set(T) & S) % 2 == rhs


# ---------------------------------------------------------------------------------------------------- degree-4 search
def test_degree4_z3_certificate_profile_random():
    cx, labels, cert = _load_s2_cert(CERT_D4, 312, {2: 24, 3: 288})
    rng = np.random.default_rng(1)
    assert verify_random(cx, cert, tucker_target, labels, 2000, rng) == 0


@pytest.mark.slow
def test_degree4_z3_certificate_profile_exhaustive():
    cx, labels, cert = _load_s2_cert(CERT_D4, 312, {2: 24, 3: 288})
    assert verify_exhaustive_F2(cx, cert, labels, 1)
