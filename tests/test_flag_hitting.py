"""Theorem 1' (SIZE_LOWER_BOUND.md): per-flag parity constraints, t(U) = minimum number of full-flag terms.

S^2: all 16 gadget domain types at the pole chain give t = 8 (>= 24 x 8 = 192 size-3 terms in every degree-3
     certificate); the label orbit of the explicit gadget alone gives 4; label orbit x realizable chain permutations
     recovers all 16 types and t = 8.  The minimum certificate (304 terms) has 12 full-flag terms per class and
     satisfies all 16 parities.
S^3: label orbit (48 gadgets) t = 6; all realizable trees with label orbits (96) still 6; realizable position
     permutations with label orbits (192 gadget types) t = 12  =>  >= 192 x 12 = 2304 size-4 terms.
"""
import os, sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'analysis'))
import flag_hitting_bound as F
from test_min_size import load_cert, CERT


def t_of(m, gadgets):
    cons = []
    for D in gadgets:
        c, dim = F.constraints_of(D)
        assert dim == 1 and len(c[0][0]) == 3 ** (m - 2) and c[0][1] == 1
        cons += c
    st, t, bound, T = F.min_parity_hitting(cons, seconds=120)
    assert st == "OPTIMAL"
    return t, cons


def test_s2_all_gadget_types_t8():
    gadgets = F.family_all_s2()
    assert len(gadgets) == 16
    t, cons = t_of(3, gadgets)
    assert t == 8 == F.min_odd_hitting_bb([S for S, _ in cons])
    # the exhaustive family is exactly label orbit x realizable chain permutations
    found, perm_gadgets = F.family_perms(3)
    assert len(found) == 4 and len(perm_gadgets) == 16
    assert {tuple(frozenset(x) for x in D) for D in gadgets} == {tuple(frozenset(x) for x in D) for D in perm_gadgets}
    _, orbit = F.family_orbit(3)
    assert len(orbit) == 8 and t_of(3, orbit)[0] == 4


def test_s2_minimum_certificate_satisfies_parities():
    from kyfan.complex import SignedComplex
    cx = SignedComplex(3)
    U = F.pole_chain(3)
    reps = [cx.rep(y) for y in U]
    cert = load_cert(CERT)
    T = set()
    for a in cert:
        if len(a) == 3 and {i for i, _ in a} == {i for i, _ in reps}:
            L = dict(a)
            T.add(tuple(s * L[i] for i, s in reps))        # vertex-label coordinates
    assert len(T) == 12
    for D in F.family_all_s2():
        (S, rhs), = F.constraints_of(D)[0]
        assert len(T & S) % 2 == rhs == 1


def test_s3_t_from_gadget_families():
    _, orbit = F.family_orbit(4)
    assert len(orbit) == 48 and t_of(4, orbit)[0] == 6
    base, trees = F.family_trees(4)
    assert len(base) == 16 and len(trees) == 96 and t_of(4, trees)[0] == 6
    found, perms = F.family_perms(4)
    assert len(found) == 32 and len(perms) == 192
    assert t_of(4, perms)[0] == 12
