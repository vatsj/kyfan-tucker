"""docs/spread_lemma.md: the coordinate-choice gadget family, the weighted spread bound, the weighted flag count
T(S) <= 2^m, and the all-degrees size bound m!/2.  Fast tier (a few seconds in total)."""
import itertools
import os
import sys
from fractions import Fraction
from math import factorial, prod

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kyfan.complex import SignedComplex, leq, neg, rank
from kyfan import spread as sp
from kyfan import linalg
from kyfan.abstract_gadget import build_dual, unsat
from kyfan.labels import label_set
from kyfan.realize_explicit import reverse_caterpillar_domains
from analysis import spread_choice_family as scf

ANALYSIS = os.path.join(os.path.dirname(__file__), '..', 'analysis')


def off_flag_vertices(m, F):
    Fset = set(F) | {neg(x) for x in F}
    return [y for y in itertools.product([-1, 0, 1], repeat=m) if any(y) and y not in Fset]


# ---------------------------------------------------------------- Lemma 1: every choice labeling is valid
@pytest.mark.parametrize("m", [3, 4, 5])
def test_lemma1_any_choice_labeling_is_nonviolating(m):
    cx = SignedComplex(m)
    rng = np.random.default_rng(m)
    for _ in range(5):
        rho = {}
        for i, v in enumerate(cx.free):
            supp = [c for c in range(1, m) if v[c - 1]]
            if supp:                                              # every vertex except e_m; no gadget condition imposed
                rho[i] = sp.choice_label(v, supp[rng.integers(len(supp))])
        assert len(rho) == cx.n_free - 1
        assert sp.is_nonviolating(cx, rho)


# ---------------------------------------------------------------- Lemma 2: the gadget criterion is exact, d = |choices|
def test_tree_domains_are_the_reverse_caterpillar():
    for n in range(2, 9):
        assert sp.tree_domains(n) == [reverse_caterpillar_domains(n)[r] for r in range(n + 1)]


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_lemma2_admissible_choices_exact(m):
    """A choice c at y keeps T+[r] inside every residual domain iff c is in supp(y) & [n] minus [a(y)];
    the number of admissible choices is d_F(y) >= 1."""
    n = m - 1
    U = sp.pole_flag(m)
    assert sp.is_flag(U)
    T = sp.tree_domains(n)
    for y in off_flag_vertices(m, U):
        al = sp.allowed_choices(y)
        assert len(al) == sp.d_flag(U, y) >= 1
        assert sp.allowed_choices(neg(y)) == al                                   # antipodal consistency
        for c in range(1, m):
            if not y[c - 1]:
                continue
            lab = sp.choice_label(y, c)
            assert sp.choice_label(neg(y), c) == -lab
            hits_tree_domain = (any((leq(U[r], y) or leq(y, U[r])) and -lab in T[r] for r in range(m)) or
                                any((leq(neg(U[r]), y) or leq(y, neg(U[r]))) and lab in T[r] for r in range(m)))
            assert hits_tree_domain == (c not in al)


@pytest.mark.parametrize("m", [3, 4, 5])
def test_lemma2_lemma3_sampled_gadgets_on_the_sphere(m):
    """From the sphere's own edge list: members of D_F (random flags F) are non-violating, label everything off F,
    and their residual CSP on F is unsatisfiable with a consistent degree-n dual (so every certificate has a term
    that is full on F: Lemma 3).  For the pole flag the domains contain the reverse-caterpillar domains."""
    cx = SignedComplex(m)
    n = m - 1
    rng = np.random.default_rng(10 + m)
    classes = sp.flag_classes(cx)
    assert len(classes) == factorial(m) * 2 ** (m - 1)
    flags = [sp.pole_flag(m)] + [classes[rng.integers(len(classes))][0] for _ in range(4)]
    for k, F in enumerate(flags):
        rho = sp.sample_gadget(cx, F, rng)
        assert len(rho) == cx.n_free - m
        assert sp.is_nonviolating(cx, rho)
        doms = sp.residual_domains(cx, rho, F)
        if k == 0:
            assert all(t <= d for t, d in zip(sp.tree_domains(n), doms))
        assert unsat(doms)
        rows, col = build_dual(doms, label_set(n), n)
        assert linalg.gf2_dense(rows, len(col)).consistent


def test_lemma3_larger_domains_keep_the_pseudo_solution():
    """Enlarging the tree domains arbitrarily never makes the degree-n dual inconsistent."""
    rng = np.random.default_rng(3)
    for n in (2, 3, 4):
        labels = label_set(n)
        for _ in range(6):
            doms = [set(t) | {l for l in labels if rng.random() < 0.4} for t in sp.tree_domains(n)]
            rows, col = build_dual(doms, labels, n)
            assert linalg.gf2_dense(rows, len(col)).consistent


# ---------------------------------------------------------------- Lemma 4: transport, d_F intrinsic
@pytest.mark.parametrize("m", [3, 4, 5])
def test_lemma4_d_is_intrinsic_and_invariant(m):
    cx = SignedComplex(m)
    rng = np.random.default_rng(20 + m)
    classes = sp.flag_classes(cx)
    U = sp.pole_flag(m)
    for k in rng.choice(len(classes), size=6, replace=False):
        F, Fidx = classes[k]
        pull = sp.to_pole_coordinates(F)
        assert tuple(pull(x) for x in F) == U
        ch = sp.gadget_choices(cx, F)
        assert set(ch) == set(range(cx.n_free)) - Fidx
        for i, labs in ch.items():
            v = cx.free[i]
            assert len(set(labs)) == len(labs) == sp.d_flag(F, v) == sp.d_flag(U, pull(v)) == sp.d_flag(F, neg(v))


@pytest.mark.parametrize("m,expected", [(3, 6), (4, 12), (5, 20), (6, 30), (7, 42)])
def test_exceptional_set_is_m_times_m_minus_1(m, expected):
    U = sp.pole_flag(m)
    ds = [sp.d_flag(U, y) for y in off_flag_vertices(m, U) if next(x for x in y if x) == 1]
    assert len(ds) == (3 ** m - 1) // 2 - m
    assert ds.count(1) == expected == m * (m - 1)


@pytest.mark.parametrize("m", [4, 5, 6])
def test_equatorial_subfamily_lies_in_G_m(m):
    """Cap-copying version: inside the family G_m of analysis/spread_mcmc.py, with exactly n(n-1) one-choice vertices."""
    movable, single, ok, _ = scf.equatorial_family(m)
    n = m - 1
    assert ok and single == n * (n - 1) and movable == (3 ** n - 1) // 2 - n


# ---------------------------------------------------------------- Lemma 5: the weighted spread bound is exact
def test_lemma5_product_distribution_m3():
    cx = SignedComplex(3)
    for F, Fidx in sp.flag_classes(cx):
        gad = list(sp.all_gadgets(cx, F))
        ch = sp.gadget_choices(cx, F)
        assert len(gad) == prod(len(v) for v in ch.values()) == 16
        A = sorted(ch)[::2]
        target = {i: gad[5][i] for i in A}
        hits = sum(all(g[i] == l for i, l in target.items()) for g in gad)
        assert Fraction(hits, len(gad)) == prod(Fraction(1, sp.d_flag(F, cx.free[i])) for i in A)


# ---------------------------------------------------------------- Lemma 6 and Theorem B: T(S) <= 2^m
def random_rich_set(cx, classes, rng):
    """A random union of a few flags plus a few random vertices (so that T(S) is not trivially 0)."""
    S = set()
    for k in rng.choice(len(classes), size=rng.integers(1, 6)):
        S |= classes[k][1]
    S |= {int(i) for i in rng.choice(cx.n_free, size=rng.integers(0, 6))}
    return frozenset(S)


@pytest.mark.parametrize("m", [3, 4, 5])
def test_lemma6_sibling_and_top_profile(m):
    cx = SignedComplex(m)
    classes = sp.flag_classes(cx)
    rng = np.random.default_rng(30 + m)
    tops = [v for v in cx.free if rank(v) == m]
    for _ in range(15):
        S = random_rich_set(cx, classes, rng)
        for F, Fidx in classes:
            if not Fidx <= S:
                continue
            prof = sp.descent_profile(cx, F, S)
            bound = Fraction(1)
            for r, ds in enumerate(prof, start=1):
                assert len(set(ds)) == len(ds) and all(1 <= d <= r for d in ds)     # distinct, in [1, r]
                bound /= factorial(len(ds))
            other = sorted(sp.d_flag(F, x) for x in tops if cx.fidx[x] in S and cx.fidx[x] not in Fidx)
            for k, d in enumerate(other, start=1):                                  # k-th smallest >= ceil(log2(k+1))
                assert 2 ** d >= k + 1
                bound /= d
            assert sp.flag_weight(cx, F, Fidx, S) <= bound
    U = sp.pole_flag(m)                                                              # top classes with d <= D: 2^D - 1
    dtop = [sp.d_flag(U, x) for x in tops if x != U[-1] and neg(x) != U[-1]]
    assert all(sum(d <= D for d in dtop) == 2 ** D - 1 for D in range(1, m))


def test_theoremB_exhaustive_m3():
    cx = SignedComplex(3)
    classes = sp.flag_classes(cx)
    best = max(sp.T_weight(cx, [i for i in range(13) if b >> i & 1], classes) for b in range(1, 2 ** 13))
    assert best == 4 <= sp.T_bound(3)


@pytest.mark.parametrize("m", [4, 5])
def test_theoremB_random_and_adversarial_sets(m):
    ts = scf.TSearch(m)
    cx, classes = ts.cx, ts.classes
    rng = np.random.default_rng(40 + m)
    for _ in range(40):
        S = random_rich_set(cx, classes, rng)
        t = sp.T_weight(cx, S, classes)
        assert t <= sp.T_bound(m)
        s = np.zeros(cx.n_free, dtype=bool)
        s[list(S)] = True
        assert abs(ts.T(s) - float(t)) < 1e-9
    val, s = ts.climb(rng, 300)
    assert 2 <= val <= sp.T_bound(m)
    assert sp.T_weight(cx, np.flatnonzero(s).tolist(), classes) <= sp.T_bound(m)
    full = np.ones(cx.n_free, dtype=bool)
    assert ts.T(full) < 1


@pytest.mark.parametrize("m", [3, 4, 5, 6])
def test_block_example_T_at_least_2_to_half_m(m):
    """S = pole flag plus the vertices created by swapping the coordinates at positions (1,2), (3,4), ...:
    2^{floor(m/2)} flags, every vertex of S has d_F <= 1 for each of them, so T(S) >= 2^{floor(m/2)}."""
    cx = SignedComplex(m)
    U = sp.pole_flag(m)
    order = [m - 1] + list(range(m - 1))                    # coordinate added at rank 1, 2, ..., m
    top = U[-1]
    flags = []
    for mask in itertools.product([0, 1], repeat=m // 2):
        o = list(order)
        for b, sw in enumerate(mask):
            if sw:
                o[2 * b], o[2 * b + 1] = o[2 * b + 1], o[2 * b]
        flags.append(tuple(tuple(top[k] if k in o[:r] else 0 for k in range(m)) for r in range(1, m + 1)))
    S = frozenset(cx.rep(x)[0] for F in flags for x in F)
    assert len(S) == m + m // 2
    for F in flags:
        assert sp.is_flag(F) and all(sp.d_flag(F, cx.free[i]) <= 1 for i in S)
    if m <= 5:
        assert sp.T_weight(cx, S) >= 2 ** (m // 2)


# ---------------------------------------------------------------- Theorem C on real certificates (S^2)
@pytest.mark.parametrize("name", ["isd_s2_d3_cert.txt", "isd_s2_d4_z3_cert.txt"])
def test_theoremC_on_stored_and_expanded_certificates(name):
    """Lemma A for all 24 x 16 gadgets, and sum_alpha sum_F P_F[alpha survives] >= #flag classes, on the stored S^2
    certificates and on certificates of degree up to ~8 obtained by expanding terms (1_alpha = sum_l 1_{alpha,v:l})."""
    cx = SignedComplex(3)
    classes = sp.flag_classes(cx)
    cert = scf.load_cert(os.path.join(ANALYSIS, name))
    rng = np.random.default_rng(7)
    for rounds in (0, 60):
        C = scf.expand_cert(cert, cx, [1, -1, 2, -2], rng, rounds)
        if rounds:
            assert max(map(len, C)) >= 5
        ok, ngad = scf.lemma_A_holds(cx, C, classes)
        assert ok and ngad == 24 * 16
        total = scf.survival_sum(cx, C, classes)
        assert total >= len(classes)
        # the chain of inequalities in the proof of Theorem C
        assert total <= sum(sp.T_weight(cx, {i for i, _ in a}, classes) for a in C) <= len(C) * sp.T_bound(3)
        assert len(C) >= sp.size_bound(3)


def test_expanded_certificate_is_a_certificate():
    """The expansion used above preserves the F_2 identity (checked on 20,000 random labelings)."""
    cx = SignedComplex(3)
    cert = scf.load_cert(os.path.join(ANALYSIS, "isd_s2_d3_cert.txt"))
    rng = np.random.default_rng(8)
    C = scf.expand_cert(cert, cx, [1, -1, 2, -2], rng, 60)
    labs = np.array([1, -1, 2, -2])
    L = labs[rng.integers(4, size=(20000, cx.n_free))]
    tot = np.zeros(20000, dtype=np.uint8)
    for a in C:
        hit = np.ones(20000, dtype=bool)
        for i, l in a:
            hit &= L[:, i] == l
        tot ^= hit
    assert tot.all()


def test_size_bound_values():
    assert [sp.size_bound(m) for m in (3, 4, 5, 8)] == [3, 12, 60, 20160]
    for m in range(3, 12):                       # #classes / 2^m
        assert sp.size_bound(m) == Fraction(factorial(m) * 2 ** (m - 1), sp.T_bound(m))
