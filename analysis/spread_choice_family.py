"""Numbers for docs/spread_lemma.md: the coordinate-choice gadget family, its exceptional set, and the weighted
flag count T(S) = sum_{F in S} prod_{y in S \\ F} 1/d_F(y).

  * exceptional set (d = 1) and entropy of D_F on the sphere, and of the equatorial sub-family inside G_m
    (the family of analysis/spread_mcmc.py), m = 3..8;
  * max_S T(S): exhaustive for m = 3, hill-climbing lower bounds for m = 4, 5, against the proved bound 2^m;
  * Lemma A and the inequality sum_alpha P[alpha survives] >= 1 on the stored S^2 certificates and on
    high-degree certificates obtained from them by expanding terms.

Usage:  PYTHONPATH=. python analysis/spread_choice_family.py            (~1 min)
"""
import sys, os, itertools, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from fractions import Fraction
from math import log2, factorial
import numpy as np

from kyfan.complex import SignedComplex, neg
from kyfan import spread as sp

HERE = os.path.dirname(__file__)


def load_cert(path):
    return [tuple((int(i), int(l)) for i, l in (tok.split(':') for tok in line.split()))
            for line in open(path) if line.strip() and not line.startswith('#')]


def expand_cert(cert, cx, labels, rng, rounds):
    """A certificate with the same sum: `rounds` times, replace a random term alpha by the 2n terms
    alpha u {v: l} for a random free vertex v outside alpha (1_alpha = sum_l 1_{alpha, v:l}); F_2 cancellation."""
    C = set(cert)
    for _ in range(rounds):
        a = sorted(C)[rng.integers(len(C))]
        used = {i for i, _ in a}
        rest = [i for i in range(cx.n_free) if i not in used]
        if not rest:
            continue
        v = rest[rng.integers(len(rest))]
        C ^= {a}
        C ^= {tuple(sorted(a + ((v, l),))) for l in labels}
    return sorted(C)


def survival_sum(cx, cert, classes=None):
    """(sum over terms alpha and flag classes F of P_{D_F}[alpha|_F full and alpha off F inside rho]) as a Fraction."""
    classes = sp.flag_classes(cx) if classes is None else classes
    choices = [(Fidx, sp.gadget_choices(cx, F)) for F, Fidx in classes]
    tot = Fraction(0)
    for a in cert:
        S = frozenset(i for i, _ in a)
        for Fidx, ch in choices:
            if Fidx <= S:
                p = Fraction(1)
                for i, l in a:
                    if i not in Fidx:
                        p *= Fraction(1 if l in ch[i] else 0, len(ch[i]))
                tot += p
    return tot


def lemma_A_holds(cx, cert, classes=None):
    """Every gadget of every D_F (m = 3: 16 per flag class) is served by some term of the certificate."""
    classes = sp.flag_classes(cx) if classes is None else classes
    sup = [(frozenset(i for i, _ in a), a) for a in cert]
    n = 0
    for F, Fidx in classes:
        cands = [a for S, a in sup if Fidx <= S]
        for rho in sp.all_gadgets(cx, F):
            n += 1
            if not any(all(i in Fidx or rho[i] == l for i, l in a) for a in cands):
                return False, n
    return True, n


class TSearch:
    """T(S) with a precomputed table d[class, vertex]; hill-climbing over S."""

    def __init__(self, m):
        self.cx = cx = SignedComplex(m)
        self.classes = sp.flag_classes(cx)
        self.member = np.zeros((len(self.classes), cx.n_free), dtype=bool)
        self.logd = np.zeros((len(self.classes), cx.n_free))
        for k, (F, Fidx) in enumerate(self.classes):
            for i, v in enumerate(cx.free):
                if i in Fidx:
                    self.member[k, i] = True
                else:
                    self.logd[k, i] = np.log(sp.d_flag(F, v))

    def T(self, s):                                  # s: boolean vector over free vertices
        inside = ~(self.member & ~s).any(axis=1)
        return float(np.exp(-(self.logd[inside][:, s]).sum(axis=1)).sum())

    def climb(self, rng, iters, start=None):
        N = self.cx.n_free
        s = np.zeros(N, dtype=bool) if start is None else start.copy()
        if start is None:
            s[list(self.classes[rng.integers(len(self.classes))][1])] = True
        cur = self.T(s)
        for _ in range(iters):
            t = s.copy()
            for i in rng.choice(N, size=rng.integers(1, 3), replace=False):
                t[i] = not t[i]
            val = self.T(t)
            if val >= cur:
                s, cur = t, val
        return cur, s


def equatorial_family(m):
    """The cap-copying sub-family inside G_m: per movable equatorial vertex of GadgetSampler(m), the labels
    u_c (n-c+1), c in supp(u) minus the initial constant run.  Returns (#movable, #with one choice, all static_ok)."""
    from analysis.spread_mcmc import GadgetSampler
    g = GadgetSampler(m)
    n = m - 1
    pos = {int(l): k for k, l in enumerate(g.labels)}
    single, ok, bits = 0, True, 0.0
    for i in g.movable:
        u = g.cxe.free[i]
        a = 0
        while a < n and u[a] == u[0] != 0:
            a += 1
        labs = [u[c - 1] * (n - c + 1) for c in range(a + 1, n + 1) if u[c - 1]]
        ok &= len(labs) >= 1 and all(g.static_ok[i, pos[l]] for l in labs) and int(g.L0[i]) in labs
        single += len(labs) == 1
        bits += log2(len(labs))
    return len(g.movable), single, ok, bits


if __name__ == "__main__":
    t0 = time.time()
    print("== exceptional set and entropy of D_F (pole flag; classes of vertices off the flag)")
    for m in range(3, 9):
        U = sp.pole_flag(m)
        Uset = set(U) | {neg(x) for x in U}
        ds = [sp.d_flag(U, y) for y in itertools.product([-1, 0, 1], repeat=m)
              if any(y) and y not in Uset and next(x for x in y if x) == 1]
        mov, single, ok, ebits = equatorial_family(m)
        print(f"m={m}: sphere: {len(ds)} classes off the flag, {ds.count(1)} with d=1 (m(m-1)={m * (m - 1)}), "
              f"entropy {sum(log2(d) for d in ds):.1f} bits;  equator (inside G_m: {ok}): {mov} movable, "
              f"{single} with one choice (n(n-1)={(m - 1) * (m - 2)}), entropy {ebits:.1f} bits")
    print("== max_S T(S)   (proved: T(S) <= 2^m for every S)")
    cx = SignedComplex(3)
    cls = sp.flag_classes(cx)
    best = max((sp.T_weight(cx, [i for i in range(13) if b >> i & 1], cls), b) for b in range(1, 2 ** 13))
    print(f"m=3: exhaustive over 2^13 sets: max T = {best[0]}  (bound 8), attained at "
          f"{[cx.free[i] for i in range(13) if best[1] >> i & 1]}")
    for m, restarts, iters in ((4, 12, 3000), (5, 6, 4000)):
        ts = TSearch(m)
        rng = np.random.default_rng(m)
        vals = sorted((ts.climb(rng, iters) for _ in range(restarts)), key=lambda z: -z[0])
        v, s = vals[0]
        print(f"m={m}: hill-climbing ({restarts} restarts x {iters} moves): best T found = {v:.4f} with |S| = {int(s.sum())} "
              f"(bound {2 ** m}); runs: {[round(x, 3) for x, _ in vals]}   [{time.time() - t0:.0f}s]")
    print("== S^2 certificates: Lemma A for every gadget of every D_F, and (1/#classes) sum_alpha sum_F P_F[alpha]  (must be >= 1)")
    labels = [1, -1, 2, -2]
    for name in ("isd_s2_d3_cert.txt", "isd_s2_d4_z3_cert.txt"):
        cert = load_cert(os.path.join(HERE, name))
        okA, ng = lemma_A_holds(cx, cert, cls)
        print(f"{name}: {len(cert)} terms, max degree {max(map(len, cert))}: Lemma A on {ng} gadgets: {okA}; "
              f"normalized survival sum = {survival_sum(cx, cert, cls) / len(cls)}")
        rng = np.random.default_rng(5)
        for rounds in (50, 400):
            big = expand_cert(cert, cx, labels, rng, rounds)
            okA, ng = lemma_A_holds(cx, big, cls)
            ss = survival_sum(cx, big, cls) / len(cls)
            print(f"   expanded x{rounds}: {len(big)} terms, max degree {max(map(len, big))}: Lemma A: {okA}; "
                  f"normalized survival sum = {ss} = {float(ss):.3f}")
    print("== proved size bound m!/2 (all degrees) versus the number of flag classes")
    for m in range(3, 13):
        print(f"m={m}: N = {3 ** m - 1}, size >= {factorial(m) // 2}   (#flag classes m!2^(m-1) = {factorial(m) * 2 ** (m - 1)})")
