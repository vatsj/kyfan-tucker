"""Explicit certificates via the random-point primal (docs/AGENT_BRIEF.md §2, "one-way sound"):
rows = evaluation of every violating size-d indicator on a random labeling, rhs = target(L).
An inconsistent sampled system proves no certificate; a consistent one yields a candidate that must be verified —
on fresh labelings (probabilistic) or exactly (`verify_exact_F2`, canonical one-hot basis; `verify_exhaustive`).
"""
import itertools
from collections import defaultdict
import numpy as np

from .violating import violating_pairs, is_violating, partial_labelings, restrict
from .labels import posalt, pos_alt_count, random_labeling
from .group import act, orbits
from . import linalg


def violating_columns(cx, labels, d, V=None, gens=None):
    """(col: alpha -> id, members) over violating size-d partial labelings, optionally quotiented by orbits."""
    V = violating_pairs(cx, labels) if V is None else V
    items = (a for a in partial_labelings(cx.n_free, labels, d) if is_violating(a, V))
    if not gens:
        col, members = {}, []
        for a in items:
            col[a] = len(members)
            members.append([a])
        return col, members
    return orbits(items, gens, lambda g, a: act(cx, g, a))


def _row_counts(L, d, col, ncols, n):
    r = np.zeros(ncols, dtype=np.int64)
    for S in itertools.combinations(range(n), d):
        c = col.get(restrict(L, S))
        if c is not None:
            r[c] += 1
    return r


def _row_bitset(L, d, col, n):
    r = 0
    for S in itertools.combinations(range(n), d):
        c = col.get(restrict(L, S))
        if c is not None:
            r ^= 1 << c
    return r


class FastRows:
    """Vectorized row evaluation (the z9a.py trick): partial labelings are encoded as base-(n*q) integers,
    looked up by searchsorted, and counted per column with bincount."""

    def __init__(self, cx, labels, d, col):
        self.n, self.q, self.d = cx.n_free, len(labels), d
        self.labels = list(labels)
        self.li = {l: k for k, l in enumerate(labels)}
        self.B = self.n * self.q
        if self.B ** d >= 2 ** 62:
            raise ValueError("key encoding overflows int64")
        self.ncols = 1 + max(col.values()) if col else 0
        items = sorted((self.key(a), c) for a, c in col.items())
        self.keys = np.array([k for k, _ in items], dtype=np.int64)
        self.oid = np.array([c for _, c in items], dtype=np.int64)
        self.tuples = np.array(list(itertools.combinations(range(self.n), d)), dtype=np.int64)

    def key(self, a):
        return sum((i * self.q + self.li[l]) * self.B ** k for k, (i, l) in enumerate(a))

    def label_indices(self, L):
        return np.array([self.li[l] for l in L], dtype=np.int64)

    def counts(self, Lidx):
        """Column counts for the labeling given by label indices."""
        k = np.zeros(len(self.tuples), dtype=np.int64)
        for pos in range(self.d):
            t = self.tuples[:, pos]
            k += (t * self.q + Lidx[t]) * self.B ** pos
        p = np.searchsorted(self.keys, k)
        p[p >= len(self.keys)] = 0
        hit = self.keys[p] == k
        return np.bincount(self.oid[p[hit]], minlength=self.ncols)

    def packed_F2_system(self, target, nrows, rng):
        """Sampled primal over F_2 as a packed uint64 matrix (rows = random labelings, last bit = target)."""
        bit_rows = []
        for _ in range(nrows):
            Lidx = rng.integers(0, self.q, self.n)
            L = [self.labels[k] for k in Lidx]
            bits = self.counts(Lidx) & 1
            bit_rows.append(np.append(bits, target(L) & 1))
        return linalg.pack_bit_rows(bit_rows, self.ncols)


def tucker_target(L):
    return 1


def kyfan_target(cx):
    """(A_+ + 1) mod 2 as a function of the labeling — the Ky Fan target over F_2."""
    return lambda L: (pos_alt_count(cx, L) + 1) & 1


def sampled_system(cx, labels, d, target, nsamp, rng, gens=None, V=None):
    """Integer sampled primal: A[r, o] = #violating size-d restrictions of labeling r in orbit o, b[r] = target."""
    n = cx.n_free
    col, members = violating_columns(cx, labels, d, V, gens)
    Ls = [random_labeling(cx, labels, rng) for _ in range(nsamp)]
    A = np.array([_row_counts(L, d, col, len(members), n) for L in Ls], dtype=np.int64)
    b = np.array([target(L) for L in Ls], dtype=np.int64)
    return A, b, col, members


def sampled_primal(cx, labels, d, target, nsamp, rng, field="F2", gens=None, V=None):
    """Try to find a degree-d certificate  target = sum_alpha c_alpha 1_alpha  from nsamp random labelings.
    Returns (cert: dict alpha -> coeff, ncols) or (None, ncols). field: 'F2' or an odd prime."""
    n = cx.n_free
    if field == "F2" and not gens:
        col, members = violating_columns(cx, labels, d, V, gens)
        ncols = len(members)
        Ls = [random_labeling(cx, labels, rng) for _ in range(nsamp)]
        x = linalg.gf2_bitset_solve([_row_bitset(L, d, col, n) for L in Ls], [target(L) for L in Ls])
        if x is None:
            return None, ncols
        coeffs = {c: 1 for c in range(ncols) if (x >> c) & 1}
    else:
        A, b, col, members = sampled_system(cx, labels, d, target, nsamp, rng, gens, V)
        ncols = len(members)
        if field == "F2":
            x = linalg.gf2_bitset_solve([int(''.join('1' if v & 1 else '0' for v in r[::-1]), 2) for r in A], list(b))
            if x is None:
                return None, ncols
            coeffs = {c: 1 for c in range(ncols) if (x >> c) & 1}
        else:
            x = linalg.fp_dense_solve(A, b, p=int(field))
            if x is None:
                return None, ncols
            coeffs = {c: int(v) for c, v in enumerate(x) if v}
    cert = {a: coeffs[c] for c, mem in enumerate(members) if c in coeffs for a in mem}
    return cert, ncols


def evaluate(cert, L, mod=2):
    """sum_alpha c_alpha 1_alpha(L) mod `mod`."""
    return sum(c for a, c in cert.items() if all(L[i] == l for i, l in a)) % mod


def verify_random(cx, cert, target, labels, nverify, rng, mod=2):
    """Number of failures of the identity on nverify fresh random labelings."""
    return sum(1 for _ in range(nverify)
               if evaluate(cert, (L := random_labeling(cx, labels, rng)), mod) != target(L) % mod)


def verify_exhaustive_F2(cx, cert, labels, target_value=1, limit=2 ** 27):
    """Check sum c_alpha 1_alpha == target_value (a constant) on ALL labelings, vectorized over base-q encodings.
    Only for q^n_free <= limit (m=3, k=2: 4^13 = 67M, ~20 s)."""
    n, q = cx.n_free, len(labels)
    N = q ** n
    if N > limit:
        raise ValueError(f"{N} labelings exceeds limit {limit}")
    li = {l: k for k, l in enumerate(labels)}
    idx = np.arange(N, dtype=np.int64)
    digits = []
    for i in range(n):
        digits.append((idx % q).astype(np.uint8))
        idx //= q
    tot = np.zeros(N, dtype=np.uint8)
    for a, c in cert.items():
        if c % 2 == 0:
            continue
        m = np.ones(N, dtype=bool)
        for i, l in a:
            m &= digits[i] == li[l]
        tot ^= m
    return bool((tot == target_value % 2).all())


def canonical_expand(alpha, labels, memo=None):
    """Expand 1_alpha in the canonical basis {1_beta : beta avoids labels[-1]} over F_2 using
    1_{beta+(v,last)} = 1_beta - sum_{l != last} 1_{beta+(v,l)}. Returns dict beta -> 1."""
    memo = {} if memo is None else memo
    last, others = labels[-1], labels[:-1]

    def rec(a):
        a = tuple(sorted(a))
        if a in memo:
            return memo[a]
        k = next((k for k, (v, l) in enumerate(a) if l == last), None)
        if k is None:
            out = {a: 1}
        else:
            v = a[k][0]
            rest = a[:k] + a[k + 1:]
            acc = defaultdict(int)
            for c, x in rec(rest).items():
                acc[c] ^= x
            for l in others:
                for c, x in rec(rest + ((v, l),)).items():
                    acc[c] ^= x
            out = {c: x for c, x in acc.items() if x}
        memo[a] = out
        return out

    return rec(alpha)


def kyfan_target_monomials(cx, labels):
    """A_+ + 1 as an F_2 combination of indicators of partial labelings (size n+1 on each top simplex, plus the empty one)."""
    out = defaultdict(int)
    out[()] ^= 1
    for c in cx.top:
        reps = [cx.rep(v) for v in c]
        for ls in itertools.product(labels, repeat=len(c)):
            if posalt([s * l for (i, s), l in zip(reps, ls)]):
                out[tuple(sorted((i, l) for (i, s), l in zip(reps, ls)))] ^= 1
    return {a: x for a, x in out.items() if x}


def verify_exact_F2(cert, target_monomials, labels):
    """Exact identity check over F_2: sum c_alpha 1_alpha == sum t_beta 1_beta as functions on all labelings,
    by comparing canonical-basis expansions (the canonical indicators are linearly independent)."""
    memo = {}
    lhs, rhs = defaultdict(int), defaultdict(int)
    for a, c in cert.items():
        if c % 2:
            for b, x in canonical_expand(a, labels, memo).items():
                lhs[b] ^= x
    for a, c in target_monomials.items():
        if c % 2:
            for b, x in canonical_expand(a, labels, memo).items():
                rhs[b] ^= x
    return {b for b, x in lhs.items() if x} == {b for b, x in rhs.items() if x}
