"""The Sherali-Adams / Nullstellensatz dual at degree d (CLAUDE.md §1, "Pseudo-solution").

Unknowns: E[alpha] for non-violating partial labelings alpha of exactly d vertices (optionally restricted to a
support and/or quotiented by orbits of a subgroup, i.e. E invariant).  E on smaller partial labelings is defined by
extension along the smallest unused vertex; the consistency equations are imposed only at level d-1 (they imply the
lower levels, see review.md), plus E[empty] = 1.

A degree-d certificate over F exists  <=>  this system is inconsistent over F.
"""
import itertools
from collections import defaultdict

from .violating import violating_pairs, is_violating, partial_labelings
from .group import act, orbits
from . import linalg


def unknowns(cx, labels, d, V, gens=None, support=None):
    """Column index of the dual system: (col: alpha -> id, members: list of orbit member lists)."""
    n = cx.n_free
    items = (a for a in partial_labelings(n, labels, d)
             if not is_violating(a, V) and (support is None or support(a)))
    if not gens:
        col, members = {}, []
        for a in items:
            col[a] = len(members)
            members.append([a])
        return col, members
    col, members = orbits(items, gens, lambda g, a: act(cx, g, a))
    if support is not None and any(not support(a) for mem in members for a in mem):
        raise ValueError("support is not invariant under the given group generators")
    return col, members


def sa_dual_system(cx, labels, d, V=None, gens=None, support=None, col=None):
    """Build the degree-d dual. Returns (rows, ncols, col) with rows = list of (dict id -> int coeff, rhs)."""
    n = cx.n_free
    V = violating_pairs(cx, labels) if V is None else V
    if col is None:
        col, _ = unknowns(cx, labels, d, V, gens, support)
    ncols = 1 + max(col.values()) if col else 0

    def term(beta, v, l):
        a = tuple(sorted(beta + ((v, l),)))
        return col.get(a)

    betas = (b for b in partial_labelings(n, labels, d - 1) if not is_violating(b, V))
    if gens:  # rows for beta and g.beta coincide on the orbit quotient: keep one representative per orbit
        _, members = orbits(betas, gens, lambda g, b: act(cx, g, b))
        betas = (mem[0] for mem in members)
    rows = []
    for beta in betas:
        used = {v for v, _ in beta}
        unused = [v for v in range(n) if v not in used]
        v0 = unused[0]
        base = defaultdict(int)
        for l in labels:
            c = term(beta, v0, l)
            if c is not None:
                base[c] += 1
        for w in unused[1:]:
            f = defaultdict(int, base)
            for l in labels:
                c = term(beta, w, l)
                if c is not None:
                    f[c] -= 1
            f = {c: x for c, x in f.items() if x}
            if f:
                rows.append((f, 0))
    f = defaultdict(int)   # E[empty] = 1 through the chain of vertices 0..d-1
    for ls in itertools.product(labels, repeat=d):
        c = col.get(tuple(zip(range(d), ls)))
        if c is not None:
            f[c] += 1
    rows.append((dict(f), 1))
    return rows, ncols, col


def solve(rows, ncols, field="F2", want_solution=False):
    """field: 'F2' or an odd prime p. Returns an object with .consistent, .rank (and .solution over F_2)."""
    if field == "F2":
        return linalg.gf2_dense(rows, ncols, want_solution=want_solution)
    consistent, rank = linalg.fp_sparse(rows, ncols, p=int(field))
    return linalg.GF2Result(consistent, rank, [], ncols, None)


def ns_degree_refutable(cx, labels, d, field="F2", gens=None, support=None):
    """True iff a degree-d certificate over `field` exists (dual inconsistent). Also returns rank and #unknowns."""
    rows, ncols, _ = sa_dual_system(cx, labels, d, gens=gens, support=support)
    res = solve(rows, ncols, field)
    return (not res.consistent), res.rank, ncols
