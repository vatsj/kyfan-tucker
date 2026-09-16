"""Task 2 helpers: structure of degree-d pseudo-solutions (orbit types, knockouts, restrictions, descriptions)."""
import itertools
from collections import Counter

from .complex import leq, neg, rank
from .group import act, orbits, full_group_generators
from .dual import unknowns, sa_dual_system, solve


def describe(cx, a):
    """Human-readable full-group-invariant-ish description of a partial labeling a = ((i,l),...).
    Geometry: ranks of the vertices; for each pair whether {±x, ±y} contains a comparable pair ('E' edge / '-' unrelated);
    whether the whole set is a simplex up to signs (a chain after choosing signs). Labels: magnitude multiset and, on
    each edge pair with equal magnitudes, whether the labels are equal on the edge (they must be, non-violating)."""
    vs = [cx.free[i] for i, _ in a]
    ls = [l for _, l in a]
    k = len(a)
    ranks = tuple(rank(v) for v in vs)
    pair = {}
    for p, q in itertools.combinations(range(k), 2):
        x, y = vs[p], vs[q]
        rel = None
        if leq(x, y) or leq(y, x):
            rel = 'E+'        # edge between the reps themselves (same sign)
        elif leq(neg(x), y) or leq(y, neg(x)):
            rel = 'E-'        # edge between rep and antipode of the other
        else:
            rel = '-'
        pair[(p, q)] = rel
    # simplex up to signs: some sign choice makes a chain
    simplex = False
    for signs in itertools.product((1, -1), repeat=k):
        ws = [tuple(s * t for t in v) for s, v in zip(signs, vs)]
        ws_sorted = sorted(ws, key=rank)
        if all(leq(ws_sorted[i], ws_sorted[i + 1]) for i in range(k - 1)):
            simplex = True
            break
    mags = tuple(sorted(abs(l) for l in ls))
    # edge-adjusted sign relation for pairs with equal magnitude
    eq = []
    for (p, q), rel in pair.items():
        if abs(ls[p]) == abs(ls[q]):
            s = 1 if rel != 'E-' else -1
            eq.append(f"{p}{q}:{'eq' if ls[p] == s * ls[q] else 'ne'}")
    edges = ''.join(('E' if pair[(p, q)] != '-' else '-') for p, q in itertools.combinations(range(k), 2))
    return f"ranks={ranks} edges={edges} {'simplex' if simplex else 'nonsimplex'} mags={mags} {' '.join(eq)}"


def full_types(cx, labels, d, V, mags):
    """Full-group orbits of non-violating size-d partial labelings: (type_of: alpha -> type id, members)."""
    return unknowns(cx, labels, d, V, gens=full_group_generators(cx.m, mags))


def knockout(cx, labels, d, V, gens, type_members, base_rows=None, col=None, verbose=False):
    """For each full-group type, solve the (gens-invariant) system with that type's columns removed.
    Returns list of (type id, consistent, rank). Uses the invariant system built once (columns filtered per type)."""
    if col is None or base_rows is None:
        col, _ = unknowns(cx, labels, d, V, gens=gens)
        base_rows, ncols, _ = sa_dual_system(cx, labels, d, V=V, gens=gens, col=col)
    ncols = 1 + max(col.values())
    out = []
    for t, mem in enumerate(type_members):
        removed = {col[a] for a in mem}
        rows = []
        for f, b in base_rows:
            g = {c: x for c, x in f.items() if c not in removed}
            if g or b:
                rows.append((g, b))
        res = solve(rows, ncols, "F2")
        out.append((t, res.consistent, res.rank))
        if verbose:
            print(f"  type {t:4d} size {len(mem):7d} removed {len(removed):6d} orbit-cols -> "
                  f"{'solution' if res.consistent else 'NO solution'} (rank {res.rank})  [{describe(cx, mem[0])}]", flush=True)
    return out
