import itertools


def violating_pairs(cx, labels):
    """Set of frozenset({(i,a),(j,b)}) over free indices i != j: partial labelings of size 2 containing a complementary edge."""
    V = set()
    for (i, si), (j, sj) in cx.free_edges():
        for a in labels:
            for b in labels:
                if si * a == -sj * b:
                    V.add(frozenset([(i, a), (j, b)]))
    return V


def is_violating(alpha, V):
    return any(frozenset(p) in V for p in itertools.combinations(alpha, 2))


def partial_labelings(n, labels, d, vertices=None):
    """All partial labelings of exactly d vertices, as sorted tuples of (vertex, label)."""
    verts = range(n) if vertices is None else vertices
    for S in itertools.combinations(verts, d):
        for ls in itertools.product(labels, repeat=d):
            yield tuple(zip(S, ls))


def restrict(L, S):
    """Restriction of a full labeling L to the vertex tuple S (sorted)."""
    return tuple((i, L[i]) for i in S)
