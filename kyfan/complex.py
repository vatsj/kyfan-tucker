import itertools
from functools import cached_property


def neg(v):
    return tuple(-x for x in v)


def leq(x, y):
    """x <= y iff every coordinate of x is 0 or equals y's."""
    return all(xi == 0 or xi == yi for xi, yi in zip(x, y))


def rank(v):
    return sum(1 for x in v if x)


class SignedComplex:
    """The signed-subset complex on [m] (barycentric subdivision of the cross-polytope boundary), a triangulation of S^{m-1}.

    verts: nonzero vectors in {-1,0,1}^m.  edges: comparable pairs (x,y) with x < y, each once.
    top: maximal chains (rank 1 < 2 < ... < m), m! 2^m of them.
    free: one representative per antipodal pair (first nonzero coordinate +1); rep(v) = (free index, sign).
    """

    def __init__(self, m):
        self.m = m
        self.n = m - 1                      # sphere dimension
        self.verts = [v for v in itertools.product([-1, 0, 1], repeat=m) if any(v)]
        self.free = [v for v in self.verts if next(x for x in v if x) == 1]
        self.fidx = {v: i for i, v in enumerate(self.free)}
        self.n_free = len(self.free)

    neg = staticmethod(neg)
    leq = staticmethod(leq)
    rank = staticmethod(rank)

    def rep(self, v):
        return (self.fidx[v], 1) if v in self.fidx else (self.fidx[neg(v)], -1)

    def label(self, v, L):
        """Label of an arbitrary vertex under the labeling L of the free vertices."""
        i, s = self.rep(v)
        return s * L[i]

    @cached_property
    def by_rank(self):
        br = {}
        for v in self.verts:
            br.setdefault(rank(v), []).append(v)
        return br

    @cached_property
    def edges(self):
        return [(x, y) for x in self.verts for y in self.verts if x != y and leq(x, y)]

    @cached_property
    def top(self):
        chains = [[v] for v in self.by_rank[1]]
        for r in range(2, self.m + 1):
            chains = [c + [w] for c in chains for w in self.by_rank[r] if leq(c[-1], w)]
        return [tuple(c) for c in chains]

    def hemisphere_top(self, j=None):
        """Top simplices of the hemisphere H = {x_j != -1} (default j = m): those whose top vertex has x_j = +1."""
        j = self.m if j is None else j
        return [c for c in self.top if c[-1][j - 1] == 1]

    def euler(self):
        # only meaningful as a sanity check for m=3 (V - E + T); higher m needs all face counts
        return len(self.verts) - len(self.edges) + len(self.top)

    def free_edges(self):
        """Edges as pairs of (free index, sign) with distinct free indices (an edge never joins x and -x)."""
        out = []
        for x, y in self.edges:
            (i, si), (j, sj) = self.rep(x), self.rep(y)
            assert i != j
            out.append(((i, si), (j, sj)))
        return out
