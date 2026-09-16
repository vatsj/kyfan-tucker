"""Task 3 — Tucker's ball form (CLAUDE.md §4).

Fix a valid labeling L_eq of the equator (the [m-1]-complex, labels ±1..±(m-1)). Variables: the cap C = {x_m = +1}
(3^{m-1} vertices; the opposite cap is determined by antipodality). Constraints: no complementary edge inside C, and for
every cap vertex y and equatorial w <= pi(y) (pi zeroes the last coordinate): lambda(y) != -L_eq(w). The latter are unary
and shrink the domain D(y). UNSAT by Tucker; any degree-d sphere certificate restricts to a degree-<= d certificate here.

The SA-dual is the same as `dual.sa_dual_system` but with per-vertex domains.
"""
import itertools
from collections import defaultdict

from .complex import SignedComplex, leq
from .labels import label_set
from . import linalg


class Ball:
    def __init__(self, m, Leq):
        self.m = m
        self.cxe = SignedComplex(m - 1)
        self.Leq = list(Leq)
        self.labels = label_set(m - 1)
        self.cap = [v + (1,) for v in itertools.product([-1, 0, 1], repeat=m - 1)]
        self.n = len(self.cap)
        self.idx = {v: i for i, v in enumerate(self.cap)}
        # domains
        self.domains = []
        for y in self.cap:
            p = y[:-1]
            forbidden = {-self.cxe.label(w, self.Leq) for w in self.cxe.verts if leq(w, p)}
            self.domains.append([l for l in self.labels if l not in forbidden])
        # binary constraints on cap edges
        self.edges = [(i, j) for i, x in enumerate(self.cap) for j, y in enumerate(self.cap) if i != j and leq(x, y)]
        self.V = set()
        for i, j in self.edges:
            for a in self.domains[i]:
                for b in self.domains[j]:
                    if a == -b:
                        self.V.add(frozenset([(i, a), (j, b)]))

    def e_m(self):
        return self.idx[(0,) * (self.m - 1) + (1,)]

    def is_violating(self, alpha):
        return any(frozenset(p) in self.V for p in itertools.combinations(alpha, 2))

    def partial_labelings(self, d):
        for S in itertools.combinations(range(self.n), d):
            for ls in itertools.product(*(self.domains[i] for i in S)):
                yield tuple(zip(S, ls))

    def sa_dual_system(self, d, support=None):
        col = {}
        for a in self.partial_labelings(d):
            if not self.is_violating(a) and (support is None or support(a)):
                col[a] = len(col)
        ncols = len(col)

        def term(beta, v, l):
            return col.get(tuple(sorted(beta + ((v, l),))))

        rows = []
        for beta in self.partial_labelings(d - 1):
            if self.is_violating(beta):
                continue
            used = {v for v, _ in beta}
            unused = [v for v in range(self.n) if v not in used]
            v0 = unused[0]
            base = defaultdict(int)
            for l in self.domains[v0]:
                c = term(beta, v0, l)
                if c is not None:
                    base[c] += 1
            for w in unused[1:]:
                f = defaultdict(int, base)
                for l in self.domains[w]:
                    c = term(beta, w, l)
                    if c is not None:
                        f[c] -= 1
                f = {c: x for c, x in f.items() if x}
                if f:
                    rows.append((f, 0))
        f = defaultdict(int)
        for ls in itertools.product(*(self.domains[i] for i in range(d))):
            c = col.get(tuple(zip(range(d), ls)))
            if c is not None:
                f[c] += 1
        rows.append((dict(f), 1))
        return rows, ncols, col

    def degree(self, dmax, field="F2", verbose=False):
        """Smallest d <= dmax with a degree-d certificate (dual inconsistent), or None. Returns (d or None, log)."""
        log = []
        for d in range(1, dmax + 1):
            rows, ncols, _ = self.sa_dual_system(d)
            if field == "F2":
                res = linalg.gf2_sparse(rows, ncols)
                consistent, rank = res.consistent, res.rank
            else:
                consistent, rank = linalg.fp_sparse(rows, ncols, p=int(field))
            log.append((d, ncols, len(rows), rank, consistent))
            if verbose:
                print(f"   d={d}: {ncols} unknowns, {len(rows)} rows, rank {rank}, {'consistent' if consistent else 'INCONSISTENT -> certificate'}", flush=True)
            if not consistent:
                return d, log
        return None, log


def valid_equatorial_labelings(m, cap=None, seed=0):
    """Valid labelings of the equator S^{m-2} with labels ±1..±(m-1) (DFS, at most `cap`)."""
    from .labels import dfs_valid_labelings
    cxe = SignedComplex(m - 1)
    return dfs_valid_labelings(cxe, label_set(m - 1), cap=cap or 10 ** 9, seed=seed)


class Residual(Ball):
    """The ball with additionally a set of cap vertices fixed to given labels (default: the pullback label
    L_eq(pi(y))). Variables = the unfixed cap vertices U. Any degree-d sphere certificate restricts to a degree-<= d
    certificate of the residual, so a consistent degree-d residual dual proves a sphere lower bound of d+1."""

    def __init__(self, m, Leq, U, fixed=None):
        super().__init__(m, Leq)
        cap_all, dom_all, idx_all = self.cap, self.domains, self.idx
        U = sorted(U)
        fixed = dict(fixed or {})
        for i, y in enumerate(cap_all):
            if i not in U and i not in fixed:
                p = y[:-1]
                assert any(p), "e_m cannot take the pullback label; put it in U or fix it explicitly"
                fixed[i] = self.cxe.label(p, self.Leq)
        # fixed part must be non-violating (unary and binary)
        for i, l in fixed.items():
            assert l in dom_all[i], f"fixed label {l} on cap vertex {cap_all[i]} violates an equatorial constraint"
        for i, j in self.edges:
            if i in fixed and j in fixed:
                assert fixed[i] != -fixed[j], f"fixed labels complementary on edge {cap_all[i]} <= {cap_all[j]}"
        self.fixed = fixed
        self.U = U
        # new variable set
        self.cap = [cap_all[i] for i in U]
        self.n = len(U)
        self.idx = {v: k for k, v in enumerate(self.cap)}
        self.domains = []
        for k, i in enumerate(U):
            forb = set()
            for a, b in self.edges:
                if a == i and b in fixed:
                    forb.add(-fixed[b])
                if b == i and a in fixed:
                    forb.add(-fixed[a])
            self.domains.append([l for l in dom_all[i] if l not in forb])
        old = {i: k for k, i in enumerate(U)}
        self.edges = [(old[a], old[b]) for a, b in self.edges if a in old and b in old]
        self.V = set()
        for a, b in self.edges:
            for x in self.domains[a]:
                for y in self.domains[b]:
                    if x == -y:
                        self.V.add(frozenset([(a, x), (b, y)]))
