"""Abstract chain gadgets: variables 0..n (0 = pole), domains D_i subset of +-[n], all pairs constrained
(no complementary labels). degree(D) = smallest d whose SA dual over F_2 is inconsistent (None if > dmax)."""
import itertools
from collections import defaultdict
from .labels import label_set
from . import linalg


def build_dual(doms, labels, d):
    """Degree-d SA dual (all levels) of the CSP with per-variable domains `doms[i]` and pairwise no-complementary
    constraints. Iterates only in-domain labels, so it stays fast even when `labels` is large."""
    n = len(doms)
    doms = [sorted(D, key=lambda l: (abs(l), l)) for D in doms]

    def comp(alpha):    # complementary pair inside a (necessarily in-domain) partial assignment
        return any(a == -b for (_, a), (_, b) in itertools.combinations(alpha, 2))
    col = {}
    for size in range(d + 1):
        for S in itertools.combinations(range(n), size):
            for ls in itertools.product(*(doms[i] for i in S)):
                alpha = tuple(zip(S, ls))
                if not comp(alpha):
                    col[alpha] = len(col)
    rows = []
    for alpha, c in col.items():
        if len(alpha) == d:
            continue
        used = {i for i, _ in alpha}
        for v in range(n):
            if v in used:
                continue
            f = defaultdict(int); f[c] += 1
            for l in doms[v]:
                beta = tuple(sorted(alpha + ((v, l),)))
                if beta in col:
                    f[col[beta]] -= 1
            rows.append((dict(f), 0))
    rows.append(({col[()]: 1}, 1))
    return rows, col


def unsat(doms):
    return not any(all(a != -b for a, b in itertools.combinations(ls, 2)) for ls in itertools.product(*doms))


def degree(doms, labels=None, dmax=None):
    """(degree or None, info) — degree d means: dual inconsistent at d, consistent at d-1."""
    n = len(doms)
    labels = label_set(max(abs(l) for D in doms for l in D)) if labels is None else labels
    dmax = n if dmax is None else dmax
    info = []
    for d in range(1, dmax + 1):
        rows, col = build_dual(doms, labels, d)
        res = linalg.gf2_dense(rows, len(col), want_solution=True)
        info.append((d, len(col), res.rank, res.consistent, int(res.solution.sum()) if res.consistent else None))
        if not res.consistent:
            return d, info
    return None, info
