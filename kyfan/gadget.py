"""Restriction lower bounds via chain gadgets (analysis/gadget2.py, verify_gadget.py).

Lemma (restriction). Let rho be a non-violating partial labeling of the free vertices of S^n and U the unfixed vertices.
A degree-d certificate 1 = sum c_alpha 1_alpha restricts (substitute rho) to a degree-<= d certificate of the residual
CSP on U: domains D(u) = labels minus {-rho(w): w adjacent to u}, constraints = no complementary pair on edges in U.
So a consistent degree-d residual dual proves: Tucker's degree on S^n is > d.

Gadgets found: rho = a valid labeling of the equator S^{n-1} with magnitudes 1..n-1 except one antipodal pair of top
simplices +-sigma labeled -+n, pulled back to the cap; U = the pole e_{n+1} plus the lift of sigma (a top simplex of S^n).
`verify(m, Leq)` re-derives everything from the sphere's own edge list and solves the full-level dual densely.
"""
import itertools
import sys
from collections import defaultdict
from .complex import SignedComplex, leq, neg
from .labels import label_set, forbidden_pairs
from .abstract_gadget import build_dual, unsat
from . import linalg

GADGETS = {   # m -> equatorial labeling (free vertices of SignedComplex(m-1)) whose chain gadget has residual degree m
    4: [-1, -2, -1, -1, 1, 2, 2, 1, 3, 2, -2, 3, 3],
    5: [2, -2, -1, 2, -3, -2, 1, -2, -3, -3, -1, -1, 2, -2, 1, 3, 1, 1, 3, 3, 3, 2, 3, 3, 3, -2, 4, 2, 3, 3, 2, 1, -2, 1, -2, 4, 2, -2, 4, 4],
    6: [-4, 4, -1, -4, -3, -3, -3, -3, 4, -3, 2, -3, -1, 2, -4, -2, -1, -4, -4, -1, -4, -4, -3, 1, -2, 2, -4, -2, 2, -4, -2, 1, -3, -2, 2, -3, -3, 2, -3, -2, 3, 2, 1, 3, 1, 3, 3, 3, 3, 4, 4, 3, 4, 1, 1, 3, 3, -2, 1, 1, -2, 1, 3, 1, 4, -2, -2, 2, 3, -2, 1, 3, 1, -2, 3, -2, 1, 4, -2, 1, 5, 3, 3, 4, 3, 2, 1, -2, 4, 4, -2, 2, 4, -2, 2, 3, 3, 3, 3, 1, 3, 3, -2, 2, 1, 3, 1, 5, 3, 3, 3, 3, 2, 2, 1, 2, 5, -2, 2, 5, 5],
    7: [5, 2, 5, 3, -3, -4, -4, -4, -4, -4, -4, 5, 3, -3, 4, 4, -3, 5, 4, 2, 5, 4, 4, -2, -2, -3, 5, -2, 4, 4, 3, -2, -3, -3, -4, -4, 5, -3, 5, 3, 1, 3, -1, -5, -5, -1, 2, 2, 2, 1, -2, -1, 3, -5, 3, -1, 3, 3, -4, 3, -2, 3, -4, -1, -1, 2, -4, -2, -2, -2, -5, 4, -1, 4, -1, 3, -5, -2, -1, -5, -5, -1, 3, -1, -1, 1, -2, -2, -4, -4, -4, -4, -4, -4, 1, 4, 4, -3, -2, -1, 4, -2, -1, -3, 4, -3, 4, -3, -1, -1, -3, -2, -3, -2, -1, -3, -3, -4, -1, -3, -2, 1, 3, -2, 1, 1, -2, 3, 1, 3, -4, -2, -4, 3, 1, -4, 1, 3, 2, 1, 3, -4, -4, 1, 3, 1, 2, 3, 4, -2, 4, 3, 1, 1, 2, 2, 5, -2, 1, -2, 3, 1, 5, 1, 4, 2, 1, -4, 1, 3, 3, -2, 3, 5, 1, 4, 4, 1, 2, 1, 1, 4, 2, 4, 1, -2, -4, 4, 5, 1, 2, 1, 4, -2, -2, 1, 1, -4, -4, 2, 1, 5, 1, -2, 1, 2, 3, 1, 3, 3, 3, 1, 1, 3, 3, 1, 3, 3, 3, 2, -4, 3, 3, 2, 1, 3, 3, 1, 1, -2, -2, 3, 1, 4, 3, 1, 5, 3, -2, -2, 1, 1, 6, 3, 3, 4, 3, -2, -4, 3, -4, -4, 3, 3, 5, 5, -2, 3, 1, 4, 5, 5, 4, 3, 4, 1, 3, -2, 1, 5, -2, 2, 3, -2, 1, 1, 1, -4, -4, -2, 1, 3, -2, 1, 3, 3, 4, 4, -2, 3, 2, 3, 1, -2, 3, 1, 4, 3, 3, 4, 2, 1, 3, -2, 3, 3, 3, 2, 3, 2, 1, 1, 3, 2, 1, -2, 3, 4, 3, 4, 1, -2, 1, 6, 4, 4, 3, 4, 1, 1, -2, 1, -4, 3, 2, 3, 3, 4, 4, 4, 4, 1, 4, 4, 4, -2, 1, 4, 4, 4, 6, 4, 4, 4, -2, 1, 1, -2, 2, 6, -2, 2, 6, 6],
    8: [-4, 4, 4, -2, -6, -4, -1, 4, -1, -3, -1, -2, -3, -5, 1, -5, 3, -5, -4, 4, 1, 6, -5, 2, -4, 4, -5, -5, 3, 4, -5, 2, -3, -5, -1, 6, -3, 4, -1, 6, 5, -4, -6, 5, 1, -6, 1, 5, 1, 5, -4, -4, 6, -4, -6, 5, 5, -2, -3, -3, -4, -1, 5, -4, 5, -2, -2, 2, 1, 2, 1, -4, -4, -2, 6, -2, 5, -4, -4, 6, -4, 5, -2, -2, -2, 2, 5, 5, -1, 6, 6, -1, 6, -3, 2, 6, 2, 1, 6, -4, 1, 1, -2, 2, 6, 6, -5, -5, 6, 3, -2, -2, -3, -4, -3, -1, -4, -3, -5, -2, 6, -6, 5, 1, -6, 3, 1, 1, 3, -2, 2, 1, -6, -6, -6, -6, 4, 1, 1, -6, -3, 2, -6, 5, 5, -2, 3, -2, 2, -6, -6, 1, 4, -2, 4, 3, -5, -5, -6, -6, -6, 1, 1, 1, -5, -5, 2, 2, 2, 4, 4, -6, -2, -2, -2, 3, 1, -6, 1, -6, -2, 3, -5, -5, -5, -5, 2, -5, 1, 1, -2, -2, -2, 2, 2, -5, 4, 4, -3, -5, 4, -3, 2, 2, 1, -6, -6, -4, 3, 5, 5, 5, 1, -4, 5, 1, -6, -2, 4, -6, 5, -4, -3, -6, -6, -3, 5, -6, -2, -5, -6, 1, 4, 1, 1, -6, -2, -2, 1, 1, 1, 4, -6, -6, 3, 4, -6, -3, -6, -3, -6, -1, -4, -6, -2, -2, -6, -4, -6, 1, 1, -5, -5, 4, 1, 2, -4, -6, 1, -6, -6, -5, -2, -6, 2, 2, -5, 2, -1, -4, -5, -2, -5, 3, -4, 2, 5, 1, 1, -2, -2, -2, 5, -4, -6, 5, -6, -6, -2, 1, 3, 2, -4, -3, 5, -4, 5, 3, 5, -2, 3, -6, 5, 5, 1, 5, 1, 1, -2, 2, 1, -4, 5, -6, 5, -2, -2, 5, -3, 5, -4, 5, 5, -2, 3, -2, -2, 1, -5, -4, 3, -4, -4, 3, 1, 1, 1, 2, -3, 1, -4, -3, 3, -5, -3, -3, -3, -4, 2, -1, -1, 3, -1, -1, 2, -3, -3, 2, 1, -3, 4, -2, -2, 1, 6, -3, 4, 4, -3, -2, 4, -2, 2, 2, 2, 2, 4, 2, 4, -2, -3, -5, 2, -3, 4, -3, 1, -5, 4, -5, 2, 2, 2, -3, 1, -3, 4, -3, -2, -3, -3, 2, 4, -3, -5, 6, 6, -2, 2, -5, -3, 2, -3, -2, -5, -2, -2, 6, -3, 2, -3, -3, -3, -2, -2, -3, -5, -3, -5, 2, -3, -2, 4, -5, -2, 1, 2, 5, 5, 5, 5, -2, 5, -2, 2, -3, -3, 1, 6, -3, 5, -2, -3, 5, 6, 2, 6, -3, -3, 6, -3, -2, 2, -3, -3, 4, 1, 6, -2, 3, 6, -3, 2, 2, 1, 4, 1, 4, 6, 6, -3, 6, -5, -3, 6, -5, 4, 4, -2, 6, 2, -5, 3, 6, -3, 4, 6, 6, 6, 2, -5, -5, -5, -3, -5, -2, -3, -5, -5, -5, -3, -5, -2, -5, -3, -3, 5, -3, -3, 2, 5, -2, 5, -2, 5, 2, 6, -3, 5, -3, -3, -3, 5, 5, -3, -3, 2, 2, 6, -3, -3, 6, -2, 6, 5, -3, 2, 6, -3, 6, 5, 6, 5, -3, 5, 2, 6, 6, 6, 5, 6, -3, 5, 2, -3, -3, -3, 5, -3, -3, 2, 2, -3, 3, -5, 6, 3, 3, 6, 6, 2, 6, 6, -5, -2, 3, -2, -3, -3, -3, -3, -3, -2, -2, 6, -3, 6, 1, -3, 1, 4, 4, 1, -2, 1, 1, 1, 2, 1, -3, -3, 1, 4, 4, 1, -3, 5, 2, 2, 4, 5, -3, 4, -2, 1, 1, 2, -5, -3, -3, 1, -3, -3, 1, 6, 2, -3, 4, -3, 1, -2, -3, -5, 2, -3, -5, 4, -3, 4, 4, -3, 1, -3, -5, -5, 1, -5, 1, 1, 1, 2, -3, 2, 4, 1, 1, -2, 6, 1, -5, -3, 2, 2, 6, -3, 4, -2, -3, 2, 2, 2, 4, 1, -3, -2, 1, -2, 2, -3, 5, 6, 5, 1, 6, 1, 6, 2, 5, -3, -3, 5, 6, 4, 5, 6, 1, 1, 2, 6, 6, 1, 4, -2, 6, 2, -3, -3, 4, 7, -3, 1, -2, -2, -3, -3, -3, 6, 4, -3, -2, -2, -3, 6, -5, -5, 1, 1, 1, 3, 1, -3, -5, -5, 6, -5, 6, -3, -2, 6, -2, -5, -3, -5, 2, 4, 6, 4, 6, 6, 1, 2, 1, 5, 5, 1, -2, -2, 5, 5, 1, 1, -3, 6, 1, -3, 5, -2, 2, -3, 5, 2, 6, 5, -3, 6, -3, 1, 6, 1, 1, 6, -2, 6, -2, -3, 2, -3, -3, 6, 6, -3, -2, -2, -3, -3, 6, 6, 2, 5, -3, 6, 5, 5, 1, 1, 6, 3, -5, 1, 1, 1, -3, 6, 2, 1, 2, 6, -2, 6, -5, -3, -3, 6, -3, 2, 6, -2, 6, -2, -3, 1, 5, 5, 4, 5, -2, -2, 1, -2, 2, 1, 1, 1, -3, -2, -2, -3, 1, -3, 5, 5, -3, 4, -2, 4, -2, -2, -5, 1, -3, 1, -2, -3, 4, -5, 1, 2, -3, 5, 1, 4, -2, 1, -2, -2, -5, -3, -3, -5, -2, -3, -5, -3, -3, 1, -3, -3, 4, -3, -2, -2, 4, 1, -5, 1, 1, 4, -5, -5, -5, -2, -3, -3, -3, -3, 4, -2, -2, -5, -5, -2, 1, 1, 1, 1, 4, -2, 5, 5, -2, 1, -3, -2, 1, 4, 5, 5, -2, 5, -3, -3, -2, -3, -3, -3, 4, 4, 5, 1, -3, 1, 4, 4, -2, 4, -2, 1, 1, -3, 5, 1, 7, -2, 1, 4, 5, 2, -3, -3, 2, 4, -2, 3, -2, 5, 2, 2, 1, 4, -3, -5, 4, -3, 1, 2, 2, 1, 4, -5, -5, -2, 4, -3, -5, -5, -3, 2, 4, -5, -5, -5, -3, 2, 5, 5, 2, -3, 5, 1, 5, 1, -3, 1, -2, 2, -3, 1, 1, -3, 1, 5, -2, -2, 2, -3, -3, 5, -2, 5, 2, 2, 5, 3, 5, 1, 3, 3, 5, 1, -3, -2, 1, 7, 5, 3, -2, -3, 2, -3, -3, 2, 5, 5, 5, -2, -3, 2, -3, -3, 2, 1, -3, 1, 1, -2, 2, 2, 1, 1, 7, -2, 3, -2, -2, -3, -3, -3, 2, 7, -3, 3, 7, 7],
}


def chain_U(m, Leq):
    """U = e_m plus the cap vertices whose equatorial projection carries label -(m-1)."""
    cxe = SignedComplex(m - 1)
    cap = [v + (1,) for v in itertools.product([-1, 0, 1], repeat=m - 1)]
    return sorted(y for y in cap if not any(y[:-1]) or cxe.label(y[:-1], Leq) == -(m - 1))


def residual(m, Leq, U):
    """(rho, domains, constraints) of the residual on U from the sphere's own edges. rho: free index -> label."""
    cx, cxe, labels = SignedComplex(m), SignedComplex(m - 1), label_set(m - 1)
    rho = {}
    for w in cxe.verts:
        i, s = cx.rep(w + (0,)); l = s * cxe.label(w, Leq)
        assert rho.setdefault(i, l) == l
    for y in (v + (1,) for v in itertools.product([-1, 0, 1], repeat=m - 1)):
        if y in U:
            continue
        assert any(y[:-1])
        i, s = cx.rep(y); assert i not in rho; rho[i] = s * cxe.label(y[:-1], Leq)
    for x, y in cx.edges:   # rho non-violating
        (i, si), (j, sj) = cx.rep(x), cx.rep(y)
        assert not (i in rho and j in rho and si * rho[i] == -sj * rho[j])
    Uf = {cx.rep(y)[0]: (cx.rep(y)[1], y) for y in U}
    dom = {y: set(labels) for y in U}
    cons = set()
    for x, y in cx.edges:
        (i, si), (j, sj) = cx.rep(x), cx.rep(y)
        for (a, sa), (b, sb) in (((i, si), (j, sj)), ((j, sj), (i, si))):
            if a in Uf and b in rho:
                s_a, ya = Uf[a]
                dom[ya].discard(s_a * (-sb * rho[b] * sa))
        if i in Uf and j in Uf:
            cons.add((Uf[i][1], Uf[j][1], si * Uf[i][0], sj * Uf[j][0]))
    return rho, dom, cons


def verify(m, Leq=None, U=None, d=None):
    """Full independent check. Returns dict with domains, unsat, dual rank/unknowns, consistent, support of E."""
    Leq = GADGETS[m] if Leq is None else Leq
    U = chain_U(m, Leq) if U is None else sorted(U)
    d = m - 1 if d is None else d
    labels = label_set(m - 1)
    rho, dom, cons = residual(m, Leq, U)

    def viol(assign):
        return any(yi in assign and yj in assign and a * assign[yi] == -b * assign[yj] for yi, yj, a, b in cons)

    unsat = not any(not viol(dict(zip(U, ls))) for ls in itertools.product(*(sorted(dom[y]) for y in U)))
    sub_sat = all(any(not viol(dict(zip(S, ls))) for ls in itertools.product(*(sorted(dom[y]) for y in S)))
                  for S in itertools.combinations(U, len(U) - 1))

    def violating(alpha):
        return any(l not in dom[y] for y, l in alpha) or viol(dict(alpha))
    col = {}
    for size in range(d + 1):
        for S in itertools.combinations(U, size):
            for ls in itertools.product(labels, repeat=size):
                alpha = tuple(zip(S, ls))
                if not violating(alpha):
                    col[alpha] = len(col)
    rows = []
    for alpha, c in col.items():
        if len(alpha) == d:
            continue
        used = {y for y, _ in alpha}
        for v in U:
            if v in used:
                continue
            f = defaultdict(int); f[c] += 1
            for l in labels:
                beta = tuple(sorted(alpha + ((v, l),)))
                if beta in col:
                    f[col[beta]] -= 1
            rows.append((dict(f), 0))
    rows.append(({col[()]: 1}, 1))
    res = linalg.gf2_dense(rows, len(col), want_solution=True)
    ok = res.consistent and all(sum(int(res.solution[c]) * x for c, x in f.items()) % 2 == b for f, b in rows)
    return dict(U=U, domains={y: sorted(dom[y], key=lambda l: (abs(l), -l)) for y in U}, unsat=unsat,
                proper_subsets_satisfiable=sub_sat, unknowns=len(col), equations=len(rows), rank=res.rank,
                consistent=res.consistent, solution_verified=ok, support=int(res.solution.sum()) if ok else None)


# ---------------------------------------------------------------------------------------------------------------
# Statement (ii): realizability. A deterministic construction of an equatorial labeling whose pole-chain residual is
# a prescribed binary conflict tree. chain_domains computes the residual domains directly from up/down-sets (no full
# S^n edge list, so it scales), validated against residual() for m <= 8.
# ---------------------------------------------------------------------------------------------------------------

def caterpillar_tree(n):
    """The caterpillar binary conflict tree on n leaves: (1; x, (2; x, (3; x, ... (n-1; x, x)))). Magnitudes 1..n-1."""
    t = 'x'
    for mag in range(n - 1, 0, -1):
        t = (mag, 'x', t)
    return t


def _tree_leaf_paths(tree):
    paths = []

    def rec(t, path):
        if t == 'x':
            paths.append(path)
            return
        mag, l, r = t
        rec(l, path + [(mag, +1)])
        rec(r, path + [(mag, -1)])

    rec(tree, [])
    return paths


def chain_domains(m, Leq):
    """Residual domains D(y_r), r = 0..n (r=0 the pole), for the pole-chain U = e_m u lift(sigma), sigma the equatorial
    top simplex w_r = -(e_1+...+e_r). Computed from equatorial down-sets and cap up-neighbours only. Returns a list of
    n+1 sets over label_set(n); reproduces `residual()` on the stored GADGETS (checked for m <= 8)."""
    n = m - 1
    cxe = SignedComplex(m - 1)
    labels = label_set(n)
    w = {r: tuple(-1 if i < r else 0 for i in range(m - 1)) for r in range(n + 1)}
    sigma_proj = {w[r] for r in range(1, n + 1)}
    out = []
    for r in range(n + 1):
        wr = w[r]
        forb = set()
        for u in cxe.verts:
            if leq(u, wr):                                  # equatorial down-neighbour (always fixed and present)
                forb.add(-cxe.label(u, Leq))
            elif leq(wr, u) and u not in sigma_proj:        # cap up-neighbour, fixed (not a U lift)
                forb.add(-cxe.label(u, Leq))
        out.append(set(labels) - forb)
    return out


def realize(m, tree=None, leaf_order=None):
    """Deterministically construct an equatorial labeling of S^{m-1} realizing `tree` (default: caterpillar) on the
    pole-chain, with leaf k (DFS order) assigned to chain rank `leaf_order[k]` (default: n, n-1, ..., 1).

    sigma is fixed to -+n; every other free vertex gets an 'allowed' set forcing the target domains, then the CSP is
    solved by greedy forward-checking (MRV, smallest label first). Returns (Leq, backtracks) or (None, reason).
    Empirically backtracks == 0 for all tested m, i.e. the construction is forced."""
    n = m - 1
    tree = caterpillar_tree(n) if tree is None else tree
    leaf_order = list(range(n, 0, -1)) if leaf_order is None else leaf_order
    cxe = SignedComplex(m - 1)
    small = label_set(n - 1)
    paths = _tree_leaf_paths(tree)
    assert len(paths) == n and len(leaf_order) == n
    target = {leaf_order[k]: {s * mag for mag, s in paths[k]} for k in range(n)}   # rank -> pick-label set
    sigma = cxe.top[0]
    w = {r: sigma[r - 1] for r in range(1, n + 1)}
    fixed = {}
    for v in sigma:
        i, s = cxe.rep(v)
        fixed[i] = -s * n
    def in_N(v, r):
        return (leq(w[r], v) or leq(v, w[r])) and v not in sigma
    allowed = {}
    for i, v in enumerate(cxe.free):
        if i in fixed:
            continue
        al = set(small)
        for r in range(1, n + 1):
            if in_N(v, r):
                al -= {-l for l in target[r]}
            if in_N(neg(v), r):
                al -= set(target[r])
        allowed[i] = sorted(al, key=lambda l: (abs(l), l))
        if not allowed[i]:
            return None, f"empty allowed set at free vertex {i}"
    F = forbidden_pairs(cxe, label_set(n))
    N = cxe.n_free
    nbr = defaultdict(dict)
    for (i, j), S in F.items():
        nbr[i][j] = S
        nbr[j][i] = {(b, a) for a, b in S}
    dom = {i: list(allowed.get(i, [])) for i in range(N)}
    L = dict(fixed)
    for i, l in L.items():
        for j, S in nbr[i].items():
            if j not in L:
                dom[j] = [b for b in dom[j] if (l, b) not in S]
    backtracks = [0]
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))

    def rec():
        free = [i for i in range(N) if i not in L]
        if not free:
            return True
        i = min(free, key=lambda v: len(dom[v]))
        first = True
        for a in dom[i]:
            if not first:
                backtracks[0] += 1
            first = False
            L[i] = a
            saved = {}
            good = True
            for j, S in nbr[i].items():
                if j not in L:
                    saved[j] = dom[j]
                    dom[j] = [b for b in dom[j] if (a, b) not in S]
                    if not dom[j]:
                        good = False
                        break
            if good and rec():
                return True
            for j, d in saved.items():
                dom[j] = d
            del L[i]
        return False

    if not rec():
        return None, "CSP infeasible"
    return [L[i] for i in range(N)], backtracks[0]


def check_realization(m, Leq, method="sparse"):
    """Fast end-to-end check that Leq's pole-chain residual is a degree-n gadget (=> Tucker F_2 degree on S^{m-1} = m):
    UNSAT, and the degree-n dual is consistent with a unique pseudo-solution. Uses chain_domains (no full sphere)."""
    n = m - 1
    doms = chain_domains(m, Leq)
    doms_list = list(doms)
    rows, col = build_dual(doms_list, label_set(n), n)
    if method == "sparse":
        res = linalg.gf2_sparse(rows, len(col))
    else:
        res = linalg.gf2_dense(rows, len(col))
    return dict(unsat=unsat(doms_list), unknowns=len(col), rank=res.rank, consistent=res.consistent,
                unique=res.consistent and res.rank == len(col))
