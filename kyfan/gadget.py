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
from collections import defaultdict
from .complex import SignedComplex
from .labels import label_set
from . import linalg

GADGETS = {   # m -> equatorial labeling (free vertices of SignedComplex(m-1)) whose chain gadget has residual degree m
    4: [-1, -2, -1, -1, 1, 2, 2, 1, 3, 2, -2, 3, 3],
    5: [2, -2, -1, 2, -3, -2, 1, -2, -3, -3, -1, -1, 2, -2, 1, 3, 1, 1, 3, 3, 3, 2, 3, 3, 3, -2, 4, 2, 3, 3, 2, 1, -2, 1, -2, 4, 2, -2, 4, 4],
    6: [-4, 4, -1, -4, -3, -3, -3, -3, 4, -3, 2, -3, -1, 2, -4, -2, -1, -4, -4, -1, -4, -4, -3, 1, -2, 2, -4, -2, 2, -4, -2, 1, -3, -2, 2, -3, -3, 2, -3, -2, 3, 2, 1, 3, 1, 3, 3, 3, 3, 4, 4, 3, 4, 1, 1, 3, 3, -2, 1, 1, -2, 1, 3, 1, 4, -2, -2, 2, 3, -2, 1, 3, 1, -2, 3, -2, 1, 4, -2, 1, 5, 3, 3, 4, 3, 2, 1, -2, 4, 4, -2, 2, 4, -2, 2, 3, 3, 3, 3, 1, 3, 3, -2, 2, 1, 3, 1, 5, 3, 3, 3, 3, 2, 2, 1, 2, 5, -2, 2, 5, 5],
    7: [5, 2, 5, 3, -3, -4, -4, -4, -4, -4, -4, 5, 3, -3, 4, 4, -3, 5, 4, 2, 5, 4, 4, -2, -2, -3, 5, -2, 4, 4, 3, -2, -3, -3, -4, -4, 5, -3, 5, 3, 1, 3, -1, -5, -5, -1, 2, 2, 2, 1, -2, -1, 3, -5, 3, -1, 3, 3, -4, 3, -2, 3, -4, -1, -1, 2, -4, -2, -2, -2, -5, 4, -1, 4, -1, 3, -5, -2, -1, -5, -5, -1, 3, -1, -1, 1, -2, -2, -4, -4, -4, -4, -4, -4, 1, 4, 4, -3, -2, -1, 4, -2, -1, -3, 4, -3, 4, -3, -1, -1, -3, -2, -3, -2, -1, -3, -3, -4, -1, -3, -2, 1, 3, -2, 1, 1, -2, 3, 1, 3, -4, -2, -4, 3, 1, -4, 1, 3, 2, 1, 3, -4, -4, 1, 3, 1, 2, 3, 4, -2, 4, 3, 1, 1, 2, 2, 5, -2, 1, -2, 3, 1, 5, 1, 4, 2, 1, -4, 1, 3, 3, -2, 3, 5, 1, 4, 4, 1, 2, 1, 1, 4, 2, 4, 1, -2, -4, 4, 5, 1, 2, 1, 4, -2, -2, 1, 1, -4, -4, 2, 1, 5, 1, -2, 1, 2, 3, 1, 3, 3, 3, 1, 1, 3, 3, 1, 3, 3, 3, 2, -4, 3, 3, 2, 1, 3, 3, 1, 1, -2, -2, 3, 1, 4, 3, 1, 5, 3, -2, -2, 1, 1, 6, 3, 3, 4, 3, -2, -4, 3, -4, -4, 3, 3, 5, 5, -2, 3, 1, 4, 5, 5, 4, 3, 4, 1, 3, -2, 1, 5, -2, 2, 3, -2, 1, 1, 1, -4, -4, -2, 1, 3, -2, 1, 3, 3, 4, 4, -2, 3, 2, 3, 1, -2, 3, 1, 4, 3, 3, 4, 2, 1, 3, -2, 3, 3, 3, 2, 3, 2, 1, 1, 3, 2, 1, -2, 3, 4, 3, 4, 1, -2, 1, 6, 4, 4, 3, 4, 1, 1, -2, 1, -4, 3, 2, 3, 3, 4, 4, 4, 4, 1, 4, 4, 4, -2, 1, 4, 4, 4, 6, 4, 4, 4, -2, 1, 1, -2, 2, 6, -2, 2, 6, 6],
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
