"""Independent verification of a chain gadget (restriction lower bound), from scratch: no Ball/Residual code reused.
usage: verify_gadget.py m pickle_index   (reads analysis/gadget2_m{m}.pkl, entry index)

Checks:
 (a) rho (equator labeling + pullback on the fixed cap vertices) is a non-violating partial labeling of the free
     vertices of S^{m-1}, using the sphere complex's own edge list;
 (b) the residual CSP on U (all cap vertices not fixed): domains = labels minus {-rho(w): w adjacent to u, w fixed},
     binary constraints = complementary labels on sphere edges inside U; compare with the pickled domains;
 (c) UNSAT by brute force; and every (|U|-1)-subset satisfiable;
 (d) the degree-(m-1) SA dual of the residual, built with consistency at EVERY level and E = 0 on every partial
     labeling that contains a forbidden label or a complementary pair, is consistent over F_2 (dense solver) — and the
     returned pseudo-solution satisfies every equation. Hence no degree-(m-1) residual certificate, hence (restriction
     lemma) no degree-(m-1) sphere certificate: Tucker's F_2 degree on S^{m-1} is >= m.
"""
import sys, itertools, pickle
from collections import defaultdict
import numpy as np
from kyfan import SignedComplex, label_set
from kyfan.complex import leq
from kyfan import linalg

m = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else 0
found = pickle.load(open(f'analysis/gadget2_m{m}.pkl', 'rb'))
gad = [g for g in found if g[2] is None]
print(f"{len(gad)} gadgets with residual degree > {m-1} in the pickle; verifying #{k}")
Leq, U_idx, _, doms_claimed = gad[k]
cx = SignedComplex(m); cxe = SignedComplex(m - 1); labels = label_set(m - 1)
assert cx.free == cx.free  # noqa
# --- rho on the sphere's free vertices ---
def rep_label(v, l):        # actual vertex v carries label l -> (free index, label on the free rep)
    i, s = cx.rep(v); return i, s * l
rho = {}
for w in cxe.verts:                                   # equator: lambda(w) = L_eq via the equatorial complex
    i, l = rep_label(w + (0,), cxe.label(w, Leq)); rho.setdefault(i, l); assert rho[i] == l
cap = [v + (1,) for v in itertools.product([-1, 0, 1], repeat=m - 1)]
U = sorted(cap[i] for i in U_idx) if all(isinstance(i, int) for i in U_idx) else sorted(U_idx)
for y in cap:
    if y in U: continue
    p = y[:-1]; assert any(p)
    i, l = rep_label(y, cxe.label(p, Leq)); assert i not in rho; rho[i] = l
print(f"(a) rho fixes {len(rho)} of {cx.n_free} free vertices; U = {U}")
# (a) non-violating: check every sphere edge
bad = 0
for x, y in cx.edges:
    (i, si), (j, sj) = cx.rep(x), cx.rep(y)
    if i in rho and j in rho and si * rho[i] == -sj * rho[j]: bad += 1
print(f"(a) complementary edges inside rho: {bad}"); assert bad == 0
# (b) residual domains and constraints from the sphere's edges
Uf = {}   # free index -> sign, for U vertices
for y in U:
    i, s = cx.rep(y); Uf[i] = s
dom = {y: set(labels) for y in U}
edgesU = set()
for x, y in cx.edges:
    (i, si), (j, sj) = cx.rep(x), cx.rep(y)
    for (a, sa), (b, sb) in (((i, si), (j, sj)), ((j, sj), (i, si))):
        if a in Uf and b in rho:
            # vertex with free index a carries label sa*lambda_free; the U-vertex y_a = free rep * Uf[a]
            # label of the U-vertex (as a cap vertex) is Uf[a] * lambda_free(a); forbidden: sa*lambda_free(a) == -sb*rho[b]
            ya = cap_y = next(y for y in U if cx.rep(y)[0] == a)
            forb_free = -sb * rho[b] * sa           # lambda_free(a) forbidden value
            dom[ya].discard(Uf[a] * forb_free)        # cap-vertex label forbidden value
    if i in Uf and j in Uf:
        yi = next(y for y in U if cx.rep(y)[0] == i); yj = next(y for y in U if cx.rep(y)[0] == j)
        edgesU.add((yi, yj, si * Uf[i], sj * Uf[j]))  # constraint: (si*Uf[i]) * lab(yi) != -(sj*Uf[j]) * lab(yj)
dom = {y: sorted(d, key=lambda l: (abs(l), -l)) for y, d in dom.items()}
print(f"(b) domains from the sphere: {dom}")
print(f"(b) claimed:                 {dict((y, sorted(d, key=lambda l: (abs(l), -l))) for y, d in doms_claimed.items())}")
assert {y: set(d) for y, d in dom.items()} == {y: set(d) for y, d in doms_claimed.items()}
def viol(assign):   # assign: dict y -> label
    return any(yi in assign and yj in assign and a * assign[yi] == -b * assign[yj] for yi, yj, a, b in edgesU)
print(f"(b) {len(edgesU)} residual edge constraints (all pairs of U comparable: {len(edgesU) == len(U) * (len(U) - 1) // 2})")
# (c) brute force
sols = [dict(zip(U, ls)) for ls in itertools.product(*(dom[y] for y in U)) if not viol(dict(zip(U, ls)))]
print(f"(c) satisfying assignments of the residual: {len(sols)} (need 0)"); assert not sols
for S in itertools.combinations(U, len(U) - 1):
    ok = any(not viol(dict(zip(S, ls))) for ls in itertools.product(*(dom[y] for y in S)))
    print(f"(c) subset {S} satisfiable: {ok}")
# (d) full dual at degree d = m-1, all levels, all labels
d = m - 1
def violating(alpha):
    return any(l not in dom[y] for y, l in alpha) or viol(dict(alpha))
col = {}
for size in range(0, d + 1):
    for S in itertools.combinations(U, size):
        for ls in itertools.product(labels, repeat=size):
            alpha = tuple(zip(S, ls))
            if not violating(alpha): col[alpha] = len(col)
rows = []
for alpha, c in col.items():
    if len(alpha) == d: continue
    used = {y for y, _ in alpha}
    for v in U:
        if v in used: continue
        f = defaultdict(int); f[c] += 1
        for l in labels:
            beta = tuple(sorted(alpha + ((v, l),)))
            if beta in col: f[col[beta]] -= 1
        rows.append((dict(f), 0))
rows.append(({col[()]: 1}, 1))
res = linalg.gf2_dense(rows, len(col), want_solution=True)
print(f"(d) degree-{d} dual: {len(col)} unknowns (all levels), {len(rows)} equations, rank {res.rank}, consistent = {res.consistent}")
assert res.consistent
E = res.solution
assert all(sum(int(E[c]) * x for c, x in f.items()) % 2 == b for f, b in rows)
print(f"(d) pseudo-solution verified on every equation; E[empty] = {E[col[()]]}; support {int(E.sum())}")
print(f"\nVERIFIED: no degree-{d} F_2 certificate for the residual  =>  Tucker F_2 degree on S^{m-1} >= {m}.")
