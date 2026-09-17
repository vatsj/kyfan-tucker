"""Theorem 1' (docs/size_lower_bound.md): per-flag parity constraints and the number t(U) of full-flag terms.

Fix the pole-chain flag U = e_m < lift(sigma) of S^{m-1} (n = m-1; ranks 0..n, rank 0 = pole). For a gadget rho of U
(a non-violating restriction whose residual on U has no degree-n certificate) let Z_rho be the space of e on the full
assignments of the residual domains with sum_{x contains beta} e(x) = 0 for every violating beta of size <= n.
Lemma B: for every degree-(n+1) certificate C, T_U = {alpha in C: supp alpha = U} satisfies
    sum_{x in T_U} e(x) = sum_x e(x)   for every e in Z_rho, every gadget rho of U.
t(U) = min |T| over T satisfying all these parities; every degree-(n+1) certificate has >= (#flag classes) * t(U)
full-flag terms (#flag classes = m! 2^m / 2, since alpha and -alpha use the same free vertices).

Gadget families (all in vertex-label coordinates: a full assignment is the tuple (lambda(y_0), ..., lambda(y_n))):
  * `--family orbit`  : the label-group orbit (signed permutations of the magnitudes, 2^n n! elements) of the explicit
                        gadget's residual (reverse caterpillar) -- computable for every n.
  * `--family trees`  : the same for every tree/leaf-order that `kyfan.gadget.realize` realizes on the pole chain.
  * `--family perms`  : every permutation of the chain positions applied to the explicit domains (and to the realized
                        trees) that is *realizable* as a restriction (CP-SAT: a valid labeling of the other free
                        vertices producing exactly those domains), then their label orbits. On S^2 this is exactly the
                        exhaustive family (label orbit x chain reversal, 16 types).
  * `--family all`    : (m = 3 only) every non-violating restriction rho at the flag, 4^10 of them -- reproduces
                        flag_hitting_bound_s2.py (16 gadget domain types, t = 8) in these coordinates.
The parity system is solved exactly by CP-SAT (min sum x subject to sum_{S} x = rhs mod 2), cross-checked by the
branch-and-bound odd-hitting-set search when every Z_rho is one-dimensional.

    python analysis/flag_hitting_bound.py --m 3 --family orbit trees all
    python analysis/flag_hitting_bound.py --m 4 --family orbit trees perms --save-realizations analysis/flag_hitting_bound_m4_realizations.txt

`--save-realizations` writes the realizing labelings found by `--family perms` (one per line: target domains, then the
labels of the free vertices off the pole chain); `load_realizations` reads them back and `check_realization` re-verifies
each one without CP-SAT (tests/test_paper_numbers.py).
"""
import argparse
import itertools
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kyfan.complex import SignedComplex, leq
from kyfan.labels import label_set
from kyfan.abstract_gadget import build_dual, unsat
from kyfan.gadget import chain_domains, realize, caterpillar_tree
from kyfan.realize_explicit import explicit_Leq
from kyfan import linalg


# ----------------------------------------------------------------------------------------------------------------
# Z_rho and the parity constraints
# ----------------------------------------------------------------------------------------------------------------
def gf2_nullspace(A):
    """Basis of {v : A v = 0} over F_2 for a dense uint8 matrix A."""
    M = A.copy() % 2
    R, C = M.shape
    piv, r = [], 0
    for c in range(C):
        if r == R:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        M[[r, p]] = M[[p, r]]
        mask = M[:, c].astype(bool); mask[r] = False
        M[mask] ^= M[r]
        piv.append(c); r += 1
    pivset = set(piv)
    basis = []
    for fc in [c for c in range(C) if c not in pivset]:
        v = np.zeros(C, dtype=np.uint8); v[fc] = 1
        for k, pc in enumerate(piv):
            if M[k, fc]:
                v[pc] = 1
        basis.append(v)
    return basis


def Zspace(doms):
    """(points, basis): points = full assignments in the domain product (tuples), basis = F_2 basis of Z_rho.
    Constraints: for every size-n beta (drop one variable) containing a complementary pair, the sum over the
    dropped variable's domain vanishes (smaller violating beta are sums of these)."""
    n1 = len(doms)
    pts = list(itertools.product(*[sorted(D, key=lambda l: (abs(l), l)) for D in doms]))
    pidx = {p: k for k, p in enumerate(pts)}
    rows = []
    for v in range(n1):
        others = [u for u in range(n1) if u != v]
        seen = set()
        for p in pts:
            beta = tuple(p[u] for u in others)
            if beta in seen:
                continue
            seen.add(beta)
            if any(a == -b for a, b in itertools.combinations(beta, 2)):
                rows.append([pidx[p[:v] + (l,) + p[v + 1:]] for l in doms[v] if p[:v] + (l,) + p[v + 1:] in pidx])
    A = np.zeros((max(len(rows), 1), len(pts)), dtype=np.uint8)
    for r, row in enumerate(rows):
        A[r, row] = 1
    return pts, gf2_nullspace(A)


def is_gadget(doms, n):
    """UNSAT residual whose degree-n SA dual is consistent."""
    if not unsat(doms):
        return False
    rows, col = build_dual(doms, label_set(n), n)
    return linalg.gf2_dense(rows, len(col)).consistent


def constraints_of(doms):
    """Parity constraints (frozenset of full assignments, rhs) from one gadget."""
    pts, basis = Zspace(doms)
    out = []
    for e in basis:
        S = frozenset(pts[k] for k in np.nonzero(e)[0])
        out.append((S, len(S) % 2))
    return out, len(basis)


# ----------------------------------------------------------------------------------------------------------------
# gadget families at the pole chain
# ----------------------------------------------------------------------------------------------------------------
def signed_perms(n):
    for pi in itertools.permutations(range(1, n + 1)):
        for eps in itertools.product((1, -1), repeat=n):
            yield {s * i: s * eps[i - 1] * pi[i - 1] for i in range(1, n + 1) for s in (1, -1)}


def label_orbit(doms, n):
    """Distinct domain tuples in the label-group orbit of `doms`."""
    seen = {}
    for g in signed_perms(n):
        key = tuple(frozenset(g[l] for l in D) for D in doms)
        seen.setdefault(key, [set(D) for D in key])
    return list(seen.values())


def family_orbit(m):
    n = m - 1
    D = chain_domains(m, explicit_Leq(m))
    assert is_gadget(D, n)
    return {'explicit': D}, label_orbit(D, n)


def family_trees(m):
    """Every tree / leaf order realizable by `realize` (n = 3: two magnitude assignments of the caterpillar shape,
    six leaf orders), then their label orbits."""
    n = m - 1
    base = {}
    shapes = set()

    def trees(mags):
        if not mags:
            return ['x']
        out = []
        for i, a in enumerate(mags):
            rest = mags[:i] + mags[i + 1:]
            for L in trees_split(rest):
                for tl in trees(L[0]):
                    for tr in trees(L[1]):
                        out.append((a, tl, tr))
        return out

    def trees_split(rest):
        outs = []
        k = len(rest)
        for mask in range(1 << k):
            outs.append(([rest[i] for i in range(k) if mask >> i & 1], [rest[i] for i in range(k) if not mask >> i & 1]))
        return outs

    for tree in trees(list(range(1, n))):
        # a binary conflict tree must have n leaves with n-1 internal nodes: every internal node has 2 children
        for order in itertools.permutations(range(1, n + 1)):
            Leq, bt = realize(m, tree=tree, leaf_order=list(order))
            if Leq is None:
                continue
            D = chain_domains(m, Leq)
            key = tuple(frozenset(x) for x in D)
            if key in shapes:
                continue
            if not is_gadget(D, n):
                continue
            shapes.add(key)
            base[(tree, order)] = D
    fam = {}
    for D in base.values():
        for D2 in label_orbit(D, n):
            fam.setdefault(tuple(frozenset(x) for x in D2), D2)
    return base, list(fam.values())


def pole_chain(m):
    """U = [e_m, y_1, ..., y_n] with y_r = lift of w_r = -(e_1 + ... + e_r) (matches `chain_domains`)."""
    return [tuple([-1] * r + [0] * (m - 1 - r) + [1]) for r in range(m)]


def realizable(m, target, seconds=60):
    """Is there a non-violating labeling of the free vertices off the pole chain U whose residual domains on U are
    exactly `target` (list of sets, rank order)?  CP-SAT feasibility. Returns the labeling (dict free index -> label)
    or None."""
    from ortools.sat.python import cp_model
    n = m - 1
    cx = SignedComplex(m)
    labels = label_set(n)
    U = pole_chain(m)
    Uset = set(U)
    Uidx = {cx.rep(y)[0] for y in U}
    others = [i for i in range(cx.n_free) if i not in Uidx]
    model = cp_model.CpModel()
    x = {(i, l): model.NewBoolVar(f"x{i}_{l}") for i in others for l in labels}
    for i in others:
        model.AddExactlyOne(x[i, l] for l in labels)
    nbrs = {y: [] for y in U}
    for a, b in cx.edges:
        (i, si), (j, sj) = cx.rep(a), cx.rep(b)
        if i in Uidx and j in Uidx:
            continue
        if i in Uidx or j in Uidx:
            for y, w in ((a, b), (b, a)):
                if y in Uset:                     # edges at -y_r are the antipodes of edges at y_r: same constraint
                    nbrs[y].append(cx.rep(w))
            continue
        for la in labels:
            for lb in labels:
                if si * la == -sj * lb:
                    model.AddBoolOr([x[i, la].Not(), x[j, lb].Not()])
    for r, y in enumerate(U):
        forb = set(labels) - set(target[r])
        allowed = {-l for l in forb}                     # lambda(w) must lie in -forb for every neighbour w of y
        for i, s in nbrs[y]:
            for l in labels:
                if s * l not in allowed:
                    model.Add(x[i, l] == 0)
        for f in forb:                                    # and every forbidden label is produced by some neighbour
            model.AddBoolOr([x[i, -s * f] for i, s in nbrs[y]])
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = 8
    st = solver.Solve(model)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        assert st == cp_model.INFEASIBLE, solver.StatusName(st)
        return None
    L = {i: next(l for l in labels if solver.Value(x[i, l])) for i in others}
    check_realization(m, target, L)      # independent of the solver
    return L


def check_realization(m, target, L):
    """Assert that the labeling L (dict: free index off the pole chain -> label) is non-violating off U and that the
    residual domains it produces on the pole chain are exactly `target`.  Pure recomputation, no CP-SAT."""
    n = m - 1
    cx = SignedComplex(m)
    labels = label_set(n)
    U = pole_chain(m)
    Uset = set(U)
    Uidx = {cx.rep(y)[0] for y in U}
    assert set(L) == {i for i in range(cx.n_free) if i not in Uidx} and all(l in labels for l in L.values())
    nbrs = {y: [] for y in U}
    for a, b in cx.edges:
        (i, si), (j, sj) = cx.rep(a), cx.rep(b)
        if i in L and j in L:
            assert si * L[i] != -sj * L[j], "complementary edge off the pole chain"
        elif i in L or j in L:
            for y, w in ((a, b), (b, a)):
                if y in Uset:
                    nbrs[y].append(cx.rep(w))
    got = [set(labels) for _ in U]
    for r, y in enumerate(U):
        for i, s in nbrs[y]:
            got[r].discard(-s * L[i])
    assert got == [set(t) for t in target], (got, target)
    return True


def save_realizations(path, found):
    with open(path, 'w') as f:
        f.write("# realizing labelings found by flag_hitting_bound.py --family perms (CP-SAT); one per line:\n")
        f.write("# <base> <perm> | <domain of rank 0> ; ... ; <domain of rank n> | i:l pairs (free index off the pole chain, label)\n")
        for key, (name, perm, L) in found.items():
            doms = ' ; '.join(' '.join(str(l) for l in sorted(k, key=lambda l: (abs(l), l))) for k in key)
            labs = ' '.join(f"{i}:{L[i]}" for i in sorted(L))
            f.write(f"{name} {','.join(map(str, perm))} | {doms} | {labs}\n")


def load_realizations(path):
    """Returns a list of (name, perm, target domains as list of sets, L)."""
    out = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        head, doms, labs = line.split('|')
        name, perm = head.strip().rsplit(None, 1)          # the base name may contain spaces (a tree literal)
        perm = tuple(int(x) for x in perm.split(','))
        target = [set(int(l) for l in d.split()) for d in doms.split(';')]
        L = {int(i): int(l) for i, l in (tok.split(':') for tok in labs.split())}
        out.append((name, perm, target, L))
    return out


def family_perms(m):
    n = m - 1
    bases = {'explicit': chain_domains(m, explicit_Leq(m))}
    tb, _ = family_trees(m)
    for k, D in tb.items():
        bases[str(k)] = D
    found = {}
    tried = set()
    for name, D in bases.items():
        for perm in itertools.permutations(range(n + 1)):
            T = [set(D[perm[r]]) for r in range(n + 1)]
            key = tuple(frozenset(t) for t in T)
            if key in tried:
                continue
            tried.add(key)
            L = realizable(m, T)
            if L is not None:
                assert is_gadget(T, n)
                found[key] = (name, perm, L)
    print(f"  realizable position permutations of the tree gadgets: {len(found)} of {len(tried)} tried:")
    for key, (name, perm, _) in found.items():
        print(f"    {name} perm {perm}: {[sorted(k, key=lambda l: (abs(l), l)) for k in key]}")
    fam = {}
    for key in found:
        for D2 in label_orbit([set(k) for k in key], n):
            fam.setdefault(tuple(frozenset(x) for x in D2), D2)
    return found, list(fam.values())


def family_all_s2():
    """m = 3: enumerate every restriction rho of the 10 free vertices off the pole chain, in vertex coordinates."""
    m, n = 3, 2
    cx = SignedComplex(m)
    labels = label_set(n)
    U = [(0, 0, 1), (-1, 0, 1), (-1, -1, 1)]           # e_3 < e_3 - e_1 < e_3 - e_1 - e_2
    Ur = [cx.rep(v) for v in U]
    Uidx = {i for i, _ in Ur}
    others = [i for i in range(cx.n_free) if i not in Uidx]
    negU = {tuple(-c for c in v) for v in U}
    Eoff = [(cx.rep(x), cx.rep(y)) for x, y in cx.edges if cx.rep(x)[0] in others and cx.rep(y)[0] in others]
    Emix = []      # (fixed rep, U rank, sign): lambda(y_r) != sign * lambda(x)
    for x, y in cx.edges:
        for a, b in ((x, y), (y, x)):
            if cx.rep(a)[0] in others and (b in U or b in negU):
                r = U.index(b) if b in U else U.index(tuple(-c for c in b))
                Emix.append((cx.rep(a), r, -1 if b in U else +1))
    types = {}
    for ls in itertools.product(labels, repeat=len(others)):
        L = dict(zip(others, ls))
        if any(si * L[i] == -sj * L[j] for (i, si), (j, sj) in Eoff):
            continue
        dom = [set(labels) for _ in U]
        for (i, si), r, sgn in Emix:
            dom[r].discard(sgn * si * L[i])
        key = tuple(frozenset(D) for D in dom)
        types[key] = types.get(key, 0) + 1
    gadgets = [[set(D) for D in key] for key in types if is_gadget([set(D) for D in key], n)]
    print(f"  m=3 all restrictions: {sum(types.values())} non-violating, {len(types)} domain types, {len(gadgets)} gadget types")
    return gadgets


# ----------------------------------------------------------------------------------------------------------------
# minimum |T|
# ----------------------------------------------------------------------------------------------------------------
def min_parity_hitting(cons, seconds=600, workers=8):
    """Exact minimum |T| with sum_{x in T cap S} 1 = rhs (mod 2) for every (S, rhs) in cons, via CP-SAT."""
    from ortools.sat.python import cp_model
    pts = sorted(set().union(*(S for S, _ in cons)))
    idx = {p: k for k, p in enumerate(pts)}
    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x{k}") for k in range(len(pts))]
    for j, (S, rhs) in enumerate(cons):
        kk = model.NewIntVar(0, len(S) // 2 + 1, f"k{j}")
        model.Add(sum(x[idx[p]] for p in S) == 2 * kk + rhs)
    model.Minimize(sum(x))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    st = solver.Solve(model)
    T = [pts[k] for k in range(len(pts)) if solver.Value(x[k])]
    return solver.StatusName(st), len(T), solver.BestObjectiveBound(), T


def min_odd_hitting_bb(supports):
    """Branch and bound for the odd-hitting number (all rhs = 1), as in flag_hitting_bound_s2.py."""
    verts = sorted(set().union(*supports))
    best = [len(verts) + 1]

    def bb(chosen, parity):
        if len(chosen) >= best[0]:
            return
        bad = [k for k, s in enumerate(supports) if parity[k] % 2 == 0]
        if not bad:
            best[0] = len(chosen); return
        for v in sorted(supports[bad[0]]):
            if v in chosen:
                continue
            chosen.add(v)
            for k, t in enumerate(supports):
                if v in t: parity[k] += 1
            bb(chosen, parity)
            for k, t in enumerate(supports):
                if v in t: parity[k] -= 1
            chosen.remove(v)
    bb(set(), [0] * len(supports))
    return best[0]


def run(m, gadgets, name, seconds, bb=True):
    n = m - 1
    t0 = time.time()
    cons, dims = [], []
    for D in gadgets:
        c, dim = constraints_of(D)
        cons += c; dims.append(dim)
    supp_sizes = sorted(len(S) for S, _ in cons)
    touched = set().union(*(S for S, _ in cons))
    maxdeg = max(sum(1 for S, _ in cons if p in S) for p in touched)
    print(f"  [{name}] {len(gadgets)} gadgets, Z-dimensions {sorted(set(dims))}, {len(cons)} parity constraints, "
          f"support sizes {min(supp_sizes)}..{max(supp_sizes)}, rhs all odd: {all(r == 1 for _, r in cons)}, "
          f"touching {len(touched)} of {(2 * n) ** (n + 1)} full labelings, max multiplicity {maxdeg}  [{time.time() - t0:.1f}s]")
    st, t, bound, T = min_parity_hitting(cons, seconds)
    print(f"  [{name}] CP-SAT: {st}, t = {t}, bound {bound:.1f}  [{time.time() - t0:.1f}s]")
    if bb and all(r == 1 for _, r in cons) and len(cons) <= 20:
        tb = min_odd_hitting_bb([S for S, _ in cons])
        print(f"  [{name}] branch-and-bound odd hitting number: {tb}")
        assert tb == t
    classes = 1
    for i in range(1, m + 1):
        classes *= 2 * i
    classes //= 2
    print(f"  [{name}] => every degree-{n + 1} certificate on S^{n} has >= {classes} x {t} = {classes * t} size-{n + 1} terms"
          + ("" if st == "OPTIMAL" else " (CP-SAT not optimal: only the bound is rigorous)"))
    return t, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, default=3)
    ap.add_argument('--family', nargs='+', default=['orbit'])
    ap.add_argument('--seconds', type=float, default=600)
    ap.add_argument('--save-realizations', default=None, help="write the labelings found by --family perms to this file")
    a = ap.parse_args()
    m, n = a.m, a.m - 1
    print(f"S^{n} (m={m}), pole-chain flag, labels +-1..+-{n}; full labelings of the flag: {(2 * n) ** (n + 1)}")
    for fam in a.family:
        if fam == 'orbit':
            base, gadgets = family_orbit(m)
            D = base['explicit']
            pts, basis = Zspace(D)
            print(f"  explicit gadget domains {[sorted(x, key=lambda l: (abs(l), l)) for x in D]}: {len(pts)} points, "
                  f"dim Z = {len(basis)}, |supp| = {int(basis[0].sum()) if basis else None} (3^(n-1) = {3 ** (n - 1)})")
            print(f"  label orbit: {len(gadgets)} distinct domain tuples of 2^n n! = {2 ** n * math.factorial(n)}")
            run(m, gadgets, 'orbit', a.seconds)
        elif fam == 'trees':
            base, gadgets = family_trees(m)
            print(f"  realizable tree gadgets at the pole chain: {len(base)} "
                  f"({sorted(set(str(k[0]) for k in base))}); with label orbits: {len(gadgets)} domain tuples")
            run(m, gadgets, 'trees', a.seconds)
        elif fam == 'perms':
            found, gadgets = family_perms(m)
            print(f"  with label orbits: {len(gadgets)} gadget domain tuples")
            if a.save_realizations:
                save_realizations(a.save_realizations, found)
            run(m, gadgets, 'perms', a.seconds, bb=False)
        elif fam == 'all':
            assert m == 3
            gadgets = family_all_s2()
            run(m, gadgets, 'all', a.seconds)


if __name__ == '__main__':
    main()
