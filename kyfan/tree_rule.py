"""The explicit degree-n pseudo-solution of the tree gadget (TREE_GADGET_THEOREM.md, statement (i)).

A binary conflict tree `tree` is a nested tuple: a leaf is `'x'`, an internal node is `(magnitude, left, right)`.
The gadget CSP has n+1 variables — the pole (index 0, domain {+n}) and the n leaves (indices 1..n in DFS order,
domain = signed ancestor magnitudes on the root->leaf path, plus the "take P" label -n) — every pair constrained to
avoid complementary labels.

- `rule_patterns(tree, S)` : the block-rule S-patterns (assignments of ancestor picks to the leaves S), one per the
  bottom-up process of the theorem; there are 3^(#full internal nodes) of them for S a proper subset of the leaves.
- `rule_E(tree, n)` : the closed-form pseudo-solution E as a set of partial assignments (frozensets of (var, label)),
  in the same encoding — pole = (0, n), a leaf pick = (leaf, +-magnitude), a leaf taking P = (leaf, -n). |supp E| = 3^n.
- `unique_E(tree, n)` : the SAME object obtained from the solver (build the degree-n SA dual and solve over F_2); the
  gadget's pseudo-solution is unique (rank = #unknowns), so this is ground truth for testing `rule_E`.
"""
import itertools

from .labels import label_set
from .abstract_gadget import build_dual
from . import linalg


def leaves_and_paths(tree):
    """Leaves numbered 1..n in DFS order, and each leaf's root->leaf path as a list of (magnitude, side sign)."""
    leaves = []
    paths = []

    def rec(t, path):
        if t == 'x':
            leaves.append(len(leaves) + 1)
            paths.append(path)
            return
        m, l, r = t
        rec(l, path + [(m, +1)])
        rec(r, path + [(m, -1)])

    rec(tree, [])
    return leaves, paths


def leaf_domains(tree, n):
    """Full gadget domains: {0: {n}} for the pole, {leaf: signed path labels | {-n}} for each leaf."""
    leaves, paths = leaves_and_paths(tree)
    assert len(leaves) == n
    doms = {0: {n}}
    for j in leaves:
        doms[j] = {s * m for m, s in paths[j - 1]} | {-n}
    return doms


def rule_patterns(tree, S):
    """All block-rule S-patterns for a leaf subset S: a set of frozensets of (leaf, signed-ancestor-label) pairs.

    Bottom-up: tokens start at the S-leaves; a set of tokens moving together is a block. At a non-full internal node
    at most one block arrives (from the full side) and takes it (its leaves pick that node) and stops. At a full node
    a block arrives from each side and exactly one of three things happens: left takes and right escapes, right takes
    and left escapes, or both escape merged. A block escaping the root is not a pattern (S = all leaves gives none)."""
    counter = [0]

    def rec(t):
        # returns list of (assignment dict leaf->label, escaping block frozenset of leaves or None)
        if t == 'x':
            counter[0] += 1
            leaf = counter[0]
            return [({}, frozenset([leaf]))] if leaf in S else [({}, None)]
        m, l, r = t
        L = rec(l)
        R = rec(r)
        out = []
        for aL, bL in L:
            for aR, bR in R:
                base = {**aL, **aR}
                if bL is None and bR is None:
                    out.append((base, None))
                elif bR is None:                        # lone block from the left takes v (label +m)
                    a = dict(base); a.update({x: +m for x in bL}); out.append((a, None))
                elif bL is None:
                    a = dict(base); a.update({x: -m for x in bR}); out.append((a, None))
                else:
                    a = dict(base); a.update({x: +m for x in bL}); out.append((a, bR))   # left takes, right escapes
                    a = dict(base); a.update({x: -m for x in bR}); out.append((a, bL))   # right takes, left escapes
                    out.append((base, bL | bR))                                          # both escape, merged
        return out

    return {frozenset(a.items()) for a, b in rec(tree) if b is None}


def rule_E(tree, n):
    """The explicit pseudo-solution E as a set of partial assignments (frozensets of (var, label))."""
    leaves, _ = leaves_and_paths(tree)
    allS = set(leaves)
    E = set()
    cache = {}

    def R(S):
        if S not in cache:
            cache[S] = rule_patterns(tree, S)
        return cache[S]

    for k in range(0, n):                               # S' a proper subset of the leaves (|S'| <= n-1)
        for Sp in itertools.combinations(leaves, k):
            Sp = frozenset(Sp)
            for pat in R(Sp):
                E.add(frozenset(pat))                    # the pattern alone
                E.add(frozenset(pat) | {(0, n)})         # + pole
                rest = allS - Sp
                if rest:                                 # + all remaining leaves taking P (size becomes n, no pole)
                    E.add(frozenset(pat) | {(x, -n) for x in rest})
    return E


def unique_E(tree, n, method="sparse"):
    """Ground truth: build the degree-n SA dual of the gadget and solve over F_2. Returns (domains, E, result).
    The pseudo-solution is unique (rank = #unknowns), so E (a set of frozensets of (var, label)) is well defined."""
    doms = leaf_domains(tree, n)
    doms_list = [doms[i] for i in range(n + 1)]
    labels = label_set(n)
    rows, col = build_dual(doms_list, labels, n)
    if method == "sparse":
        res = linalg.gf2_sparse(rows, len(col), want_solution=True)
    else:
        res = linalg.gf2_dense(rows, len(col), want_solution=True)
    assert res.consistent and res.rank == len(col), "gadget pseudo-solution is not unique"
    E = {frozenset(alpha) for alpha, idx in col.items() if res.solution[idx]}
    return doms, E, res
