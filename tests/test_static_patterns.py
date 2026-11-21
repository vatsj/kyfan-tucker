"""The two ingredients of the written proof of the extension lemma (docs/tree_gadget.md), checked without the solver.

Lemma 1: the static conditions (C1)-(C3) describe exactly the block-rule S-patterns.
Lemma 2 (path lemma): for k : [1..j] -> [1..j+1] with k(i) >= i, the admissible t are {j+1} if k is the identity and
otherwise either none or exactly {a, k(a)} with a the least i such that k(i) > i.
"""
import itertools

import pytest

from kyfan.tree_rule import rule_patterns, all_tree_shapes


def _build(tree):
    nodes = {}
    leafnodes = []

    def rec(t, parent):
        nid = len(nodes)
        nodes[nid] = dict(parent=parent, kids=[], leaf=None, mag=None, under=set())
        if t == 'x':
            leafnodes.append(nid)
            nodes[nid]['leaf'] = len(leafnodes)
            nodes[nid]['under'] = {len(leafnodes)}
        else:
            m, l, r = t
            nodes[nid]['mag'] = m
            a = rec(l, nid)
            b = rec(r, nid)
            nodes[nid]['kids'] = [a, b]
            nodes[nid]['under'] = nodes[a]['under'] | nodes[b]['under']
        return nid

    rec(tree, None)
    return nodes, leafnodes


def static_patterns(tree, S):
    """All pick assignments on S satisfying (C1)-(C3), in the encoding of rule_patterns."""
    nodes, leafnodes = _build(tree)
    S = set(S)
    full = {v: nodes[v]['under'] <= S for v in nodes}
    node_of_leaf = {nodes[v]['leaf']: v for v in leafnodes}

    def ancestors(v):
        out = []
        u = nodes[v]['parent']
        while u is not None:
            out.append(u)
            u = nodes[u]['parent']
        return out

    anc = {x: ancestors(v) for x, v in node_of_leaf.items()}

    def at_or_above(p, u):          # p == u or p is an ancestor of u
        return p == u or p in ancestors(u)

    top = {}
    for x in S:
        v = node_of_leaf[x]
        while nodes[v]['parent'] is not None and full[nodes[v]['parent']]:
            v = nodes[v]['parent']
        top[x] = nodes[v]['parent']

    Sl = sorted(S)
    out = set()
    for picks in itertools.product(*[anc[x] for x in Sl]):
        pi = dict(zip(Sl, picks))
        if any(top[x] is None or not at_or_above(top[x], pi[x]) for x in Sl):                 # (C1)
            continue
        ok = True
        for u in nodes:
            if not nodes[u]['kids']:
                continue
            takes = []
            for c in nodes[u]['kids']:
                if not full[c]:
                    continue
                arrive = [x for x in nodes[c]['under'] if at_or_above(pi[x], u)]              # B(c)
                take = [x for x in arrive if pi[x] == u]
                if take and len(take) != len(arrive):                                         # (C2)
                    ok = False
                    break
                takes.append(bool(arrive) and len(take) == len(arrive))
            if not ok or (len(takes) == 2 and all(takes)):                                    # (C3)
                ok = False
                break
        if ok:
            def label(x):
                u = pi[x]
                return nodes[u]['mag'] * (+1 if x in nodes[nodes[u]['kids'][0]]['under'] else -1)
            out.add(frozenset((x, label(x)) for x in Sl))
    return out


@pytest.mark.parametrize("n", [2, 3, 4, 5, 6])
def test_static_characterization_matches_rule(n):
    for tree in all_tree_shapes(n):
        for r in range(n):
            for S in itertools.combinations(range(1, n + 1), r):
                assert static_patterns(tree, S) == rule_patterns(tree, frozenset(S)), (tree, S)


def _admissible(t, k, j):
    c = t
    for i in range(1, j + 1):
        if c == i:
            if k[i] == i:
                return False
            c = k[i]
        elif k[i] not in (i, c):
            return False
    return c == j + 1


@pytest.mark.parametrize("j", range(0, 8))
def test_path_lemma(j):
    for ks in itertools.product(*[range(i, j + 2) for i in range(1, j + 1)]):
        k = {i + 1: v for i, v in enumerate(ks)}
        adm = [t for t in range(1, j + 2) if _admissible(t, k, j)]
        N = [i for i in k if k[i] != i]
        if not N:
            assert adm == [j + 1]
        else:
            a = min(N)
            assert adm in ([], sorted({a, k[a]}))
