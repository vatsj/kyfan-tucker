import itertools
import random
from collections import defaultdict


def label_set(k):
    """[+1,-1,+2,-2,...,+k,-k]"""
    return [s * i for i in range(1, k + 1) for s in (1, -1)]


def posalt(ls):
    """Positively alternating: distinct magnitudes and, sorted by magnitude, signs +,-,+,-,...  (posalt(()) is True)."""
    s = sorted(ls, key=abs)
    return len({abs(l) for l in s}) == len(s) and all((l > 0) == (i % 2 == 0) for i, l in enumerate(s))


def negalt(ls):
    return posalt([-l for l in ls])


def complementary(ls):
    return any(a == -b for a, b in itertools.combinations(ls, 2))


def pos_alt_count(cx, L):
    """A_+(L): number of positively alternating top simplices."""
    return sum(posalt([cx.label(v, L) for v in c]) for c in cx.top)


def complementary_edges(cx, L):
    return sum(1 for x, y in cx.edges if cx.label(x, L) == -cx.label(y, L))


def is_valid(cx, L):
    return complementary_edges(cx, L) == 0


def random_labeling(cx, labels, rng):
    return [rng.choice(labels) for _ in range(cx.n_free)]


def forbidden_pairs(cx, labels):
    """For free i<j: set of (a,b) with a on i and b on j making a complementary edge."""
    F = defaultdict(set)
    for (i, si), (j, sj) in cx.free_edges():
        for a in labels:
            for b in labels:
                if si * a == -sj * b:
                    F[(min(i, j), max(i, j))].add((a, b) if i < j else (b, a))
    return F


def dfs_valid_labelings(cx, labels, cap=200000, seed=0):
    """DFS enumeration of valid labelings (no complementary edge), at most `cap`, in a seeded random order."""
    F = forbidden_pairs(cx, labels)
    n = cx.n_free
    L = [None] * n
    out = []
    rng = random.Random(seed)

    def ok(i, a):
        return all((L[j], a) not in F.get((j, i), ()) for j in range(i))

    def rec(i):
        if len(out) >= cap:
            return
        if i == n:
            out.append(tuple(L))
            return
        ls = labels[:]
        rng.shuffle(ls)
        for a in ls:
            if ok(i, a):
                L[i] = a
                rec(i + 1)
                L[i] = None

    rec(0)
    return out


def octahedron_pos_alt(l1, l2, l3):
    """A_+ on the UNSUBDIVIDED octahedron with labels l_i on e_i (and -l_i on -e_i)."""
    c = 0
    for s in itertools.product((1, -1), repeat=3):
        if posalt((s[0] * l1, s[1] * l2, s[2] * l3)):
            c += 1
    return c


def random_valid_labeling(cx, labels, fixed=None, seed=0, node_budget=20000, restarts=50, allowed=None):
    """A random valid labeling (no complementary edge) extending `fixed` (dict free index -> label), by randomized
    backtracking with forward checking and smallest-domain-first, restarting after node_budget nodes. None if not found.
    `allowed`: label list for the unfixed vertices (default: all of `labels`)."""
    F = forbidden_pairs(cx, labels)
    n = cx.n_free
    nbr = defaultdict(dict)   # nbr[i][j] = set of (a on i, b on j) forbidden
    for (i, j), S in F.items():
        nbr[i][j] = S
        nbr[j][i] = {(b, a) for a, b in S}
    rng = random.Random(seed)
    for _ in range(restarts):
        dom = {i: list(allowed if allowed is not None else labels) for i in range(n)}
        L = {}
        ok = True
        for i, l in (fixed or {}).items():
            L[i] = l
        # propagate fixed
        for i, l in L.items():
            for j, S in nbr[i].items():
                if j not in L:
                    dom[j] = [b for b in dom[j] if (l, b) not in S]
                    if not dom[j]:
                        ok = False
        if not ok:
            return None
        nodes = [0]

        def rec():
            nodes[0] += 1
            if nodes[0] > node_budget:
                return False
            free = [i for i in range(n) if i not in L]
            if not free:
                return True
            i = min(free, key=lambda v: (len(dom[v]), rng.random()))
            vals = dom[i][:]
            rng.shuffle(vals)
            for a in vals:
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
                if nodes[0] > node_budget:
                    return False
            return False

        if rec():
            return [L[i] for i in range(n)]
    return None
