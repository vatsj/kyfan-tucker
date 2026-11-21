"""Structure of the estimation instance on S^n with k labels, over antipodal classes (free vertices)."""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from kyfan.complex import SignedComplex, leq, neg, rank

def instance(m, k):
    cx = SignedComplex(m); free = cx.free; N = len(free)
    labels = [s * a for a in range(1, k + 1) for s in (1, -1)]           # index -> label
    edges = {}
    for x, y in cx.edges:
        (i, si), (j, sj) = cx.rep(x), cx.rep(y)
        if i != j:
            a, b = (i, j) if i < j else (j, i)
            edges[(a, b)] = si * sj          # constraint: L[a] != -(si*sj) * L[b]
    # flag classes: maximal chains up to antipode, as tuples of (free index, sign)
    byrank = {}
    for v in cx.verts: byrank.setdefault(rank(v), []).append(v)
    chains = [[v] for v in byrank[1]]
    for r in range(2, m + 1):
        chains = [c + [w] for c in chains for w in byrank[r] if leq(c[-1], w)]
    flags = {}
    for c in chains:
        reps = tuple(cx.rep(v) for v in c)
        key = tuple(sorted(i for i, _ in reps))
        if key not in flags: flags[key] = reps      # one representative per class; the antipode has all signs flipped
    return cx, N, labels, edges, list(flags.values())

def posalt(ls):
    s = sorted(ls, key=abs)
    return len({abs(l) for l in s}) == len(s) and all((l > 0) == (i % 2 == 0) for i, l in enumerate(s))
