"""Realize a target tree gadget on S^{m-1}: target domains D_r (r = rank of the lift, 1..n) for the chain through the pole.
The fixed neighbours of lift(w_r) are N_r = (star(w_r) u down(w_r)) \ sigma, so the equatorial labeling must use on N_r
only labels l with -l not in D_r (and cover them). usage: gadget3.py m ntries tree leaf_order
tree e.g. "(1;(2;x,x),(3;x,(4;x,x)))", leaf_order = permutation of ranks for the leaves in left-to-right order, e.g. "5,4,3,2,1"."""
import sys, time, random, itertools, ast
from collections import Counter, defaultdict
from kyfan import SignedComplex, label_set
from kyfan.labels import forbidden_pairs, complementary_edges
from kyfan.complex import leq, neg
from kyfan.gadget import verify, chain_U

m = int(sys.argv[1]); ntries = int(sys.argv[2]); tree_s = sys.argv[3]; order = [int(x) for x in sys.argv[4].split(",")]
n = m - 1; P = n
def parse(s):
    s = s.strip()
    if s == 'x': return 'x'
    assert s[0] == '(' and s[-1] == ')'
    body = s[1:-1]; mag, rest = body.split(';', 1)
    depth = 0
    for k, ch in enumerate(rest):
        depth += (ch == '(') - (ch == ')')
        if ch == ',' and depth == 0: return (int(mag), parse(rest[:k]), parse(rest[k + 1:]))
tree = parse(tree_s)
def leaves(t, path=()):
    if t == 'x': return [set(path)]
    return leaves(t[1], path + (t[0],)) + leaves(t[2], path + (-t[0],))
paths = leaves(tree); assert len(paths) == n == len(order)
target = {r: paths[k] for k, r in enumerate(order)}          # rank -> path label set (domain = path + {-P})
print(f"m={m}: target domains by rank: { {r: sorted(D) + [-P] for r, D in target.items()} }, pole {{{P}}}")
cxe = SignedComplex(m - 1); small = label_set(n - 1)
sigma = cxe.top[0]; w = {r: sigma[r - 1] for r in range(1, n + 1)}
fixed = {}
for v in sigma:
    i, s = cxe.rep(v); fixed[i] = -s * P
def in_N(v, r):
    return (leq(w[r], v) or leq(v, w[r])) and v not in sigma
allowed = {}
for i, v in enumerate(cxe.free):
    if i in fixed: continue
    al = set(small)
    for r in range(1, n + 1):
        if in_N(v, r): al -= {-l for l in target[r]}
        if in_N(neg(v), r): al -= set(target[r])
    allowed[i] = sorted(al, key=lambda l: (abs(l), -l))
print("allowed-set sizes:", Counter(len(a) for a in allowed.values()), " empty:", sum(1 for a in allowed.values() if not a))
F = forbidden_pairs(cxe, label_set(n)); N = cxe.n_free
nbr = defaultdict(dict)
for (i, j), S in F.items():
    nbr[i][j] = S; nbr[j][i] = {(b, a) for a, b in S}
def sample(seed, budget=30000):
    rng = random.Random(seed)
    for _ in range(30):
        dom = {i: list(allowed.get(i, [])) for i in range(N)}; L = dict(fixed)
        for i, l in L.items():
            for j, S in nbr[i].items():
                if j not in L: dom[j] = [b for b in dom[j] if (l, b) not in S]
        nodes = [0]
        def rec():
            nodes[0] += 1
            if nodes[0] > budget: return False
            free = [i for i in range(N) if i not in L]
            if not free: return True
            i = min(free, key=lambda v: (len(dom[v]), rng.random()))
            vals = dom[i][:]; rng.shuffle(vals)
            for a in vals:
                L[i] = a; saved = {}; good = True
                for j, S in nbr[i].items():
                    if j not in L:
                        saved[j] = dom[j]; dom[j] = [b for b in dom[j] if (a, b) not in S]
                        if not dom[j]: good = False; break
                if good and rec(): return True
                for j, d in saved.items(): dom[j] = d
                del L[i]
                if nodes[0] > budget: return False
            return False
        if rec(): return [L[i] for i in range(N)]
    return None
t = time.time(); hist = Counter(); found = []
for seed in range(ntries):
    L = sample(seed)
    if L is None: hist['no labeling'] += 1; continue
    assert complementary_edges(cxe, L) == 0
    r = verify(m, L)
    doms = r["domains"]; U = r["U"]
    key = 'target' if all(set(doms[y]) == (target[sum(1 for t in y[:-1] if t)] | {-P}) for y in U if any(y[:-1])) else 'other'
    hist[(key, 'consistent' if r["consistent"] else 'inconsistent')] += 1
    if r["consistent"]:
        found.append((L, r))
    print(f"seed {seed}: domains {[doms[y] for y in U]} -> degree-{n} dual {'CONSISTENT (sphere degree >= %d)' % m if r['consistent'] else 'inconsistent'} rank {r['rank']}/{r['unknowns']}", flush=True)
print(f"\nhistogram: {dict(hist)}  [{time.time()-t:.0f}s]")
import pickle; pickle.dump(found, open(f'analysis/gadget3_m{m}.pkl', 'wb'))
if found: print("FOUND: L_eq =", found[0][0])
