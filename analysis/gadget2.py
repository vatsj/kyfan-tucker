"""Chain gadget at m with designed equatorial labelings: sigma (top simplex of the equator) carries -(m-1); the star of
sigma's rank-r vertex (vertices >= w_r) is labeled avoiding the values in AVOID[r] (so the lift of w_r keeps those
values' negatives... i.e. keeps -x for x in AVOID[r] available). usage: gadget2.py m ntries avoid_spec
avoid_spec like "1:+1;2:+1,+2" meaning star(w_1) avoids +1, star(w_2) avoids +1 and +2."""
import sys, time, random
from collections import Counter
from kyfan import SignedComplex, label_set
from kyfan.labels import forbidden_pairs, complementary_edges, random_valid_labeling
from kyfan.ball import Ball, Residual
from kyfan.complex import leq, rank

m = int(sys.argv[1]); ntries = int(sys.argv[2]); spec = sys.argv[3] if len(sys.argv) > 3 else ""
cxe = SignedComplex(m - 1); small = label_set(m - 2); big = m - 1
sigma = cxe.top[0]
avoid = {}
for part in filter(None, spec.split(";")):
    r, vals = part.split(":"); avoid[int(r)] = {int(x) for x in vals.split(",")}
fixed = {}
for v in sigma:
    i, s = cxe.rep(v); fixed[i] = -s * big
sig_free = set(fixed)
# per-vertex allowed sets: a free vertex i (rep v) is in star(w_r) if v >= w_r or -v >= w_r (then its label is -lambda(-v))
def allowed_for(i):
    v = cxe.free[i]; al = set(small)
    for r, vals in avoid.items():
        w = sigma[r - 1]
        if leq(w, v): al -= vals                    # lambda(v) must avoid vals
        if leq(w, cxe.neg(v)): al -= {-x for x in vals}   # lambda(-v) = -lambda(v) must avoid vals
    return sorted(al, key=lambda l: (abs(l), -l))
allowed = {i: allowed_for(i) for i in range(cxe.n_free) if i not in sig_free}
print(f"m={m} sigma={sigma} avoid={avoid}; allowed-set sizes: {Counter(len(a) for a in allowed.values())}")

def sample(seed):
    # random_valid_labeling with per-vertex allowed sets: emulate by fixing nothing extra and filtering domains
    from kyfan.labels import forbidden_pairs
    from collections import defaultdict
    F = forbidden_pairs(cxe, label_set(m - 1)); n = cxe.n_free
    nbr = defaultdict(dict)
    for (i, j), S in F.items():
        nbr[i][j] = S; nbr[j][i] = {(b, a) for a, b in S}
    rng = random.Random(seed)
    for _ in range(50):
        dom = {i: list(allowed.get(i, [])) for i in range(n)}; L = dict(fixed)
        for i, l in L.items():
            for j, S in nbr[i].items():
                if j not in L: dom[j] = [b for b in dom[j] if (l, b) not in S]
        nodes = [0]
        def rec():
            nodes[0] += 1
            if nodes[0] > 20000: return False
            free = [i for i in range(n) if i not in L]
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
                if nodes[0] > 20000: return False
            return False
        if rec(): return [L[i] for i in range(n)]
    return None
hist = Counter(); found = []
t = time.time()
for seed in range(ntries):
    L = sample(seed)
    if L is None: hist['none'] += 1; continue
    assert complementary_edges(cxe, L) == 0
    B = Ball(m, L); em = B.e_m()
    U = [em] + [i for i, y in enumerate(B.cap) if any(y[:-1]) and cxe.label(y[:-1], L) == -big]
    R = Residual(m, L, U)
    d, log = R.degree(m - 1)
    hist[d if d else m] += 1
    if d is None or d == m - 1:
        found.append((L, U, d, {R.cap[k]: R.domains[k] for k in range(R.n)}))
print(f"residual degree histogram over {ntries} labelings: {dict(hist)}  [{time.time()-t:.0f}s]")
for L, U, d, doms in found[:6]:
    print(f"  degree {d if d else '>%d (SPHERE DEGREE >= %d)' % (m-1, m)}: domains {doms}\n     L_eq={L}")
import pickle; pickle.dump(found, open(f'analysis/gadget2_m{m}.pkl', 'wb'))
