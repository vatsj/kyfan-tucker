"""Minimum-weight (fewest-monomial) Tucker certificates on S^2: information-set decoding + an exact CP-SAT model.

A degree-<=d certificate is a solution x of  A x = 1  over F_2: columns of A = violating partial labelings of sizes
2..d (one monomial each), rows = labelings.  Its *size* is wt(x); the minimum is syndrome decoding.

  * rows = R random labelings (random-point primal, one-way sound: every certificate found here is re-verified
    exhaustively on all 4^13 labelings with `kyfan.primal.verify_exhaustive_F2`).
  * Gauss-Jordan once.  The reduced system is very sparse (a size-2 violating monomial is the sum of its four
    size-3 extensions at any third vertex, etc.), so the free columns have weight 4..17 and plain Stern/Dumer on
    random rows degenerates (almost every key is 0).  We therefore run ISD in its sparse-friendly form:
      - Canteaut-Chabaud information-set updates (one pivot swap = one row broadcast) between rounds,
      - at each information set, with the free columns as python ints over the pivot rows: Prange (free vars 0),
        Dumer p=1 (every single free column), Dumer p=1+1 birthday matching on l rows chosen inside the current
        rhs support, and a first-improvement descent in the coset x0 + span(null basis) using the sparse null basis
        (e_c + col_c), i.e. p arbitrary but greedy.
  * `--cpsat`: the exact model.  With pivot rows i and free columns c, x_piv[i] = rhs[i] + sum_c M[i,c] x_c, so
    size = sum_c x_c + sum_i y_i with y_i = rhs[i] XOR (XOR_c M[i,c] x_c): a pure XOR-constrained pseudo-boolean
    minimisation that CP-SAT (AddBoolXor) solves with a proof of optimality when it finishes.

Z_3-invariant reduction (degree 4): unknowns = orbits of the coordinate 3-cycle g (odd order, valid over F_2).  Let w*
be the minimum orbit weight.  (upper) an orbit solution of weight w* expands to <= 3w* monomials.  (lower) if x is any
certificate of size s, then x + gx + g^2x is an invariant certificate (a sum of three certificates) with <= 3s
monomials, hence orbit weight <= s + (#fixed orbits) = s + 8; so  s >= w* - 8.  The true degree-4 minimum therefore lies
in [w* - 8, 3w*] -- the "factor-3 slack".
"""
import sys, os, time, argparse, pickle
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from collections import Counter

from kyfan.complex import SignedComplex
from kyfan.labels import label_set
from kyfan.violating import violating_pairs, is_violating, partial_labelings
from kyfan.group import coordinate_cycle, orbits, act, assert_odd_order
from kyfan.primal import FastRows, verify_exhaustive_F2
from kyfan import linalg

ONE = np.uint64(1)


# ---------------------------------------------------------------- system
def columns(cx, labels, d, z3=False):
    """Violating partial labelings of sizes 2..d as columns (optionally Z_3-orbits). Returns (col: alpha->id, members)."""
    V = violating_pairs(cx, labels)
    items = (a for s in range(2, d + 1) for a in partial_labelings(cx.n_free, labels, s) if is_violating(a, V))
    if not z3:
        col, members = {}, []
        for a in items:
            col[a] = len(members)
            members.append([a])
        return col, members
    gens = [coordinate_cycle(cx.m, list(range(1, cx.m)), coords=(0, 1, 2))]
    assert_odd_order(gens, cx.m, list(range(1, cx.m)))
    return orbits(items, gens, lambda g, a: act(cx, g, a))


def build(m, k, d, nrows, seed=0, z3=False):
    """Packed uint64 system (rows = random labelings, bit ncols = rhs 1). Returns (M, members, labels, cx)."""
    cx = SignedComplex(m)
    labels = label_set(k)
    col, members = columns(cx, labels, d, z3)
    C = len(members)
    frs = []
    for s in range(2, d + 1):
        cs = {a: c for a, c in col.items() if len(a) == s}
        fr = FastRows(cx, labels, s, cs)
        fr.ncols = C
        frs.append(fr)
    rng = np.random.default_rng(seed)
    q = len(labels)
    W = (C + 1 + 63) // 64
    M = np.zeros((nrows, W), dtype=np.uint64)
    t = time.time()
    for r in range(nrows):
        Lidx = rng.integers(0, q, cx.n_free)
        bits = np.zeros(C, dtype=np.int64)
        for fr in frs:
            bits += fr.counts(Lidx)
        bits = np.append(bits & 1, 1).astype(np.uint8)
        packed = np.packbits(bits, bitorder='little').tobytes().ljust(W * 8, b'\0')
        M[r] = np.frombuffer(packed, dtype=np.uint64)
    print(f"  system: {C} columns ({'Z_3 orbits' if z3 else 'monomials'}), {nrows} rows  [{time.time() - t:.0f}s]", flush=True)
    return M, members, labels, cx


# ---------------------------------------------------------------- reduced form
def getbit(M, c):
    return (M[:, c // 64] >> np.uint64(c % 64)) & ONE


def reduce(M, ncols):
    """Gauss-Jordan in place. Returns (rank, pivots: row i has pivot column pivots[i])."""
    t = time.time()
    res = linalg.gf2_dense_packed(M, ncols)
    assert res.consistent, "sampled system inconsistent: no certificate of this degree"
    print(f"  rank {res.rank} (of {ncols} columns)  [{time.time() - t:.0f}s]", flush=True)
    return res.rank, list(res.pivots)


def col_int(M, c):
    """Column c of the packed matrix as a python int over the rows."""
    return int.from_bytes(np.packbits(getbit(M, c).astype(np.uint8), bitorder='little').tobytes(), 'little')


class ISD:
    """Keeps M[:rank] in reduced row-echelon form w.r.t. the current pivot set (row i <-> column piv[i]);
    the rhs is bit `ncols`.  Tracks the best solution found (column-weight, vector)."""

    def __init__(self, M, ncols, rank, pivots, rng):
        self.M, self.C, self.r, self.rng = M[:rank].copy(), ncols, rank, rng
        self.piv = list(pivots)
        self.is_piv = np.zeros(ncols, dtype=bool)
        self.is_piv[self.piv] = True
        self.free = np.nonzero(~self.is_piv)[0]
        self.best, self.best_x = None, None

    # -- information set
    def swap(self, k=1):
        """k Canteaut-Chabaud steps: pick a random free column c and a random row i with M[i,c]=1, pivot on (i,c)."""
        M = self.M
        for _ in range(k):
            for _ in range(50):
                c = int(self.rng.choice(self.free))
                nz = np.nonzero(getbit(M, c))[0]
                if len(nz):
                    break
            else:
                return
            i = int(self.rng.choice(nz))
            idx = nz[nz != i]
            if len(idx):
                M[idx] ^= M[i]
            self.is_piv[self.piv[i]] = False
            self.is_piv[c] = True
            self.piv[i] = c
        self.free = np.nonzero(~self.is_piv)[0]

    def record(self, F, xp, w):
        """F: set of free columns used; xp: python int over pivot rows."""
        if self.best is not None and w >= self.best:
            return False
        x = np.zeros(self.C, dtype=np.uint8)
        x[list(F)] = 1
        pv = np.array(self.piv)
        bits = np.unpackbits(np.frombuffer(xp.to_bytes((self.r + 7) // 8, 'little'), dtype=np.uint8), bitorder='little')[:self.r]
        x[pv[bits.astype(bool)]] = 1
        assert int(x.sum()) == w
        self.best, self.best_x = w, x
        return True

    # -- one round at the current information set
    def round(self, ell=4, descent=True):
        M, r = self.M, self.r
        rhs = col_int(M, self.C)
        fc = [int(c) for c in self.free]
        cols = {c: col_int(M, c) for c in fc}
        # Prange and Dumer p=1
        w0 = rhs.bit_count()
        self.record((), rhs, w0)
        for c in fc:
            v = rhs ^ cols[c]
            self.record((c,), v, 1 + v.bit_count())
        # Dumer p=1+1 on ell rows chosen inside the rhs support (target = all ones there)
        supp = [i for i in range(r) if (rhs >> i) & 1]
        if len(supp) >= ell:
            rows = self.rng.choice(supp, size=ell, replace=False)
            key = lambda v: sum(((v >> int(i)) & 1) << j for j, i in enumerate(rows))
            target = (1 << ell) - 1
            half = len(fc) // 2
            buckets = {}
            for c in fc[:half]:
                buckets.setdefault(key(cols[c]) ^ target, []).append(c)
            for c2 in fc[half:]:
                for c1 in buckets.get(key(cols[c2]), ()):
                    v = rhs ^ cols[c1] ^ cols[c2]
                    self.record((c1, c2), v, 2 + v.bit_count())
        # first-improvement descent in the coset using the sparse null basis
        if descent:
            F, cur = set(), rhs
            w = cur.bit_count()
            improved = True
            while improved:
                improved = False
                order = self.rng.permutation(len(fc))
                for j in order:
                    c = fc[j]
                    v = cur ^ cols[c]
                    nw = (len(F) + (-1 if c in F else 1)) + v.bit_count()
                    if nw < w:
                        F ^= {c}
                        cur, w = v, nw
                        improved = True
            self.record(F, cur, w)
        return w0

    def run(self, seconds, ell=4, swaps=1, log=None):
        t = time.time()
        it = 0
        while time.time() - t < seconds:
            self.swap(swaps)
            self.round(ell)
            it += 1
            if log and it % log == 0:
                print(f"    round {it}: best {self.best}  [{time.time() - t:.0f}s]", flush=True)
        return it


# ---------------------------------------------------------------- exact model
def cpsat_min_weight(M, ncols, rank, pivots, seconds, workers=8, hint=None, log=True):
    """Minimise wt(x) subject to the reduced system, with CP-SAT.  Returns (status, weight, x, bound)."""
    from ortools.sat.python import cp_model
    C, r = ncols, rank
    is_piv = np.zeros(C, dtype=bool)
    is_piv[pivots] = True
    free = np.nonzero(~is_piv)[0]
    bits = np.unpackbits(M[:r].view(np.uint8), axis=1, bitorder='little')
    rhs = bits[:, C].astype(bool)
    B = bits[:, free].astype(bool)                      # r x nfree
    mdl = cp_model.CpModel()
    xf = {int(c): mdl.NewBoolVar(f"x{c}") for c in free}
    y = [mdl.NewBoolVar(f"y{i}") for i in range(r)]
    for i in range(r):
        terms = [xf[int(free[j])] for j in np.nonzero(B[i])[0]]
        # y_i = rhs_i XOR (XOR terms)  <=>  XOR(y_i, terms) = rhs_i  <=>  XOR(y_i, terms, [1 if rhs_i==0]) = 1
        lits = [y[i]] + terms
        if not rhs[i]:
            lits.append(mdl.NewConstant(1))
        if len(lits) == 1:
            mdl.Add(lits[0] == 1)
        else:
            mdl.AddBoolXOr(lits)
    mdl.Minimize(sum(xf.values()) + sum(y))
    if hint is not None:
        for c, v in xf.items():
            mdl.AddHint(v, int(hint[c]))
        pv = np.array(pivots)
        for i in range(r):
            mdl.AddHint(y[i], int(hint[pv[i]]))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_workers = workers
    solver.parameters.log_search_progress = log
    st = solver.Solve(mdl)
    name = solver.StatusName(st)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return name, None, None, solver.BestObjectiveBound()
    x = np.zeros(C, dtype=np.uint8)
    for c, v in xf.items():
        x[c] = solver.Value(v)
    pv = np.array(pivots)
    for i in range(r):
        x[pv[i]] = solver.Value(y[i])
    return name, int(x.sum()), x, solver.BestObjectiveBound()


# ---------------------------------------------------------------- verification
def expand(x, members):
    """Solution vector over columns -> certificate dict alpha -> 1 (orbit columns are expanded)."""
    cert = {}
    for c in np.nonzero(x)[0]:
        for a in members[c]:
            cert[a] = cert.get(a, 0) ^ 1
    return {a: 1 for a, v in cert.items() if v}


def check_solution(M0, x, ncols):
    """A x == rhs on the sampled rows (M0 = the unreduced packed system)."""
    xs = np.packbits(np.append(x, 0).astype(np.uint8), bitorder='little').tobytes().ljust(M0.shape[1] * 8, b'\0')
    xv = np.frombuffer(xs, dtype=np.uint64)
    prod = np.ascontiguousarray(M0 & xv)
    par = np.unpackbits(prod.view(np.uint8), axis=1).sum(axis=1)
    return bool(((par & 1) == getbit(M0, ncols).astype(np.int64)).all())


def report(cx, labels, members, x, M0, C, tag):
    cert = expand(x, members)
    print(f"  [{tag}] column weight {int(x.sum())}; monomials after orbit expansion {len(cert)}; "
          f"sampled rows satisfied: {check_solution(M0, x, C)}", flush=True)
    t = time.time()
    ok = verify_exhaustive_F2(cx, cert, labels, 1)
    print(f"  [{tag}] EXHAUSTIVE verification on all {len(labels)}^{cx.n_free} labelings: {ok}  [{time.time() - t:.0f}s]")
    print(f"  [{tag}] size profile: {dict(sorted(Counter(len(a_) for a_ in cert).items()))}", flush=True)
    return cert, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, default=3)
    ap.add_argument("--rows", type=int, default=9500)
    ap.add_argument("--z3", action="store_true")
    ap.add_argument("--seconds", type=float, default=230)
    ap.add_argument("--ell", type=int, default=4)
    ap.add_argument("--swaps", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cpsat", type=float, default=0, help="seconds for the exact CP-SAT model (0 = skip)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--save", default=None)
    a = ap.parse_args()
    m, k = 3, 2
    print(f"S^2 (m=3, labels +-1,+-2), degree <= {a.d}, z3={a.z3}")
    M, members, labels, cx = build(m, k, a.d, a.rows, a.seed, a.z3)
    C = len(members)
    M0 = M.copy()
    rank, piv = reduce(M, C)
    out = dict(d=a.d, z3=a.z3, rank=rank, ncols=C)
    isd = ISD(M, C, rank, piv, np.random.default_rng(a.seed + 1))
    w0 = isd.round(a.ell, descent=False)
    print(f"  first information set: Prange {w0}, after Dumer p<=2: {isd.best}", flush=True)
    it = isd.run(a.seconds, a.ell, a.swaps, log=50)
    print(f"  ISD: {it} information sets in {a.seconds:.0f}s; best {isd.best}")
    cert, ok = report(cx, labels, members, isd.best_x, M0, C, "ISD")
    out.update(isd_weight=isd.best, isd_cert=cert, isd_verified=ok)
    if a.cpsat > 0:
        print(f"  CP-SAT exact model, {a.cpsat:.0f}s, {a.workers} workers, hinted with the ISD solution", flush=True)
        t = time.time()
        st, w, x, bound = cpsat_min_weight(M, C, rank, piv, a.cpsat, a.workers, hint=isd.best_x, log=False)
        print(f"  CP-SAT: status {st}, weight {w}, lower bound {bound:.1f}  [{time.time() - t:.0f}s]", flush=True)
        out.update(cpsat_status=st, cpsat_weight=w, cpsat_bound=bound)
        if x is not None:
            cert, ok = report(cx, labels, members, x, M0, C, "CP-SAT")
            out.update(cpsat_cert=cert, cpsat_verified=ok)
    if a.save:
        with open(a.save, 'wb') as f:
            pickle.dump(out, f)


if __name__ == "__main__":
    main()
