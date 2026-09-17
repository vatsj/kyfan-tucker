"""Spread-lemma experiments (SIZE_LOWER_BOUND.md §4) with MCMC diagnostics.

The set G_m of *gadget-compatible* equatorial labelings for the pole flag of S^{m-1}: valid labelings L of the
equator (the [m-1]-complex, labels +-[n], n = m-1) that agree with the explicit labeling on the chain
sigma = {w_r = -(e_1+...+e_r)} and whose pole-chain residual domains CONTAIN the reverse-caterpillar domains
(D(0) contains +n; D(r) contains T[r] = {+(n-r)} u {-(n-r+1),...,-n}).  Every such L is a gadget (tree theorem,
label-agnostic), and the uniform distribution on G_m is the natural candidate for the D_F of the spread lemma.

Sampler: Glauber dynamics (heat bath) on G_m -- pick a uniformly random movable free vertex, resample its label
uniformly among the labels allowed by the static (domain) constraints and the current neighbours.  Stationary
distribution = uniform on the connected component of the explicit gadget (ergodicity NOT proved; see diagnostics).

Diagnostics (new here; `gadget_spread_experiments.py` had none):
  * several independent chains (different seeds; all start at the explicit gadget, the only known member of G_m),
  * per scalar statistic: mean, integrated autocorrelation time tau (Sokal window), ESS, split-R-hat across chains,
  * Hamming distance from the start, and the fraction of moves with <= 1 allowed label (frozen vertices),
  * match probabilities  P[ L'|_A = L|_A ]  estimated between two INDEPENDENT chains (no self-matching), with
    standard errors over random sets A, and the per-vertex decay rate beta fitted from log P vs |A|.

Usage: python analysis/spread_mcmc.py --m 7 --chains 4 --sweeps 3000
"""
import sys, os, time, argparse, pickle
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from collections import Counter

from kyfan.complex import SignedComplex, leq, neg, rank
from kyfan.realize_explicit import explicit_Leq


class GadgetSampler:
    def __init__(self, m, seed=0):
        self.m, self.n = m, m - 1
        n = self.n
        cxe = SignedComplex(n)
        self.cxe = cxe
        free = cxe.free
        N = self.N = len(free)
        self.rng = np.random.default_rng(seed)
        self.labels = np.array([s * k for k in range(1, n + 1) for s in (1, -1)])
        self.L0 = np.array(explicit_Leq(m), dtype=np.int64)
        self.L = self.L0.copy()
        # neighbours among free reps, with the sign product: constraint  si*L[i] != -sj*L[j]  <=>  L[i] != -(si*sj)*L[j]
        nb = [[] for _ in range(N)]
        for x, y in cxe.edges:
            (i, si), (j, sj) = cxe.rep(x), cxe.rep(y)
            if i != j:
                nb[i].append((j, si * sj))
                nb[j].append((i, si * sj))
        self.nbr_idx = [np.array([j for j, _ in l], dtype=np.int64) for l in nb]
        self.nbr_sgn = [np.array([-s for _, s in l], dtype=np.int64) for l in nb]   # forbidden = sgn * L[j]
        # static forbidden labels from the gadget domain conditions
        w = {r: tuple(-1 if k < r else 0 for k in range(n)) for r in range(1, n + 1)}
        sigma = set(w.values())
        T = {1: {n - 1}}
        for r in range(2, n + 1):
            T[r] = ({n - r} if r < n else set()) | {-k for k in range(n - r + 1, n)}
        self.static_ok = np.ones((N, len(self.labels)), dtype=bool)
        fixed = set()
        lab_pos = {int(l): k for k, l in enumerate(self.labels)}
        for i, v in enumerate(free):
            for vv, s in ((v, 1), (neg(v), -1)):             # lambda(vv) = s * L[i]
                if vv in sigma:
                    fixed.add(i)
                    continue
                forb = {s * (-n)}                             # pole (+n) is adjacent to the lift of vv: lambda(vv) != -n
                for r in range(1, n + 1):
                    if leq(vv, w[r]) or leq(w[r], vv):        # vv is a fixed neighbour of leaf r
                        forb |= {s * (-l) for l in T[r] | {-n}}
                for l in forb:
                    self.static_ok[i, lab_pos[l]] = False
        self.fixed = fixed
        self.movable = np.array([i for i in range(N) if i not in fixed], dtype=np.int64)
        assert all(self.L0[i] in self.allowed(i) for i in self.movable), "explicit gadget violates its own constraints"
        # flag-density bookkeeping: vertices by rank, and for each vertex the vertices one rank below (zero one coordinate)
        self.byrank = {}
        for v in cxe.verts:
            self.byrank.setdefault(rank(v), []).append(v)
        self.below = {v: [v[:k] + (0,) + v[k + 1:] for k in range(n) if v[k]] for v in cxe.verts}
        self.nflags = int(np.prod(range(1, n + 1))) * 2 ** n
        self.frozen = 0
        self.moves = 0

    def allowed(self, i):
        dyn = self.nbr_sgn[i] * self.L[self.nbr_idx[i]]
        ok = self.static_ok[i] & ~np.isin(self.labels, dyn)
        return self.labels[ok]

    def sweep(self, k=1):
        rng, mv = self.rng, self.movable
        picks = rng.choice(mv, size=k * self.N)
        for i in picks:
            al = self.allowed(i)
            self.moves += 1
            if len(al) <= 1:
                self.frozen += 1
            self.L[i] = al[rng.integers(len(al))]

    # -- statistics of the current state
    def label_of(self, v):
        return self.cxe.label(v, self.L)

    def flag_density(self):
        """max over magnitudes k of  #maximal chains of the equator entirely of magnitude k / #maximal chains."""
        best = 0.0
        n = self.n
        for k in range(1, n + 1):
            cnt = {v: 1 if abs(self.label_of(v)) == k else 0 for v in self.byrank[1]}
            for r in range(2, n + 1):
                cnt = {v: (sum(cnt[u] for u in self.below[v]) if abs(self.label_of(v)) == k else 0)
                       for v in self.byrank[r]}
            best = max(best, sum(cnt.values()) / self.nflags)
        return best

    def stats(self):
        mags = Counter(np.abs(self.L).tolist())
        labs = Counter(self.L.tolist())
        return dict(top_mag=max(mags.values()) / self.N,
                    top_label=max(labs.values()) / self.N,
                    n_mags=len(mags),
                    flag_density=self.flag_density(),
                    hamming=float((self.L != self.L0).mean()))


# ---------------------------------------------------------------- diagnostics
def autocorr_time(x):
    """Integrated autocorrelation time with Sokal's automatic window (c = 5)."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    T = len(x)
    if T < 8 or x.var() == 0:
        return 1.0
    f = np.fft.rfft(x, n=2 * T)
    acf = np.fft.irfft(f * np.conj(f))[:T] / (x.var() * np.arange(T, 0, -1))
    tau = 1.0
    for M in range(1, T):
        tau = 1 + 2 * acf[1:M + 1].sum()
        if M >= 5 * tau:
            break
    return max(tau, 1.0)


def split_rhat(chains):
    """Gelman-Rubin split-R-hat on a list of equal-length 1-d traces."""
    xs = [np.asarray(c, dtype=float) for c in chains]
    T = min(len(c) for c in xs) // 2
    halves = [c[:T] for c in xs] + [c[T:2 * T] for c in xs]
    X = np.array(halves)
    Mn, Tn = X.shape
    W = X.var(axis=1, ddof=1).mean()
    B = Tn * X.mean(axis=1).var(ddof=1)
    if W == 0:
        return 1.0
    return float(np.sqrt(((Tn - 1) / Tn * W + B / Tn) / W))


def run_chains(m, chains, sweeps, burn_frac=0.25, record_every=1, seed0=0, log=True):
    """Run `chains` independent Glauber chains from the explicit gadget; record the scalar statistics every
    `record_every` sweeps after burn-in.  Returns (traces: stat -> list of per-chain arrays, samplers, times)."""
    traces = {}
    samplers = []
    t0 = time.time()
    burn = int(sweeps * burn_frac)
    for c in range(chains):
        S = GadgetSampler(m, seed=seed0 + 1000 * (c + 1))
        rec = {k: [] for k in ("top_mag", "top_label", "n_mags", "flag_density", "hamming")}
        for s in range(sweeps):
            S.sweep()
            if s >= burn and (s - burn) % record_every == 0:
                for k, v in S.stats().items():
                    rec[k].append(v)
        for k in rec:
            traces.setdefault(k, []).append(np.array(rec[k]))
        samplers.append(S)
        if log:
            print(f"    chain {c}: {sweeps} sweeps, frozen-move rate {S.frozen / S.moves:.3f}, "
                  f"final hamming {rec['hamming'][-1]:.2f}  [{time.time() - t0:.0f}s]", flush=True)
    return traces, samplers


def summarize(traces, record_every):
    out = {}
    for k, ch in traces.items():
        allx = np.concatenate(ch)
        taus = [autocorr_time(c) for c in ch]
        tau = float(np.mean(taus)) * record_every           # in sweeps
        ess = sum(len(c) / t for c, t in zip(ch, taus))
        out[k] = dict(mean=float(allx.mean()), sd=float(allx.std()), tau_sweeps=tau, ess=float(ess),
                      rhat=split_rhat(ch), chain_means=[float(c.mean()) for c in ch])
    return out


# ---------------------------------------------------------------- match probabilities
def sample_sets(S, t, kind, rng):
    """A random set of t movable free-vertex indices: 'random' (uniform), 'chain' (from a random maximal chain of
    the equator), 'ball' (from the down-set of a random top vertex)."""
    cxe, n, mv = S.cxe, S.n, set(S.movable.tolist())
    if kind == "random":
        return rng.choice(S.movable, size=min(t, len(S.movable)), replace=False)
    if kind == "chain":
        x = [0] * n
        x[rng.integers(n)] = int(rng.choice([1, -1]))
        ch = [tuple(x)]
        while sum(1 for c in x if c) < n:
            zeros = [k for k in range(n) if x[k] == 0]
            x[zeros[rng.integers(len(zeros))]] = int(rng.choice([1, -1]))
            ch.append(tuple(x))
        idx = [cxe.rep(c)[0] for c in ch]
        idx = [i for i in idx if i in mv]
        return np.array(idx[:t])
    if kind == "ball":
        top = tuple(int(s) for s in rng.choice([1, -1], size=n))
        ds = [v for v in cxe.verts if leq(v, top)]
        rng.shuffle(ds)
        idx, seen = [], set()
        for v in ds:
            i = cxe.rep(v)[0]
            if i in mv and i not in seen:
                seen.add(i)
                idx.append(i)
        return np.array(idx[:t])
    raise ValueError(kind)


def match_probabilities(m, nsamp, burn, thin, sizes=(1, 2, 3, 4, 6, 8, 12), trials=40, ntargets=200, seed=0, log=True):
    """Two independent chains A, B.  For each set type and size t: over `trials` random sets A and `ntargets`
    B-samples b, the fraction of A-samples a with a|_A = b|_A.  Returns dict kind -> (sizes, means, stderrs)."""
    t0 = time.time()
    SA, SB = GadgetSampler(m, seed=seed + 11), GadgetSampler(m, seed=seed + 22)
    SA.sweep(burn)
    SB.sweep(burn)
    A = np.empty((nsamp, SA.N), dtype=np.int64)
    B = np.empty((nsamp, SB.N), dtype=np.int64)
    for s in range(nsamp):
        SA.sweep(thin)
        SB.sweep(thin)
        A[s], B[s] = SA.L, SB.L
    rng = np.random.default_rng(seed + 33)
    out = {}
    for kind in ("random", "chain", "ball"):
        means, errs = [], []
        for t in sizes:
            ps = []
            for _ in range(trials):
                idx = sample_sets(SA, t, kind, rng)
                if len(idx) < t:
                    continue
                tg = B[rng.choice(nsamp, size=ntargets, replace=False)][:, idx]
                sub = A[:, idx]
                hits = (sub[None, :, :] == tg[:, None, :]).all(axis=2).mean(axis=1)   # per target
                ps.append(hits.mean())
            if ps:
                means.append(float(np.mean(ps)))
                errs.append(float(np.std(ps) / np.sqrt(len(ps))))
            else:
                means.append(float("nan"))
                errs.append(float("nan"))
        out[kind] = (list(sizes), means, errs)
    if log:
        print(f"  match probabilities m={m}: {nsamp} samples/chain, burn {burn}, thin {thin}  [{time.time() - t0:.0f}s]")
        for kind, (sz, mu, se) in out.items():
            valid = [(t, p) for t, p in zip(sz, mu) if p == p and p > 0]
            beta = float(np.exp(np.polyfit([t for t, _ in valid], [np.log(p) for _, p in valid], 1)[0])) if len(valid) > 1 else float("nan")
            print(f"    {kind:6s}: " + "  ".join(f"|A|={t}: {p:.4f}±{e:.4f}" for t, p, e in zip(sz, mu, se)) + f"   fitted beta/vertex {beta:.3f}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, nargs="+", default=[4, 5, 6, 7])
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--sweeps", type=int, default=2000)
    ap.add_argument("--every", type=int, default=2)
    ap.add_argument("--match", action="store_true")
    ap.add_argument("--nsamp", type=int, default=2000)
    ap.add_argument("--thin", type=int, default=3)
    ap.add_argument("--save", default=None)
    a = ap.parse_args()
    results = {}
    for m in a.m:
        S = GadgetSampler(m)
        e = S.stats()
        print(f"m={m} (equator S^{m-2}, {S.N} free vertices, {len(S.movable)} movable): explicit gadget: "
              f"top_mag {e['top_mag']:.3f}, top_label {e['top_label']:.3f}, flag_density {e['flag_density']:.4f} (=1/(2n)={1 / (2 * S.n):.4f})", flush=True)
        traces, _ = run_chains(m, a.chains, a.sweeps, record_every=a.every)
        summ = summarize(traces, a.every)
        for k, v in summ.items():
            print(f"    {k:13s}: mean {v['mean']:.4f} (sd {v['sd']:.4f})  tau {v['tau_sweeps']:.1f} sweeps  ESS {v['ess']:.0f}  "
                  f"R-hat {v['rhat']:.3f}  chain means {['%.4f' % x for x in v['chain_means']]}", flush=True)
        results[m] = dict(explicit=e, summary=summ)
        if a.match:
            results[m]["match"] = match_probabilities(m, a.nsamp, burn=max(200, a.sweeps // 4), thin=a.thin)
    if a.save:
        with open(a.save, "wb") as f:
            pickle.dump(results, f)


if __name__ == "__main__":
    main()
