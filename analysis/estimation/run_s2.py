import sys, itertools, time, random
import numpy as np
from setup import instance, posalt
from lowdeg import estimate

def build(m, k):
    cx, N, labels, edges, flags = instance(m, k)
    q = len(labels)
    ef = []
    for (a, b), s in edges.items():
        t = np.array([[0.0 if labels[x] == -s * labels[y] else 1.0 for y in range(q)] for x in range(q)])
        ef.append(((a, b), t))
    sf = []
    for f in flags:
        order = sorted(range(m), key=lambda i: f[i][0]); axes = tuple(f[i][0] for i in order); sg = [f[i][1] for i in order]
        t = np.ones((q,) * m)
        for idx in itertools.product(range(q), repeat=m):
            ls = [sg[i] * labels[idx[i]] for i in range(m)]
            if posalt(ls) or posalt([-l for l in ls]): t[idx] = -1.0
        sf.append((axes, t))
    return N, q, ef, sf

if __name__ == "__main__":
    # --- sanity: S^1 with 3 magnitudes, d = N must be exact
    N, q, ef, sf = build(2, 3)
    tot_Z = tot_N = 0.0
    for idx in itertools.product(range(q), repeat=N):
        v = 1.0
        for ax, t in ef: v *= t[tuple(idx[a] for a in ax)]
        s = 1.0
        for ax, t in sf: s *= t[tuple(idx[a] for a in ax)]
        tot_Z += v; tot_N += v * s
    print(f"S^1 sanity (N={N} vars): brute Z={tot_Z/q**N:.6f} N={tot_N/q**N:.6f};  d={N}: Z={estimate(ef,N):.6f} N={estimate(ef+sf,N):.6f}")
    # --- S^2, k = 3
    N, q, ef, sf = build(3, 3)
    Zexact = 22193664 / q**N
    print(f"S^2: {N} variables, {q} labels, {len(ef)} edge factors, {len(sf)} simplex factors;  exact P[valid] = {Zexact:.6e}, exact E[(-1)^A; valid] = {-Zexact:.6e}")
    rng = random.Random(0); mixed = ef + sf; rng.shuffle(mixed)
    for d in map(int, sys.argv[1:]):
        t0 = time.time()
        Zd = estimate(ef, d); Nd = estimate(ef + sf, d); Ns = estimate(sf + ef, d); Nm = estimate(mixed, d)
        print(f"  d={d}: P[valid] est {Zd:.4e} (x{Zd/Zexact:.2f});  numerator est, edges-then-simplices {Nd:+.4e} (ratio to Zest {Nd/Zd:+.3f}), "
              f"simplices-first {Ns:+.4e}, shuffled {Nm:+.4e};  truth ratio -1   [{time.time()-t0:.0f}s]")
        sys.stdout.flush()
