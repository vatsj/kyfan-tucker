import sys, time
from run_s2 import build
from lowdeg import estimate
N, q, ef, sf = build(4, 4)
print(f"S^3, 4 magnitudes: {N} variables, {q} labels, {len(ef)} edge factors, {len(sf)} simplex factors; truth: ratio = -1 (Ky Fan)", flush=True)
for d in map(int, sys.argv[1:]):
    t0 = time.time(); Zd = estimate(ef, d); Nd = estimate(ef + sf, d)
    print(f"  d={d}: P[valid] est {Zd:.4e};  numerator est {Nd:+.4e};  ratio {Nd/Zd:+.4f}   [{time.time()-t0:.0f}s]", flush=True)
