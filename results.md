# Results log

All degrees are Sherali–Adams level in the one-hot quotient (see CLAUDE.md §1). All runs from `scripts/`.
Runtimes are from a modest sandbox; a laptop will be faster.

## Settled

| # | statement | field | degree | evidence | command | runtime |
|---|---|---|---|---|---|---|
| 1 | Tucker, S^2 (m=3, labels ±1,±2) | F_2 | = 3 | dual consistent at d=2, inconsistent at d=3; explicit certificate (458 monomials) verified on 20,000 fresh labelings | `python3 nsdeg.py`, `python3 primal.py` | 10 s, 60 s |
| 2 | Tucker, S^2 | Q | = 5 | dual consistent at d=3,4 (F_p, p=1000003) and float-rank; full-group symmetric certificate at d=5 verified on 3,000 fresh labelings (max residual 2e-14) | `primal.py` (d=3,4); d=5 via symmetric orbit system (657→5482 orbits) — see note | 10 s; 5 min |
| 3 | Tucker, S^3 (m=4, labels ±1..±3) | F_2 | = 4 | d=3: Z_3×Z_3-averaged unrestricted system inconsistent (rank 27,584 of 33,312 orbit-unknowns); d=4: tower | `python3 z9a.py && python3 z9b.py` | 1 min + 15 min |
| 4 | Ky Fan, S^2, labels ±1..±3 | F_2 | = 3 | explicit certificate (1,286 monomials) verified on 20,000 fresh | `primal.py` | 30 s |
| 5 | Ky Fan, S^n | F_2 | = n+1 | tower theorem; local lemma exhaustive for n≤3, k≤5; identity checked on 400 random labelings each for m=3,4,5 | `python3 tower.py` | 2 min |

Note on #2, d=5: the script that produced it (orbit enumeration + float rank + fresh verification) was run in stages; reconstruct from `primal.py`'s `primal_symmetric_Fp` with d=5 and ~1.6×#orbits sample rows, then lstsq + fresh-sample check. Expected: 924,480 violating unknowns in 5,482 orbits; rank 1,976; residual ~1e-14.

## Structural facts (degree-2 pseudo-solution on S^2)

| fact | command |
|---|---|
| solution space: 1,104 unknowns, rank 507, 597 free parameters | `python3 pseudo2b.py` |
| no full-group-invariant solution (20 orbits) | `python3 pseudo2.py` |
| no label-group-invariant, no hemisphere-stabilizer-invariant, no distinct-magnitudes-only solution | `python3 pseudo2b.py` |
| Z_3-invariant solution exists (368 orbits, 201 free) | `python3 pseudo2b.py` |
| necessary pair types: distinct-magnitude edges at ranks (1,2),(1,3),(2,3); equal labels on two unrelated rank-1 vertices | `python3 pseudo2c.py` |

## Negative / methodological

- Unsubdivided octahedron: A_+ ≡ 1 on all 192 valid labelings (degenerate). `build.py`.
- Fully-symmetric F_2 certificates do not exist on S^2 at d=3 or d=4, although an asymmetric d=3 one does. `s3.py` (set M=3). Never conclude from symmetric F_2 searches.
- Full-group symmetric search on S^3 at d=3: none over F_2 (uninformative), none over Q (real: residual 7e-2). `s3.py`.
- Pullback of an equatorial pseudo-solution to S^n loses one degree (hand argument, CLAUDE.md §2).

## Open

- Tucker F_2 lower bound n+1 for n ≥ 4.
- Structure of the degree-3 pseudo-solution on S^3 (Task 2).
- Ball-version degrees (Task 3).
