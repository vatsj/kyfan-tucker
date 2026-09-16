# Results log

All degrees are Sherali–Adams level in the one-hot quotient (see CLAUDE.md §1).
Every number below is an assertion in `tests/` (the tests *are* this table); `./reproduce.sh` runs the fast subset (~2 min),
`./reproduce.sh all` everything (~15 min). Runtimes: M-series laptop, 10 cores, 16 GB, Python 3.13, numpy 2.5 (2026-09-16).
The original `scripts/` (exec-chained, m=3 hard-coded) were reviewed in `review.md`, reproduced number-for-number, and then
replaced by the package `kyfan/` (git history has them).

## Settled

| # | statement | field | degree | evidence | test | runtime |
|---|---|---|---|---|---|---|
| 1 | Tucker, S^2 (m=3, labels ±1,±2) | F_2 | = 3 | dual: 1,104 unknowns, consistent at d=2 (rank 507, 597 free); 12,688 unknowns, inconsistent at d=3 (rank 8,400). Explicit certificate (5,616 violating size-3 unknowns; ~458 monomials) verified exactly: canonical one-hot basis expansion, and exhaustively on all 4^13 labelings | `test_nsdeg.py`; `test_primal.py::test_tucker_s2_F2_degree3_certificate_exact`, `…_exhaustive` (slow) | 2 s; 10 s, 30 s |
| 2 | Tucker, S^2 | Q | = 5 | dual consistent at d=3 over F_p (p=1000003, rank 8,401); full-group orbit primal (valid over Q/F_p, |G|=384): inconsistent at d=3 (55 orbits) and d=4 (657 orbits) — one-way sound; at d=5: 924,480 violating unknowns in 5,482 orbits, float rank 1,976, relative residual 3e-15, max error 1e-14 on fresh labelings, and an exact F_p certificate | `test_nsdeg.py`; `test_primal.py::test_tucker_s2_Fp_symmetric_degree3_4_none`; `…::test_tucker_s2_Q_degree5_symmetric_certificate` (veryslow) | 5 s; 5 s; 5 min |
| 3 | Tucker, S^3 (m=4, labels ±1..±3) | F_2 | = 4 | d=3, via the dual: the Z_3×Z_3-invariant degree-3 SA-dual (1,834,800 non-violating size-3 partial labelings in 204,048 orbits; 111,233 rows) is consistent, rank 94,064 — an explicit invariant pseudo-solution (support 8,036 orbits), checked against all 987,265 unreduced consistency equations. Averaging is valid over F_2 since |Z_3×Z_3| = 9 is odd. The unrestricted dual (987,457 × 1,834,800) is also consistent, rank 827,876. Cross-check via the sampled primal (z9): inconsistent, 33,312 orbit-unknowns, rank 27,584. d=4: tower | `test_sparse.py::test_s3_degree3_pseudo_solution_via_dual` (30 s); `test_s3.py::test_s3_degree3_F2_no_certificate_z9` (veryslow, 7 min) | 30 s; 7 min |
| 4 | Ky Fan, S^2, labels ±1..±3 | F_2 | = 3 | ≤: explicit certificate (13,176 violating unknowns; ~1,286 monomials) verified exactly (canonical basis: Σc_α1_α ≡ A_+ + 1). ≥: row 5 | `test_primal.py::test_kyfan_s2_F2_degree3_certificate_exact` | 60 s |
| 5 | Ky Fan, S^n, k ≥ n+1 | F_2 | = n+1 | ≤: tower theorem — local lemma g_n = 0 on non-complementary tuples, exhaustive n≤3, k≤5; telescoped identity on random labelings m=3,4,5 (29 / 221 / 2,141 local-lemma instances). ≥: separating functional (CLAUDE.md §2): E(A_+)=1, E(1)=0, E(1_α)=0 on size-n α (exhaustive m=3, sampled m=4,5); the negation-closed box gives E(A_+)=0 | `test_tower.py` | 15 s |
| 6 | Tucker, S^n | F_2 | ≤ n+1; = n+1 for n=2,3 | upper: tower (row 5, A_+ ≡ 0 when k ≤ n); lower: rows 1, 3. Open for n ≥ 4 — the Ky Fan functional kills constants and cannot be reused (a Tucker lower bound needs E[∅]=1) | — | — |

## Task 1 — sparse GF(2) solver (`gf2solve/`, Rust; `linalg.gf2_sparse`, default in `dual.solve`)

Sparse Gaussian elimination with Markowitz-style pivoting (min-count column, shortest row in it; rows as sorted index
lists, lazy column index, min-heap with stale entries). Validated against the dense solver on the m=3 systems
(507; 8400, inconsistent), on random systems, and by checking returned solutions / null-space vectors row by row
(`test_sparse.py`). Timings (this laptop):

| system | rows × cols | nnz | rank | consistent | time |
|---|---|---|---|---|---|
| S^2 d=3 dual | 11,041 × 12,688 | 65,640 | 8,400 | no | 0.05 s (dense: 0.95 s) |
| S^3 d=3 dual, Z_3×Z_3-invariant | 111,233 × 204,048 | 1.13 M | 94,064 | yes | 0.8 s |
| S^3 d=3 dual, unrestricted | 987,457 × 1,834,800 | 10.06 M | 827,876 | yes | 21 s (longest row 750) |

Not attempted: the z9 *primal* (38k × 33k, ~25% dense) is the wrong shape for a sparse solver; its statement is now
covered by the dual above. Sage/M4RI and block Wiedemann were unnecessary.

## Structural facts (degree-2 pseudo-solution on S^2) — `test_pseudo2.py`, 5 s

| fact |
|---|
| solution space: 1,104 unknowns, rank 507, 597 free parameters |
| no full-group-invariant solution (20 orbits) — a statement about invariant solutions only, |G| even |
| no label-group-invariant (198 orbits), no hemisphere-stabilizer-invariant (44), no distinct-magnitudes-only (624 unknowns) solution; nor with both restrictions (78; 19) |
| coordinate-Z_3-invariant solution exists (368 orbits, 201 free) |
| knockout by full-group orbit type: exactly 4 of 20 types are necessary — distinct-magnitude edges at rank pairs (1,2), (1,3), (2,3); equal labels on two unrelated rank-1 vertices |

## Negative / methodological

- Unsubdivided octahedron: A_+ ≡ 1 on all 192 valid labelings with labels ±1..±4 (degenerate). `test_complex.py::test_octahedron_degenerate`.
- Fully-symmetric F_2 certificates do not exist on S^2 at d=3 (55 orbits) or d=4 (657), although an asymmetric d=3 one does. Never conclude from symmetric F_2 searches. `test_primal.py::test_tucker_s2_F2_full_group_symmetric_none_d3_d4` (slow, 3 s).
- Full-group symmetric search on S^3 at d=3 (299,280 violating in 232 orbits): none over F_2 (uninformative), none over R (residual 7e-2 ≫ 1e-3). `test_s3.py::test_s3_full_group_symmetric_degree3_none` (slow, 20 s).
- Pullback of an equatorial pseudo-solution to S^n loses one degree (hand argument, CLAUDE.md §2).
- Consistency equations at level d−1 suffice: building only those gives the same ranks (507; 8400/8401) as building all levels. The package builds only level d−1.

## Open

- Tucker F_2 lower bound n+1 for n ≥ 4 (row 6).
- Structure of the degree-3 pseudo-solution on S^3 (Task 2).
- Ball-version degrees (Task 3).
