# docs/ — index: statement → proof text → mechanical check → tier → runtime

Numbering follows the paper. Tiers: **fast** = `make test-fast` (`-m "not slow and not veryslow"`, ~1.5 min total),
**slow** = `make test-slow` (~6 min total), **veryslow** = `make test-veryslow` (~20 min total). Runtimes are wall-clock on an
Apple-silicon laptop (Python 3.13.7, numpy 2.5.3, Rust 1.91); the sparse GF(2) solver `gf2solve/` is used wherever "sparse" is
noted. Setup and the paper-cited commands: top-level `README.md`.

## Proofs (all unconditional unless marked)

| statement | proof text | what is checked mechanically | test | tier | runtime |
|---|---|---|---|---|---|
| Tower certificate: A_+ + 1 = Σ_j Σ_{σ ∈ H^{(j)}} g_{j−1}(λ|_σ), degree m = n+1 (Tucker and Ky Fan upper bounds) | `tower.md` | local lemma g_n = 0 on every non-complementary (n+1)-tuple, exhaustively for n ≤ 3, k ≤ 5; the telescoped identity on 200 random labelings for m = 3, 4, 5 (29 / 221 / 2,141 local-lemma instances) | `tests/test_tower.py::test_local_lemma_exhaustive`, `::test_telescoped_identity` | fast | 0.1 s + 0.3 / 0.8 / 5.6 s |
| Ky Fan over F_2 has degree exactly n+1 (separating functional) | `kyfan_lower_bound.md` | E(A_+) = 1, E(1) = 0, E(1_α) = 0 on 1,000 sampled size-n α, m = 3, 4, 5; the negation-closed box gives E(A_+) = 0 | `tests/test_tower.py::test_kyfan_lower_functional`, `::test_kyfan_lower_needs_collision` | fast | < 5 s |
| Statement (i): the tree gadget has F_2 degree n+1 for every binary conflict tree — explicit pseudo-solution (block rule) + extension lemma (★) | `tree_gadget.md` | Lemma (★) exhaustively for every unordered tree shape with n ≤ 6 (12 shapes) and all 11 shapes with n = 7; the closed-form E equals the solver's unique pseudo-solution, support 3^n, same shapes (n = 7: 36,409 unknowns each, sparse) | `tests/test_extension_lemma.py`, `tests/test_rule_vs_solver.py` | fast (n ≤ 6); slow (n = 7) | 1.5 s + 1.5 s; 24 s + 11 s (all 11 shapes) |
| Restriction lemma | `tree_gadget.md` §"The restriction lemma" | (used through statement (ii); the searched gadgets `tests/test_gadget.py` and the direct S^2/S^3 duals are independent confirmations) | `tests/test_gadget.py`, `tests/test_nsdeg.py`, `tests/test_sparse.py` | fast / slow / veryslow | see HISTORY.md |
| Statement (ii): the explicit labeling λ(u) = s·(n−q+1) is valid and realizes the reverse caterpillar; hence Tucker's F_2 degree on S^n = n+1 for all n | `realization.md` | equatorial validity + exact residual domains + consistent, unique degree-n residual dual, m = 3..8 (m = 8: 95,901 unknowns, sparse); validity + exact domains alone for m = 9, 10; full-sphere `verify` m = 4, 5, 6 | `tests/test_realize_explicit.py` | fast (m ≤ 8); slow (m = 9); veryslow (m = 10) | 4 s total for m ≤ 8; 21 s; 3.5 min |
| Theorem 1 (size ≥ #flags for degree-(n+1) certificates), Lemma A, Lemma B, Theorem 1′ (per-flag parities), Theorem 2, Theorem 3 | `size_lower_bound.md` §1–§3 | see "Paper-cited computations" below for the numbers | | | |
| Theorem 4 (size ≥ m!2^m/K^m for all degrees) — conditional on the Spread Lemma (conjecture) | `size_lower_bound.md` §4 | empirical support only: `spread_evidence.md` | `tests/test_spread_mcmc.py` | fast | 4 s |

## Paper-cited computations

| claim | script and stored output | test (re-verifies the stored artifact) | tier | runtime |
|---|---|---|---|---|
| minimum certificate size at degree ≤ 3 on S^2 = **304** (16 size-2 + 288 size-3 monomials) | `analysis/isd.py --d 3 …` → `analysis/isd_s2_d3_cert.txt` (2.5 min) | `tests/test_paper_numbers.py::test_min_certificate_304_exhaustive` (identity on all 4^13 labelings; `_random`: 2,000 labelings, fast); lower bound: `tests/test_min_size.py::test_s2_degree3_min_size_304` (CP-SAT optimality on the sampled system) | slow (fast variant); slow | 26 s (0.5 s); 70–110 s |
| Theorem 1′: **t = 8 on S^2** (all 16 gadget domain types at the pole flag) | `analysis/flag_hitting_bound.py --m 3 --family orbit trees all` → `flag_hitting_bound_m3.out`; independent coordinates: `analysis/flag_hitting_bound_s2.py` → `flag_hitting_bound_s2.out` | `tests/test_paper_numbers.py::test_t8_on_s2_branch_and_bound` (exhaustive 3,568 restrictions → 16 types, branch and bound); `tests/test_flag_hitting.py::test_s2_all_gadget_types_t8` (CP-SAT cross-check, families) | fast | 2 s; 4 s |
| Theorem 1′: **t ≥ 12 on S^3** (192 gadget types from 32 realizable chain permutations × label orbits) | `analysis/flag_hitting_bound.py --m 4 --family orbit trees perms --save-realizations …` → `flag_hitting_bound_m4.out`, `flag_hitting_bound_m4_realizations.txt` | `tests/test_paper_numbers.py::test_t12_on_s3_from_stored_realizations` (re-verifies the 32 stored labelings without CP-SAT, expands to 192, solves the parity system: OPTIMAL 12, proven bound 12); `tests/test_flag_hitting.py::test_s3_t_from_gadget_families` (re-searches the realizations with CP-SAT) | fast | 1 s; 3 s |
| degree-≤ 4 Z_3-invariant search on S^2: best certificate 312 monomials, profile **{2: 24, 3: 288}, no size-4 term** (CP-SAT FEASIBLE 104 orbits, proven bound 24 — not decisive) | `analysis/isd.py --d 4 --z3 …` → `isd_s2_d4_z3.out`, `isd_s2_d4_z3_cert.txt` (42 min; not rerun) | `tests/test_paper_numbers.py::test_degree4_z3_certificate_profile_exhaustive` (profile + identity on all 4^13 labelings; `_random` fast) | slow (fast variant) | 26 s (0.5 s) |
| Spread Lemma evidence, m = 4..8 (flag densities, match probabilities, fitted β) | `analysis/spread_mcmc.py --m 4 5 6 7 8 --chains 4 --sweeps 3000 --match --nsamp 2000` → `spread_mcmc.out` (10 min; reproduces byte-for-byte modulo `[Ns]` stamps) | `tests/test_spread_mcmc.py` (m = 4, 5 short chains; qualitative assertions) | fast | 4 s |

## Other documents

- `spread_evidence.md` — what `spread_mcmc.py` samples, every diagnostic, every caveat, and the stored numbers.
- `AGENT_BRIEF.md` — the original working brief (definitions §1, settled table §2, file map §3), kept as a historical document.
- `../HISTORY.md` — the development log: every settled number with the command that produced it, the negative results and
  the exploratory drivers (`../analysis/legacy/`).
