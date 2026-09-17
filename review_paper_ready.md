# Paper-readiness review (Phase 0)

Read: `README.md`, `results.md`, `CLAUDE.md`, `TREE_GADGET_THEOREM.md`, `REALIZATION_THEOREM.md`, `SIZE_LOWER_BOUND.md`,
`review.md`, `analysis/README.md`, every file in `tests/`, `analysis/spread_mcmc.py`, the two stored certificates.
Nothing was recomputed or edited for this review.

## A. Numerical claims and where each is produced

### The three computations the paper cites

| claim | stated in | produced by | tier / runtime |
|---|---|---|---|
| min certificate size at degree ≤ 3 on S^2 = **304** (16 size-2 + 288 size-3) | SIZE_LOWER_BOUND §1′, results.md "Certificate size", analysis/README | `analysis/isd.py --d 3 …` (2.5 min); stored `analysis/isd_s2_d3_cert.txt`; `tests/test_min_size.py::test_s2_degree3_min_size_304` (exhaustive check of the stored certificate + CP-SAT optimality on the sampled system) | slow, ~70 s |
| degree-4 Z_3 search: best 312 monomials = 104 orbits, profile **{2: 24, 3: 288}** (no size-4 term), CP-SAT FEASIBLE 104 / proven bound 24 | results.md "Certificate size" | `analysis/isd.py --d 4 --z3 …` (42 min, NOT to be rerun); stored `analysis/isd_s2_d4_z3_cert.txt` (profile confirmed by counting lines: 24 + 288) | no test yet |
| **t = 8** on S^2 (16 gadget types, 3,568 restrictions; label orbit alone 4) | SIZE_LOWER_BOUND §1′, results.md, tests docstring | `analysis/flag_hitting_bound.py --m 3 --family orbit trees all`, `analysis/flag_hitting_bound_s2.py`; `tests/test_flag_hitting.py::test_s2_all_gadget_types_t8` (CP-SAT + branch and bound) | fast, ~4 s total |
| **t ≥ 12** on S^3 (label orbit 48 → 6; trees 96 → 6; realizable permutations 32 of 96, 192 gadgets → 12) | same | `analysis/flag_hitting_bound.py --m 4 --family orbit trees perms`; `tests/test_flag_hitting.py::test_s3_t_from_gadget_families` | fast |
| min certificate has 12 full-flag terms per class, all 16 parities odd | same | `tests/test_flag_hitting.py::test_s2_minimum_certificate_satisfies_parities` | fast |
| Spread-lemma MCMC, m = 4..8: top-magnitude fraction 0.461/0.371/0.287/0.226/0.183; flag density 0.094/0.024/0.0051/0.0009/0.0001; fitted β (random) 0.65/0.48/0.36/0.29/0.22; match-probability tables | results.md "Spread-lemma MCMC" | `analysis/spread_mcmc.py --m 4 5 6 7 8 --chains 4 --sweeps 3000 --match --nsamp 2000 --save …` (~10 min; seeds fixed: chain c uses 1000(c+1), match 11/22/33); output `analysis/spread_mcmc.out` | no test yet; `.out` lines carry wall-clock `[Ns]` stamps that will not reproduce byte-for-byte |

### The five proof checks

| claim | stated in | produced by | tier / runtime |
|---|---|---|---|
| Extension lemma (★) exhaustively, all tree shapes n ≤ 6; n = 7 one shape | TREE_GADGET_THEOREM, CLAUDE.md, results.md | `tests/test_extension_lemma.py` (8 shapes n = 3..6 fast; `test_extension_lemma_n7` slow, one shape) | fast + slow |
| block rule = solver's unique pseudo-solution, support 3^n, n ≤ 6 (+ n = 7 sparse, 36,409 unknowns) | same | `tests/test_rule_vs_solver.py` | fast + slow |
| explicit labeling valid, exact reverse-caterpillar domains, m ≤ 10; residual dual consistent + unique m ≤ 8 (m = 8: 95,901 unknowns) | REALIZATION_THEOREM, CLAUDE.md, results.md | `tests/test_realize_explicit.py` — **covers m = 3..8 only**; the m = 9, 10 validity/domain check was an inline run and is not a test (see B.6) | fast (m ≤ 8 with solve, ~3 s at m = 8) |
| tower: local lemma g_n = 0 on non-complementary tuples, n ≤ 3, k ≤ 5; telescoped identity m = 3, 4, 5 (200 random labelings; 29 / 221 / 2141 local-lemma instances) | CLAUDE.md §2, results.md row 5 | `tests/test_tower.py::test_local_lemma_exhaustive`, `::test_telescoped_identity` | fast |
| Ky Fan separating functional: E(A_+) = 1, E(1) = 0, m = 3, 4, 5; negation-closed box gives 0 | CLAUDE.md §2, results.md | `tests/test_tower.py::test_kyfan_lower_functional`, `::test_kyfan_lower_needs_collision` | fast |

### Other settled numbers (degree table; stay in HISTORY.md, not cited by the paper)

S^2 Tucker: F_2 degree 3 (1104 unknowns rank 507 / 12688 rank 8400 inconsistent) — `test_nsdeg.py`; Q degree 5 — `test_primal.py` (veryslow).
S^3 Tucker F_2 degree 4 (Z_3×Z_3 dual, rank 27,584 inconsistent) — `test_s3.py` (veryslow, 7 min) and `test_sparse.py`.
Chain gadgets m = 4..8 (support 3^{m−1}, unique dual) — `test_gadget.py`, `test_realize.py`.
Degree-2 pseudo-solution structure on S^2 (597 free parameters, necessary pair types) — `test_pseudo2.py`.
Octahedron degeneracy, symmetric-search negatives — `test_complex.py`, `test_primal.py`, `test_s3.py`.

## B. Places that state the same number differently

1. **Explicit S^2 certificate size.** results.md settled-table row 1 says "~458 monomials" (and row 4 "~1,286" for Ky Fan);
   the minimum is 304 (results.md "Certificate size" already footnotes this). 458 / 1,286 are sizes of the first sampled-primal
   certificates, not minima. Resolution: HISTORY.md keeps the row verbatim (it is a log); README/docs state only 304.
2. **Number of flags in Theorem 1.** SIZE_LOWER_BOUND §0/§1/§5 and results.md say size ≥ #flags = m!·2^m (48 on S^2), and the
   min-size table's "lower bound from theory" column says 48 (Thm 1) and 12 (Thm 2 = 48/C(4,3)). §1′'s closing note observes that
   the "distinct flags ⇒ distinct terms" argument only separates antipodal *classes*, m!·2^{m−1} (24 on S^2), and that
   #flags is recovered because t ≥ 2. So Theorem 1's *statement* is correct only via Theorem 1′; Theorem 2's bound as written
   (#flags / C(D, m)) inherits the same factor-2 gap and is not repaired by 1′ (it applies to any degree). **This is a mathematical
   content issue, so I am not touching it**: in `docs/size_lower_bound.md` I will carry the text over unchanged, including the
   existing note, and flag the Theorem 2 factor of 2 in the same note (as a remark, not a change of statement). Please confirm
   whether the paper states Theorem 1/2 with m!·2^m or m!·2^{m−1}.
3. **Spread-lemma evidence numbers.** SIZE_LOWER_BOUND §3 and §4 cite `gadget_spread_experiments.py` (not in the repo) with
   m = 4..7 top-class fractions 0.44/0.36/0.29/0.23, flag densities 0.104/0.027/0.005/0.001, and m = 6 match probabilities
   0.39/0.19/0.06/0.013/0.004/0.001 (random), chains 0.90/0.47/0.14/0.07, balls 0.37/0.17/0.09/0.03/0.009/0.003. The repo's
   `analysis/spread_mcmc.out` (the reproducible run) gives 0.461/0.371/0.287/0.226 (m = 8: 0.183), 0.094/0.024/0.0051/0.0009
   (m = 8: 0.0001), and m = 6 random 0.39/0.16/0.046/0.014/0.0033/0.0002. Same qualitative picture, different run (earlier,
   shorter, different seeds). §5 says "empirical support m ≤ 7" while the MCMC covers m ≤ 8. Resolution: `docs/size_lower_bound.md`
   §4 Evidence will point to `docs/spread_evidence.md` and quote only `spread_mcmc.out` numbers. Since SIZE_LOWER_BOUND.md
   is being *moved*, the old paragraph would disappear; I will keep it verbatim in HISTORY.md as an appendix ("earlier run,
   script not retained") so no number is lost.
4. **Runtimes of the test suite.** README.md: fast ~2 min / all ~15 min; results.md header: fast ~3 min / all ~30 min;
   reproduce.sh comments: ~2 min / ~10 min. Will be re-measured in Phase 3 and stated once (docs/README.md).
5. **Tower identity sample count.** review.md says the original script checked 400 random labelings; the test uses 200. Cosmetic
   (the identity is a theorem); docs will say 200 = what the test does.
6. **Realization "m ≤ 10".** REALIZATION_THEOREM, CLAUDE.md and results.md all claim validity + exact domains checked for
   m = 3..10, but `tests/test_realize_explicit.py` parametrizes m = 3..8 (with solve). The m = 9, 10 check was run inline and
   never made a test. Phase 2 item 3 asks for exactly this test (validity + exact domains m ≤ 10, no solve); adding it does not
   change any number.
7. **Extension lemma "n ≤ 7".** Exhaustive for all 8 tree shapes at n ≤ 6 (that is every shape — 1, 2, 3, 2 shapes for n = 3..6);
   n = 7 is a single shape (`TREE_N7`), slow tier. The paper's "n ≤ 7" should read "all shapes n ≤ 6, one shape n = 7" —
   the docs will say so explicitly.
8. **Theorem 1′ S^2 bound "24 · 8 = 192"** appears consistently (SIZE_LOWER_BOUND, results.md, test docstring, .out). Fine.
   Min-size table's "lower bound from theory 48 (Thm 1)" is superseded by 192 (Thm 1′) — HISTORY keeps it; docs use 192.

## C. Session / collaborator / task-number / date references (to be removed or neutralized in Phase 1)

- `CLAUDE.md` — is a session brief throughout ("You are continuing a research computation", Phase 0, Task 1–4, "DONE 2026-09-16",
  "Working style"). → moved to `docs/AGENT_BRIEF.md`; top-level `CLAUDE.md` becomes a pointer.
- `README.md` — "the full brief … Read first" pointing at CLAUDE.md; runtime claims. → rewritten.
- `results.md` — "Status (2026-09-16) — what has happened since CLAUDE.md was written"; "the author supplied the
  separating-functional proof"; "collaborator's write-up"; "Task 1/2/3" headings; "Certificate size (2026-09-17)";
  "(NEW, settles the conjecture at n=4)" style remarks. → renamed HISTORY.md, content unchanged, 2-line header.
- `review.md` — Phase-0 script review, dated, refers to deleted `scripts/`. → moved to `analysis/legacy/review_scripts.md`
  (it documents commit 5a4dd50's scripts; nothing else references it).
- `SIZE_LOWER_BOUND.md` — "Computation (2026-09-17, …)"; §3/§4 cite the non-repo `gadget_spread_experiments.py`. → docs/size_lower_bound.md.
- `TREE_GADGET_THEOREM.md` — "Settles statement (i) of `results.md` ('What a proof for all n now needs')"; "now **proved**";
  "The 'never backtracks' question is retired". → docs/tree_gadget.md, narrative tense removed.
- `REALIZATION_THEOREM.md` — clean; only file references change.
- `analysis/README.md` — "Tasks 2, 3", "collaborator's original S^2 computation", "results.md points at". → reindexed.
- Test docstrings: `test_gadget.py` "(NEW, settles the conjecture at n=4)", `test_complex.py`/`test_nsdeg.py`/… "build.py:",
  "nsdeg.py:", "results.md #1" references to deleted scripts. → docstrings reworded to point at docs/ and HISTORY.md
  (no assertion changes).
- No file says "the other Claude". No `scripts/` directory exists (deleted at commit c2d8578 after `review.md`; originals at 5a4dd50),
  so Phase 1 item 4 reduces to moving `review.md`.

## D. Environment facts for Phase 3

Python 3.13.7; numpy 2.5.3; pytest 9.1.1; ortools 9.15.6755; cargo/rustc 1.91.1 (gf2solve builds on first use).
No `LICENSE`, `CITATION.cff`, `Makefile`, or `.github/` exist yet.

## E. Things I will not do without confirmation

- Change the m!·2^m ↔ m!·2^{m−1} statement of Theorems 1/2 (B.2).
- Rerun `isd.py --d 4` (42 min) — the degree-4 profile test will read the stored certificate and re-verify it exhaustively (18 s).
- Regenerate `spread_mcmc.out` beyond confirming the fixed-seed run reproduces its numbers (the m = 8 stage alone is ~8 min;
  I will run the full command once in the background and diff modulo the `[Ns]` timing stamps).
