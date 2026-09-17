# kyfan-tucker

Companion repository for the paper on the Nullstellensatz / Sherali–Adams degree of Tucker's lemma and Ky Fan's lemma on the
subdivided cross-polytope over F_2. Everything the paper cites from the repository — proofs, mechanical checks, and the exact
small numbers — is indexed in `docs/README.md`.

## What is proved here

Degrees are Sherali–Adams levels in the one-hot quotient (a certificate of degree d uses partial labelings of ≤ d vertices).
On the sphere S^n (m = n+1 coordinates):

- **Tucker over F_2 has degree exactly n+1** for every n ≥ 2. Upper bound: the tower certificate (`docs/tower.md`). Lower bound:
  the restriction lemma, the tree gadget with its explicit pseudo-solution (statement (i), `docs/tree_gadget.md`), and a closed-form
  equatorial labeling realizing the reverse caterpillar (statement (ii), `docs/realization.md`).
- **Ky Fan over F_2 has degree exactly n+1** for k ≥ n+1 labels: the tower and a separating functional (`docs/kyfan_lower_bound.md`).
- **Certificate size** (`docs/size_lower_bound.md`): Theorems 1, 1′, 2, 3 (unconditional), Theorem 4 conditional on the
  Spread Lemma, which is a conjecture with empirical support (`docs/spread_evidence.md`).

## The three paper-cited computations

Run from the repository root after the install below (`PYTHONPATH=.` is set by `make`; set it yourself for direct calls).

1. **Minimum certificate size on S^2 at degree ≤ 3 = 304.** Stored certificate `analysis/isd_s2_d3_cert.txt`, produced by
   `PYTHONPATH=. .venv/bin/python analysis/isd.py --d 3 --seconds 60 --cpsat 600 --save analysis/isd_s2_d3.pkl` (2.5 min).
   Re-verified by `tests/test_paper_numbers.py` (identity on all 4^13 labelings) and `tests/test_min_size.py` (CP-SAT optimality).
2. **Per-flag parities, Theorem 1′: t = 8 on S^2, t ≥ 12 on S^3.**
   `PYTHONPATH=. .venv/bin/python analysis/flag_hitting_bound.py --m 3 --family orbit trees all` and
   `… --m 4 --family orbit trees perms --save-realizations analysis/flag_hitting_bound_m4_realizations.txt` (2 s each);
   outputs `analysis/flag_hitting_bound_m3.out`, `_m4.out`. Re-verified by `tests/test_paper_numbers.py`, `tests/test_flag_hitting.py`.
   The degree-4 Z_3-invariant search (`analysis/isd_s2_d4_z3.out`, 42 min) found no certificate using a size-4 monomial; its
   certificate `analysis/isd_s2_d4_z3_cert.txt` (profile {2: 24, 3: 288}) is re-verified by the same test file.
3. **Spread Lemma evidence.** `PYTHONPATH=. .venv/bin/python analysis/spread_mcmc.py --m 4 5 6 7 8 --chains 4 --sweeps 3000 --match --nsamp 2000 --save analysis/spread_mcmc.pkl`
   (10 min) regenerates `analysis/spread_mcmc.out` byte-for-byte except for the `[Ns]` wall-clock stamps. Diagnostics and caveats:
   `docs/spread_evidence.md`; qualitative fast test `tests/test_spread_mcmc.py`.

The five mechanical proof checks (tower identity m ≤ 5, extension lemma all shapes n ≤ 7, rule = solver n ≤ 7, realization m ≤ 10,
separating functional m ≤ 5) are in `tests/test_tower.py`, `tests/test_extension_lemma.py`, `tests/test_rule_vs_solver.py`,
`tests/test_realize_explicit.py`; tiers and runtimes in `docs/README.md`.

## Install

Python 3.13 (3.10+ should work; the stored outputs were produced with 3.13.7) and Rust (cargo; the outputs used rustc 1.91.1).

    python3 -m venv .venv && .venv/bin/pip install -r requirements.txt     # numpy, pytest, ortools (pinned)
    make gf2solve                                                          # cargo build --release in gf2solve/
    make test-fast

`ortools` (CP-SAT) is needed only by `analysis/isd.py`, `analysis/flag_hitting_bound.py` and the tests that re-verify their numbers
(`test_min_size.py`, `test_paper_numbers.py`, `test_flag_hitting.py`); everything else is numpy + the Rust solver. The Rust solver
is built automatically on first use if `cargo` is on the path.

## Tests

| tier | command | contents | time |
|---|---|---|---|
| fast | `make test-fast` | every proof check that fits, all paper numbers with random-labeling certificate checks, MCMC sanity | ~1.5 min |
| slow | `make test-slow` | exhaustive certificate checks, CP-SAT optimality of 304, all 11 tree shapes at n = 7, S^3 sparse dual, realization m = 9 | ~6 min |
| veryslow | `make test-veryslow` | S^7 chain-gadget bound (16 min, pure-Python edge enumeration), S^3 degree-3 dense GF(2) lower bound (12 min), S^2 Q degree-5 certificate (7 min), realization m = 10 (3 min), m = 8 caterpillar realization | ~40 min |
| paper | `make reproduce-paper` | the three cited computations and the five proof checks (fast + slow of those files) | ~5 min |

`./reproduce.sh` = fast + slow (~8 min); `./reproduce.sh all` = everything (~50 min). CI (`.github/workflows/test-fast.yml`) runs the fast tier.

## Layout

- `docs/` — one file per proof, paper numbering; `docs/README.md` is the index (statement → file → test → tier → runtime).
- `kyfan/` — the package: `SignedComplex(m)`, labelings, violating pairs, group actions, SA-dual builder, primal, solvers,
  `tower.py`, `lower.py`, `tree_rule.py`, `realize_explicit.py`, `gadget.py`, `ball.py`; `gf2solve/` is the Rust sparse GF(2) solver.
- `tests/` — the settled numbers as assertions (`HISTORY.md` maps each number to its test).
- `analysis/` — the paper-cited scripts and their stored outputs (`analysis/README.md`); `analysis/legacy/` the exploratory drivers.
- `HISTORY.md` — development log (all settled numbers, negative results, runtimes). `CITATION.cff`, `LICENSE` (MIT).
