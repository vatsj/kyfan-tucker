# analysis/ — the computations cited by the paper

Run from the repository root with `PYTHONPATH=. .venv/bin/python analysis/<script> ...`. Each `.out` file is the captured output of
the command listed. `.pkl` files are gitignored (regenerate by rerunning). The numbers these scripts produce are asserted by
`tests/test_paper_numbers.py`, `tests/test_min_size.py`, `tests/test_flag_hitting.py`, `tests/test_spread_mcmc.py`.

| script | command | output | what it shows |
|---|---|---|---|
| `isd.py` | `analysis/isd.py --d 3 --seconds 60 --cpsat 600 --save analysis/isd_s2_d3.pkl` (2.5 min) | `isd_s2_d3_cert.txt` | **minimum certificate size on S^2 at degree ≤ 3 = 304** (16 size-2 + 288 size-3 monomials): sparse-aware information-set decoding (Canteaut–Chabaud + Prange/Dumer + coset descent) finds 304, CP-SAT proves optimality on the reduced sampled system, the certificate is verified on all 4^13 labelings (`docs/size_lower_bound.md`) |
| `isd.py` | `analysis/isd.py --d 4 --z3 --rows 45000 --seconds 300 --cpsat 1800 --save analysis/isd_s2_d4_z3.pkl` (42 min) | `isd_s2_d4_z3.out`, `isd_s2_d4_z3_cert.txt` | degree ≤ 4, Z_3-invariant search: best certificate 312 monomials = 104 orbits, profile {2: 24, 3: 288} — **no size-4 monomial**; CP-SAT FEASIBLE 104, proven bound 24 (not decisive) |
| `flag_hitting_bound.py` | `analysis/flag_hitting_bound.py --m 3 --family orbit trees all` (2 s); `--m 4 --family orbit trees perms --save-realizations analysis/flag_hitting_bound_m4_realizations.txt` (2 s) | `flag_hitting_bound_m3.out`, `flag_hitting_bound_m4.out`, `flag_hitting_bound_m4_realizations.txt` (the 32 realizing labelings, re-verified without CP-SAT by `tests/test_paper_numbers.py`) | **Theorem 1′**: per-flag parity constraints from gadget pseudo-solutions, t(U) = minimum number of full-flag terms (CP-SAT exact, branch-and-bound cross-check); gadget families = label orbit / realized trees / realizable chain permutations (CP-SAT realizer); **S^2: t = 8, S^3: t ≥ 12** |
| `flag_hitting_bound_s2.py` | `analysis/flag_hitting_bound_s2.py` (1 s) | `flag_hitting_bound_s2.out` | independent S^2 computation in vertex coordinates: exhaustive 3,568 restrictions, 16 gadget types, t = 8; the minimum certificate's 12 full-flag terms satisfy all 16 parities |
| `spread_mcmc.py` | `analysis/spread_mcmc.py --m 4 5 6 7 8 --chains 4 --sweeps 3000 --match --nsamp 2000 --save analysis/spread_mcmc.pkl` (~10 min) | `spread_mcmc.out` | **Spread Lemma evidence**: Glauber chains over gadget-compatible equatorial labelings with τ / ESS / split-R-hat, flag densities, two-chain match probabilities and fitted β (`docs/spread_evidence.md`) |

`legacy/` holds the exploratory drivers from the degree computation (pseudo-solution structure on S^3, the ball form, the
searched chain gadgets) and the review of the original scripts; see `legacy/README.md`.
