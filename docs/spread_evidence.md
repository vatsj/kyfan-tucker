# Empirical support for the Spread Lemma (`docs/size_lower_bound.md` §4)

Script: `analysis/spread_mcmc.py`. Output: `analysis/spread_mcmc.out`, produced by

    PYTHONPATH=. .venv/bin/python analysis/spread_mcmc.py --m 4 5 6 7 8 --chains 4 --sweeps 3000 --match --nsamp 2000 \
        --save analysis/spread_mcmc.pkl        # ~10 min; the .pkl is gitignored

All randomness is seeded (see "Reproducibility"); rerunning the command regenerates `spread_mcmc.out` byte-for-byte except for the
wall-clock stamps `[Ns]`. Fast test: `tests/test_spread_mcmc.py` (m = 4, 5, short chains; qualitative assertions only).

This document says exactly what was sampled, what each printed number means, and what the experiment does *not* establish. Nothing
here is a proof; the Spread Lemma is a conjecture.

## 1. What is sampled

Fix m and the pole flag of S^n, n = m − 1, whose equatorial part is the chain σ = { w_r = −(e_1 + … + e_r) : 1 ≤ r ≤ n } of the
[n]-complex. The *gadget-compatible* set G_m consists of the labelings L of the free vertices of the equator (labels ±[n]) such
that

1. L is valid (no complementary comparable pair);
2. L agrees with the explicit labeling of `docs/realization.md` on σ (L(w_r) = −n);
3. the residual domains on the pole chain *contain* the reverse-caterpillar domains: +n ∈ D(0) (so no equatorial vertex carries
   −n except σ, i.e. every vertex outside ±σ avoids the label −n, its lift being adjacent to the pole), and
   T[r] ∪ {−n} ⊆ D(r) with T[1] = {+(n−1)}, T[r] = {+(n−r)} ∪ {−(n−r+1), …, −(n−1)} for 2 ≤ r < n, T[n] = {−1, …, −(n−1)}
   (so every neighbour u of w_r avoids the labels −T[r] ∪ {+n}).

Condition 3 is encoded as a per-vertex *static* forbidden set (`GadgetSampler.static_ok`); conditions 1–2 as the neighbour
constraints and as the set of *fixed* vertices (the reps of σ, which are never resampled). By the tree-gadget theorem (statement (i),
which is label-agnostic — `docs/tree_gadget.md`) every L ∈ G_m realizes a gadget on the pole flag: its residual has a
degree-n pseudo-solution supported on the caterpillar sub-domains. The uniform distribution on G_m is the candidate D_F of the
Spread Lemma for this flag; by symmetry (the signed permutations of the coordinates act transitively on flags) it induces a
candidate for every flag.

The printed header for each m gives the number of free vertices N = (3^n − 1)/2 and the number of movable ones N − n.

## 2. The Markov chain

Glauber dynamics (single-site heat bath): repeat, pick a uniformly random movable free vertex i, compute the set of labels allowed
at i by the static forbidden set and by the current labels of i's neighbours (the constraint on an edge {x, y} with reps (i, s_i),
(j, s_j) is s_i L[i] ≠ −s_j L[j]), and set L[i] to a uniformly random allowed label. One *sweep* is N such single-site updates
(N = number of free vertices, not the number of movable ones). The label of i itself is always allowed, so the chain never gets
stuck; a move at a vertex with at most one allowed label is a *frozen move* and the fraction of frozen moves is reported per chain.

Heat-bath updates are reversible with respect to the uniform distribution on G_m, so the chain's stationary distribution is uniform
**on the connected component of its start state**. All chains start at the explicit gadget, the only member of G_m known a priori.

**Ergodicity is not proved.** If the single-site move graph on G_m is disconnected, every number below describes the component of the
explicit gadget, not G_m. The Hamming-distance trace (fraction of free vertices whose label differs from the explicit gadget) shows
that the chains move far from the start (mean 0.23, 0.43, 0.58, 0.68, 0.74 for m = 4..8, with sd ≈ 0.01–0.09), and the frozen-move
rate falls with m (0.47, 0.24, 0.12, 0.06, 0.03), but neither is evidence of connectivity.

## 3. Diagnostics (`run_chains`, `summarize`)

Per m: `--chains 4` independent chains with seeds 1000·(c + 1), c = 0..3, each run for `--sweeps 3000` sweeps; the first 25 % (750
sweeps) are burn-in, after which the scalar statistics below are recorded every `--every 2` sweeps (1,125 records per chain, 4,500
in total). **All four chains start at the same state** (the explicit gadget), so agreement between chains tests mixing from a common
start, not convergence from over-dispersed starts as in the textbook Gelman–Rubin protocol; a between-chain diagnostic cannot detect
a component of G_m that no chain reaches.

Statistics of the current state (`GadgetSampler.stats`):

- `top_mag`: the largest fraction of free vertices carrying a single magnitude;
- `top_label`: the largest fraction of free vertices carrying a single signed label;
- `n_mags`: the number of distinct magnitudes in use (equal to n at every recorded state of the stored run);
- `flag_density`: max over magnitudes k of (number of maximal chains of the equator all of whose n vertices have |label| = k) /
  (n! 2^n = number of maximal chains). This is the quantity Theorem 3's "monochromatic region" terms exploit: for the explicit
  gadget it equals 1/(2n) (printed as a check on the first line of each block);
- `hamming`: the fraction of free vertices whose label differs from the explicit gadget's.

For each statistic the output lists: the mean and standard deviation over all 4,500 records; the integrated autocorrelation time
τ in sweeps (Sokal's automatic windowing with c = 5 on the per-chain traces of records, averaged over chains and multiplied by the
recording interval); the effective sample size ESS = Σ_chains (records / τ_chain); the split-R̂ (each chain cut in halves, 8
half-chains, Gelman–Rubin statistic); and the per-chain means. In the stored run τ ≤ 5.4 sweeps for every statistic, ESS ≥ 1,688,
and every R̂ ≤ 1.002; the printed values are what the run produced and were not tuned.

The standard deviations are of the *statistic across samples*, not standard errors of the mean; a standard error can be formed as
sd/√ESS.

## 4. Match probabilities (`match_probabilities`)

The Spread Lemma bounds P_{ρ∼D_F}[ρ|_A = a] ≤ β^{|A| − cm}. The experiment estimates a proxy: the *collision probability* of two
independent draws from (the sampled approximation of) D_F on a set A,

    P_A = P_{ρ, ρ' independent} [ ρ|_A = ρ'|_A ] = Σ_a P[ρ|_A = a]²,

which satisfies max_a P[ρ|_A = a] ≥ P_A ≥ (max_a P[ρ|_A = a])² — so a geometric decay of P_A in |A| is necessary for the lemma and
implies a geometric decay of the maximum with rate √β at worst; it does not by itself give the uniform bound the lemma asks for.

Procedure: two independent chains A and B (seeds 11 and 22), each burned in for max(200, sweeps/4) = 750 sweeps, then `--nsamp 2000`
states recorded from each after every `--thin 3` sweeps. For each set type and each size t ∈ {1, 2, 3, 4, 6, 8, 12}: 40 random sets A
of t movable vertices are drawn (rng seed 33); for each, 200 B-states are chosen uniformly as targets and the reported number is the
fraction of the 2,000 A-states that agree with the target on A, averaged over targets and then over the 40 sets. The `±` value is the
standard deviation over the 40 sets divided by √40 — the spread across *sets*; it does not account for the autocorrelation of the
2,000 states within a chain (thinning by 3 sweeps ≈ τ leaves residual correlation). `nan` means no set of that size and type exists
(a chain has only n vertices; small down-sets have few movable vertices).

Set types (`sample_sets`): `random` = t uniformly random movable vertices; `chain` = the movable vertices of a uniformly random
maximal chain of the equator (built coordinate by coordinate with random signs), truncated to t; `ball` = t movable vertices from the
down-set of a uniformly random top vertex (rank n), in random order. Chains are the sets on which labels are most constrained (a
chain is a clique of the comparability graph, and consecutive labels can never be complementary), so `chain` is the adversarial
type for a spread bound; `ball` sets are the local neighbourhoods on which a certificate term would be supported.

The *fitted β per vertex* is exp of the slope of a least-squares line through (t, log P_t) over the sizes with a positive estimate;
it is a summary of the decay rate, not an estimate of the lemma's β (which must dominate the *maximum* over a, uniformly in A and
m, and may carry an additive c·m in the exponent). Estimates below ≈ 1/2000 are dominated by the sample size (2,000 states per
chain) and print as 0.0000; they were excluded from the fit when exactly zero.

## 5. Results of the stored run (m = 4..8)

| m | n | movable | top_mag (explicit) | flag_density (explicit 1/(2n)) | τ_max | ESS_min | R̂_max |
|---|---|---|---|---|---|---|---|
| 4 | 3 | 10 | 0.461 (0.462) | 0.094 (0.167) | 2.6 | 3,517 | 1.002 |
| 5 | 4 | 36 | 0.371 (0.450) | 0.024 (0.125) | 3.9 | 2,364 | 1.001 |
| 6 | 5 | 116 | 0.287 (0.446) | 0.0051 (0.100) | 4.7 | 2,170 | 1.001 |
| 7 | 6 | 358 | 0.226 (0.445) | 0.0009 (0.083) | 4.3 | 2,125 | 1.002 |
| 8 | 7 | 1,086 | 0.183 (0.445) | 0.0001 (0.071) | 5.4 | 1,688 | 1.001 |

Match probabilities, `random` sets (mean over 40 sets; standard errors in the `.out` file):

| m | t=1 | 2 | 3 | 4 | 6 | 8 | 12 | fitted β |
|---|---|---|---|---|---|---|---|---|
| 4 | 0.671 | 0.485 | 0.270 | 0.181 | 0.076 | 0.033 | — | 0.645 |
| 5 | 0.473 | 0.287 | 0.082 | 0.045 | 0.0079 | 0.0025 | 0.0002 | 0.479 |
| 6 | 0.391 | 0.161 | 0.046 | 0.014 | 0.0033 | 0.0002 | 0.0000 | 0.355 |
| 7 | 0.264 | 0.099 | 0.033 | 0.0049 | 0.0004 | 0.0000 | 0.0000 | 0.294 |
| 8 | 0.233 | 0.066 | 0.020 | 0.0023 | 0.0002 | 0.0000 | 0.0000 | 0.220 |

Fitted β for `chain` sets: 0.587, 0.440, 0.449, 0.320, 0.270; for `ball` sets: 0.650, 0.515, 0.466, 0.406, 0.335 (m = 4..8). On
chain sets P_1 ≈ 0.8–0.95: the movable vertex of lowest rank on a chain through σ's neighbourhood is nearly forced, which is the
rigid zone the additive c·m term in the lemma is meant to absorb.

What these numbers support:

- the flag density of the largest magnitude class under the (sampled) uniform distribution on G_m decays much faster than the
  explicit gadget's 1/(2n) — roughly by a factor 4–5 per step in m — so the monochromatic-region terms that defeat Theorem 3
  against the explicit gadget are rare against a random gadget;
- collision probabilities decay geometrically in |A| for all three set types, with a per-vertex rate that decreases with m.

What they do not establish: ergodicity of the chain (§2), the uniform-in-a bound (§4), the behaviour as m → ∞ (m ≤ 8 only), and
any value of β or c. A `nan` or a 0.0000 entry is a sample-size limit, not a measurement.

## 6. Reproducibility

Every random choice derives from `numpy.random.default_rng` with a fixed seed: chain c of `run_chains` uses `seed0 + 1000·(c+1)`
with `seed0 = 0`; `match_probabilities` uses seeds 11, 22 (the two chains) and 33 (set and target sampling), all with `seed = 0`.
The only floating-point computations are the statistics (means, FFT-based autocorrelation, the polyfit for β); on the machine that
produced the stored output (Python 3.13, numpy 2.5) a rerun of the command above reproduced `analysis/spread_mcmc.out` exactly,
modulo the `[Ns]` timing stamps. On another platform the last printed digit of τ, ESS, R̂ or β could in principle differ through
FFT/BLAS rounding; the sampled states, and hence the means, flag densities and match fractions, are integer computations on seeded
integer draws and should not.

`tests/test_spread_mcmc.py` (fast tier) runs m = 4 and m = 5 with short chains and asserts only qualitative facts with generous
tolerances: the mean flag density is below the explicit gadget's 1/(2n) by a factor ≥ 1.5, every gadget-compatibility constraint
holds at every recorded state, and the random-set match probability decreases at least geometrically over t = 1, 2, 3, 4 with a fitted
β < 0.9. It does not assert any digit of the stored output.
