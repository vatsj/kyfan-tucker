# Results log

All degrees are Sherali–Adams level in the one-hot quotient (see CLAUDE.md §1).
Every number below is an assertion in `tests/` (the tests *are* this table); `./reproduce.sh` runs the fast subset (~3 min),
`./reproduce.sh all` everything (~30 min). Runtimes: M-series laptop, 10 cores, 16 GB, Python 3.13, numpy 2.5 (2026-09-16).

## Status (2026-09-16) — what has happened since CLAUDE.md was written

Everything in CLAUDE.md §3a and §4 Tasks 1–4 is done, in one session; the brief is kept as the record of what was asked.
- **Phase 0**: the original `scripts/` were reviewed (`review.md`), reproduced number-for-number, then replaced by the package `kyfan/`
  and `tests/` (git history has the scripts at commit 5a4dd50). Two findings: `primal.py`'s "group preserves violating pairs: False"
  was a broken check, not a bug; and the Ky Fan "exactly n+1" row lacked a lower bound — the author supplied the separating-functional
  proof now in CLAUDE.md §2 (`kyfan/lower.py`).
- **Task 1**: Rust sparse GF(2) solver (`gf2solve/`), 1.8M unknowns in 20 s. **Task 2**: structure of the degree-3 pseudo-solution on S^3.
  **Task 3**: the ball is tight over F_2 at m=3,4.
- **Main result (theorem): Tucker's F_2 degree on S^n is exactly n+1 for every n.** A *restriction lemma* reduces the lower
  bound to a **tree gadget** (one top simplex through the pole, residual domains = a binary conflict tree). Two parts, both
  now proved: **(i)** the tree gadget has degree n+1 for every tree (explicit pseudo-solution + extension lemma,
  `TREE_GADGET_THEOREM.md`); **(ii)** an explicit closed-form equatorial labeling λ(u) = s·(n−q+1) realizes it on S^n for
  all n (`REALIZATION_THEOREM.md`). The earlier searched `GADGETS` (n=4..7) and the deterministic `realize` (n≤8) are now
  superseded by (ii) but kept as independent checks; the "never backtracks" question is moot.

Reading order for a new session: this Status; row 6 of the table; `TREE_GADGET_THEOREM.md` (i) and `REALIZATION_THEOREM.md`
(ii); `kyfan/tree_rule.py`, `kyfan/realize_explicit.py`, `kyfan/gadget.py` (the restriction lemma and `verify`); then
`analysis/README.md` for the exploratory drivers. To reproduce the theorem in seconds:
`.venv/bin/python -m pytest tests/test_realize_explicit.py tests/test_rule_vs_solver.py tests/test_extension_lemma.py -q`.

Layout: `kyfan/` package (see CLAUDE.md §3), `tests/` (fast by default; `slow` runs by default, `veryslow` needs `-m veryslow`),
`gf2solve/` (Rust; built automatically on first use, needs `cargo`), `analysis/` (exploratory scripts + their `.out` outputs, indexed in `analysis/README.md`),
`review.md` (Phase-0 review of the deleted scripts). Setup: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`.

## Settled

| # | statement | field | degree | evidence | test | runtime |
|---|---|---|---|---|---|---|
| 1 | Tucker, S^2 (m=3, labels ±1,±2) | F_2 | = 3 | dual: 1,104 unknowns, consistent at d=2 (rank 507, 597 free); 12,688 unknowns, inconsistent at d=3 (rank 8,400). Explicit certificate (5,616 violating size-3 unknowns; ~458 monomials — the minimum is 304, see "Certificate size") verified exactly: canonical one-hot basis expansion, and exhaustively on all 4^13 labelings | `test_nsdeg.py`; `test_primal.py::test_tucker_s2_F2_degree3_certificate_exact`, `…_exhaustive` (slow) | 2 s; 10 s, 30 s |
| 2 | Tucker, S^2 | Q | = 5 | dual consistent at d=3 over F_p (p=1000003, rank 8,401); full-group orbit primal (valid over Q/F_p, |G|=384): inconsistent at d=3 (55 orbits) and d=4 (657 orbits) — one-way sound; at d=5: 924,480 violating unknowns in 5,482 orbits, float rank 1,976, relative residual 3e-15, max error 1e-14 on fresh labelings, and an exact F_p certificate | `test_nsdeg.py`; `test_primal.py::test_tucker_s2_Fp_symmetric_degree3_4_none`; `…::test_tucker_s2_Q_degree5_symmetric_certificate` (veryslow) | 5 s; 5 s; 5 min |
| 3 | Tucker, S^3 (m=4, labels ±1..±3) | F_2 | = 4 | d=3, via the dual: the Z_3×Z_3-invariant degree-3 SA-dual (1,834,800 non-violating size-3 partial labelings in 204,048 orbits; 111,233 rows) is consistent, rank 94,064 — an explicit invariant pseudo-solution (support 8,036 orbits), checked against all 987,265 unreduced consistency equations. Averaging is valid over F_2 since |Z_3×Z_3| = 9 is odd. The unrestricted dual (987,457 × 1,834,800) is also consistent, rank 827,876. Cross-check via the sampled primal (z9): inconsistent, 33,312 orbit-unknowns, rank 27,584. d=4: tower | `test_sparse.py::test_s3_degree3_pseudo_solution_via_dual` (30 s); `test_s3.py::test_s3_degree3_F2_no_certificate_z9` (veryslow, 7 min) | 30 s; 7 min |
| 4 | Ky Fan, S^2, labels ±1..±3 | F_2 | = 3 | ≤: explicit certificate (13,176 violating unknowns; ~1,286 monomials) verified exactly (canonical basis: Σc_α1_α ≡ A_+ + 1). ≥: row 5 | `test_primal.py::test_kyfan_s2_F2_degree3_certificate_exact` | 60 s |
| 5 | Ky Fan, S^n, k ≥ n+1 | F_2 | = n+1 | ≤: tower theorem — local lemma g_n = 0 on non-complementary tuples, exhaustive n≤3, k≤5; telescoped identity on random labelings m=3,4,5 (29 / 221 / 2,141 local-lemma instances). ≥: separating functional (CLAUDE.md §2): E(A_+)=1, E(1)=0, E(1_α)=0 on size-n α (exhaustive m=3, sampled m=4,5); the negation-closed box gives E(A_+)=0 | `test_tower.py` | 15 s |
| 6 | Tucker, S^n | F_2 | **= n+1 for all n (theorem)** | upper: tower (row 5, A_+ ≡ 0 when k ≤ n). Lower: restriction lemma + tree gadget. **(i)** the tree gadget has degree n+1 for every binary conflict tree (`TREE_GADGET_THEOREM.md`, `kyfan/tree_rule.py`); **(ii)** the explicit equatorial labeling λ(u) = s·(n−q+1), q = start of u's final same-sign run, realizes the reverse-caterpillar tree as the pole-chain residual — validity and domains proved in `REALIZATION_THEOREM.md`, no search. Mechanically checked: (i) all tree shapes n≤6 (+n=7 sparse); (ii) validity+domains m≤10 and residual degree-n dual consistent+unique m≤8. Also verified for n=2,3 by direct dual (rows 1,3) and by the searched `GADGETS` m=4..8 (`kyfan.gadget.verify`) | `test_realize_explicit.py`, `test_rule_vs_solver.py`, `test_extension_lemma.py`, `test_gadget.py`, `test_realize.py` | (i) 2 s; (ii) m≤8: 4 s |

## Task 3 + restriction gadgets — the ball, and the S^4, S^5, S^6, S^7 lower bounds

**Ball (CLAUDE.md §4 Task 3; `kyfan/ball.py`).** Fix a valid equatorial labeling L_eq; variables = the cap {x_m = +1};
constraints = no complementary cap edge, plus unary domain restrictions λ(y) ≠ −L_eq(w) for equatorial w ≤ π(y).
m=3: F_2 degree **3 for all 80** valid L_eq (= sphere degree; tight) and F_p degree 3 (sphere over Q: 5 — the ball is strictly easier over Q).
m=4: F_2 degree > 3 for **40/40** sampled L_eq, hence 4 by restriction of the tower certificate (= sphere degree; tight). F_p rank equals F_2 rank at every level tested.
The full m=5 ball at degree 4 has ~5×10^8 unknowns (not attempted) — but a *smaller* restriction suffices, see below. `analysis/ball_m3.out`, `ball_m4.out`.

**Restriction lemma.** Let ρ be a non-violating partial labeling of the free vertices of S^n, U the unfixed vertices.
Substituting ρ into a degree-d sphere certificate 1 = Σ c_α 1_α gives 1 = Σ_{α: α|_off-U ⊆ ρ} c_α 1_{α|_U} on all labelings of U;
each surviving α|_U contains a complementary pair inside U or a label forbidden by a fixed neighbour (a unary violation),
so this is a degree-≤d certificate of the residual CSP (domains D(u) = labels ∖ {−ρ(w): w ~ u fixed}, pairwise non-complementary on edges).
Hence: **a consistent degree-d residual dual ⇒ sphere degree > d.** (`kyfan/gadget.py` docstring; the ball is the special case U = cap.)

**Chain gadgets (`analysis/residual_*.py`, `gadget*.py`; `kyfan/gadget.py`).** ρ = a valid labeling of the equator S^{n−1} with magnitudes 1..n−1 except one antipodal
pair of top simplices ±σ labeled ∓n, pulled back to the cap; U = e_{n+1} ∪ lift(σ), a top simplex of S^n through the pole, n+1 variables.
Then D(e) = {n}, every other domain contains −n, and the remaining labels form a **binary conflict tree**: root ±1 splits the chain into a +1 group and a −1 group,
each group splits recursively on a fresh magnitude, a leaf's domain = its root-to-leaf sign path. Found by search at m=4,5 (`gadget2.py`: forbid one value on the star of σ's rank-1 vertex),
then realized to order at m=6,7,8 by deriving per-vertex allowed label sets from the target domains (`gadget3.py`; success 6/6, 3/3, 2/2 seeds; ~1 s, ~20 s, ~12 min per seed — the m=8 time is the pure-Python edge enumeration of S^7).

| m (S^{m−1}) | U domains (rank n..1, pole) | degree-n dual | conclusion |
|---|---|---|---|
| 4 | {1,−2,−3} {1,2,−3} {−1,−3} · {3} | 49 unknowns, rank 49, consistent, support 27 | S^3 ≥ 4 (known) |
| 5 | {1,−2,−4} {1,2,−4} {−1,3,−4} {−1,−3,−4} · {4} | 215, rank 215, support 81 | **S^4 ≥ 5** |
| 6 | {1,2,−5} {1,−2,−5} {−1,3,−5} {−1,−3,4,−5} {−1,−3,−4,−5} · {5} | 1,131, rank 1,131, support 243 | **S^5 ≥ 6** |
| 7 | {1,2,−6} {1,−2,−6} {−1,3,4,−6} {−1,3,−4,−6} {−1,−3,5,−6} {−1,−3,−5,−6} · {6} | 5,915, rank 5,915, support 729 | **S^6 ≥ 7** |
| 8 | {1,2,3,−7} {1,2,−3,−7} {1,−2,−7} {−1,4,5,−7} {−1,4,−5,−7} {−1,−4,6,−7} {−1,−4,−6,−7} · {7} | 30,423, rank 30,423 | **S^7 ≥ 8** (2/2 seeds, 12 min each) |

In every case the residual is UNSAT with every proper subset satisfiable, and the degree-n pseudo-solution is **unique** with support exactly 3^n.
The equatorial labelings are stored in `kyfan.gadget.GADGETS`; a uniform *deterministic* alternative is `kyfan.gadget.realize(m)` (statement (ii), below).

**Abstract tree gadgets (`kyfan/abstract_gadget.py`).** Forgetting the sphere: variables = pole {n} plus the leaves of any binary conflict tree on n leaves
(magnitudes 1..n−1 on internal nodes, leaf domain = path ∪ {−n}), all pairs constrained. Every tree tried (n = 3,4,5, all shapes incl. caterpillars)
has degree n+1 with a unique degree-n pseudo-solution of support 3^n (`test_gadget.py::test_abstract_tree_gadgets`, `analysis/abstract_trees.out`).
Non-tree designs (e.g. {−1,3},{−1,−3,4},{−1,−4}) have lower degree.

**Proof for all n — both parts done.** (i) The tree gadget has F_2 degree exactly n+1 for every binary conflict tree
(`TREE_GADGET_THEOREM.md`): explicit pseudo-solution E (the block rule, `kyfan/tree_rule.py`), consistency via the extension
Lemma (★). Checked against the solver for all tree shapes n ≤ 6 and against the sparse solver at n = 7
(`test_rule_vs_solver.py`, `test_extension_lemma.py`). (ii) **Explicit realization** (`REALIZATION_THEOREM.md`,
`kyfan/realize_explicit.py`): the closed-form labeling λ(u) = s·(n−q+1), where q is the start of u's final maximal
same-sign run, is valid (no complementary comparable pair — one-line proof) and its pole-chain residual is exactly the
reverse caterpillar D(r) = {+(n−r)} ∪ {−(n−r+1),…,−n}, pole {+n}. No search. Validity + domains checked m ≤ 10, residual
degree-n dual consistent+unique m ≤ 8 (`test_realize_explicit.py`). Together with the restriction lemma and the tower upper
bound this is a theorem: **Tucker's F_2 degree on S^n is n+1 for all n.** The searched `GADGETS` (n=4..7) and the greedy
`realize` (n≤8, zero backtracks) remain as independent constructions; the "never backtracks" question is retired.

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

## Task 2 — the degree-3 pseudo-solution on S^3 (labels ±1..±3), with the n=2 comparison

All exact: support restrictions and knockouts are full-group-invariant sets, so solving the Z_3×Z_3-invariant system
(odd order) decides existence of *any* pseudo-solution with that support. Scripts in `analysis/` (outputs `*.out` alongside).

**Object.** 1,834,800 non-violating size-3 partial labelings in 808 full-group orbit types and 204,048 Z_3×Z_3 orbits;
invariant system 111,233 × 204,048, rank 94,064, 109,984 free parameters (`test_sparse.py`). Unrestricted: rank 827,876, 1,006,924 free.

**Type-level knockout** (`s3_degree3_types.py`, 16 min): exactly 4 of 808 types are necessary — the 2-simplices (chains) of
rank patterns (1,2,3), (1,2,4), (1,3,4), (2,3,4) with all-distinct magnitudes. (n=2: distinct-magnitude edges at the 3 rank pairs, plus "two unrelated rank-1 vertices, equal labels".)

**Marginal knockout** (`s3_degree3_marginals.py`, 3 min): forcing E = 0 on a full-group type of size-2 partial labelings —
necessary: the distinct-magnitude edges at all 6 rank pairs; the n=2 type "unrelated rank-1 pair, equal labels" is *not* necessary as a marginal at n=3. Size-1: all 4 rank types (trivially).

**Class-level knockout / minimal families** (`class_knockout.py`, invariant classes = (simplex?, #edge pairs, magnitude pattern AAA/AAB/ABC)):

| | n=2, degree 2 (4 classes) | n=3, degree 3 (12 classes) |
|---|---|---|
| necessary | unrelated-AA, edge-AA, edge-AB | AAB with 1 edge, AAB with 2 edges, simplex-AAB, simplex-ABC |
| minimal sufficient family | = the necessary classes (99 free params) | necessary + unrelated-AAB (115,823 free); no other class substitutes |
| dispensable (each alone) | unrelated-AB | all AAA (constant magnitude) classes, all non-simplex ABC classes |

Persisting pattern: a degree-n pseudo-solution must be supported on (i) the (n−1)-simplices with all-distinct magnitudes
and (ii) "exactly two equal magnitudes" configurations in every geometry (with and without edges), and it can drop
all-distinct non-simplices entirely. (At n=2 "exactly two equal" = AA.)

**Restrictions** (`s3_degree3_restrictions.py`; all NO, as at n=2): distinct magnitudes only (52,768 unknowns, exact);
label-group invariant (85,248 orbits); hemisphere-stabilizer invariant (2,597); both with distinct magnitudes (9,880; 360).
Also NO at both n=2,3 (`simplex_support.py`, `magnitude_support.py`): simplices only; no simplices; contains an edge; pairwise unrelated only;
unrelated pairs → distinct magnitudes; edge pairs → distinct magnitudes. "Unrelated pairs → equal magnitude" works at n=2 (99 free) but NOT at n=3.

**Base labeling + corrections** (`base_labeling.py`): L_0 = pullback of a valid equatorial labeling with a single ±1 pair,
+1 on e_m (exactly one antipodal pair of complementary edges); support = partial labelings disagreeing with L_0 in ≤ h positions.
NO solution for every h < d, for 5 different L_0 at n=2 and 4 at n=3 (h=2 system: 766k unknowns, 13 s). Random L_0: ~12% admit
an h=1 solution at n=2 (2–7 free parameters), 0 of 12 admit h=2 at n=3. The pseudo-solution is not "a labeling plus local corrections".

**Greedy type-level family** (`s3_degree3_greedy.py`, 20 min): 379 of 808 types (509k partial labelings), order-dependent; no sparse type structure.

## Structural facts (degree-2 pseudo-solution on S^2) — `test_pseudo2.py`, 5 s

| fact |
|---|
| solution space: 1,104 unknowns, rank 507, 597 free parameters |
| no full-group-invariant solution (20 orbits) — a statement about invariant solutions only, |G| even |
| no label-group-invariant (198 orbits), no hemisphere-stabilizer-invariant (44), no distinct-magnitudes-only (624 unknowns) solution; nor with both restrictions (78; 19) |
| coordinate-Z_3-invariant solution exists (368 orbits, 201 free) |
| knockout by full-group orbit type: exactly 4 of 20 types are necessary — distinct-magnitude edges at rank pairs (1,2), (1,3), (2,3); equal labels on two unrelated rank-1 vertices |

## Certificate size (2026-09-17) — `SIZE_LOWER_BOUND.md`, `analysis/isd.py`, `analysis/spread_mcmc.py`

Size = number of monomials (violating partial labelings) in a certificate 1 = Σ 1_α; #flags = m!·2^m (48 on S^2).

**Status of the theory (`SIZE_LOWER_BOUND.md`, collaborator's write-up, all ingredients proved in this repo):**
- **Proved**: every certificate of degree exactly n+1 has size ≥ #flags (Thm 1); every certificate of degree ≤ D has size
  ≥ #flags / C(D, m), superpolynomial (N^{Θ(log log N)}) for D = O(n) (Thm 2). Both use only Lemma A (every flag forces a
  full-flag term) = restriction lemma + tree theorem (i) + realization (ii) + B_m-transport of the explicit gadget to every flag.
- **Reduces to the spread lemma for all degrees** (Thm 4): if for every flag there is a distribution on gadget labelings with
  P[ρ|_A = a] ≤ β^{|A| − cm} for all vertex sets A off the flag, then every certificate of any degree has size ≥ #flags/K^m.
  Thm 3 (proved) is the special case where every term uses ≥ μ magnitudes off the flags it serves.
- **Empirical support for the lemma**: the MCMC numbers below (the explicit gadget itself violates the lemma — its
  magnitude-1 class contains a 1/(2n) fraction of all flags — so D_F must be random; uniform on gadget-compatible labelings is the candidate).

**Minimum certificate size on S^2 (m=3, labels ±1,±2)** — `analysis/isd.py`; unknowns = all violating partial labelings of
sizes 2..d, rows = random labelings (one-way sound; every certificate below was re-verified on all 4^13 labelings).
The reduced system is very sparse (a size-2 violating monomial is the sum of its 4 size-3 extensions at any third vertex;
free columns have weight 4–17), so textbook Stern/Dumer on random rows degenerates (almost all keys 0); the ISD used is
Canteaut–Chabaud information-set updates + Prange + Dumer p=1, p=1+1 (matching on rows inside the rhs support) + a
first-improvement descent in the coset with the sparse null basis. Exactness comes from CP-SAT (`ortools`) on the reduced
system: x_piv = rhs + M_free x_free, minimise Σx_free + Σ y_i with y_i = rhs_i ⊕ (⊕ M_ij x_j) as XOR constraints.
True certificates are a subset of the sampled system's solutions, so the CP-SAT optimum is a lower bound; a verified
certificate attaining it makes the number exact.

| degree | unknowns | rank | first Prange | ISD best (time) | CP-SAT | lower bound from theory | verified |
|---|---|---|---|---|---|---|---|
| ≤ 3 | 5,760 monomials | 4,176 | 360 | **304** (≈40 s, 3.5k information sets) | **OPTIMAL 304** (31 s, 8 workers) | 48 (Thm 1) | exhaustive, 17 s; `tests/test_min_size.py` (slow, 70 s) |
| ≤ 4, Z_3-invariant | 33,512 orbits of 100,520 monomials | 16,902 | 136 (Dumer: 133) | **104 orbits = 312 monomials** (300 s, 435 information sets) | FEASIBLE 104, proven bound only **24** at the 1800 s cap | 12 (Thm 2: 48/C(4,3)) | exhaustive, 18 s; `analysis/isd_s2_d4_z3_cert.txt` |

- Degree 3: the minimum is exactly **304 = 16 size-2 + 288 size-3 monomials**, 6.3× the flag bound. (The ~458 quoted in row 1 was
  a first Prange solution, not a minimum.) Command: `python analysis/isd.py --d 3 --seconds 60 --cpsat 600 --save analysis/isd_s2_d3.pkl` (2.5 min).
- Degree 4 (Z_3 = coordinate 3-cycle on coords 0,1,2, odd order, so orbit sums are legitimate F_2 certificates): the best Z_3-invariant certificate found has orbit weight w* = 104 (312 monomials), and it uses **no size-4 monomial
  at all** (profile {2: 24, 3: 288}) — it is a Z_3-symmetrised degree-3 certificate, barely larger than the degree-3 optimum 304.
  CP-SAT (1800 s, 8 workers, hinted) did not improve on 104 and only proved w* ≥ 24, so the invariant minimum is in [24, 104].
  For the true (non-invariant) degree-4 minimum s the reduction gives s ∈ [w* − 8, 3 w*] (x + gx + g²x is invariant; 8 fixed orbits),
  i.e. only **s ∈ [16, 312]** from what is proved, and s ≤ 304 from degree 3. **Not decisive**: neither "degree 4 brings the size down
  a lot" nor "it doesn't" is established. The elimination (349 s dense on 45,000 × 33,512) and the CP-SAT relaxation are the
  bottlenecks; a longer CP-SAT run, or a direct non-invariant degree-4 model (100,520 columns), is the next step if this matters.
  Command: `python analysis/isd.py --d 4 --z3 --rows 45000 --seconds 300 --cpsat 1800 --save analysis/isd_s2_d4_z3.pkl` (≈ 6 min elimination + 5 min ISD + 30 min CP-SAT + 2 × 18 s verification ≈ 42 min).

**Spread-lemma MCMC (`analysis/spread_mcmc.py`, output `analysis/spread_mcmc.out`, pickle `analysis/spread_mcmc.pkl`).**
Glauber dynamics (heat bath, one random movable free vertex per move, uniform over its allowed labels) on the set G_m of
valid equatorial labelings that agree with the explicit gadget on the chain σ and whose pole-chain residual domains contain
the reverse-caterpillar domains (every element is a gadget by the label-agnostic tree theorem). Stationary distribution =
uniform on the component of the explicit gadget. Diagnostics: 4 chains × 3000 sweeps (25% burn-in), integrated
autocorrelation time τ (Sokal window), ESS, split-R-hat; match probabilities between two *independent* chains
(2000 samples each, thin 3), ± = standard error over 40 random sets A. Ergodicity of the dynamics is not proved.

| m (equator S^{m−2}, movable vertices) | frozen-move rate | τ (sweeps), R-hat | largest-magnitude fraction (explicit) | flag density of largest class (explicit 1/(2n)) | Hamming from start | runtime |
|---|---|---|---|---|---|---|
| 4 (10) | 0.47 | 2.0–2.6, 1.00 | 0.461 (0.462) | 0.094 (0.167) | 0.23 | 2 s |
| 5 (36) | 0.24 | 2.0–3.9, 1.00 | 0.371 (0.450) | 0.024 (0.125) | 0.43 | 7 s |
| 6 (116) | 0.12 | 2.0–4.7, 1.00 | 0.287 (0.446) | 0.0051 (0.100) | 0.58 | 22 s |
| 7 (358) | 0.062 | 2.0–4.3, 1.00 | 0.226 (0.445) | 0.0009 (0.083) | 0.68 | 77 s |
| 8 (1086) | 0.031 | 2.0–5.4, 1.00 | 0.183 (0.445) | 0.0001 (0.071) | 0.74 | 240 s + 240 s match |

Match probability P[L'|_A = L|_A] for two independent gadget samples, by |A| (random / chain / ball sets), and the fitted
per-vertex decay β (least squares on log P):

| m | |A|=1 | 2 | 3 | 4 | 6 | 8 | 12 | β (random / chain / ball) |
|---|---|---|---|---|---|---|---|---|
| 4 random | 0.67 | 0.49 | 0.27 | 0.18 | 0.076 | 0.033 | – | 0.65 / 0.59 / 0.65 |
| 5 random | 0.47 | 0.29 | 0.082 | 0.045 | 0.0079 | 0.0025 | 0.0002 | 0.48 / 0.44 / 0.52 |
| 6 random | 0.39 | 0.16 | 0.046 | 0.014 | 0.0033 | 0.0002 | <1e-4 | 0.36 / 0.45 / 0.47 |
| 6 ball | 0.47 | 0.21 | 0.14 | 0.048 | 0.011 | 0.0019 | 0.0001 | |
| 7 random | 0.26 | 0.099 | 0.033 | 0.0049 | 0.0004 | <1e-4 | <1e-4 | 0.29 / 0.32 / 0.41 |
| 7 chain | 0.95 | 0.54 | 0.15 | 0.040 | 0.0039 | – | – | |
| 7 ball | 0.45 | 0.18 | 0.047 | 0.030 | 0.0031 | 0.0013 | <1e-4 | |
| 8 random | 0.23 | 0.066 | 0.020 | 0.0023 | 0.0002 | <1e-4 | <1e-4 | 0.22 / 0.27 / 0.34 |
| 8 chain | 0.95 | 0.49 | 0.17 | 0.038 | 0.0016 | – | – | |
| 8 ball | 0.31 | 0.12 | 0.038 | 0.014 | 0.0020 | 0.0002 | <1e-4 | |

What this shows: all statistics have ESS in the thousands and R-hat = 1.00 (same-start chains, so R-hat is a weak check;
the small τ and the stable Hamming distance from the start are the stronger evidence of mixing). The largest magnitude
class shrinks toward the uniform share 1/n and its flag density decays much faster than the explicit gadget's 1/(2n);
match probabilities decay geometrically in |A| on all three set types, and the fitted β itself *decreases* with m
(random sets: 0.65, 0.48, 0.36, 0.29, 0.22 for m = 4..8, roughly 1/n as for independent uniform magnitudes; chains 0.59→0.27;
balls 0.65→0.34). Chains start high (P ≈ 0.95 at |A|=1: the rank-1 vertex of a random chain is often in the rigid zone)
but decay as fast as random sets from |A| = 2 on. This is consistent with the spread lemma with a fixed β < 1 and a small c. The residual frozen-move rate (vertices with ≤ 1 allowed
label) falls like ~1/n — the "rigid zone" near σ allowed by the c·m term.

## Negative / methodological

- Unsubdivided octahedron: A_+ ≡ 1 on all 192 valid labelings with labels ±1..±4 (degenerate). `test_complex.py::test_octahedron_degenerate`.
- Fully-symmetric F_2 certificates do not exist on S^2 at d=3 (55 orbits) or d=4 (657), although an asymmetric d=3 one does. Never conclude from symmetric F_2 searches. `test_primal.py::test_tucker_s2_F2_full_group_symmetric_none_d3_d4` (slow, 3 s).
- Full-group symmetric search on S^3 at d=3 (299,280 violating in 232 orbits): none over F_2 (uninformative), none over R (residual 7e-2 ≫ 1e-3). `test_s3.py::test_s3_full_group_symmetric_degree3_none` (slow, 20 s).
- Pullback of an equatorial pseudo-solution to S^n loses one degree (hand argument, CLAUDE.md §2).
- Consistency equations at level d−1 suffice: building only those gives the same ranks (507; 8400/8401) as building all levels. The package builds only level d−1.

## Open

- Tucker's F_2 degree on S^n is **settled** (= n+1 for all n, row 6): (i) `TREE_GADGET_THEOREM.md` + (ii) `REALIZATION_THEOREM.md`.
- Task 2 follow-ups (not needed for the theorem, but of interest): explain the AAB/simplex-ABC necessity pattern of the S^3
  degree-3 pseudo-solution for general n; explain why constant-magnitude configurations are dispensable at n=3 but not n=2.
- Independent corroboration recorded but not a test: the deterministic `realize` reached S^8 = 9 (876,809-unknown dual) before
  its exploratory run was stopped.
