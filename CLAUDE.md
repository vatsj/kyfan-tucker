# Brief for Claude Code: Nullstellensatz degree of Tucker / Ky Fan over F_2

You are continuing a research computation. Everything you need is in this repository. `results.md` is the ledger of every settled number (read it first — its top section says what has been done since this brief was written); `kyfan/` + `tests/` are the computations that produce them. Read this file fully before running anything. Prior results are settled — do not re-derive them; extend them.

## 1. The objects (precise definitions)

**Complex.** For an integer m, the vertices are the nonzero vectors in {-1,0,+1}^m ("signed subsets of [m]"). x <= y iff for every coordinate, x_i = 0 or x_i = y_i. Edges are comparable pairs; simplices are chains. This triangulates the sphere S^n with n = m-1. Vertex count 3^m - 1; top simplices m! * 2^m (maximal chains). The antipode of x is -x.

**Free vertices.** One representative per antipodal pair: the vertices whose first nonzero coordinate is +1. There are (3^m - 1)/2 of them. A labeling is specified on free vertices and extended by lambda(-x) = -lambda(x).

**Labels.** lambda(x) in {+-1, ..., +-k}. An edge {x,y} is *complementary* if lambda(x) = -lambda(y). A labeling is *valid* if it has no complementary edge.

**Tucker:** with k = n labels (k = m-1), no valid labeling exists (UNSAT).
**Ky Fan:** with any k, for every valid labeling the number of *positively alternating* top simplices is odd. A simplex is positively alternating if its n+1 labels have distinct magnitudes and, sorted by magnitude, the signs read +,-,+,-,... Count only positively alternating ones on the sphere (alternating simplices come in antipodal pairs, so "all alternating" is always even).

**Hemisphere / equator.** H = {x : x_m != -1}. Its equator {x_m = 0} is exactly the complex for [m-1]. Every (n-1)-simplex in H lies in two top simplices of H if interior, one if equatorial.

**Degree convention (important).** We work in the one-hot quotient: variables y_{v,l} = [lambda(v) = l], with "each vertex has exactly one label" taken for free. A monomial is a partial labeling. *Degree d* = the certificate only uses partial labelings of <= d vertices. This equals Sherali-Adams level d. All degrees below are in this convention.

**Certificate (Nullstellensatz over a field F).** A degree-d certificate of Tucker is an identity, as functions on ALL labelings of the free vertices,
    1 = sum_{alpha violating, |alpha| <= d} c_alpha * 1_alpha
where alpha ranges over partial labelings containing a complementary pair and 1_alpha is its indicator. Size-d violating alpha suffice (smaller ones are sums of larger). For Ky Fan replace 1 by A_+ + 1 where A_+ = number of positively alternating top simplices (a sum of size-(n+1) indicators); this only makes sense over F_2.

**Pseudo-solution (the dual).** No degree-d certificate over F exists iff there is a linear functional E on partial labelings of size <= d with E[empty] = 1, E[alpha] = 0 on violating alpha, and consistent under extension: for every beta with |beta| < d and every vertex v not in beta, sum_l E[beta + (v->l)] = E[beta]. Equivalently, a global function mu on full labelings with odd total mass whose size-d violating partial labelings are all covered an even number of times. **Consistency equations are only needed at level d-1** (they imply the lower levels).

## 2. Settled results (do not recompute except as a sanity check)

| statement | field | degree |
|---|---|---|
| Tucker, S^2 (m=3, labels +-1,+-2) | F_2 | exactly 3 |
| Tucker, S^2 | Q | exactly 5 (fails at 3 and 4) |
| Tucker, S^3 (m=4, labels +-1..+-3) | F_2 | exactly 4 (fails at 3; tower gives 4) |
| Tucker, S^n | F_2 | = n+1 for all n (theorem: tower upper bound + restriction lemma + tree gadget, below) |
| Ky Fan, S^2, labels +-1..+-3 | F_2 | exactly 3 |
| Ky Fan, S^n, k >= n+1 | F_2 | exactly n+1 (theorem, below: tower + separating functional) |

**Theorem (tower).** Define the local lemma on an (n+1)-tuple of labels:
    g_n(l) = sum_i [l minus l_i is positively alternating] + [l positively alternating] + [l negatively alternating]   (mod 2).
Then g_n = 0 on every tuple with no complementary pair (verified exhaustively for n <= 3, k <= 5; proof is a case analysis). Handshake on H and telescoping to A^{(0)} = 1 gives
    A^{(m)} + 1 = sum_{j=1..m} sum_{sigma in H^{(j)}} g_{j-1}(lambda|_sigma)
as an identity on all labelings, where H^{(j)} is the hemisphere of the [j]-complex (padded with zeros). Right side has degree m = n+1. Verified numerically for m = 3,4,5 (`kyfan/tower.py`, `tests/test_tower.py`). Tucker follows because A_+ = 0 syntactically when k <= n.

**Theorem (Ky Fan degree).** Over F_2, the degree is exactly n+1 for all n >= 1, k >= n+1. Upper bound: the tower. Lower bound: the functional E(f) = sum_{lambda(x_i) in S_i} f(lambda) over one top simplex sigma_0 = {x_1 < ... < x_{n+1}} with S_1 = {+1,+2}, S_i = {+i,-i} for i = 2..n+1, all other vertices fixed. Even |S_i| makes E vanish on every function of <= n vertices (including constants, so E(1) = 0); only +-sigma_0 contribute to E(A_+) (a chain other than +-sigma_0 cannot use exactly the reps of sigma_0, since consecutive chain elements cannot be flipped independently), and the collision at magnitude 2 leaves exactly one alternating tuple in the box, so E(A_+ + 1) = 1. Hence A_+ + 1 is not in the span of size-<= n indicators. Code: `kyfan/lower.py`, `tests/test_tower.py` (verified m = 3,4,5; the naive negation-closed box S_1 = {+1,-1} gives E(A_+) = 0, which is why the collision is needed).

**Theorem (restriction lemma).** If rho is a non-violating partial labeling of the free vertices of S^n with unfixed set U, every degree-d certificate restricts (substitute rho) to a degree-<= d certificate of the residual CSP on U (domains D(u) = labels minus {-rho(w) : w fixed, adjacent to u}; pairwise non-complementary on edges). So a consistent degree-d residual dual proves the sphere degree is > d.

**Theorem (Tucker degree).** Tucker's F_2 degree on S^n is exactly n+1 for every n >= 2. Upper bound: the tower (A_+ = 0 when k <= n). Lower bound: the restriction lemma with U = one top simplex through the pole e_{n+1}. Two parts, both proved:
- **(i)** For every binary conflict tree the tree gadget (pole {+n}, leaves with the signed root-to-leaf paths ∪ {-n}) has F_2 degree n+1 — explicit pseudo-solution (the block rule) and consistency via the extension lemma. `TREE_GADGET_THEOREM.md`, `kyfan/tree_rule.py`, `tests/test_rule_vs_solver.py`, `tests/test_extension_lemma.py`.
- **(ii)** The explicit equatorial labeling lambda(u) = s*(n-q+1), where p = last nonzero index of u, s = u_p, and q = start of the final maximal same-sign run ending at p, is valid (no complementary comparable pair) and its pole-chain residual is exactly the reverse-caterpillar tree. No search. `REALIZATION_THEOREM.md`, `kyfan/realize_explicit.py`, `tests/test_realize_explicit.py`.

Mechanically checked: (i) all tree shapes n <= 6 (and n = 7 sparse); (ii) validity + residual domains m <= 10 and the residual degree-n dual consistent+unique m <= 8. Also verified for n = 2, 3 by the direct dual (`nsdeg`, the S^3 Z_3xZ_3 dual) and for n = 4..7 by the searched `GADGETS` (`kyfan.gadget.verify`) and the deterministic `realize` (n <= 8). The Ky Fan functional cannot be reused here (it kills constants; Tucker's certificate is for the constant 1) — the tree gadget is what supplies a pseudo-solution with E[empty] = 1.

**Structural facts about the degree-2 pseudo-solution on S^2 (`tests/test_pseudo2.py`):**
- No pseudo-solution is invariant under the full symmetry group, nor under the label group (signed permutations of magnitudes), nor under the stabilizer of a hemisphere.
- No pseudo-solution is supported only on distinct-magnitude pairs.
- Necessary pair types (removing any one kills all solutions): edges with two different magnitudes, at every rank pair; and two unrelated rank-1 vertices (singletons e_i, e_j) with equal labels.
- The unrestricted solution space has 597 free parameters (1104 unknowns, rank 507).

**Negative/methodological facts you must respect:**
- The unsubdivided octahedron (cross-polytope boundary without subdivision) is degenerate: A_+ is identically 1. Only the subdivided complex above tests anything.
- Averaging a certificate over a symmetry group is valid only when the group order is invertible in the field. Over F_2, use only odd-order subgroups (e.g. Z_3 x Z_3 from a coordinate 3-cycle and a magnitude 3-cycle). Over Q the full group is fine. Fully-symmetric F_2 certificates do NOT exist even one degree above the minimum, so a "no symmetric certificate" over F_2 means nothing.
- The random-point primal (rows = evaluation on random labelings) is one-way sound: an inconsistent sampled system proves no certificate; a consistent one must be verified on fresh labelings.
- Floating-point rank over R is an acceptable proxy for Q only with a residual check (we used relative residual < 1e-12 as "certificate", > 1e-3 as "none").
- Pullback of a pseudo-solution from the equator to S^n loses exactly one degree (a violating n-window can project to n distinct equatorial vertices). Do not expect a naive induction to work.

## 3. Files

Package `kyfan/` (parametrized by m; `from kyfan import SignedComplex, label_set, ...`):
- `complex.py` — `SignedComplex(m)`: verts, free reps (`rep(v)` = (index, sign)), `label(v, L)`, edges, top chains, `hemisphere_top(j)`.
- `labels.py` — `label_set(k)`, `posalt`/`negalt`, `pos_alt_count` (A_+), DFS over valid labelings, octahedron demo.
- `violating.py` — `violating_pairs`, `is_violating`, `partial_labelings(n, labels, d)` (sorted tuples of (free index, label)).
- `group.py` — elements (pi, eps, lpi, leps), `act`, `compose`, `generated_group`, `assert_odd_order`; generators for the full group, hemisphere stabilizer, label group, and `z3z3_generators` (the only F_2-valid averaging); `orbits`.
- `dual.py` — `sa_dual_system(cx, labels, d, gens=, support=)`: unknowns = non-violating size-d partial labelings (orbit-quotiented / support-restricted), consistency rows at level d-1 only, plus E[empty]=1; `solve(rows, ncols, 'F2' | p)`.
- `primal.py` — random-point primal (`sampled_primal`, `FastRows` for large vectorized systems), `verify_random`, `verify_exhaustive_F2`, `verify_exact_F2` (canonical one-hot basis), `kyfan_target_monomials`.
- `linalg.py` — `gf2_sparse` (wrapper for the Rust solver in `gf2solve/`, Task 1: 1.8M unknowns in ~20 s; returns rank, consistency, particular solution, optional null-space basis), dense bit-packed GF(2) Gauss-Jordan (`gf2_dense`, `gf2_dense_packed`; for dense primal systems), sparse mod-p elimination, python-int bitset solver, dense mod-p, real lstsq residual.
- `gf2solve/` — Rust crate (`cargo build --release`; the wrapper builds it on first use). Input/output format in `src/main.rs`.
- `tower.py` — local lemma `g`, `local_lemma_violations`, telescoped identity `check_identity`.
- `lower.py` — the Ky Fan separating functional (`box`, `E`, `check`).
- `ball.py` — Tucker's ball form (`Ball`) and general restrictions (`Residual`), domain-aware SA dual.
- `gadget.py` — restriction lemma + chain gadgets: `GADGETS[m]` (searched equatorial labelings, m = 4..8), `chain_U`, `verify(m)` (independent end-to-end check); statement (ii): `chain_domains` (fast, no full sphere), `caterpillar_tree`, `realize(m)` (deterministic construction), `check_realization`.
- `abstract_gadget.py` — tree gadgets without the sphere: `build_dual`, `degree(domains)`.
- `tree_rule.py` — statement (i): the explicit block-rule pseudo-solution `rule_E`, the block rule `rule_patterns`, and the solver ground truth `unique_E`.
- `realize_explicit.py` — statement (ii): the closed-form labeling `explicit_label`/`explicit_Leq` and `check(m)` (validity, reverse-caterpillar domains, residual dual).
- `analysis.py` — Task 2 helpers (orbit types, knockouts, descriptions). Exploratory drivers and their outputs live in `analysis/`.

`tests/` are the settled table (see results.md for the test-to-number map); `./reproduce.sh` (fast, ~2 min) or `./reproduce.sh all` (~15 min).
The original `scripts/` were reviewed in `review.md`, reproduced, and deleted after the tests passed (git history: commit 5a4dd50).


## 3a. Phase 0 — review and reproduce before anything else

**Status: completed 2026-09-16** (see review.md, results.md). Kept for the record.

Do this first, and do not start Task 1 until every check passes.

1. **Read every script in `scripts/` end to end.** For each, write one line in `review.md`: what it computes, which of the settled results it supports, and any fragility (the `exec()` chains, hard-coded `m=3`, the `[::-1]` bit-packing in `s3.py`, sampled-primal assumptions). Flag anything that looks like a bug — a suspected bug in these scripts is a finding, not an inconvenience; report it before "fixing" it, because a fix could change a number in the table.
2. **Run `./reproduce.sh`** and diff its output against `results.md`. Then run `primal.py`, then `z9a.py && z9b.py`. Every number in the settled table must reproduce. If one does not, stop and report: which number, what you got, and your best guess why. Do not proceed on a table you cannot reproduce.
3. **Check the conventions against the code**, not just the prose: confirm the degree is partial-labeling size (count vertices in a monomial, not polynomial degree), confirm "positively alternating" is what `pos_alt` computes, confirm free-vertex representatives are "first nonzero coordinate is +1", confirm consistency equations are built only where needed.
4. **Refactor into a package** `kyfan/` (complex(m), labelings, violating pairs, group actions with an explicit odd-order-subgroup helper, SA-dual builder, solvers) with a `tests/` directory whose tests *are* the settled table. Keep `scripts/` untouched as the reference until the tests pass, then delete the duplicated logic.
5. Only then, Task 1.

## 4. Tasks, in priority order

### Task 1 — A real sparse GF(2) solver  (DONE 2026-09-16: `gf2solve/`, see results.md)
The bottleneck everywhere is linear algebra over F_2 at 10^5–10^6 unknowns. The SA-dual system is sparse (~2*(2k) nonzeros per consistency equation). Options, try in this order:
1. Sparse Gaussian elimination with Markowitz pivoting, written in Rust or C++ (bitset rows are wrong here — use sorted index lists or hash sets; fill-in is the enemy, so pivot on low-degree columns/rows first).
2. SageMath if installed (`matrix(GF(2), ..., sparse=True)`), or M4RI/M4RIE bindings for anything that fits dense in memory (dense is fine up to ~100k x 100k on 32 GB).
3. Block Wiedemann only if 1–2 fail; probably unnecessary.
Validate on the m=3 systems (must reproduce: F_2 degree 2 no, 3 yes) and on the z9 system (must reproduce: inconsistent).

### Task 2 — Extract and analyze the degree-3 pseudo-solution on S^3  (DONE 2026-09-16: see results.md "Task 2", analysis/)
This is the object the conjecture's proof has to generalize from. Build the degree-3 SA-dual for m=4, labels +-1..+-3: ~1.8M non-violating size-3 partial labelings before symmetry; average over the odd subgroup Z_3 x Z_3 (coordinate 3-cycle on coords 0,1,2; magnitude 3-cycle 1->2->3->1) to cut ~9x. Solve. Then:
- Repeat the knockout analysis of pseudo2c.py at the level of full-group orbit *types* of size-3 partial labelings: which types are necessary?
- Test the same restrictions as at n=2 (`tests/test_pseudo2.py`): label-symmetric (expect impossible), hemisphere-stabilizer (expect impossible), distinct-magnitudes-only (expect impossible).
- Look for a sparse or structured solution (e.g. minimize support greedily, or impose support on "base labeling + corrections": pick a labeling L_0 with exactly one antipodal pair of complementary edges — pull back a valid equatorial labeling that uses label -1 exactly once, and put +1 on e_m — and search for a solution supported on labelings within small Hamming distance of L_0).
Report whatever pattern is or isn't there. A pattern that persists from n=2 to n=3 is a conjecture to prove; its absence is evidence against the n+1 conjecture.

### Task 3 — The ball version  (DONE 2026-09-16: tight over F_2 for m = 3, 4; and the restriction idea gives S^4, S^5, S^6 — see results.md)
Fix a valid labeling L_eq on the equator (m-1 complex, labels +-1..+-(m-1)) and consider only the cap {x_m = +1} (and its antipode). Constraints: no complementary edge inside the cap, and none between a cap vertex y and an equatorial vertex w <= pi(y) where pi zeroes the last coordinate. This is Tucker's ball form; it is UNSAT, and any degree-d sphere certificate restricts to a degree-<=d certificate here. Compute its F_2 degree for m = 3, 4 (fewer variables than the sphere). If it equals the sphere degree, the reduction is tight and the ball is the right object for a proof; note that with the pullback labeling every cap vertex except e_m is fine and e_m is adjacent to the entire cap.

### Task 4 — results.md
Every number with the exact command and runtime that produced it, plus a short "what this shows" line. This file is half the value of the repo.

### Optional (only after 1–4)
- Sherali-Adams LP (real, nonnegative) feasibility at degree n for Tucker on S^n, m = 3, 4 — does an actual local distribution exist, or only a signed one?
- Kneser transfer: substitute a color variable for each label variable via S_x = lex-first r-subset inside x^+ or x^- (this makes the substitution degree 1) and verify the certificate on random colorings of KG(6,2) (m=6, r=2; elementary case, plumbing test) — deprioritized; do not build first.

## 5. Working style
- Every claimed degree comes with both directions: a certificate verified on fresh random labelings, or an inconsistent (sub)system.
- Never conclude from a full-group-symmetric F_2 search.
- Log runtimes and memory; if a solve exceeds ~1 hour, stop and report the size rather than waiting.
- When a construction fails, record *why* in results.md — the failure modes so far (pullback loses a degree; no symmetric solutions) are as useful as the numbers.
