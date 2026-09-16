# Brief for Claude Code: Nullstellensatz degree of Tucker / Ky Fan over F_2

You are continuing a research computation. Everything you need is in this repository. `results.md` is the ledger of every settled number; `scripts/` are the computations that produced them. Read this file fully before running anything. Prior results are settled — do not re-derive them; extend them.

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
| Ky Fan, S^2, labels +-1..+-3 | F_2 | exactly 3 |
| Ky Fan, any n | F_2 | exactly n+1 (theorem, below) |

**Theorem (tower).** Define the local lemma on an (n+1)-tuple of labels:
    g_n(l) = sum_i [l minus l_i is positively alternating] + [l positively alternating] + [l negatively alternating]   (mod 2).
Then g_n = 0 on every tuple with no complementary pair (verified exhaustively for n <= 3, k <= 5; proof is a case analysis). Handshake on H and telescoping to A^{(0)} = 1 gives
    A^{(m)} + 1 = sum_{j=1..m} sum_{sigma in H^{(j)}} g_{j-1}(lambda|_sigma)
as an identity on all labelings, where H^{(j)} is the hemisphere of the [j]-complex (padded with zeros). Right side has degree m = n+1. Verified numerically for m = 3,4,5 (tower.py). Tucker follows because A_+ = 0 syntactically when k <= n.

**Conjecture (open).** Tucker's F_2 degree on S^n is exactly n+1 for all n. Upper bound is the theorem; lower bound is known only for n = 2, 3.

**Structural facts about the degree-2 pseudo-solution on S^2 (pseudo2b.py, pseudo2c.py):**
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

- `scripts/build.py` — complex for m=3, brute-force Tucker/Ky Fan checks.
- `scripts/nsdeg.py` — SA-dual system, NS degree over F_2 and F_p, m=3 (dense bit-packed GF(2) elimination in numpy; fine to ~40k unknowns, useless beyond).
- `scripts/primal.py` — explicit certificates via random-point primal; full-group orbit reduction (F_p only).
- `scripts/s3.py`, `z9a.py`, `z9b.py` — m=4: symmetric search (uninformative over F_2) and the Z_3xZ_3-averaged unrestricted F_2 search (this is what proved degree >= 4 on S^3). z9b is checkpointed elimination from a sandbox with a 300 s limit; on a laptop just run it straight.
- `scripts/tower.py` — the theorem: exhaustive local-lemma check and the telescoped identity for m = 3,4,5.
- `scripts/pseudo2.py`, `pseudo2b.py`, `pseudo2c.py` — degree-2 pseudo-solution structure on S^2. These use `exec()` on other files; refactor before extending.

(Refactoring is step 4 of Phase 0 below.) The plan: refactor everything into one module `kyfan/` (complex construction parametrized by m, labeling utilities, violating pairs, group action with odd-order subgroups, SA-dual system builder, sparse GF(2)/GF(p) solvers) with a test that reproduces every number in the table above. Keep the exec-based scripts as reference until the tests pass.


## 3a. Phase 0 — review and reproduce before anything else

Do this first, and do not start Task 1 until every check passes.

1. **Read every script in `scripts/` end to end.** For each, write one line in `review.md`: what it computes, which of the settled results it supports, and any fragility (the `exec()` chains, hard-coded `m=3`, the `[::-1]` bit-packing in `s3.py`, sampled-primal assumptions). Flag anything that looks like a bug — a suspected bug in these scripts is a finding, not an inconvenience; report it before "fixing" it, because a fix could change a number in the table.
2. **Run `./reproduce.sh`** and diff its output against `results.md`. Then run `primal.py`, then `z9a.py && z9b.py`. Every number in the settled table must reproduce. If one does not, stop and report: which number, what you got, and your best guess why. Do not proceed on a table you cannot reproduce.
3. **Check the conventions against the code**, not just the prose: confirm the degree is partial-labeling size (count vertices in a monomial, not polynomial degree), confirm "positively alternating" is what `pos_alt` computes, confirm free-vertex representatives are "first nonzero coordinate is +1", confirm consistency equations are built only where needed.
4. **Refactor into a package** `kyfan/` (complex(m), labelings, violating pairs, group actions with an explicit odd-order-subgroup helper, SA-dual builder, solvers) with a `tests/` directory whose tests *are* the settled table. Keep `scripts/` untouched as the reference until the tests pass, then delete the duplicated logic.
5. Only then, Task 1.

## 4. Tasks, in priority order

### Task 1 — A real sparse GF(2) solver
The bottleneck everywhere is linear algebra over F_2 at 10^5–10^6 unknowns. The SA-dual system is sparse (~2*(2k) nonzeros per consistency equation). Options, try in this order:
1. Sparse Gaussian elimination with Markowitz pivoting, written in Rust or C++ (bitset rows are wrong here — use sorted index lists or hash sets; fill-in is the enemy, so pivot on low-degree columns/rows first).
2. SageMath if installed (`matrix(GF(2), ..., sparse=True)`), or M4RI/M4RIE bindings for anything that fits dense in memory (dense is fine up to ~100k x 100k on 32 GB).
3. Block Wiedemann only if 1–2 fail; probably unnecessary.
Validate on the m=3 systems (must reproduce: F_2 degree 2 no, 3 yes) and on the z9 system (must reproduce: inconsistent).

### Task 2 — Extract and analyze the degree-3 pseudo-solution on S^3
This is the object the conjecture's proof has to generalize from. Build the degree-3 SA-dual for m=4, labels +-1..+-3: ~1.8M non-violating size-3 partial labelings before symmetry; average over the odd subgroup Z_3 x Z_3 (coordinate 3-cycle on coords 0,1,2; magnitude 3-cycle 1->2->3->1) to cut ~9x. Solve. Then:
- Repeat the knockout analysis of pseudo2c.py at the level of full-group orbit *types* of size-3 partial labelings: which types are necessary?
- Test the same restrictions as pseudo2b.py: label-symmetric (expect impossible), hemisphere-stabilizer (expect impossible), distinct-magnitudes-only (expect impossible).
- Look for a sparse or structured solution (e.g. minimize support greedily, or impose support on "base labeling + corrections": pick a labeling L_0 with exactly one antipodal pair of complementary edges — pull back a valid equatorial labeling that uses label -1 exactly once, and put +1 on e_m — and search for a solution supported on labelings within small Hamming distance of L_0).
Report whatever pattern is or isn't there. A pattern that persists from n=2 to n=3 is a conjecture to prove; its absence is evidence against the n+1 conjecture.

### Task 3 — The ball version
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
