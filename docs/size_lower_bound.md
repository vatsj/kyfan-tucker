# Size lower bounds for Nullstellensatz certificates of Tucker on S^n

Companion to `docs/tree_gadget.md` (statement (i)) and `docs/realization.md` (statement (ii)). Computations:
`analysis/isd.py` (minimum certificate size, `tests/test_min_size.py`, `tests/test_paper_numbers.py`),
`analysis/flag_hitting_bound.py` (Theorem 1′, `tests/test_flag_hitting.py`, `tests/test_paper_numbers.py`),
`analysis/spread_mcmc.py` (Spread Lemma evidence, `docs/spread_evidence.md`, `tests/test_spread_mcmc.py`).
Everything here is over F_2 in the one-hot quotient: a certificate is an identity of functions on all labelings of the free vertices,

    1 = Σ_{α ∈ C} 1_α ,

C a set of violating partial labelings (each containing a complementary pair), and **size(C) = |C|**, the number of
terms. **Degree** = max |α|. A *flag* is a top simplex (maximal chain, ranks 1..m = n+1); #flags = m!·2^m. With
N = 3^m − 1 vertices, m!2^m = N^{Θ(log log N)}. "full(U)" means a labeling of every vertex of the flag U.

## 0. Ingredients (proved elsewhere in this repo)

- **Restriction lemma** (`docs/tree_gadget.md`; code `kyfan/gadget.py`). If ρ is a non-violating partial labeling of all free vertices outside a
  flag U, then substituting ρ into a certificate C gives a certificate of the residual CSP on U whose terms are
  { α|_U : α ∈ C, α|_{off U} ⊆ ρ }.
- **Tree theorem** (i). The residual CSP has F_2 degree n+1 whenever its domains contain the domains of a binary
  conflict tree gadget; this is label-agnostic (any injective assignment of magnitudes to internal nodes, any sign
  conventions), because the block-rule proof uses only distinctness of magnitudes and the sign split at each node.
- **Realization** (ii). The explicit labeling λ(u) = s·(n−q+1) is a gadget ρ_{U_0} for the pole flag U_0.
- **Transport.** The hyperoctahedral group B_m acts simply transitively on flags (they are the chambers of the type-B
  Coxeter complex; |B_m| = m!2^m = #flags) and preserves the order relation, so g·ρ_{U_0} is a gadget for g(U_0).
  Hence **every flag U has a gadget ρ_U**.

**Lemma A (every flag forces a term).** For every flag U and gadget ρ_U, every certificate C contains a term α with
α|_U = full(U) and α|_{off U} ⊆ ρ_U.
*Proof.* The residual certificate has degree ≥ n+1 by (i)+(ii), so it contains a term of size n+1 = |U|. ∎

## 1. Theorem 1 (proved): minimal-degree certificates

Every certificate of degree exactly n+1 has size ≥ #flags = m!·2^m.
*Proof.* The term α given by Lemma A has |α| ≤ n+1 and α|_U = full(U), so α is exactly a full labeling of U.
Distinct flags have distinct vertex sets, so the assigned terms are distinct. ∎

## 1′. Theorem 1′ (proved): per-flag parity constraints

Fix a flag U (identified with its antipode: both use the same free vertices, so a *flag class* is an antipodal pair;
there are m!·2^{m−1} classes). A *gadget* for U is any non-violating partial labeling ρ of the free vertices outside U
whose residual CSP on U has no degree-n certificate; let Z_ρ be the space of vectors e on the residual's full
assignments with Σ_{x ⊇ β} e(x) = 0 for every violating β of size ≤ n (pseudo-solutions are the e ∈ Z_ρ with
Σ_x e(x) = 1; when the residual is a tree gadget, Z_ρ is one-dimensional, spanned by μ_ρ).

**Lemma B.** Let C be a certificate of degree exactly n+1 and T_U = {α ∈ C : supp(α) = U}. Then for every gadget ρ of U
and every e ∈ Z_ρ,   Σ_{x ∈ T_U} e(x) = Σ_x e(x).   In particular |T_U ∩ supp μ_ρ| is odd for every tree-gadget ρ.
*Proof.* Every α ∈ T_U has empty off-U part, so it survives the restriction by ρ unchanged; every other surviving term
has size ≤ n and is violating in the residual. Hence Σ_{α ∈ T_U} 1_α + 1 lies in the span of the violating
low-degree indicators, and pairing with e ∈ Z_ρ (which kills that span) gives the identity. ∎

**Theorem 1′.** size(C) ≥ Σ_{flag classes U} t(U), where t(U) = min{ |T| : T ⊆ full labelings of U satisfying Lemma B
for every gadget ρ of U }. On S^2 (m = 3): at one flag there are 3,568 non-violating restrictions, 16 gadget domain types,
each with Z_ρ one-dimensional and |supp μ_ρ| = 3; the 16 supports touch 32 of the 64 full labelings with maximum
multiplicity 2, and the exact minimum odd-hitting set has size **t = 8** (branch and bound, `analysis/flag_hitting_bound_s2.py`).
Hence every degree-3 certificate of Tucker on S^2 has ≥ 24 · 8 = **192** size-3 terms. (The true minimum, 304 = 288 + 16,
has exactly 12 full-flag terms on every class; those 12 satisfy all 16 parities; the gap 8 → 12 is cross-flag coupling
through the shared size-2 terms, invisible to a per-flag bound.)

Remark. Restricting to the label-group orbit of the explicit gadget gives a computable sub-family at any flag, hence
a valid lower bound on t(U) for every n: 2^n n! tree gadgets with supports of size 3^{n−1} inside (2n)^{n+1} labelings.

**Computation (`analysis/flag_hitting_bound.py`, outputs `analysis/flag_hitting_bound_m{3,4}.out`;
`tests/test_flag_hitting.py`).** Since every gadget is transported to every flag class by the hyperoctahedral group,
t(U) is the same for all classes; the pole chain is used. Gadget families (each a subset of the true family, so each
gives a valid lower bound on t; the parity system is solved exactly by CP-SAT, cross-checked by branch and bound):

| S^n | family | gadget domain tuples | t (exact for the family) | full-flag terms ≥ classes × t |
|---|---|---|---|---|
| S^2 | label orbit of the explicit gadget | 8 | 4 | 24 × 4 = 96 |
| S^2 | label orbit × realizable chain permutations | 16 = all gadget types | **8** | 24 × 8 = **192** |
| S^2 | all 3,568 restrictions (exhaustive, independent coordinates) | 16 | 8 | 192 |
| S^3 | label orbit of the explicit gadget | 48 (supports of size 9 in 1296) | 6 | 192 × 6 = 1,152 |
| S^3 | all 4 tree shapes × 6 leaf orders realized by `realize` (16 base gadgets), label orbits | 96 | 6 | 1,152 |
| S^3 | label orbit × realizable chain permutations | 192 (32 realizable of 96 permutations) | **12** | 192 × 12 = **2,304** |

"Realizable chain permutations": the residual CSP on a chain is symmetric under permuting the n+1 positions, so a
permuted domain tuple is again a (tree) gadget *if* some restriction produces it; realizability is decided by CP-SAT
(a valid labeling of the other free vertices whose forbidden sets are exactly the complements of the target domains;
the returned labeling is re-checked independently). On S^2 the realizable permutations of the reverse caterpillar are
exactly the four keeping the singleton domain at an end of the chain (identity, swap of the other two, and their
reversals), and their label orbits are exactly the 16 exhaustive gadget types. On S^3, 32 of the 96 position
permutations of the four tree shapes are realizable (8 per shape, all with the singleton domain at an end of the
chain), and they double t from 6 to 12.
So on S^3 every degree-4 certificate has ≥ 2,304 full-flag terms, versus Theorem 1's count of 384 (= #flags); the true
t(U) may be larger still (the exhaustive family is out of reach: 6^36 restrictions).

Note on Theorem 1's count: a full labeling of U is the same term as a full labeling of −U (same free vertices), so the
"distinct flags have distinct terms" argument gives #flags/2 = m!·2^{m−1} classes, and #flags = m!·2^m needs t(U) ≥ 2,
which Theorem 1′ supplies with room to spare (t ≥ 8 on S^2, ≥ 12 on S^3). The same antipodal identification affects
Theorem 2 as stated below: counting flag classes instead of flags, its unconditional form is size ≥ m!·2^{m−1} / C(D, m)
(the asymptotic statement N^{Θ(log log N)} is unaffected).

## 2. Theorem 2 (proved): certificates of degree O(n)

Every certificate of degree ≤ D has size ≥ m!·2^m / C(D, m). In particular, for D = cm (c ≥ 1 constant),

    size ≥ m! · (2/(e·c))^m = N^{Θ(log log N)}.

*Proof.* Assign to each flag U a term α_U from Lemma A. A term α is assigned to at most the number of flags whose
vertex set lies inside supp(α), which is ≤ C(|supp α|, m) ≤ C(D, m) since a flag is determined by its m vertices.
So #terms ≥ #flags / C(D, m); and C(cm, m) ≤ (ec)^m. ∎

Remark: this uses only *existence* of a gadget at every flag. It is unconditional now that (ii) is proved.

## 3. Theorem 3 (proved, but conditional on a property of the certificate): label randomization

Let G_L be the group of signed permutations of the magnitudes [n], |G_L| = 2^n n!. For a flag F let D_F be the uniform
distribution on { π·ρ_F : π ∈ G_L }; each element is a gadget (tree theorem, label-agnostic). For a term α and a flag F,
say α *serves* F if α|_F = full(F) and α|_{off F} ⊆ ρ for some ρ in supp D_F.

(a) *Every complementary pair of α meets every flag it serves* (else α|_{off F} would be violating), so α serves at
most 2·(m−1)!2^{m−1} = #flags/m flags (flags through a fixed vertex of rank r number r!(m−r)!2^{m−r} ≤ (m−1)!2^{m−1}).

(b) If α uses μ distinct magnitudes on supp(α) \ F, then P_{ρ∼D_F}[α|_{off F} ⊆ ρ] ≤ 2^{n−μ}(n−μ)! / (2^n n!)
(π is forced on the μ magnitudes that must be matched; the rest is free).

**Theorem 3.** If every term of C uses at least μ distinct magnitudes outside each flag it serves, then

    size(C) ≥ m!·2^{m−1} / (2^{n−μ} (n−μ)!) .

*Proof.* Pick U uniformly among flags and ρ ∼ D_U. By Lemma A some α ∈ C has α|_U = full(U), α|_{off U} ⊆ ρ, so
Σ_{α∈C} P[α survives] ≥ 1, and P[α survives] = (1/#flags) Σ_{F served by α} P_{D_F}[α|_{off F} ⊆ ρ] ≤
(1/#flags)(#flags/m)·2^{n−μ}(n−μ)!/(2^n n!). Invert, using m·2^n n! = m!·2^{m−1}. ∎

This is superpolynomial as soon as μ → ∞ with m. It fails exactly for terms that are large but use O(1) magnitudes
outside the flags they serve — the "monochromatic region" terms. Against the *explicit* gadget those terms are
dangerous: the explicit labeling has a 1/(2n) fraction of all flags inside its largest-magnitude class
(`analysis/spread_mcmc.py`, `flag_density`; see `docs/spread_evidence.md`), so a single monochromatic term serves
~#flags/(2m) flags.

## 4. The spread lemma (conjecture), and Theorem 4 (proved conditional on it)

**Spread Lemma (conjectured).** There exist constants β < 1 and c ≥ 0, and for every flag F a distribution D_F on
gadget labelings for F, such that for every set A of vertices disjoint from F and every labeling a of A,

    P_{ρ ∼ D_F} [ ρ|_A = a ] ≤ β^{|A| − c·m} .

(The additive c·m allows a rigid zone of O(m) vertices near F — the vertices comparable to σ — where the labels are
forced by the domain conditions.)

**Theorem 4.** Under the Spread Lemma, every certificate (of any degree) has size ≥ m!·2^m / K^m for a constant
K = K(β, c), i.e. size = N^{Θ(log log N)}.
*Proof.* As in Theorem 3, Σ_α P[α survives] ≥ 1 with
P[α survives] ≤ (1/#flags) Σ_{F ⊆ supp α} β^{|supp α| − m − cm} ≤ (1/#flags)·C(s, m)·β^{s−m}·β^{−cm}, s = |supp α|.
Writing s = m + r and x = r/m: C(m+r, m) β^r ≤ (e(1+x)β^x)^m ≤ K_0(β)^m, where K_0 = max_{x≥0} e(1+x)β^x < ∞.
So P[α survives] ≤ K_0^m β^{−cm}/#flags and size ≥ #flags/(K_0 β^{−c})^m. ∎

**Evidence** (`analysis/spread_mcmc.py`, output `analysis/spread_mcmc.out`; Glauber dynamics on the valid equatorial
labelings that agree with the explicit gadget on the pole chain and whose residual domains contain the reverse-caterpillar
domains, started from the explicit gadget; ergodicity of the dynamics is not proved. Full description, diagnostics and caveats
in `docs/spread_evidence.md`):
- largest-magnitude class fraction: 0.461, 0.371, 0.287, 0.226, 0.183 for m = 4..8 (explicit gadget: ≈ 0.45 for all m);
- flag density of the largest magnitude class: 0.094, 0.024, 0.0051, 0.0009, 0.0001 for m = 4..8 (explicit: 1/(2n) = 0.167,
  0.125, 0.100, 0.083, 0.071);
- match probability of two independent gadget samples on a set A, m = 6, random A with |A| = 1, 2, 3, 4, 6, 8:
  0.39, 0.16, 0.046, 0.014, 0.0033, 0.0002; ball-shaped A: 0.47, 0.21, 0.14, 0.048, 0.011, 0.0019.
  Geometric decay on all three set types (random, chain, ball), with fitted per-vertex β decreasing in m
  (random sets: 0.65, 0.48, 0.36, 0.29, 0.22 for m = 4..8), consistent with the lemma with a small c.
- The explicit gadget violates the lemma (polynomial flag density), so D_F must be genuinely random; the uniform
  distribution over gadgets is the natural candidate. No proof is known; the obvious constructed distributions
  (symmetry transports) carry only O(n log n) bits and cannot suffice.

## 5. Summary

| statement | status |
|---|---|
| size ≥ m!2^m for degree-(n+1) certificates | proved (Thm 1) |
| size ≥ Σ_classes t(U) for degree-(n+1) certificates; t = 8 on S^2 (exact), ≥ 12 on S^3 | proved (Thm 1′), computed |
| size ≥ m!2^m / C(D,m) for degree ≤ D; superpoly for D = O(n) | proved (Thm 2) |
| size ≥ m!2^{m−1}/(2^{n−μ}(n−μ)!) if terms use ≥ μ magnitudes off served flags | proved (Thm 3) |
| size ≥ m!2^m/K^m for all degrees | conditional on the Spread Lemma (Thm 4) |
| Spread Lemma | conjectured; empirical support m ≤ 8 (`docs/spread_evidence.md`) |

Upper bound for comparison: the tower certificate has Σ_j j!2^{j−1} ≈ m!2^{m−1} local-lemma instances, each expanding
to ≤ (2k)^m monomials — N^{Θ(log log N)} as well. So Theorem 1 is tight up to the per-lemma factor for minimal-degree
certificates, and the open question is whether *high* degree can buy a polynomial-size certificate; the Spread Lemma
says it cannot.
