# Size lower bounds for Nullstellensatz certificates of Tucker on S^n

Companion to `TREE_GADGET_THEOREM.md` (statement (i)) and `REALIZATION_THEOREM.md` (statement (ii)). Everything here
is over F_2 in the one-hot quotient: a certificate is an identity of functions on all labelings of the free vertices,

    1 = Σ_{α ∈ C} 1_α ,

C a set of violating partial labelings (each containing a complementary pair), and **size(C) = |C|**, the number of
terms. **Degree** = max |α|. A *flag* is a top simplex (maximal chain, ranks 1..m = n+1); #flags = m!·2^m. With
N = 3^m − 1 vertices, m!2^m = N^{Θ(log log N)}. "full(U)" means a labeling of every vertex of the flag U.

## 0. Ingredients (proved elsewhere in this repo)

- **Restriction lemma** (`kyfan/gadget.py`). If ρ is a non-violating partial labeling of all free vertices outside a
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
dangerous: the explicit labeling has a 1/(2n) fraction of all flags inside its magnitude-1 class
(`gadget_spread_experiments.py`, `flag_density`), so a single monochromatic term serves ~#flags/(2m) flags.

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

**Evidence** (`gadget_spread_experiments.py`, Glauber dynamics over valid labelings satisfying the tree-domain
constraints, started from the explicit gadget; MCMC ergodicity not verified, samples correlated):
- largest-magnitude class fraction: 0.44, 0.36, 0.29, 0.23 for m = 4..7 (explicit gadget: ≈ 0.45 for all m);
- flag density of the largest magnitude class: 0.104, 0.027, 0.005, 0.001 for m = 4..7 (explicit: 1/(2n));
- match probability of a fresh gadget against a fixed labeling on A, m = 6: |A| = 1,2,3,4,6,8 →
  random sets 0.39, 0.19, 0.06, 0.013, 0.004, 0.001; chains 0.90, 0.47, 0.14, 0.07; balls 0.37, 0.17, 0.09, 0.03, 0.009, 0.003.
  Geometric decay with β ≈ 0.5–0.6 per vertex on all three set types, consistent with the lemma with a small c.
- The explicit gadget violates the lemma (polynomial flag density), so D_F must be genuinely random; the uniform
  distribution over gadgets is the natural candidate. No proof is known; the obvious constructed distributions
  (symmetry transports) carry only O(n log n) bits and cannot suffice.

## 5. Summary

| statement | status |
|---|---|
| size ≥ m!2^m for degree-(n+1) certificates | proved (Thm 1) |
| size ≥ m!2^m / C(D,m) for degree ≤ D; superpoly for D = O(n) | proved (Thm 2) |
| size ≥ m!2^{m−1}/(2^{n−μ}(n−μ)!) if terms use ≥ μ magnitudes off served flags | proved (Thm 3) |
| size ≥ m!2^m/K^m for all degrees | conditional on the Spread Lemma (Thm 4) |
| Spread Lemma | conjectured; empirical support m ≤ 7 |

Upper bound for comparison: the tower certificate has Σ_j j!2^{j−1} ≈ m!2^{m−1} local-lemma instances, each expanding
to ≤ (2k)^m monomials — N^{Θ(log log N)} as well. So Theorem 1 is tight up to the per-lemma factor for minimal-degree
certificates, and the open question is whether *high* degree can buy a polynomial-size certificate; the Spread Lemma
says it cannot.
