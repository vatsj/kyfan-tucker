# A weighted spread lemma and an all-degrees size lower bound: size ≥ m!/2

Companion to `docs/size_lower_bound.md` §4. Code: `kyfan/spread.py` (reference implementation of every object below),
`analysis/spread_choice_family.py` (numbers; output `analysis/spread_choice_family.out`). Tests: `tests/test_spread_lemma.py`
(fast tier, 5 s; every lemma below names its test).

**Outcome.** The Spread Lemma *as conjectured* (P[ρ|_A = a] ≤ β^{|A| − cm}) is **not** proved, and it is false for the family
constructed here (that family has m(m−1) deterministic vertices per flag, more than any c·m). What is proved, for all n, is a
*weighted* spread statement for an explicit product family of gadgets (Lemma 5), a counting theorem that absorbs the
Θ(m²) exceptional vertices (Theorem B), and from them the unconditional replacement for Theorem 4:

> **Theorem C.** For every n ≥ 2 (m = n+1), every F_2 Nullstellensatz certificate of Tucker on S^n in the one-hot encoding,
> **of any degree**, has at least **m!/2** terms. With N = 3^m − 1 vertices this is ≥ (N+1)^{log_3 m − 1/ln 3}/2 = N^{Ω(log log N)}.

Status labels used below: **PROVED** (for all n, proof text here), **VERIFIED** (exhaustive or sampled computation, range stated),
**TRUSTED** (proved elsewhere in this repository and used as stated).

Conventions: m = n+1; vertices are the nonzero y ∈ {−1,0,1}^m; x ≤ y iff y agrees with x on supp(x); rank = |supp|; labels ±[n],
antipodal, so a *variable* is an antipodal class [y] = {y, −y} (there are (3^m − 1)/2 of them) and a labeling of a set of classes
is a labeling of its vertices with λ(−y) = −λ(y). A *flag* F = (x_1 < … < x_m) is a maximal chain; F and −F have the same variables
and form a *flag class*; there are m!·2^{m−1} classes. For a flag F and a vertex y put

    d_F(y) = rank(y) − max{ rank(x) : x ∈ F ∪ −F, x ≤ y }        (the maximum of the empty set is 0).

d_F(y) = 0 iff y ∈ ±F; d_F(−y) = d_F(y) = d_{−F}(y); and d_{gF}(gy) = d_F(y) for every signed permutation g of the coordinates
(the definition uses only the order and the antipodal map).

## 1. Choice labelings

Pole flag: U_0 = (v_0 < v_1 < … < v_n), v_0 = e_m, v_r = e_m − (e_1 + … + e_r). (This is the flag `e_{n+1} ∪ lift(σ)` of
`docs/realization.md`.) A **choice function** assigns to every vertex y ≠ ±e_m a coordinate c(y) ∈ supp(y) ∩ [n] with
c(−y) = c(y); its **choice labeling** is

    ρ_c(y) = y_{c(y)} · (n − c(y) + 1) ∈ ±[n].

The explicit labeling of `docs/realization.md` is the case c(y) = start of the final constant run of y|_{[n]}.

**Lemma 1 (PROVED).** Every choice labeling is antipodal and has no complementary comparable pair.

*Proof.* Antipodality is c(−y) = c(y). If y ≤ y′ and ρ_c(y) = −ρ_c(y′), the magnitudes agree, so c(y) = c(y′) =: c, and the signs
differ, y_c = −y′_c. But y_c ≠ 0 and y ≤ y′ give y′_c = y_c. ∎
*Test:* `test_lemma1_any_choice_labeling_is_nonviolating` (random unconstrained choice functions on the full edge list, m = 3, 4, 5).

In coordinates, the label s·(n−c+1) is the pair (c, s). The reverse-caterpillar domains of `docs/realization.md`
(with −n included) read

    T⁺[0] = {(1,+)},    T⁺[r] = {(c,−) : 1 ≤ c ≤ r} ∪ {(r+1,+)}  (1 ≤ r < n),    T⁺[n] = {(c,−) : 1 ≤ c ≤ n}.

(*Test:* `test_tree_domains_are_the_reverse_caterpillar`, n = 2..8, against `kyfan.realize_explicit.reverse_caterpillar_domains`.)

For y ∉ ±U_0 let a(y) = max{a : v_a ≤ y or −v_a ≤ y}, and a(y) = 0 if no ±v_a lies below y. Concretely: a(y) = 0 if y_m = 0;
if y_m = s ≠ 0, a(y) is the length of the initial run of −s in y_1, y_2, ….

**Lemma 2 (PROVED; gadget criterion).** Let c be a choice function and ρ_c its labeling restricted to the vertices outside ±U_0.
The residual domain of v_r contains T⁺[r] for every r = 0..n **iff** c(y) ∉ [a(y)] for every y ∉ ±U_0. The set of admissible
choices at y is A(y) = (supp(y) ∩ [n]) ∖ [a(y)], and |A(y)| = d_{U_0}(y) ≥ 1.

*Proof.* The residual domain of v_r is D(r) = {ℓ : no y ∉ ±U_0 comparable to v_r has ρ_c(y) = −ℓ}. So the condition is: for every
y ∉ ±U_0 comparable to v_r, with c = c(y), the pair (c, −y_c) is not in T⁺[r]; that is, **not** (c ≤ r and y_c = +1) and **not**
(c = r+1 and y_c = −1).

- y < v_r (r ≥ 1). Then y_i ∈ {0, −1} for i ≤ r and y_i = 0 for r < i ≤ n. So c ≤ r and y_c = −1: both clauses are false.
  Vertices below the flag impose nothing.
- y > v_r (r ≥ 0). Then y_m = +1 and y_1 = … = y_r = −1, so the first clause is false, and the condition is: c ≠ r+1 whenever
  y_{r+1} = −1. Ranging over all r with v_r ≤ y, i.e. r ≤ a(y): c ∉ {r+1 : r ≤ a(y), y_{r+1} = −1} = [a(y)] (y_{a+1} ≠ −1 by
  maximality of a).

A vertex comparable to −v_r is the antipode of one comparable to v_r, and ρ_c, a, A are antipodally symmetric, so the condition at y
for all of ±U_0 is exactly c(y) ∉ [a(y)] when some ±v_a lies below y, and is empty otherwise (if y_m = 0 no vertex of ±U_0 lies
below y, so only the first bullet can apply). Finally |A(y)| = |supp(y) ∩ [n]| − a(y): if y_m = 0 this is rank(y) = d_{U_0}(y); if y_m ≠ 0 the largest vertex
of ±U_0 below y is ±v_{a(y)} of rank a(y)+1 and |supp(y) ∩ [n]| = rank(y) − 1, so again |A(y)| = d_{U_0}(y). It is ≥ 1 because
y ∉ ±U_0. ∎
*Tests:* `test_lemma2_admissible_choices_exact` — for every y ∉ ±U_0 and every c ∈ supp(y) ∩ [n], "(c, y_c) is excluded by some
T⁺[r]" ⇔ c ∉ A(y), and |A(y)| = d_{U_0}(y) ≥ 1; exhaustive for m = 3..6.
`test_lemma2_lemma3_sampled_gadgets_on_the_sphere` recomputes the residual domains from the sphere's edge list (m = 3, 4, 5).

Remarks. (1) Cap vertices choose independently of their equatorial projections: Lemma 1 needs no relation between c(u, ±1) and
c(u, 0). (2) Equatorial vertices (y_m = 0) are completely unconstrained, including the vertices (w_r, 0) that the cap-copying
construction of `docs/realization.md` ties to the label −n. This answers "issue 2" of the brief: on the sphere there is one
independent choice per antipodal class, not one per three vertices.

**Lemma 3 (PROVED from TRUSTED ingredients; Lemma A for containing domains).** Let ρ be a non-violating labeling of all classes
outside a flag U whose residual domains satisfy D(x_r) ⊇ T′[r], where (T′[r])_r are the domains of a binary conflict tree gadget
in the sense of `docs/tree_gadget.md` (for the families below: the reverse caterpillar T⁺).
Then every certificate C contains a term α with supp(α) ⊇ U and α|_{off U} ⊆ ρ.

*Proof.* Let E be the degree-n pseudo-solution of the tree gadget (TRUSTED: `docs/tree_gadget.md`), defined on partial labelings of
≤ n vertices of U with labels in T′, and extend it by E[β] = 0 whenever β uses a label outside T′. Then E[∅] = 1; E vanishes on
every β that is violating for the residual CSP with the larger domains D (a complementary pair is a violation of the tree gadget
too; an out-of-D label is in particular out of T′); and for |β| < n, v ∉ β, Σ_ℓ E[β ∪ {v ↦ ℓ}] = E[β], because the new labels
ℓ ∉ T′(v) contribute 0 and, if β itself uses a label outside T′, both sides are 0. So E is a degree-n pseudo-solution of the
residual CSP with domains D. By the restriction lemma (TRUSTED: `docs/tree_gadget.md`), substituting ρ into C gives a certificate
of that residual CSP with terms {α|_U : α ∈ C, α|_{off U} ⊆ ρ}; it cannot have degree ≤ n, so some surviving term has
|α|_U| = n+1 = |U|. ∎
*Tests:* `test_lemma3_larger_domains_keep_the_pseudo_solution` (random enlargements of the tree domains, n = 2, 3, 4: degree-n dual
stays consistent); `test_lemma2_lemma3_sampled_gadgets_on_the_sphere` (degree-n dual of the actual residual of sampled family
members at random flags is consistent, m = 3, 4, 5); `test_theoremC_on_stored_and_expanded_certificates` (the conclusion itself,
for all 24 × 16 gadgets of the family on S^2, on the two stored certificates and on two expanded ones with terms of size > m).

**Lemma 4 (PROVED; transport).** For every flag F there is a family 𝒢_F of labelings of all classes outside F, each satisfying the
hypothesis of Lemma 3, such that 𝒢_F is a product: every class [y] ∉ F has a set Λ_F(y) of exactly d_F(y) labels, and 𝒢_F is the
set of all labelings ρ with ρ(y) ∈ Λ_F(y) for every y.

*Proof.* For U_0 take the choice labelings with c(y) ∈ A(y) (Lemmas 1, 2): Λ(y) = {y_c (n−c+1) : c ∈ A(y)}, d_{U_0}(y) distinct
labels (distinct magnitudes). The hyperoctahedral group B_m acts simply transitively on flags and preserves order and
antipodes; for F = g U_0 put 𝒢_F = {ρ ∘ g^{-1}}. Non-violation and the residual domains are transported by g, and
d_{gU_0}(gy) = d_{U_0}(y). ∎
*Test:* `test_lemma4_d_is_intrinsic_and_invariant` (random flags, m = 3, 4, 5: `to_pole_coordinates(F)` maps F to U_0, and the number
of distinct admissible labels of every class equals d_F).

## 2. The weighted spread lemma, and why the conjectured form fails for this family

Let D_F be the uniform distribution on 𝒢_F: independent uniform labels ρ(y) ∈ Λ_F(y), one per class.

**Lemma 5 (PROVED; weighted spread).** For every flag F, every set A of classes outside F and every labeling a of A,

    P_{ρ∼D_F}[ρ|_A = a] = Π_{y∈A} 1/d_F(y)  if a(y) ∈ Λ_F(y) for all y ∈ A,  and 0 otherwise;   in particular ≤ Π_{y∈A} 1/d_F(y).

*Proof.* Product distribution with uniform marginals on sets of size d_F(y). ∎
*Test:* `test_lemma5_product_distribution_m3` (m = 3: all 24 flag classes, |𝒢_F| = 16, exact marginal count).

**Proposition 1 (PROVED; the exceptional set).** The classes off F with d_F = 1 (deterministic label under D_F) number exactly
m(m−1): the n classes of rank 1 other than x_1's, and for each a = 0..n−1 the 2(n−a−1)+1 classes of x_{a+1} plus one further
coordinate, other than x_{a+2}. In the cap-copying version on the equator (c(u,±1) = c(u,0), the setting of the family G_m of
`docs/spread_evidence.md`) the admissible choices of an equatorial u are supp(u) minus its initial constant run, this sub-family
lies in G_m, and exactly n(n−1) of the movable vertices have a single choice.

*Proof.* d_F(y) = 1 iff y covers its largest ±F-vertex x_j (or has rank 1 with nothing below, j = 0) and y ∉ ±F. For F = U_0:
y_m = 0 forces rank 1 (n classes); y ≥ v_a maximal with one more coordinate i ∈ [n]∖[a]: i = a+1 must carry +1 (1 choice), i > a+1
either sign (2(n−a−1) choices); sum over a = 0..n−1 is n + n(n−1) = n², total n + n² = m(m−1). The equatorial count is the same
kind of count: u_1 = 0 forces a single nonzero coordinate j ≥ 2 (n−1 classes); an initial constant run of length a ≥ 1 followed
by one coordinate gives 1 class at position a+1 (opposite sign) and 2(n−a−1) at later positions, and summing over a = 1..n−1 gives
(n−1) + (n−1)(n−2); total n(n−1). Membership in G_m is the static condition of `GadgetSampler`, checked by the test. ∎
*Tests:* `test_exceptional_set_is_m_times_m_minus_1` (m = 3..7), `test_equatorial_subfamily_lies_in_G_m` (m = 4, 5, 6; the analysis
script prints m ≤ 8: 6, 12, 20, 30, 42 one-choice vertices among 10, 36, 116, 358, 1086 movable ones).

Consequences for the questions in `docs/spread_evidence.md` (Addendum):

- *Closed form for the alternative labels.* alt(u) = u_c·(n−c+1) for any second admissible coordinate c; no pairwise conditions
  at all (Lemma 1), and more than two labels per vertex wherever d ≥ 3. The exceptional set of this closed form on the equator is
  E_n = {u : exactly one nonzero coordinate beyond the initial constant run}, |E_n| = n(n−1) = 6, 12, 20, 30, 42 for n = 3..7.
  The greedy families of the Addendum have smaller exceptional sets (5, 9, 13, 22, 30) because they also use labels that are not
  of choice type; the truly forced set of G_m has n vertices. None of this matters for the size bound.
- *The conjectured Spread Lemma fails for D_F:* with A = the exceptional set, P[ρ|_A = a] = 1 and |A| = m(m−1) > c·m. Whether some
  other family satisfies the conjectured form is open and no longer needed.
- *Issue 1 of the brief* (a slack of Θ(m²) vertices kills the count C(s,m)β^{s−m}β^{−cm}) is real for this family and is resolved by
  Theorem B: the count has to be done with the weights 1/d_F(y), flag by flag, not with |supp α|.

## 3. The weighted flag count

For a set S of classes and a flag class F ⊆ S put W_F(S) = Π_{y∈S∖F} 1/d_F(y) and T(S) = Σ_{F ⊆ S} W_F(S) (sum over flag classes).

**Lemma 6 (PROVED).** Let F = (x_1 < … < x_m) be a flag with all classes in S, and let W_r = supp(x_r).
(a) *Siblings.* For 1 ≤ r ≤ m−1 let u_{r,i} = x_{r+1} with coordinate i zeroed, for i ∈ W_{r+1}, and let
K_r = {i : u_{r,i} ≠ x_r, [u_{r,i}] ∈ S}, k_r = |K_r|. The classes [u_{r,i}] (all r, all i ∈ K_r) are pairwise distinct and lie in
S∖F; for fixed r the numbers d_F(u_{r,i}), i ∈ K_r, are pairwise distinct elements of {1,…,r}.
(b) *Tops.* For D = 1..m−1 exactly 2^D − 1 classes [x′] ≠ [x_m] of rank m have d_F(x′) ≤ D.
(c) Consequently, if S contains t classes of rank m,

    W_F(S) ≤ g(t) · Π_{r=1}^{m−1} 1/k_r! ,      g(t) = Π_{i=1}^{t−1} 1/⌈log_2(i+1)⌉ .

*Proof.* (a) u := u_{r,i} has rank r and the signs of x_m on its support. It is not x_r by assumption, hence not in F; it is not in
−F, whose vertices have the opposite signs. Two such vertices in the same class are equal (the antipode has opposite signs), and
equal vertices have the same rank r and the same missing coordinate i. No vertex of −F lies below u (signs). x_j ≤ u iff
W_j ⊆ W_{r+1}∖{i} iff j ≤ r and i ∉ W_j. Since i ∈ K_r means i ∈ W_r, let j(i) = min{j : i ∈ W_j} ∈ {1,…,r}; then the largest
F-vertex below u has rank j(i) − 1 and d_F(u) = r − j(i) + 1 ∈ {1,…,r}. The map i ↦ j(i) is injective (one coordinate enters at each
rank).
(b) x′ ∈ {±1}^m has x_j ≤ x′ iff x′ agrees with x_m on W_j, and −x_j ≤ x′ iff it agrees with −x_m on W_j; so d_F(x′) ≤ D iff x′
agrees with x_m or with −x_m on W_{m−D}. There are 2·2^D such vectors (W_{m−D} ≠ ∅), i.e. 2^D classes, one of which is [x_m].
(c) All the classes in (a) and the t−1 other top classes are distinct members of S∖F, and every omitted factor of W_F(S) is ≤ 1.
For fixed r the k_r sibling factors are reciprocals of distinct positive integers, product ≤ 1/k_r!. By (b) the i-th smallest value
of d_F among the other tops is ≥ ⌈log_2(i+1)⌉. ∎
*Test:* `test_lemma6_sibling_and_top_profile` (m = 3, 4, 5: all flags inside 15 random sets S each; distinctness and range of the
sibling values, the 2^D − 1 count, and W_F(S) ≤ the bound of (c) computed from the actual values).

**Theorem B (PROVED).** For every set S of classes, T(S) ≤ 2^m.

*Proof.* Let x^{(1)},…,x^{(t)} be representatives of the rank-m classes in S. Every flag class has a unique representative F whose
top is one of the x^{(j)}, and the flags with top x and all classes in S are exactly the root-to-depth-(m−1) paths of the following
tree 𝒯_x: the root is W_m = [m]; a node W with |W| = r+1 ≥ 2 has the children W∖{i} for those i ∈ W with [x|_{W∖i}] ∈ S. Along the
path of F the node W_{r+1} has exactly k_r + 1 children (the notation of Lemma 6). In any rooted tree,
Σ_{paths to depth m−1} Π_{nodes on the path, above the leaf} 1/(number of children) ≤ 1 (it is the probability that the uniform random
descent reaches depth m−1). Hence, by Lemma 6(c),

    Σ_{F with top x} W_F(S) ≤ g(t) Σ_F Π_r (k_r+1)/k_r! · Π_r 1/(k_r+1) ≤ g(t) · max_k((k+1)/k!)^{m−1} = g(t)·2^{m−1},

since (k+1)/k! takes the values 1, 2, 3/2, 2/3, … . Summing over the t tops, T(S) ≤ t·g(t)·2^{m−1}, and t·g(t) ≤ 2:
it equals 1, 2, 3/2, 1 for t = 1, 2, 3, 4 and decreases from t = 2 on, since (t+1)g(t+1)/(t g(t)) = (t+1)/(t⌈log_2(t+1)⌉) < 1
for t ≥ 2. ∎
*Tests:* `test_theoremB_exhaustive_m3` (all 2^13 sets on S^2: max T(S) = 4 ≤ 8), `test_theoremB_random_and_adversarial_sets`
(m = 4, 5: random unions of flags and a hill-climb).

Remark (the weights matter; 2^m is the right order). Take S = U_0 together with the ⌊m/2⌋ vertices obtained by swapping the
coordinates that enter at ranks (1,2), (3,4), …. S has m + ⌊m/2⌋ classes and contains 2^{⌊m/2⌋} flag classes, and every class of S
has d_F ≤ 1 for each of them, so T(S) ≥ 2^{⌊m/2⌋} (PROVED by inspection; `test_block_example_T_at_least_2_to_half_m`, m = 3..6).
So sup_S T(S) is between 2^{⌊m/2⌋} and 2^m, and no argument through this family can give more than
m!·2^{m−1}/2^{⌊m/2⌋}. VERIFIED values: sup T = 4 for m = 3 (exhaustive); hill-climbing finds 6 (m = 4) and 9 (m = 5), which are
lower bounds on the supremum only.

## 4. The size bound

**Theorem C (PROVED, given the TRUSTED tree theorem and restriction lemma).** Let n ≥ 2, m = n+1. Every F_2 certificate
1 = Σ_{α∈C} 1_α of Tucker on S^n (α violating partial labelings; no bound on |α|) satisfies

    |C| ≥ m!·2^{m−1} / 2^m = m!/2.

*Proof.* Choose a flag class F uniformly at random (probability 1/(m!2^{m−1}) each) and ρ ∼ D_F. By Lemmas 4 and 3 there is always
a term α ∈ C with supp(α) ⊇ F and α|_{off F} ⊆ ρ. Taking expectations and using Lemma 5,

    1 ≤ Σ_{α∈C} (1/(m!2^{m−1})) Σ_{F ⊆ supp α} P_{D_F}[α|_{off F} ⊆ ρ] ≤ (1/(m!2^{m−1})) Σ_{α∈C} T(supp α) ≤ |C|·2^m/(m!2^{m−1})

by Theorem B. ∎

Since m! ≥ (m/e)^m and 3^m = N+1, m!/2 ≥ (N+1)^{log_3 m − 1/ln 3}/2, superpolynomial in N; the tower certificate
(`docs/tower.md`) has N^{O(log log N)} terms, so the minimum size is N^{Θ(log log N)} at every degree. For comparison, Theorem 2
gives m!2^{m−1}/C(D,m) for degree ≤ D, which is larger than m!/2 exactly when C(D, m) < 2^m.

Values: m!/2 = 3, 12, 60, 360, 2520, 20160 for m = 3..8 (N = 26, 80, 242, 728, 2186, 6560). On S^2 the true minimum at degree 3 is
304 (`docs/size_lower_bound.md`), so the constant is far from sharp at small m; the bound is about growth.

*Tests:* `test_theoremC_on_stored_and_expanded_certificates` — on the two stored S^2 certificates and on certificates with terms of
size ≥ 5 obtained from them by the identity 1_α = Σ_ℓ 1_{α ∪ {v ↦ ℓ}} (`test_expanded_certificate_is_a_certificate`):
Lemma A holds for each of the 24 × 16 gadgets, and the three inequalities of the proof hold with exact rational arithmetic
(the normalized survival sum Σ_α Σ_F P_F[α]/#classes is exactly 12 on both stored certificates — their 12 full-flag terms per class —
and 12.8–15.5 on the expanded ones; `analysis/spread_choice_family.out`). `test_size_bound_values`.

## 5. Status

| statement | status |
|---|---|
| Lemma 1 (choice labelings are valid), Lemma 2 (gadget criterion, d_F = number of admissible choices), Lemma 4 (transport) | PROVED for all n; VERIFIED m ≤ 6 (static, exhaustive), m ≤ 5 (sphere edge list, sampled) |
| Lemma 3 (Lemma A when the residual domains contain tree domains) | PROVED from the tree theorem and the restriction lemma (TRUSTED, `docs/tree_gadget.md`); VERIFIED on S^2 certificates and on sampled residual duals m ≤ 5 |
| Lemma 5 (weighted spread: P[ρ|_A = a] ≤ Π 1/d_F(y)), Proposition 1 (exceptional set = m(m−1) classes) | PROVED; VERIFIED m ≤ 7 |
| Lemma 6, Theorem B (T(S) ≤ 2^m) | PROVED; VERIFIED exhaustively for m = 3, sampled for m = 4, 5 |
| Theorem C (size ≥ m!/2 at every degree) | PROVED (given the TRUSTED ingredients) |
| Spread Lemma in the conjectured form β^{|A| − cm} | false for the family D_F of this note; open for other families; not needed |
| sup_S T(S) | between 2^{⌊m/2⌋} and 2^m (PROVED); = 4 for m = 3 (VERIFIED); exact growth rate open |

What remains open: the right constant K in size ≥ m!·2^{m−1}/K^m (here K = 2; at degree exactly n+1 Theorem 1 gives the full count), i.e.
whether high degree can save a factor c^m at all; this needs either a family with fewer low-d vertices (the non-choice labels found by
the greedy search are a start) or a sharper count than Lemma 6(c), which ignores every vertex of S that is not a sibling or a top.
