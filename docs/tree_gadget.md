# Statement (i): the tree gadget has F_2 degree exactly n+1, for every n

Code: `kyfan/tree_rule.py` (`rule_E`, `rule_patterns`, `unique_E`, `all_tree_shapes`), `kyfan/abstract_gadget.py`. Tests:
`tests/test_rule_vs_solver.py` (the closed-form E equals the solver's unique pseudo-solution and has support 3^n, for every
unordered binary tree shape with n ≤ 6 leaves — 1, 2, 3, 6 shapes; fast tier — and all 11 shapes with n = 7, 36,409 unknowns
each, sparse solver — slow tier, 11 s) and `tests/test_extension_lemma.py` (Lemma (★) exhaustively for the same shapes: n ≤ 6
fast, n = 7 slow, 24 s), `tests/test_static_patterns.py` (Lemmas 1 and 2 of the written proof; fast). Since the rule and the lemma are label-agnostic and mirror-symmetric, one representative per unordered
shape (magnitudes in preorder) covers every binary conflict tree with n ≤ 7 leaves.

Together with statement (ii) (`docs/realization.md`) and the restriction lemma below this gives Tucker's F_2 degree on S^n = n+1
for all n.

## Setup

A *binary conflict tree* T has n leaves and n−1 internal nodes carrying distinct magnitudes from {1,…,n−1}. The gadget CSP has n+1
variables: the pole p with domain {+n}, and each leaf ℓ with domain
  D(ℓ) = { s·m(v) : v an ancestor of ℓ, s = +1 if ℓ is in v's left subtree, −1 if right } ∪ {−n}.
Every pair of variables is constrained: no complementary labels. "Leaf ℓ picks v" means λ(ℓ) = ±m(v) with the sign of ℓ's side;
"ℓ takes P" means λ(ℓ) = −n.
Degree = Sherali–Adams level over F_2 (size of the partial assignments used). UNSAT and the upper bound n+1 are immediate (n+1
variables). The claim is the lower bound: a degree-n pseudo-solution E exists.

The construction and proof use only that the magnitudes on internal nodes are distinct and that each node splits the signs of its two
subtrees; they do not depend on which magnitude sits on which node or on the sign convention. This label-agnosticism is what
`docs/size_lower_bound.md` uses to transport gadgets.

## The pseudo-solution (explicit)

For S ⊆ leaves call an internal node v *full* if every leaf below v lies in S. (v full ⇒ all descendants full.) A *maximal full node*
is a full node whose parent is not full; an S-leaf whose parent is not full counts as a maximal full node too.

**Pattern rule.** An assignment of picks to the leaves of S ⊊ leaves is an *S-pattern* iff it is produced by this bottom-up process,
and R_S is the indicator of S-patterns:
- tokens start at the S-leaves; a set of tokens moving together is a *block*;
- at a non-full internal node, at most one block arrives (from the side whose child is full); if one arrives it takes that node (all
  its leaves pick it) and stops;
- at a full internal node v, one block arrives from each side (L from the left, R from the right — a block always escapes a full
  node, see below), and exactly one of three things happens: L takes v and R escapes upward; R takes v and L escapes; both escape
  and merge into one block.

Because a full node always emits, the root is full iff S = all leaves, in which case the two blocks at the root have nowhere to go:
R_{all} = ∅. For S ⊊ leaves the number of S-patterns is 3^{#full internal nodes}. Each S-pattern is non-violating (the only possible
complementary pair is "both blocks take v", which is excluded).

**Definition of E** on partial assignments α of ≤ n variables. Let S' = leaves of α with a pick, S_P = leaves of α with label −n.
- E[α] = 0 if α is violating (domain violation, complementary pair, or −n together with the pole);
- E[α] = 0 if S_P ≠ ∅ and |α| < n;
- otherwise E[α] = R_{S'}(α|_{S'}), the pole being ignored. (When |α| = n and p ∉ α, this drops the P-leaves; when S' = all leaves
  it is 0.)
So E[∅] = 1, E vanishes on violating α, and |supp E| = 3^n.

## Patterns, statically

The proofs below use a static description of the S-patterns. Notation: L(v) is the set of leaves at or below a node v. Fix
S ⊊ leaves; a node v is *full* if L(v) ⊆ S (a leaf is full iff it lies in S). The root is not full. For x ∈ S let F(x) be the
highest full node on the path from x to the root and top(x) the parent of F(x). A *pick assignment* on S is a map π sending each
x ∈ S to a proper ancestor π(x); it determines the labels (x gets ±m(π(x)), with the sign of x's side of π(x)), and two leaves carry
complementary labels iff they pick the same node from opposite sides. For a full node c with parent u put

    B(c) = { x ∈ L(c) : π(x) is u or an ancestor of u }      ("the leaves that arrive at u from c").

**Lemma 1 (static form of the rule).** A pick assignment π on S is an S-pattern iff
- (C1) for every x ∈ S, π(x) lies on the path from x to top(x), top(x) included;
- (C2) for every full node c with parent u: either π(x) = u for all x ∈ B(c) ("B(c) takes u"), or π(x) is a proper ancestor of u
  for all x ∈ B(c) ("B(c) passes u");
- (C3) for every internal node u both of whose children are full, the two children's sets do not both take u.

Moreover (C2) and (C3) imply B(c) ≠ ∅ for every full c.

*Proof.* Non-emptiness, by induction on c: for a leaf c ∈ S, B(c) = {c} since π(c) is a proper ancestor of c. For a full internal c
with children a, b (both full), B(a) and B(b) are non-empty by induction, each takes or passes c by (C2), and by (C3) one of them
passes c; its elements pick proper ancestors of c, i.e. lie in B(c).

(⇒) Let π be produced by the process. By induction on the tree, a node emits a block iff it is full: an S-leaf emits its token; a
node with two emitting children emits in each of the three options; a node with exactly one emitting child absorbs that block; a
node with none emits nothing. Since the root emits nothing in a pattern, every token is eventually assigned, and the block emitted
by a full node c is the set of leaves below c that are assigned above c, which is B(c). At the parent u this block either takes u
or escapes upward as a whole, which is (C2); the three options at a node with two full children exclude "both take", which is (C3).
The token of x moves up through the full nodes above x and is absorbed at the latest at the first non-full node, top(x); this is (C1).

(⇐) Let π satisfy (C1)–(C3). Run the process, choosing at each internal node u with two full children a, b the option dictated by
π: by (C2) each of B(a), B(b) takes or passes u, by (C3) not both take, so exactly one of "B(a) takes, B(b) passes", "B(b) takes,
B(a) passes", "both pass" holds, and these are the three options. By induction the block emitted by a full node c is B(c) and the
labels assigned at or below c agree with π: at u as above, the leaves of B(a) ∪ B(b) that pass u are exactly B(u). At a non-full
node u with exactly one full child c, every x ∈ B(c) has F(x) = c and top(x) = u, so π(x) = u by (C1) and the definition of B(c),
which is what the process assigns. The root is not full and emits nothing. So the process outputs π. ∎

Each S-pattern is non-violating: two leaves on opposite sides of u that both pick u would, by (C1), both lie below full children of
u (a leaf below a non-full child c' of u has top(x) at or below c'), contradicting (C3).

## Consistency

We need, for every β with |β| < n and every variable v ∉ β: Σ_x E[β ∪ {v ↦ x}] = E[β] over F_2, the sum over all labels x (labels
outside v's domain are violating and contribute 0). If β is violating, so is every extension and both sides are 0. Let β be
non-violating, S its pick-leaves and Q its P-leaves; write β_S for β restricted to S. Note that Q ≠ ∅ forces p ∉ β.

*v = pole.* The only in-domain label is +n. If Q ≠ ∅ then E[β] = 0 (as |β| < n) and the extension is violating. If Q = ∅ then
E[β ∪ {p ↦ n}] = R_S(β_S) = E[β], since the pole is ignored and S ⊊ leaves.

*v = a leaf ℓ.* The in-domain labels are −n and the picks of ℓ. A pick extension is non-violating iff it creates no complementary
pair among the pick-leaves, and R_{S∪ℓ} vanishes on violating assignments, so the pick extensions contribute
#{ picks t of ℓ : β_S ∪ {ℓ ↦ t} is an (S ∪ ℓ)-pattern } whenever E is given by the third clause of its definition. Three cases:

1. Q ≠ ∅ and |β| ≤ n−2. Then E[β] = 0, and every extension has a P-leaf and size < n, so contributes 0.
2. Q ≠ ∅ and |β| = n−1. Then E[β] = 0, p ∉ β, and S ∪ Q ∪ {ℓ} is the set of all leaves, so |S| ≤ n−2. The extension ℓ ↦ −n
   has size n and contributes R_S(β_S); the pick extensions have size n and contribute the count above. By Lemma (★) the count is
   1 if β_S is an S-pattern and 0 or 2 otherwise, so the total is 0 mod 2 in both cases.
3. Q = ∅, so E[β] = R_S(β_S).
   - If |S| = n−1 (hence p ∉ β and S ∪ ℓ = all leaves): ℓ ↦ −n has size n and contributes R_S(β_S); the pick extensions have
     pick-set all leaves and contribute R_all = 0. Total R_S(β_S) = E[β].
   - If |S| ≤ n−2: ℓ ↦ −n contributes 0 (it has size < n, or else |β| = n−1 forces p ∈ β and the extension is violating). The pick
     extensions contribute the count above, which by Lemma (★) is ≡ R_S(β_S) = E[β] mod 2.

So consistency reduces to:

**Lemma (★) (extension lemma).** For S ⊊ leaves with |S| ≤ n−2, ℓ ∉ S, and any pick assignment β on S:
  #{ picks t of ℓ : β ∪ {ℓ ↦ t} is an (S ∪ ℓ)-pattern } = 1 if β is an S-pattern, and ∈ {0, 2} otherwise.

The proof has a combinatorial core that does not mention trees.

**Lemma 2 (path lemma).** Let j ≥ 0 and let k : {1,…,j} → {1,…,j+1} satisfy k(i) ≥ i. Call t ∈ {1,…,j+1} *admissible* if the
following run succeeds: c_0 = t, and for i = 1,…,j,
- if c_{i−1} = i: require k(i) > i, and set c_i = k(i);
- if c_{i−1} > i: require k(i) ∈ {i, c_{i−1}}, and set c_i = c_{i−1};

and finally c_j = j+1. (Inductively c_{i−1} ≥ i, so the two cases are exhaustive.) Let N = { i : k(i) > i }. If N = ∅ the only
admissible t is j+1. If N ≠ ∅, with a = min N, the set of admissible t is either empty or {a, k(a)}, and a ≠ k(a).

*Proof.* In a successful run, c_i = k(i) for i ∈ N (if c_{i−1} = i this is the first case; otherwise the second case forces
k(i) = c_{i−1} = c_i) and c_i = c_{i−1} for i ∉ N (the first case is excluded since k(i) = i). So define c by c_0 = t,
c_i = k(i) for i ∈ N, c_i = c_{i−1} for i ∉ N. Then t is admissible iff
- (R1) c_{i−1} ≠ i for every i ∉ N,
- (R2) c_{i−1} ∈ {i, k(i)} for every i ∈ N,
- (R3) c_j = j+1:

these are necessary by the two cases, and conversely they make every step of the run legal with the values c_i just defined.
If N = ∅ then c_i = t for all i, (R1) says t ∉ {1,…,j}, so t = j+1, and (R3) holds. If N ≠ ∅ and a = min N, then c_i = t for
i < a, while c_a = k(a) and every later c_i do not depend on t. The conditions that involve t are (R1) for i < a, i.e. t ≥ a, and
(R2) at i = a, i.e. t ∈ {a, k(a)}; both a and k(a) are ≥ a. All remaining conditions are independent of t. Hence either both a and
k(a) are admissible or no t is. ∎

*Proof of Lemma (★).* Let S' = S ∪ {ℓ} ≠ all leaves. Let ℓ = u_0, u_1, …, u_r = root be the path from ℓ to the root and w_i the
sibling of u_{i−1} (a child of u_i). Let j ≥ 0 be largest such that w_1, …, w_j are all full (with respect to S). Then, with
respect to S', the nodes u_0, …, u_j are full and u_{j+1}, …, u_r are not; j < r because the root is not full for S'; and every
node off the path has the same fullness for S and S'. Consequently:
- for x ∈ S ∩ L(w_i) with i ≤ j: top(x) = u_i with respect to S, and top(x) = u_{j+1} with respect to S';
- for every other x ∈ S, top(x) is the same for S and S'.

For 1 ≤ i ≤ j let A_i = B(w_i) = { x ∈ L(w_i) : β(x) ∈ {u_i, …, u_r} }; it is non-empty whenever (C2), (C3) hold below w_i.
Let (K) be the conjunction of: (C1) for all x ∈ S not below any w_i with i ≤ j; (C2) for all full c other than w_1, …, w_j;
(C3) for all u with two full children — fullness taken with respect to S. These conditions read the same for S and for S' and do
not involve ℓ. By Lemma 1:

(a) β is an S-pattern iff (K) holds and β(x) = u_i for all x ∈ A_i and all i ≤ j. (For x ∈ A_i, (C1) says β(x) is at or below
u_i, hence equal to u_i; then (C2) for c = w_i holds; leaves of L(w_i) outside A_i satisfy (C1) trivially.)

(b) β' = β ∪ {ℓ ↦ u_t} is an S'-pattern iff (K) holds and, writing Z_{i−1} = B'(u_{i−1}) for the leaves of S' ∩ L(u_{i−1})
that β' sends to u_i or above (so Z_0 = {ℓ}):
- (P1) t ≤ j+1, and β(x) is at or below u_{j+1} for x ∈ A_i, i ≤ j     [(C1)];
- (P2) each A_i (i ≤ j) takes or passes u_i, and each Z_{i−1} (1 ≤ i ≤ j+1) takes or passes u_i     [(C2)];
- (P3) for i ≤ j, A_i and Z_{i−1} do not both take u_i     [(C3)].

Suppose (P1), (P2) hold. Then every Z_i (0 ≤ i ≤ j) is *constant*: all its elements pick the same node u_{c_i}, i < c_i ≤ j+1. For
i = j this is (P2) at u_{j+1} together with (P1). For i < j, Z_i takes u_{i+1} (constant) or passes it, in which case
Z_i ⊆ Z_{i+1}, constant by induction. Likewise every A_i is constant, with value u_{k(i)}, i ≤ k(i) ≤ j+1: it takes u_i, or
passes it and is contained in Z_i.

So if (K) fails, or some A_i is not constant with a value in {u_i, …, u_{j+1}}, then β is not an S-pattern by (a) and no t gives
an S'-pattern by (b): the count is 0. Otherwise k is defined, β is an S-pattern iff k(i) = i for all i by (a), and we claim that
β' is an S'-pattern iff t is admissible for k in the sense of Lemma 2, with c_i the common pick index of Z_i. Indeed Z_0 = {ℓ}
has c_0 = t. At u_i (1 ≤ i ≤ j), Z_i consists of the elements of Z_{i−1} ∪ A_i that pass u_i. If c_{i−1} = i, then Z_{i−1}
takes u_i, (P3) requires A_i to pass, i.e. k(i) > i, and then Z_i = A_i, c_i = k(i). If c_{i−1} > i, (P3) holds; if k(i) = i then
Z_i = Z_{i−1} and c_i = c_{i−1}; if k(i) > i then Z_i = Z_{i−1} ∪ A_i, which is constant iff k(i) = c_{i−1}. Finally (P1), (P2)
at u_{j+1} say c_j = j+1. Conversely, if t is admissible then every Z_i is constant with value in {u_{i+1}, …, u_{j+1}}, which
gives (P1)–(P3). Lemma 2 now gives the count: 1 if k is the identity, and 0 or 2 otherwise. ∎

Consistency follows, so E is a degree-n pseudo-solution and the gadget's F_2 degree is exactly n+1 for every binary conflict tree.

Checks of the proof's two ingredients, independent of the solver: `tests/test_static_patterns.py` verifies Lemma 1 (the static
conditions (C1)–(C3) give exactly the block-rule patterns, every tree shape with n ≤ 6 and every S) and Lemma 2 (all maps k with
j ≤ 7).

## The restriction lemma

**Lemma (restriction).** Let ρ be a non-violating partial labeling of the free vertices of S^n with unfixed set U, and let the
*residual CSP* on U have domains D(u) = ±[k] ∖ { −ρ(w) : w fixed, w adjacent to u } and the constraint "no complementary pair on an
edge inside U". Every degree-d F_2 certificate of Tucker on S^n restricts to a degree-≤ d certificate of the residual: an identity
1 = Σ c_α 1_{α|_U} over partial labelings of U that are violating for the residual (an out-of-domain label, or a complementary pair).
Consequently, if the residual has a consistent degree-d pseudo-solution (E[∅] = 1, E = 0 on residual-violating partial labelings,
consistent under extension), the sphere has no degree-d certificate.

*Proof.* Substitute ρ into 1 = Σ c_α 1_α. A term 1_α becomes 1_{α|_U} if α agrees with ρ off U and 0 otherwise, so the identity
becomes 1 = Σ_{α ⊇-compatible} c_α 1_{α|_U} on labelings of U. Each surviving α contains a complementary pair; it is not inside ρ
(ρ is non-violating), so either it lies inside U — a complementary edge of the residual — or it joins u ∈ U to a fixed w with
α(u) = −ρ(w) ∉ D(u). Either way α|_U is residual-violating, and |α|_U| ≤ |α| ≤ d. ∎

So Tucker's F_2 degree on S^n is n+1 whenever some equatorial labeling of S^{n−1} realizes a binary conflict tree as the residual on a
top simplex through the pole. This is statement (ii), proved by an explicit closed-form labeling in `docs/realization.md`.

## Independent checks

Two earlier constructions of realizing labelings remain in the repository as independent checks, both superseded by the closed form:
the searched `GADGETS` (n = 4..7, `kyfan.gadget.verify`, `tests/test_gadget.py`) and the deterministic greedy realizer
`kyfan/gadget.py:realize(m)`, which fills the caterpillar target domains by forward-checking and runs with zero backtracks for every
tested m (3..8) (`tests/test_realize.py`). `chain_domains(m, Leq)` reduces the residual check to the pole chain only (no
O(3^m·3^m) sphere edge list); the explicit labeling needs no CSP at all, and its residual dual is the only cost past m ≈ 8.
