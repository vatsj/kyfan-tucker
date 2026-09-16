# The tree gadget has F_2 degree exactly n+1, for every n

Settles statement (i) of `results.md` ("What a proof for all n now needs"). The remaining piece for Tucker on S^n is (ii) realizability.

## Setup

A *binary conflict tree* T has n leaves and n−1 internal nodes carrying distinct magnitudes from {1,…,n−1}. The gadget CSP has n+1 variables: the pole p with domain {+n}, and each leaf ℓ with domain
  D(ℓ) = { s·m(v) : v an ancestor of ℓ, s = +1 if ℓ is in v's left subtree, −1 if right } ∪ {−n}.
Every pair of variables is constrained: no complementary labels. "Leaf ℓ picks v" means λ(ℓ) = ±m(v) with the sign of ℓ's side; "ℓ takes P" means λ(ℓ) = −n.
Degree = Sherali–Adams level over F_2 (size of the partial assignments used). UNSAT and the upper bound n+1 are immediate (n+1 variables). The claim is the lower bound: a degree-n pseudo-solution E exists.

## The pseudo-solution (explicit)

For S ⊆ leaves call an internal node v *full* if every leaf below v lies in S. (v full ⇒ all descendants full.) A *maximal full node* is a full node whose parent is not full; an S-leaf whose parent is not full counts as a maximal full node too.

**Pattern rule.** An assignment of picks to the leaves of S ⊊ leaves is an *S-pattern* iff it is produced by this bottom-up process, and R_S is the indicator of S-patterns:
- tokens start at the S-leaves; a set of tokens moving together is a *block*;
- at a non-full internal node, at most one block arrives (from the side whose child is full); if one arrives it takes that node (all its leaves pick it) and stops;
- at a full internal node v, one block arrives from each side (L from the left, R from the right — a block always escapes a full node, see below), and exactly one of three things happens: L takes v and R escapes upward; R takes v and L escapes; both escape and merge into one block.

Because a full node always emits, the root is full iff S = all leaves, in which case the two blocks at the root have nowhere to go: R_{all} = ∅. For S ⊊ leaves the number of S-patterns is 3^{#full internal nodes}. Each S-pattern is non-violating (the only possible complementary pair is "both blocks take v", which is excluded).

**Definition of E** on partial assignments α of ≤ n variables. Let S' = leaves of α with a pick, S_P = leaves of α with label −n.
- E[α] = 0 if α is violating (domain violation, complementary pair, or −n together with the pole);
- E[α] = 0 if S_P ≠ ∅ and |α| < n;
- otherwise E[α] = R_{S'}(α|_{S'}), the pole being ignored. (When |α| = n and p ∉ α, this drops the P-leaves; when S' = all leaves it is 0.)
So E[∅] = 1, E vanishes on violating α, and |supp E| = 3^n. (Checked against the unique solver solution for every tree shape with n ≤ 6: `test_rule_vs_solver.py`.)

## Consistency

We need, for every non-violating β with |β| < n and every variable v ∉ β: Σ_x E[β ∪ {v ↦ x}] = E[β].
- v = pole: the only extension is p ↦ n; it is violating iff β has a P-leaf, in which case E[β] = 0 too; otherwise both sides equal R_{S'}(β).
- v = leaf ℓ. Write β' = β restricted to its pick-leaves S. The extension ℓ ↦ −n contributes E[β'] when |β| = n−1 and p ∉ β, and 0 otherwise. The pick extensions contribute Σ_x R_{S∪ℓ}(β' ∪ ℓ ↦ x) (0 if S ∪ ℓ = all leaves). Checking the cases (|S| = n−1 with p ∉ β; |S| ≤ n−2 with or without P-leaves) reduces everything to:

**Lemma (★).** For S ⊊ leaves with |S| ≤ n−2, ℓ ∉ S, and any assignment β on S:
  #{ x : β ∪ ℓ ↦ x is an (S ∪ ℓ)-pattern } = 1 if β is an S-pattern, and ∈ {0, 2} otherwise.
(Verified exhaustively for all tree shapes with n ≤ 6: `test_extension_lemma.py`.)

**Proof.** Let ℓ = u_0, u_1, …, u_k = root be ℓ's ancestors and w_i the sibling of u_{i−1} at u_i. Adding ℓ to S changes fullness only on this path: u_1,…,u_j become full, where j is the largest index with w_1,…,w_j all full in S; u_{j+1} exists because S ∪ ℓ ≠ all leaves. Off the path both processes are identical, so any failure there gives 0 extensions. On the path, in the S-process each w_i is a maximal full node whose escaping block B_i must take u_i. So β is an S-pattern iff (off-path OK and) each B_i takes u_i. Write k_i for the index of the node B_i takes in β; validity of β anywhere forces k_i ∈ [i, j+1] (a label below u_i is not an ancestor of w_i's leaves; a label above u_{j+1} cannot be realized in the (S∪ℓ)-process either, since the block reaching u_{j+1} is alone there and takes it). Let N = {i : k_i ≠ i}.

In the (S∪ℓ)-process, ℓ's block travels up the path; at each u_i (i ≤ j) it meets B_i; the choice at u_i is forced by β except for ℓ's own pick x = u_t. Tracking the common target c of the traveling block: at u_i, if c = i the block takes u_i and B_i must escape (k_i > i), after which c := k_i; if c > i then k_i ∈ {i, c} (B_i takes u_i, or merges), and c is unchanged; finally c must equal j+1. Unwinding: t = u_t is valid iff the chain t = t_0, t_1 = k_{t_0}, t_2 = k_{t_1}, … reaches j+1 without hitting a fixed point, and every i ∈ N points to the next chain element above it. Equivalently the fibers of k on N are consecutive intervals F_1 < F_2 < … < F_q with target(F_a) = min F_{a+1} and target(F_q) = j+1, and then the valid t are exactly t = min F_1 and t = min F_2 (or t = j+1 when q = 1). Hence: N = ∅ gives exactly one valid t (t = j+1); N ≠ ∅ gives exactly two if the fiber structure holds and zero otherwise. ∎

Consistency follows, so E is a degree-n pseudo-solution and the gadget's F_2 degree is exactly n+1 for every binary conflict tree.

The closed-form E and the extension lemma are checked against the solver for every tree shape with n ≤ 6 and (against the
sparse solver) for n = 7: `tests/test_rule_vs_solver.py`, `tests/test_extension_lemma.py`; the rule is `kyfan/tree_rule.py`.

## Consequence and statement (ii)

By the restriction lemma (`kyfan/gadget.py`), Tucker's F_2 degree on S^n is n+1 for every n for which some equatorial
labeling of S^{n−1} realizes a binary conflict tree as the residual on a top simplex through the pole. This is statement (ii).

**Progress on (ii): a deterministic construction.** `kyfan/gadget.py:realize(m)` builds an explicit equatorial labeling
realizing the *caterpillar* tree (magnitudes 1..n−1, leaf k → chain rank n+1−k): fix the equatorial top simplex σ to ∓n;
each other free vertex v gets the label set forced by the target domains (for each chain rank r, if v or −v lies in the
star/down-set of σ's rank-r vertex, forbid the corresponding path labels); then fill by greedy forward-checking (smallest
domain, then smallest label). **This runs with zero backtracks for every tested m (3..8)** — i.e. the construction is
forced, not a search — and the resulting residual is a degree-n gadget (checked end-to-end: `residual` on the full sphere
for m ≤ 6, and `chain_domains` + the sparse degree-n dual up to m = 8). `tests/test_realize.py`.

So (ii) is reduced to a clean statement with strong computational evidence: **the greedy caterpillar construction never
backtracks (equivalently, the explicit allowed-set CSP on S^{n−1} is satisfiable) for all n.** A proof of that — e.g.
exhibiting the labeling in closed form, or showing the allowed-set CSP is arc-consistent hence greedily solvable —
finishes Tucker's F_2 degree = n+1 for all n. `chain_domains(m, Leq)` reduces the check to the pole-chain only (no
O(3^m·3^m) sphere edge list), so it scales; the equatorial CSP's O(V^2) `forbidden_pairs` is the remaining bottleneck
past m ≈ 8.
