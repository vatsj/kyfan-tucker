# Ky Fan over F_2 has degree exactly n+1: the separating functional

Code: `kyfan/lower.py` (`box`, `E`, `check`). Tests: `tests/test_tower.py::test_kyfan_lower_functional` (m = 3, 4, 5:
E(A_+) = 1, E(1) = 0, and E(1_α) = 0 on 1,000 sampled size-n partial labelings; Claim 1 below proves the last for all α) and
`::test_kyfan_lower_needs_collision` (the negation-closed box gives E(A_+) = 0). Fast tier, < 5 s.

**Theorem.** For every n ≥ 1 and k ≥ n+1 the F_2 Nullstellensatz degree (Sherali–Adams level, in the one-hot quotient) of Ky Fan
on S^n — the identity A_+ + 1 = Σ c_α 1_α over violating partial labelings α — is exactly n+1.

The upper bound is the tower certificate (`docs/tower.md`). For the lower bound we exhibit a linear functional E on functions of
the labeling that vanishes on every function of at most n vertices but has E(A_+ + 1) = 1; then A_+ + 1 is not in the span of the
violating indicators of size ≤ n.

## The box

Fix a top simplex σ_0 = { x_1 < x_2 < … < x_{n+1} } (x_i of rank i). The n+1 vertices have distinct free representatives
(x_i and x_j are comparable, and a vertex is never comparable to the antipode of another vertex of the same chain: x_i ≤ −x_j
would force x_i's nonzero coordinates to agree with both x_j and −x_j). Let

    S_1 = {+1, +2},   S_i = {+i, −i} for i = 2, …, n+1,

and let B be the set of labelings λ with λ(x_i) ∈ S_i for all i and every other free vertex carrying a fixed label (1 in the
code). Define E(f) = Σ_{λ ∈ B} f(λ) mod 2.

**Claim 1.** E(f) = 0 for every function f that depends on the labels of at most n free vertices (in particular E(1) = 0 and
E(1_α) = 0 for every partial labeling α with |α| ≤ n).

*Proof.* f misses at least one of the n+1 representatives of σ_0, say x_i. Group the elements of B into pairs that differ only in
λ(x_i) (|S_i| = 2); f takes the same value on both elements of a pair, so the sum is even. ∎

**Claim 2.** E(A_+) = 1.

*Proof.* A_+ = Σ_σ [λ|_σ posalt] over top simplices σ. If the free-representative set of σ is not exactly that of σ_0, then
[λ|_σ posalt] depends on at most n of the box variables and E kills it by Claim 1. If it is, then σ = { ε_1 x_1, …, ε_{n+1} x_{n+1} }
with signs ε_i; consecutive elements must be comparable, and ε_i x_i ≤ ε_{i+1} x_{i+1} forces ε_i x_i and ε_{i+1} x_{i+1} to agree
on supp(x_i), i.e. ε_i = ε_{i+1}. So σ = ±σ_0, and since −σ_0 is posalt iff σ_0 is negalt,

    E(A_+) = #{ λ ∈ B : λ|_{σ_0} posalt } + #{ λ ∈ B : λ|_{σ_0} negalt }   (mod 2).

In the box, λ(x_1) ∈ {1, 2} and |λ(x_i)| = i for i ≥ 2, so the magnitudes are distinct only when λ(x_1) = 1 (the value 2 collides
with |λ(x_2)|). Then, sorted by magnitude, the tuple is (1, ±2, …, ±(n+1)); it is posalt for exactly one choice of signs
(λ(x_i) = (−1)^{i−1} i) and never negalt (the smallest label is +1). Hence E(A_+) = 1. ∎

So E(A_+ + 1) = 1 while E vanishes on the span of all violating indicators of size ≤ n: no degree-n certificate exists. ∎

## Remark: the collision is essential

With the negation-closed box S_1 = {+1, −1} (all S_i = {±i}), the magnitude tuple is always (1, 2, …, n+1) and the box contains
exactly one posalt and exactly one negalt labeling of σ_0, so E(A_+) = 0 and the functional proves nothing. The collision at
magnitude 2 removes the negalt contribution and leaves an odd count. `test_kyfan_lower_needs_collision` checks both values on S^3.

## Relation to Tucker

The same functional cannot bound Tucker's degree: Tucker's certificate is for the constant 1, and E kills constants. The Tucker
lower bound is obtained instead from the restriction lemma and the tree gadget (`docs/tree_gadget.md`, `docs/realization.md`),
which supply a pseudo-solution with E[∅] = 1.
