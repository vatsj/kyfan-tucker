# The tower certificate: degree n+1 upper bound for Ky Fan and Tucker over F_2

Code: `kyfan/tower.py` (`g`, `local_lemma_violations`, `check_identity`). Tests: `tests/test_tower.py::test_local_lemma_exhaustive`
(g_n = 0 on every non-complementary tuple, exhaustively for n ≤ 3, k ≤ 5) and `::test_telescoped_identity` (the identity below on
200 random labelings of S^2, S^3, S^4 with k = m+1; 29 / 221 / 2,141 local-lemma instances). Both fast tier, < 15 s.

## Setting

Vertices of the complex K_m are the nonzero vectors of {−1, 0, +1}^m, x ≤ y iff x_i ∈ {0, y_i} for all i; simplices are chains;
K_m triangulates S^n, n = m−1. A labeling λ: V → ±[k] is antipodal, λ(−x) = −λ(x). A tuple of labels is *positively alternating*
(posalt) if its magnitudes are distinct and, sorted by magnitude, the signs read +, −, +, −, …; *negatively alternating* (negalt) if
they read −, +, −, …. The empty tuple is both. A_+^{(m)}(λ) = number of top simplices (maximal chains, m vertices) of K_m whose
label tuple is posalt. Everything is mod 2.

## The local lemma

For an (n+1)-tuple of labels l = (l_0, …, l_n) put

    g_n(l) = Σ_{i=0}^{n} [l ∖ l_i is posalt] + [l is posalt] + [l is negalt]   (mod 2),

where l ∖ l_i is the n-tuple with the i-th entry deleted. (Whether a tuple is posalt depends only on its multiset, so the order of
the entries is irrelevant.)

**Lemma (local lemma).** If l contains no complementary pair (no i ≠ j with l_i = −l_j), then g_n(l) = 0.

*Proof.* Two cases according to whether the magnitudes |l_0|, …, |l_n| are distinct.

*Case A: all n+1 magnitudes distinct.* Sort l by magnitude and let s_0, …, s_n be the signs in that order; let a_j = (−1)^j be
the sign pattern of a posalt tuple and P = { j : s_j = a_j } the set of positions where l agrees with it. Then l is posalt iff
P = {0, …, n} and negalt iff P = ∅. Deleting the entry at position i shifts the entries after i one place down, so l ∖ l_i is
posalt iff s_j = a_j for j < i and s_j = a_{j−1} = −a_j for j > i, i.e. iff P ∩ [0, i) = [0, i) and P ∩ (i, n] = ∅; the
position i itself is unconstrained. Hence the number of i with l ∖ l_i posalt is

  - 0 if P is not an initial segment {0, …, p−1} for some p ∈ {0, …, n+1};
  - otherwise [p ≤ n] + [p ≥ 1] (the choices i = p, for which P = [0, i), and i = p−1, for which P = [0, i]).

So: if P is not an initial segment, all three kinds of term vanish and g_n(l) = 0. If P = ∅ (l negalt): one deletion (i = 0) plus
the negalt term, g_n(l) = 2 = 0. If P = {0, …, n} (l posalt): one deletion (i = n) plus the posalt term, g_n(l) = 0. If
P = {0, …, p−1} with 1 ≤ p ≤ n: two deletions, l neither posalt nor negalt, g_n(l) = 0. (This case does not use the hypothesis.)

*Case B: some magnitude repeats.* Then l is neither posalt nor negalt, and l ∖ l_i can have distinct magnitudes only if exactly
one magnitude μ is repeated, exactly twice, at positions i and i', and the deleted entry is one of these two; otherwise every term
is 0. In that situation the hypothesis gives l_i = l_{i'} (same magnitude, not complementary, so the same sign), hence
l ∖ l_i and l ∖ l_{i'} are the same multiset and contribute equal terms: g_n(l) = 2·[l ∖ l_i posalt] = 0. ∎

The hypothesis is used only in Case B: if l_i = −l_{i'} is the unique repeated magnitude, l ∖ l_i and l ∖ l_{i'} differ in the
sign at μ and at most one of them is posalt, so g_n(l) can be 1. That is exactly where the certificate's terms come from.

Note g_0 ≡ 0 (the empty tuple is posalt, and a single label is either posalt or negalt), which is consistent with the
convention A_+^{(0)} = 1 (K_0 has the empty chain as its only top simplex).

## The hemisphere handshake

Let H = H^{(m)} = { x ∈ K_m : x_m ≠ −1 } be the closed upper hemisphere; its boundary is the equator { x_m = 0 } = K_{m−1}.
A maximal chain of K_m has top element with all coordinates nonzero, and every element of the chain has x_m ∈ {0, x_m(top)}; so
the top simplices of K_m are partitioned into those in H (top's x_m = +1) and their antipodes in −H. Since negating all labels
swaps posalt and negalt,

    A_+^{(m)} = Σ_{σ ∈ top(H)} ( [λ|_σ posalt] + [λ|_σ negalt] ).                                   (1)

Every chain τ of m−1 vertices in H lies in exactly two top simplices of H if some vertex of τ has x_m = +1 (the two completions
at the missing rank take their coordinates from the element above, so they stay in H; at the top rank both u_{m−1} ± e_i are in H
because x_m is already +1), and in exactly one if τ is equatorial (τ is a maximal chain u_1 < … < u_{m−1} of K_{m−1}, and its only
completion in H is by u_{m−1} + e_m). Therefore, mod 2,

    Σ_{σ ∈ top(H)} Σ_{i} [λ|_{σ ∖ σ_i} posalt] = Σ_{τ equatorial} [λ|_τ posalt] = A_+^{(m−1)}(λ|_{K_{m−1}}).            (2)

Adding (1) and (2):

    Σ_{σ ∈ top(H^{(m)})} g_{m−1}(λ|_σ) = A_+^{(m)} + A_+^{(m−1)}   (mod 2).                                            (3)

## Telescoping and the certificate

Identify K_j with the sub-complex { x_{j+1} = … = x_m = 0 } of K_m (pad with zeros); then λ restricts to a labeling of every K_j,
and (3) holds at every level j = 1, …, m. Summing and using A_+^{(0)} = 1:

    A_+^{(m)} + 1 = Σ_{j=1}^{m} Σ_{σ ∈ top(H^{(j)})} g_{j−1}(λ|_σ)   (mod 2), as functions of the labeling λ.

(The j = 1 term is identically zero since g_0 ≡ 0; A_+^{(1)} = 1 because exactly one of ±e_1 carries a positive label.)

Each summand g_{j−1}(λ|_σ) is a function of the j ≤ m labels on σ, so in the one-hot quotient it is a sum of indicators
1_{σ ↦ l} over tuples l with g_{j−1}(l) = 1; by the local lemma every such l contains a complementary pair, so every indicator is
a violating partial labeling of size j ≤ m = n+1. Hence:

**Theorem (tower).** Over F_2, A_+ + 1 is a sum of indicators of violating partial labelings of size ≤ n+1: Ky Fan on S^n has a
Nullstellensatz certificate of degree n+1 for every k. When k ≤ n no top simplex can be alternating, A_+ ≡ 0 syntactically, and the
same identity reads 1 = Σ (violating indicators): Tucker on S^n has a degree-(n+1) certificate.

The number of local-lemma instances is Σ_{j=1}^{m} |top(H^{(j)})| = Σ_j j!·2^{j−1} (29, 221, 2,141 for m = 3, 4, 5), each expanding
to at most (2k)^j monomials; see `docs/size_lower_bound.md` for the comparison with the lower bounds.

## What the tests check

- `test_local_lemma_exhaustive`: for n = 0..3 and k = 1..5, g_n(l) = 0 for every non-complementary (n+1)-tuple over ±[k]
  (the lemma is proved above for all n, k; the test is a sanity check of the code's `posalt`/`negalt` conventions).
- `test_telescoped_identity`: the identity holds on 200 random labelings of S^2, S^3, S^4 (k = m+1), and the instance counts are
  29 / 221 / 2,141.
