# Statement (ii): an explicit realization, and Tucker's F_2 degree on S^n = n+1 for all n

Code: `kyfan/realize_explicit.py` (`explicit_label`, `explicit_Leq`, `reverse_caterpillar_domains`, `check`). Tests:
`tests/test_realize_explicit.py` — equatorial validity, exact reverse-caterpillar residual domains, and a consistent, unique
degree-(m−1) residual dual for m = 3..8 (m = 8: 95,901 unknowns, sparse solver, ~3 s; fast tier), plus validity and exact
domains alone for m = 9, 10 (`test_explicit_valid_and_domains_m9_m10`: 21 s, slow tier; 3.5 min, veryslow tier); the
full-sphere `verify` for m = 4, 5, 6.

## The labeling

For a signed subset u ∈ {-1,0,+1}^n \ {0} of [n], let p be the last nonzero index and let q ≤ p be the start of the
final maximal run of equal signs (u_q = u_{q+1} = ... = u_p = s ≠ 0, and either q = 1 or u_{q-1} ≠ s). Define

    λ(u) = s · (n − q + 1).

λ is antipodal (λ(−u) = −λ(u)), takes values in ±[n], and equals ∓n exactly on the chain σ = {w_r = −(e_1+…+e_r)}
and its antipode (q = 1 with a single run means u is an initial segment of constant sign).

## Validity (no complementary comparable pair)

If u ≤ v and λ(u) = −λ(v), then q(u) = q(v) =: q and the signs of the final runs are opposite. But v ⊇ u, so
v_q = u_q = s(u); while v's final run starts at q with sign s(v) = −s(u), i.e. v_q = −s(u). Contradiction. ∎

Hence the pullback of λ to S^n minus the pole chain U = e_{n+1} ∪ lift(σ) is a non-violating partial labeling ρ:
cap edges are lifted equatorial edges, cap–equator edges y ≥ u are lifted comparabilities, and the two caps are
incomparable.

## The residual domains are the reverse caterpillar

Forb(r) = { −λ(u) : u ≤ w_r } ∪ { −λ(u) : u ≥ w_r, u ∉ σ }, D(r) = ±[n] \ Forb(r).

- down(w_r) = { −e_J : ∅ ≠ J ⊆ [r] }: the final run of −e_J starts at some q ∈ [1, r], and every q ∈ [1, r] occurs
  (J = {q,…,p}). So λ(down(w_r)) = { −n, −(n−1), …, −(n−r+1) }, forbidding +n, …, +(n−r+1).
- up(w_r) \ σ = { u : u_1 = … = u_r = −1, u ∉ σ }: let r' ≥ r be the maximal initial run of −1's; u ∉ σ means some
  nonzero coordinate lies beyond r'. If the last coordinate is +1, the final run starts at q ∈ [r+1, n], all occurring:
  labels +(n−r), …, +1. If it is −1, the final run is separated from the initial run by a break, so q ∈ [r+2, n]:
  labels −(n−r−1), …, −1. Forbidden: −1, …, −(n−r) and +1, …, +(n−r−1).
- Together: D(r) = { +(n−r) } ∪ { −(n−r+1), …, −n }   (for r = n: { −1, …, −n }).
- Pole: its fixed neighbours are all cap vertices u + e_{n+1}, u ∉ σ, carrying λ(u); −n occurs only on σ, +n on −σ,
  and every label in ±[n−1] occurs (q ∈ [2, n], both signs). So D(0) = {+n}.

Dropping −n from the leaf domains, the pick-sets D(1) = {+(n−1)}, D(2) = {+(n−2), −(n−1)}, …,
D(n−1) = {+1, −2, …, −(n−1)}, D(n) = {−1, …, −(n−1)} are precisely the leaf paths of the binary conflict tree
(n−1; x, (n−2; x, (… (1; x, x)))) — the caterpillar with magnitudes in reverse order — with the pole {+n}.

## Conclusion

By the tree-gadget theorem (statement (i), `docs/tree_gadget.md`), the residual has a degree-n pseudo-solution over
F_2, so by the restriction lemma no degree-n certificate of Tucker on S^n exists. With the tower upper bound
(`docs/tower.md`):

    Tucker's F_2 Nullstellensatz degree (Sherali–Adams level) on S^n is exactly n+1, for every n ≥ 2.

Checked mechanically: validity and the exact residual domains for m = 3..10 (S^2..S^9), and consistency + uniqueness of
the residual degree-n dual for m ≤ 8 (sparse solver). For n = 2, 3 the degree is also confirmed by the direct SA dual on the
sphere (`tests/test_nsdeg.py`, `tests/test_sparse.py`), and for n = 4..7 by the searched gadgets (`tests/test_gadget.py`).
