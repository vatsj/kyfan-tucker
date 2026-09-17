"""Statement (ii), explicit: a closed-form equatorial labeling whose pole-chain residual is a binary conflict tree,
completing Tucker's F_2 degree = n+1 for all n (docs/realization.md).

For a signed subset u in {-1,0,+1}^n \\ {0}, let p be the last nonzero index and q <= p the start of the final maximal
run of equal signs (u_q = ... = u_p = s, and q = 1 or u_{q-1} != s). Then

    lambda(u) = s * (n - q + 1).

lambda is antipodal, takes values in +-[n], equals -+n exactly on the chain sigma = {-(e_1+...+e_r)} and its antipode,
and has no complementary comparable pair (proof in docs/realization.md). Pulled back to S^n it fixes everything off
the pole chain U = e_{n+1} u lift(sigma); the residual domains on U are the reverse caterpillar
    D(0) = {+n} (pole),  D(r) = {+(n-r)} u {-(n-r+1), ..., -n}  (r = 1..n; D(n) = {-1,...,-n}),
i.e. the leaf paths of (n-1; x, (n-2; x, (... (1; x, x)))) plus the pole. No search, no CSP: this is the construction.
"""
from .complex import SignedComplex, leq
from .labels import label_set
from .gadget import chain_domains
from .abstract_gadget import build_dual, unsat
from . import linalg


def explicit_label(u):
    """lambda(u) for a signed subset u given as a tuple in {-1,0,1}^n (not all zero)."""
    n = len(u)
    p = max(i for i in range(n) if u[i])
    s = u[p]
    q = p
    while q > 0 and u[q - 1] == s:
        q -= 1
    return s * (n - q)                       # q is 0-based here: n - (q+1) + 1 = n - q


def explicit_Leq(m):
    """The explicit labeling on the free vertices of the equator S^{m-2} (= SignedComplex(m-1))."""
    cxe = SignedComplex(m - 1)
    return [explicit_label(v) for v in cxe.free]


def reverse_caterpillar_domains(n):
    """Target residual domains D(0..n) produced by the explicit labeling (pole first)."""
    target = {0: {n}}
    for r in range(1, n + 1):
        picks = {n - r} if r < n else set()
        target[r] = picks | {-k for k in range(n - r + 1, n)} | {-n}
    return target


def check(m, solve=True, method="sparse"):
    """Verify the explicit construction for S^{m-1}: equatorial validity, the exact residual domains, and (if solve)
    that the residual degree-(m-1) dual is consistent with a unique pseudo-solution. Returns a dict of results."""
    n = m - 1
    cxe = SignedComplex(m - 1)
    Leq = explicit_Leq(m)
    violations = sum(1 for x, y in cxe.edges if cxe.label(x, Leq) == -cxe.label(y, Leq))
    D = chain_domains(m, Leq)
    target = reverse_caterpillar_domains(n)
    domains_ok = all(D[r] == target[r] for r in range(n + 1))
    out = dict(m=m, equatorial_violations=violations, domains_match=domains_ok, domains=D)
    if solve:
        doms_list = [set(D[r]) for r in range(n + 1)]
        rows, col = build_dual(doms_list, label_set(n), n)
        res = (linalg.gf2_sparse if method == "sparse" else linalg.gf2_dense)(rows, len(col))
        out.update(unsat=unsat(doms_list), unknowns=len(col), rank=res.rank,
                   consistent=res.consistent, unique=res.consistent and res.rank == len(col))
    return out
