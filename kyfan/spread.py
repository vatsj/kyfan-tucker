"""The coordinate-choice gadget family and the weighted flag count (docs/spread_lemma.md).

Sphere S^n = SignedComplex(m), m = n+1, labels +-[n].  Pole flag U_0 = (v_0 < v_1 < ... < v_n) with v_0 = e_m and
v_r = e_m - (e_1 + ... + e_r).

A *choice labeling* picks for every vertex y != +-e_m a coordinate c(y) in supp(y) & [n], with c(-y) = c(y), and sets

    rho_c(y) = y_c * (n - c + 1).

Every choice labeling is free of complementary edges (two comparable vertices agree on every coordinate of the
smaller one).  It is a gadget for U_0 -- its residual domains on U_0 contain the reverse-caterpillar domains -- iff
c(y) avoids [a(y)] for every y outside +-U_0, where a(y) = max{a : v_a <= y or -v_a <= y} (0 if there is none).
The number of admissible choices is

    d_F(y) = rank(y) - rank(largest vertex of +-F below y)        (rank of "none" = 0),

an expression that makes sense for every flag F and is invariant under the hyperoctahedral group.

Functions here are deliberately naive (they are the reference implementation the tests check the lemmas against).
"""
import itertools
from fractions import Fraction
from math import factorial

from .complex import SignedComplex, leq, neg, rank


# ------------------------------------------------------------------ flags
def pole_flag(m):
    """U_0 = (e_m, e_m - e_1, e_m - e_1 - e_2, ...) as a tuple of vertices, ranks 1..m."""
    return tuple(tuple([-1] * r + [0] * (m - 1 - r) + [1]) for r in range(m))


def is_flag(F):
    return all(rank(x) == r + 1 for r, x in enumerate(F)) and all(leq(F[r], F[r + 1]) for r in range(len(F) - 1))


def to_pole_coordinates(F):
    """The signed permutation g^{-1} with g(U_0) = F, as a function on vertices: y -> y' with y in gadget position
    relative to F iff y' is in that position relative to U_0."""
    m = len(F)
    order, sign = [], []
    prev = (0,) * m
    for x in F:
        (i,) = [k for k in range(m) if x[k] and not prev[k]]
        order.append(i)
        sign.append(x[i])
        prev = x

    def pull(y):
        out = [0] * m
        out[m - 1] = sign[0] * y[order[0]]
        for r in range(1, m):
            out[r - 1] = -sign[r] * y[order[r]]
        return tuple(out)

    return pull


def flag_classes(cx):
    """One flag per antipodal pair {F, -F}, with its set of free-vertex indices."""
    seen, out = set(), []
    for F in cx.top:
        key = frozenset(cx.rep(x)[0] for x in F)
        if key not in seen:
            seen.add(key)
            out.append((F, key))
    return out


# ------------------------------------------------------------------ the family
def a_pole(y):
    """max{a >= 0 : v_a <= y or -v_a <= y} for the pole flag, or None if no +-v_a lies below y."""
    m = len(y)
    s = y[m - 1]
    if s == 0:
        return None
    a = 0
    while a < m - 1 and y[a] == -s:
        a += 1
    return a


def allowed_choices(y):
    """Admissible coordinates c (1-based, in [n]) for a vertex y outside +-U_0 of the pole flag."""
    m = len(y)
    a = a_pole(y) or 0
    return [c for c in range(a + 1, m) if y[c - 1]]


def choice_label(y, c):
    return y[c - 1] * (len(y) - c)            # n - c + 1 with n = len(y) - 1


def d_flag(F, y):
    """d_F(y) = rank(y) - rank(largest vertex of +-F below y); 0 iff y in +-F."""
    below = [rank(x) for x in F if leq(x, y) or leq(neg(x), y)]
    return rank(y) - (max(below) if below else 0)


def gadget_choices(cx, F):
    """For the flag F: dict free index -> list of admissible labels of that free vertex (off +-F)."""
    pull = to_pole_coordinates(F)
    Fidx = {cx.rep(x)[0] for x in F}
    out = {}
    for i, v in enumerate(cx.free):
        if i in Fidx:
            continue
        yp = pull(v)
        out[i] = [choice_label(yp, c) for c in allowed_choices(yp)]
    return out


def sample_gadget(cx, F, rng):
    """A uniform member of D_F: dict free index -> label, defined off F."""
    return {i: int(ls[rng.integers(len(ls))]) for i, ls in gadget_choices(cx, F).items()}


def all_gadgets(cx, F):
    """Every member of supp D_F (only sensible for m = 3)."""
    ch = gadget_choices(cx, F)
    keys = sorted(ch)
    for ls in itertools.product(*(ch[k] for k in keys)):
        yield dict(zip(keys, ls))


# ------------------------------------------------------------------ independent checks on the sphere
def is_nonviolating(cx, rho):
    for x, y in cx.edges:
        (i, si), (j, sj) = cx.rep(x), cx.rep(y)
        if i in rho and j in rho and si * rho[i] == -sj * rho[j]:
            return False
    return True


def residual_domains(cx, rho, F):
    """Residual domains of the vertices of F (in chain order), from the sphere's own edge list:
    D(x) = labels l such that no labelled neighbour of x carries -l."""
    n = cx.m - 1
    labels = {s * k for k in range(1, n + 1) for s in (1, -1)}
    dom = {x: set(labels) for x in F}
    for x, y in cx.edges:
        for a, b in ((x, y), (y, x)):
            if a in dom:
                j, sj = cx.rep(b)
                if j in rho:
                    dom[a].discard(-sj * rho[j])
    return [dom[x] for x in F]


def tree_domains(n):
    """T+[r] = reverse-caterpillar domain of v_r, r = 0..n, as label sets (docs/realization.md)."""
    out = [{n}]
    for r in range(1, n + 1):
        out.append(({n - r} if r < n else set()) | {-k for k in range(n - r + 1, n + 1)})
    return out


# ------------------------------------------------------------------ the weighted flag count
def flag_weight(cx, F, Fidx, S):
    """W_F(S) = prod over classes y in S \\ F of 1/d_F(y)   (S a set of free indices containing Fidx)."""
    w = Fraction(1)
    for i in S:
        if i not in Fidx:
            w /= d_flag(F, cx.free[i])
    return w


def T_weight(cx, S, classes=None):
    """T(S) = sum over flag classes F inside S of W_F(S)."""
    S = frozenset(S)
    classes = flag_classes(cx) if classes is None else classes
    return sum((flag_weight(cx, F, Fidx, S) for F, Fidx in classes if Fidx <= S), Fraction(0))


def T_bound(m):
    """Proved upper bound on T(S) for every S (Theorem B): 2^m."""
    return 2 ** m


def size_bound(m):
    """Proved lower bound on the number of terms of any certificate on S^{m-1} (Theorem C): m!/2."""
    return Fraction(factorial(m), 2)


def descent_profile(cx, F, S):
    """For a flag F inside S (free indices): the numbers k_r, r = 1..m-1, of 'siblings' u = x_{r+1} minus one
    coordinate, u != x_r, whose class lies in S, together with their d_F values (proof of Lemma 4)."""
    m = cx.m
    out = []
    for r in range(1, m):
        top = F[r]                                    # rank r+1
        ds = []
        for k in range(m):
            if top[k]:
                u = top[:k] + (0,) + top[k + 1:]
                if u != F[r - 1] and cx.rep(u)[0] in S:
                    ds.append(d_flag(F, u))
        out.append(ds)
    return out
