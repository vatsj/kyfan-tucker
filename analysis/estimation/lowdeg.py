"""Low-degree deduction-projection estimator for E_x[ prod_a F_a(x) ], x = uniform labeling.

g is kept as a sum of ANOVA components g_S(x_S), |S| <= d (each with zero mean along every axis).  One step:
g <- P_{<=d}( g * F ) for the next local factor F, where P_{<=d} is the orthogonal (L^2) projection onto functions of
degree <= d.  The estimate is the constant component at the end.  With d = number of variables it is exact.
"""
import itertools, sys, time
import numpy as np

def anova(table, axes):
    """ANOVA components of a table over the given axes: dict subset(tuple, sorted) -> table on that subset."""
    comps = {}
    n = len(axes)
    for r in range(n + 1):
        for R in itertools.combinations(range(n), r):
            t = table.mean(axis=tuple(i for i in range(n) if i not in R)) if r < n else table.copy()
            for ax in range(r):
                t = t - t.mean(axis=ax, keepdims=True)
            comps[tuple(axes[i] for i in R)] = t
    return comps

def multiply_project(g, Fc, d, tol=0.0):
    """g: dict S -> table (axes sorted by vertex id); Fc: ANOVA components of the factor. Returns P_{<=d}(g*F)."""
    out = {}
    for S, gs in g.items():
        Sset = set(S)
        for R, fr in Fc.items():
            O = [v for v in R if v in Sset]
            base = sorted(Sset.symmetric_difference(R))
            if len(base) > d:
                continue
            U = sorted(Sset | set(R))
            h = np.einsum(gs, list(S), fr, list(R), U) if (S or R) else np.array(gs * fr)
            # ANOVA over the overlap variables only
            for r in range(len(O) + 1):
                if len(base) + r > d:
                    break
                for Q in itertools.combinations(O, r):
                    drop = tuple(U.index(v) for v in O if v not in Q)
                    t = h.mean(axis=drop) if drop else h
                    axes = [v for v in U if v not in O or v in Q]
                    for v in Q:
                        t = t - t.mean(axis=axes.index(v), keepdims=True)
                    key = tuple(axes)
                    if key in out:
                        out[key] = out[key] + t
                    else:
                        out[key] = np.array(t, dtype=float)
    return out

def estimate(factors, d):
    """factors: list of (axes tuple sorted, table). Returns the degree-d deduction-projection estimate of E[prod]."""
    g = {(): np.array(1.0)}
    for axes, table in factors:
        g = multiply_project(g, anova(table, axes), d)
    return float(g[()])
