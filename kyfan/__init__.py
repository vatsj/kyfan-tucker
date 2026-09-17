"""Nullstellensatz / Sherali-Adams degree of Tucker and Ky Fan on the signed-subset complex.

Conventions (docs/AGENT_BRIEF.md §1): vertices are nonzero vectors in {-1,0,1}^m; free vertices have first
nonzero coordinate +1; a labeling is a list L of labels on free vertices, extended by lambda(-x) = -lambda(x);
a *partial labeling* alpha is a sorted tuple of (free index, label); degree d = partial labelings of <= d vertices.
"""
from .complex import SignedComplex
from .labels import label_set, posalt, negalt, complementary, pos_alt_count, complementary_edges, random_labeling
from .violating import violating_pairs, is_violating, partial_labelings, restrict
from . import group, dual, primal, linalg, tower, lower
